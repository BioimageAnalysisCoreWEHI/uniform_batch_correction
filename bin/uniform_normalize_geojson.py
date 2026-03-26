#!/usr/bin/env python3
"""
UniFORM-style normalization for GeoJSON measurements and OME-TIFF/TIFF images.
"""

import argparse
import csv
import json
import math
import re
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import numpy as np
import tifffile
from scipy import sparse
from scipy.signal import correlate


def get_anndata_module():
    try:
        import anndata as ad
    except ModuleNotFoundError as error:
        raise ModuleNotFoundError(
            "AnnData mode requires the 'anndata' package. "
            "Install it in the runtime environment or run non-adata modes without adata dependencies."
        ) from error
    return ad


def parse_args():
    parser = argparse.ArgumentParser(description="Normalize GeoJSON measurements, OME-TIFF intensities, or AnnData features")
    parser.add_argument("--inputs", nargs="+", required=True, help="Input files")
    parser.add_argument("--mode", choices=["geojson", "ome_tiff", "adata"], default="geojson")
    parser.add_argument("--num-bins", type=int, default=1024)
    parser.add_argument("--min-value", type=float, default=1.0)
    parser.add_argument("--exclude-pattern", default=r"^(kronos_|emb_)")
    parser.add_argument("--output-suffix", default="_uniform")
    parser.add_argument("--pixel-output-suffix", default="_unifrom")
    parser.add_argument("--pixel-sample-size", type=int, default=200000)
    parser.add_argument("--adata-group-by", default="image")
    parser.add_argument("--adata-sample-size", type=int, default=200000)
    parser.add_argument("--adata-target", choices=["all", "cell_mean"], default="all")
    parser.add_argument("--adata-filter-column", default="")
    parser.add_argument("--adata-filter-regex", default="")
    parser.add_argument("--qc-dir", default="qc")
    parser.add_argument("--generate-plots", default="true")
    parser.add_argument("--qc-top-n-keys", type=int, default=12)
    parser.add_argument("--qc-max-heatmap-keys", type=int, default=40)
    return parser.parse_args()


def parse_bool(value):
    value_lower = str(value).strip().lower()
    if value_lower in {"1", "true", "yes", "y"}:
        return True
    if value_lower in {"0", "false", "no", "n"}:
        return False
    raise ValueError(f"Invalid boolean value: {value}")


def is_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and np.isfinite(value)


def log_transform(values, min_value=1.0):
    arr = np.asarray(values, dtype=float)
    valid = arr >= min_value
    out = np.full(arr.shape, np.nan, dtype=float)
    out[valid] = np.log(arr[valid])
    return out


def choose_reference(hist_matrix):
    mean_hist = np.mean(hist_matrix, axis=0)
    distances = np.linalg.norm(hist_matrix - mean_hist, axis=1)
    return int(np.argmin(distances))


def compute_fft_shifts(reference_hist, hist_list):
    shifts = []
    for hist in hist_list:
        corr_fft = correlate(hist.flatten(), reference_hist.flatten(), mode="full", method="fft")
        shift_fft = np.argmax(corr_fft) - (len(reference_hist) - 1)
        shifts.append(int(shift_fft))
    return shifts


def output_path_for_image(path, suffix):
    path = Path(path)
    name = path.name
    if name.endswith(".ome.tif"):
        base = name[: -len(".ome.tif")]
        return path.with_name(f"{base}{suffix}.ome.tif")
    if name.endswith(".ome.tiff"):
        base = name[: -len(".ome.tiff")]
        return path.with_name(f"{base}{suffix}.ome.tiff")
    if name.endswith(".tiff"):
        base = name[: -len(".tiff")]
        return path.with_name(f"{base}{suffix}.tiff")
    if name.endswith(".tif"):
        base = name[: -len(".tif")]
        return path.with_name(f"{base}{suffix}.tif")
    return path.with_name(f"{path.stem}{suffix}{path.suffix}")


def output_path_for_adata(path, suffix):
    path = Path(path)
    name = path.name
    if name.endswith(".h5ad"):
        base = name[: -len(".h5ad")]
        return path.with_name(f"{base}{suffix}.h5ad")
    return path.with_name(f"{path.stem}{suffix}{path.suffix}")


# ------------------------------ GEOJSON MODE ------------------------------

def load_geojsons(paths):
    records = []
    for path in paths:
        with open(path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
        sample_id = Path(path).stem
        cell_features = [
            (idx, feat)
            for idx, feat in enumerate(data.get("features", []))
            if feat.get("properties", {}).get("objectType") == "cell"
        ]
        records.append(
            {
                "path": path,
                "sample_id": sample_id,
                "data": data,
                "cell_features": cell_features,
            }
        )
    return records


def collect_measurement_keys(records, exclude_pattern):
    regex = re.compile(exclude_pattern)
    key_counts = defaultdict(int)
    sample_count = len(records)

    for rec in records:
        keys_here = set()
        for _, feat in rec["cell_features"]:
            measurements = feat.get("properties", {}).get("measurements", {})
            for key, value in measurements.items():
                if regex.search(key):
                    continue
                if is_number(value):
                    keys_here.add(key)
        for key in keys_here:
            key_counts[key] += 1

    return sorted([key for key, count in key_counts.items() if count == sample_count])


def normalize_geojson(records, num_bins, min_value, exclude_pattern, output_suffix):
    sample_ids = [rec["sample_id"] for rec in records]
    records_before = load_geojsons([rec["path"] for rec in records])
    keys = collect_measurement_keys(records, exclude_pattern)

    print(f"[geojson] Loaded {len(records)} sample(s)")
    print(f"[geojson] Found {len(keys)} shared numeric key(s) after exclusion pattern: {exclude_pattern}")

    key_diagnostics = {}
    ranked_keys = []

    n_samples = len(records)
    if keys and n_samples >= 2:
        scales_by_key = {}
        for key_idx, key in enumerate(keys, start=1):
            if key_idx == 1 or key_idx % 25 == 0 or key_idx == len(keys):
                print(f"[geojson] Processing key {key_idx}/{len(keys)}: {key}")
            sample_values = []
            global_min = float("inf")
            global_max = float("-inf")

            for rec in records:
                values = []
                for _, feat in rec["cell_features"]:
                    measurements = feat.get("properties", {}).get("measurements", {})
                    value = measurements.get(key)
                    if is_number(value):
                        values.append(float(value))

                log_vals = log_transform(values, min_value=min_value)
                log_vals = log_vals[np.isfinite(log_vals)]
                if log_vals.size == 0:
                    log_vals = np.array([0.0], dtype=float)

                sample_values.append(log_vals)
                global_min = min(global_min, float(np.min(log_vals)))
                global_max = max(global_max, float(np.max(log_vals)))

            if not np.isfinite(global_min) or not np.isfinite(global_max) or global_max <= global_min:
                scales_by_key[key] = [1.0] * n_samples
                continue

            hist_list = []
            for log_vals in sample_values:
                counts, _ = np.histogram(log_vals, bins=num_bins, range=(global_min, global_max))
                hist_list.append(counts.astype(float))

            hist_matrix = np.stack(hist_list, axis=0)
            ref_idx = choose_reference(hist_matrix)
            shifts = compute_fft_shifts(hist_matrix[ref_idx], hist_list)

            increment = (global_max - global_min) / max(1, num_bins - 1)
            scales = [math.exp(-shift * increment) for shift in shifts]
            scales_by_key[key] = scales
            key_diagnostics[key] = {
                "reference_sample": sample_ids[ref_idx],
                "ref_idx": ref_idx,
                "shifts": shifts,
                "scales": scales,
                "global_min": global_min,
                "global_max": global_max,
            }

        for sample_idx, rec in enumerate(records):
            for _, feat in rec["cell_features"]:
                measurements = feat.get("properties", {}).get("measurements", {})
                for key in keys:
                    value = measurements.get(key)
                    if is_number(value):
                        measurements[key] = float(value) * scales_by_key[key][sample_idx]

        ranked_keys = sorted(
            keys,
            key=lambda key: float(np.std(np.log(np.asarray(scales_by_key[key], dtype=float) + 1e-12))),
            reverse=True,
        )

    outputs = []
    for rec in records:
        in_path = Path(rec["path"])
        out_path = in_path.with_name(f"{in_path.stem}{output_suffix}.geojson")
        with out_path.open("w", encoding="utf-8") as handle:
            json.dump(rec["data"], handle)
        outputs.append(str(out_path))

    if ranked_keys:
        print(f"[geojson] Top changed keys: {', '.join(ranked_keys[:10])}")

    return records_before, records, sample_ids, key_diagnostics, ranked_keys, outputs


# ------------------------------ PIXEL MODE ------------------------------

def load_image_as_cyx(path):
    image = tifffile.imread(path)
    image = np.asarray(image)

    if image.ndim == 2:
        return image[np.newaxis, ...]

    if image.ndim != 3:
        raise ValueError(f"Unsupported image dimensions for {path}: {image.shape}")

    # Heuristic for (C, Y, X) vs (Y, X, C)
    if image.shape[0] <= 64 and image.shape[1] > 64 and image.shape[2] > 64:
        return image
    if image.shape[2] <= 64 and image.shape[0] > 64 and image.shape[1] > 64:
        return np.transpose(image, (2, 0, 1))

    # fallback assume channels first
    return image


def sample_pixels(channel_arr, sample_size, min_value):
    flat = np.ravel(channel_arr)
    valid = flat[flat >= min_value]
    if valid.size == 0:
        return np.array([0.0], dtype=float)
    if valid.size <= sample_size:
        return np.asarray(valid, dtype=float)
    rng = np.random.default_rng(0)
    idx = rng.choice(valid.size, size=sample_size, replace=False)
    return np.asarray(valid[idx], dtype=float)


def extract_ome_channel_names(path, n_channels):
    default_names = [f"channel_{idx}" for idx in range(n_channels)]

    try:
        with tifffile.TiffFile(path) as handle:
            ome_xml = handle.ome_metadata
    except Exception:
        return default_names

    if not ome_xml:
        return default_names

    try:
        root = ET.fromstring(ome_xml)
        channel_nodes = [node for node in root.iter() if node.tag.endswith("Channel")]
    except Exception:
        return default_names

    names = []
    for idx, node in enumerate(channel_nodes[:n_channels]):
        name = (node.get("Name") or "").strip()
        names.append(name if name else f"channel_{idx}")

    if len(names) < n_channels:
        names.extend([f"channel_{idx}" for idx in range(len(names), n_channels)])

    return names


def normalize_pixel(paths, num_bins, min_value, pixel_output_suffix, pixel_sample_size):
    sample_ids = [Path(path).stem.replace('.ome', '') for path in paths]
    n_samples = len(paths)

    print(f"[pixel] Loaded {n_samples} image sample(s)")

    n_channels = None
    for path in paths:
        image = load_image_as_cyx(path)
        channels_here = int(image.shape[0])
        n_channels = channels_here if n_channels is None else min(n_channels, channels_here)

    if n_channels is None or n_channels <= 0:
        raise ValueError("No valid image channels found for pixel normalization")

    channel_names = extract_ome_channel_names(paths[0], n_channels)
    print(f"[pixel] Detected {n_channels} channel(s)")

    scales_by_channel = {}
    ch_diagnostics = {}

    if n_samples >= 2:
        for channel in range(n_channels):
            if channel == 0 or (channel + 1) % 5 == 0 or channel == n_channels - 1:
                print(f"[pixel] Processing channel {channel + 1}/{n_channels}: {channel_names[channel]}")
            sample_values = []
            global_min = float("inf")
            global_max = float("-inf")

            for path in paths:
                image = load_image_as_cyx(path)
                sampled = sample_pixels(image[channel], pixel_sample_size, min_value=min_value)
                log_vals = log_transform(sampled, min_value=min_value)
                log_vals = log_vals[np.isfinite(log_vals)]
                if log_vals.size == 0:
                    log_vals = np.array([0.0], dtype=float)
                sample_values.append(log_vals)
                global_min = min(global_min, float(np.min(log_vals)))
                global_max = max(global_max, float(np.max(log_vals)))

            if not np.isfinite(global_min) or not np.isfinite(global_max) or global_max <= global_min:
                scales_by_channel[channel] = [1.0] * n_samples
                continue

            hist_list = []
            for log_vals in sample_values:
                counts, _ = np.histogram(log_vals, bins=num_bins, range=(global_min, global_max))
                hist_list.append(counts.astype(float))

            hist_matrix = np.stack(hist_list, axis=0)
            ref_idx = choose_reference(hist_matrix)
            shifts = compute_fft_shifts(hist_matrix[ref_idx], hist_list)
            increment = (global_max - global_min) / max(1, num_bins - 1)
            scales = [math.exp(-shift * increment) for shift in shifts]

            scales_by_channel[channel] = scales
            ch_diagnostics[channel] = {
                "reference_sample": sample_ids[ref_idx],
                "scales": scales,
                "shifts": shifts,
                "global_min": global_min,
                "global_max": global_max,
            }
    else:
        for channel in range(n_channels):
            scales_by_channel[channel] = [1.0]
            ch_diagnostics[channel] = {
                "reference_sample": sample_ids[0],
                "scales": [1.0],
                "shifts": [0],
                "global_min": 0.0,
                "global_max": 0.0,
            }

    outputs = []
    for sample_idx, path in enumerate(paths):
        image = load_image_as_cyx(path)
        corrected = np.empty_like(image)
        orig_dtype = image.dtype

        for channel in range(n_channels):
            scale = float(scales_by_channel[channel][sample_idx])
            channel_scaled = image[channel].astype(np.float32, copy=False) * scale
            if np.issubdtype(orig_dtype, np.integer):
                info = np.iinfo(orig_dtype)
                channel_scaled = np.clip(channel_scaled, info.min, info.max)
            corrected[channel] = channel_scaled.astype(orig_dtype, copy=False)

        if image.shape[0] > n_channels:
            corrected[n_channels:] = image[n_channels:]

        out_path = output_path_for_image(path, pixel_output_suffix)
        tifffile.imwrite(out_path, corrected)
        outputs.append(str(out_path))
        print(f"[pixel] Wrote normalized image {sample_idx + 1}/{n_samples}: {out_path}")

    ranked_channels = sorted(
        range(n_channels),
        key=lambda ch: float(np.std(np.log(np.asarray(scales_by_channel[ch], dtype=float) + 1e-12))),
        reverse=True,
    )

    if ranked_channels:
        top_channels = ranked_channels[:10]
        top_desc = [f"{ch}:{channel_names[ch]}" for ch in top_channels]
        print(f"[pixel] Top changed channels: {', '.join(top_desc)}")

    return sample_ids, channel_names, ch_diagnostics, ranked_channels, outputs


# ------------------------------ ADATA MODE ------------------------------

def sample_adata_group_feature_values(adata_obj, obs_idx, feature_idx, sample_size, min_value, rng):
    if obs_idx.size == 0:
        return np.array([0.0], dtype=float)

    if obs_idx.size > sample_size:
        chosen_obs = np.asarray(rng.choice(obs_idx, size=sample_size, replace=False), dtype=np.int64)
    else:
        chosen_obs = np.asarray(obs_idx, dtype=np.int64)

    matrix = adata_obj.X
    if sparse.issparse(matrix):
        values = np.asarray(matrix[chosen_obs, feature_idx].toarray(), dtype=float).ravel()
    else:
        values = np.asarray(matrix[chosen_obs, feature_idx], dtype=float).ravel()

    valid = values[values >= min_value]
    if valid.size == 0:
        return np.array([0.0], dtype=float)
    return valid


def infer_adata_feature_names(adata_obj):
    preferred_cols = [col for col in ["feature_type", "marker", "compartment", "statistic"] if col in adata_obj.var.columns]
    if preferred_cols:
        names = []
        for row_idx in range(int(adata_obj.n_vars)):
            parts = []
            for col in preferred_cols:
                value = adata_obj.var.iloc[row_idx][col]
                if value is None:
                    continue
                text = str(value).strip()
                if text and text.lower() != "nan":
                    parts.append(text)
            names.append("|".join(parts) if parts else str(adata_obj.var_names[row_idx]))
        return names
    if "marker_name" in adata_obj.var.columns:
        return [str(v) for v in adata_obj.var["marker_name"].tolist()]
    if "feature_name" in adata_obj.var.columns:
        return [str(v) for v in adata_obj.var["feature_name"].tolist()]
    return [str(v) for v in adata_obj.var_names.tolist()]


def select_adata_feature_indices(adata_obj, feature_names, adata_target, filter_column, filter_regex):
    if adata_target == "cell_mean":
        base_indices = []

        if "statistic" in adata_obj.var.columns:
            stat_values = [str(v).strip().lower() for v in adata_obj.var["statistic"].tolist()]

            if "compartment" in adata_obj.var.columns:
                compartment_values = [str(v).strip().lower() for v in adata_obj.var["compartment"].tolist()]
                base_indices = [
                    idx
                    for idx, (stat, compartment) in enumerate(zip(stat_values, compartment_values))
                    if stat == "mean" and compartment == "cell"
                ]
            else:
                base_indices = [idx for idx, stat in enumerate(stat_values) if stat == "mean"]

        if not base_indices:
            var_names = [str(v) for v in adata_obj.var_names.tolist()]
            base_indices = [idx for idx, name in enumerate(var_names) if re.search(r"(?i)_cell_mean$", name)]

        if not base_indices:
            base_indices = [idx for idx, name in enumerate(feature_names) if re.search(r"(?i)_cell_mean$", str(name))]
    else:
        base_indices = list(range(len(feature_names)))

    if not filter_regex:
        return base_indices

    regex = re.compile(filter_regex)

    if filter_column:
        if filter_column not in adata_obj.var.columns:
            available = ", ".join([str(col) for col in adata_obj.var.columns[:30]])
            raise ValueError(
                f"AnnData var column '{filter_column}' was not found for adata filtering. Available columns (first 30): {available}"
            )
        values = [str(v) for v in adata_obj.var[filter_column].tolist()]
    else:
        values = feature_names

    selected = [idx for idx in base_indices if regex.search(values[idx])]
    return selected


def build_adata_group_records(adatas, paths, group_by):
    records = []
    for adata_idx, adata_obj in enumerate(adatas):
        if group_by not in adata_obj.obs.columns:
            available = ", ".join([str(col) for col in adata_obj.obs.columns[:20]])
            raise ValueError(
                f"AnnData file '{paths[adata_idx]}' is missing obs column '{group_by}'. Available columns (first 20): {available}"
            )

        grouping = adata_obj.obs[group_by].astype(str)
        grouped_indices = grouping.groupby(grouping).indices

        for group_value, indices in grouped_indices.items():
            obs_idx = np.asarray(indices, dtype=np.int64)
            if obs_idx.size == 0:
                continue
            label = group_value if len(adatas) == 1 else f"{Path(paths[adata_idx]).stem}:{group_value}"
            records.append(
                {
                    "adata_idx": adata_idx,
                    "group": group_value,
                    "label": label,
                    "obs_idx": obs_idx,
                }
            )

    return records


def normalize_adata(
    paths,
    num_bins,
    min_value,
    output_suffix,
    group_by,
    adata_sample_size,
    adata_target,
    adata_filter_column,
    adata_filter_regex,
):
    ad = get_anndata_module()
    n_files = len(paths)
    print(f"[adata] Loaded {n_files} AnnData file(s)")

    adatas = [ad.read_h5ad(path) for path in paths]
    if not adatas:
        raise ValueError("No AnnData inputs found")

    group_records = build_adata_group_records(adatas, paths, group_by)
    sample_ids = [record["label"] for record in group_records]
    n_samples = len(sample_ids)
    print(f"[adata] Grouping by obs['{group_by}'] yielded {n_samples} sample group(s)")

    if n_samples <= 0:
        raise ValueError(f"No non-empty groups were found using obs column '{group_by}'")

    n_features = min(int(adata_obj.n_vars) for adata_obj in adatas)
    if n_features <= 0:
        raise ValueError("No features found in AnnData inputs")

    feature_names = infer_adata_feature_names(adatas[0])[:n_features]
    print(f"[adata] Detected {n_features} feature(s)")

    selected_features = select_adata_feature_indices(
        adatas[0],
        feature_names,
        adata_target=adata_target,
        filter_column=adata_filter_column,
        filter_regex=adata_filter_regex,
    )

    if adata_target == "cell_mean":
        print(f"[adata] Target preset 'cell_mean' selected: {len(selected_features)}/{n_features} features")

    if adata_filter_regex:
        source_desc = adata_filter_column if adata_filter_column else "feature names"
        print(
            f"[adata] Filtered features using regex '{adata_filter_regex}' on {source_desc}: "
            f"{len(selected_features)}/{n_features} selected"
        )
    elif adata_target == "all":
        print(f"[adata] No adata feature filter set; normalizing all {n_features} features")

    if len(selected_features) == 0:
        print("[adata] Warning: zero features matched filter. AnnData outputs will be written unchanged.")

    scales_by_feature = {}
    diagnostics = {}
    rng = np.random.default_rng(0)

    if n_samples >= 2 and selected_features:
        for loop_idx, feature_idx in enumerate(selected_features, start=1):
            if feature_idx == 0 or (feature_idx + 1) % 50 == 0 or feature_idx == n_features - 1:
                print(f"[adata] Processing feature {feature_idx + 1}/{n_features}: {feature_names[feature_idx]}")
            elif loop_idx == 1 or loop_idx % 50 == 0 or loop_idx == len(selected_features):
                print(f"[adata] Processing selected feature {loop_idx}/{len(selected_features)}: {feature_names[feature_idx]}")

            sample_values = []
            global_min = float("inf")
            global_max = float("-inf")

            for record in group_records:
                adata_obj = adatas[record["adata_idx"]]
                values = sample_adata_group_feature_values(
                    adata_obj,
                    record["obs_idx"],
                    feature_idx,
                    sample_size=adata_sample_size,
                    min_value=min_value,
                    rng=rng,
                )
                log_vals = log_transform(values, min_value=min_value)
                log_vals = log_vals[np.isfinite(log_vals)]
                if log_vals.size == 0:
                    log_vals = np.array([0.0], dtype=float)

                sample_values.append(log_vals)
                global_min = min(global_min, float(np.min(log_vals)))
                global_max = max(global_max, float(np.max(log_vals)))

            if not np.isfinite(global_min) or not np.isfinite(global_max) or global_max <= global_min:
                scales_by_feature[feature_idx] = [1.0] * n_samples
                diagnostics[feature_idx] = {
                    "reference_sample": sample_ids[0],
                    "scales": [1.0] * n_samples,
                    "shifts": [0] * n_samples,
                    "global_min": 0.0,
                    "global_max": 0.0,
                }
                continue

            hist_list = []
            for log_vals in sample_values:
                counts, _ = np.histogram(log_vals, bins=num_bins, range=(global_min, global_max))
                hist_list.append(counts.astype(float))

            hist_matrix = np.stack(hist_list, axis=0)
            ref_idx = choose_reference(hist_matrix)
            shifts = compute_fft_shifts(hist_matrix[ref_idx], hist_list)

            increment = (global_max - global_min) / max(1, num_bins - 1)
            scales = [math.exp(-shift * increment) for shift in shifts]

            scales_by_feature[feature_idx] = scales
            diagnostics[feature_idx] = {
                "reference_sample": sample_ids[ref_idx],
                "scales": scales,
                "shifts": shifts,
                "global_min": global_min,
                "global_max": global_max,
            }
    else:
        for feature_idx in selected_features:
            scales_by_feature[feature_idx] = [1.0]
            diagnostics[feature_idx] = {
                "reference_sample": sample_ids[0],
                "scales": [1.0],
                "shifts": [0],
                "global_min": 0.0,
                "global_max": 0.0,
            }

    per_sample_scale_vectors = []
    for sample_idx in range(n_samples):
        vec = np.ones(n_features, dtype=float)
        for feature_idx in selected_features:
            vec[feature_idx] = float(scales_by_feature[feature_idx][sample_idx])
        per_sample_scale_vectors.append(vec)

    for adata_idx, adata_obj in enumerate(adatas):
        matrix = adata_obj.X
        if sparse.issparse(matrix):
            matrix = matrix.tocsr(copy=True)
        else:
            matrix = np.asarray(matrix, dtype=float)

        for sample_idx, record in enumerate(group_records):
            if record["adata_idx"] != adata_idx:
                continue
            obs_idx = record["obs_idx"]
            scale_vector = per_sample_scale_vectors[sample_idx]
            if sparse.issparse(matrix):
                matrix[obs_idx, :] = matrix[obs_idx, :].multiply(scale_vector)
            else:
                matrix[obs_idx, :] = matrix[obs_idx, :] * scale_vector

        adata_obj.X = matrix

    outputs = []
    for path, adata_obj in zip(paths, adatas):
        out_path = output_path_for_adata(path, output_suffix)
        adata_obj.write_h5ad(out_path)
        outputs.append(str(out_path))
        print(f"[adata] Wrote normalized AnnData: {out_path}")

    ranked_features = sorted(
        selected_features,
        key=lambda idx: float(np.std(np.log(np.asarray(scales_by_feature[idx], dtype=float) + 1e-12))),
        reverse=True,
    )
    feature_name_map = {idx: feature_names[idx] for idx in range(n_features)}

    return sample_ids, group_records, feature_name_map, diagnostics, ranked_features, outputs


# ------------------------------ QC ------------------------------

def collect_per_sample_values(records, key):
    per_sample = []
    for rec in records:
        vals = []
        for _, feat in rec["cell_features"]:
            measurements = feat.get("properties", {}).get("measurements", {})
            value = measurements.get(key)
            if is_number(value):
                vals.append(float(value))
        per_sample.append(np.asarray(vals, dtype=float))
    return per_sample


def write_qc_summary(qc_dir, diagnostics, ranked_items, item_prefix="key", item_name_map=None):
    qc_dir.mkdir(parents=True, exist_ok=True)

    summary_path = qc_dir / f"uniform_{item_prefix}_summary.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        header = [item_prefix]
        if item_name_map is not None:
            header.append(f"{item_prefix}_name")
        header.extend(
            [
                "reference_sample",
                "scale_min",
                "scale_max",
                "scale_median",
                "scale_std",
                "shift_min",
                "shift_max",
                "global_min_log",
                "global_max_log",
            ]
        )
        writer.writerow(header)

        for item in ranked_items:
            diag = diagnostics[item]
            scales = np.asarray(diag["scales"], dtype=float)
            shifts = np.asarray(diag["shifts"], dtype=float)
            row = [item]
            if item_name_map is not None:
                row.append(item_name_map.get(item, ""))
            row.extend(
                [
                    diag["reference_sample"],
                    float(np.min(scales)),
                    float(np.max(scales)),
                    float(np.median(scales)),
                    float(np.std(scales)),
                    float(np.min(shifts)),
                    float(np.max(shifts)),
                    float(diag["global_min"]),
                    float(diag["global_max"]),
                ]
            )
            writer.writerow(row)

    run_summary_path = qc_dir / f"uniform_{item_prefix}_run_summary.json"
    with run_summary_path.open("w", encoding="utf-8") as handle:
        json.dump(
            {
                f"n_{item_prefix}s": len(ranked_items),
                f"top_changed_{item_prefix}s": [str(i) for i in ranked_items[:20]],
            },
            handle,
            indent=2,
        )

    return summary_path, run_summary_path


def plot_geojson_qc(qc_dir, records_before, records_after, sample_ids, diagnostics, ranked_keys, min_value, top_n_keys, max_heatmap_keys):
    try:
        import matplotlib.pyplot as plt
    except Exception as error:
        print(f"Plotting skipped: matplotlib unavailable ({error})")
        return []

    saved = []
    selected_keys = ranked_keys[: max(0, top_n_keys)]
    for key in selected_keys:
        before_vals = collect_per_sample_values(records_before, key)
        after_vals = collect_per_sample_values(records_after, key)
        diag = diagnostics[key]
        global_min = diag["global_min"]
        global_max = diag["global_max"]

        fig, axes = plt.subplots(1, 2, figsize=(14, 5), dpi=140)
        for sid, arr in zip(sample_ids, before_vals):
            vals = log_transform(arr, min_value=min_value)
            vals = vals[np.isfinite(vals)]
            if vals.size:
                axes[0].hist(vals, bins=80, range=(global_min, global_max), histtype="step", linewidth=1.2, label=sid)
        for sid, arr in zip(sample_ids, after_vals):
            vals = log_transform(arr, min_value=min_value)
            vals = vals[np.isfinite(vals)]
            if vals.size:
                axes[1].hist(vals, bins=80, range=(global_min, global_max), histtype="step", linewidth=1.2, label=sid)

        axes[0].set_title(f"Before normalization\n{key}")
        axes[1].set_title(f"After normalization\n{key}")
        axes[0].set_xlabel("log(value)")
        axes[1].set_xlabel("log(value)")
        axes[0].set_ylabel("Count")
        axes[1].set_ylabel("Count")

        handles, labels = axes[1].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc="center right", frameon=False)
        fig.tight_layout(rect=[0, 0, 0.88, 1])

        out_path = qc_dir / f"hist_before_after_{re.sub(r'[^A-Za-z0-9_.-]+', '_', key)}.png"
        fig.savefig(out_path)
        plt.close(fig)
        saved.append(str(out_path))

    if ranked_keys:
        heatmap_keys = ranked_keys[: max(1, max_heatmap_keys)]
        scale_matrix = np.array([diagnostics[key]["scales"] for key in heatmap_keys], dtype=float)
        log_scale = np.log(scale_matrix + 1e-12)

        fig, ax = plt.subplots(figsize=(max(7, len(sample_ids) * 1.2), max(6, len(heatmap_keys) * 0.25)), dpi=140)
        image = ax.imshow(log_scale, aspect="auto", cmap="coolwarm")
        ax.set_title("Log scale factors by key and sample")
        ax.set_xlabel("Sample")
        ax.set_ylabel("Measurement key")
        ax.set_xticks(range(len(sample_ids)))
        ax.set_xticklabels(sample_ids, rotation=45, ha="right")
        ax.set_yticks(range(len(heatmap_keys)))
        ax.set_yticklabels(heatmap_keys)
        fig.colorbar(image, ax=ax, label="log(scale)")
        fig.tight_layout()

        out_path = qc_dir / "scale_factor_heatmap.png"
        fig.savefig(out_path)
        plt.close(fig)
        saved.append(str(out_path))

    return saved


def plot_pixel_qc(qc_dir, sample_ids, channel_names, diagnostics, ranked_channels, max_heatmap_keys):
    try:
        import matplotlib.pyplot as plt
    except Exception as error:
        print(f"Plotting skipped: matplotlib unavailable ({error})")
        return []

    saved = []
    if ranked_channels:
        channels = ranked_channels[: max(1, max_heatmap_keys)]
        scale_matrix = np.array([diagnostics[ch]["scales"] for ch in channels], dtype=float)
        log_scale = np.log(scale_matrix + 1e-12)

        fig, ax = plt.subplots(figsize=(max(7, len(sample_ids) * 1.2), max(6, len(channels) * 0.25)), dpi=140)
        image = ax.imshow(log_scale, aspect="auto", cmap="coolwarm")
        ax.set_title("Log scale factors by channel and sample (pixel mode)")
        ax.set_xlabel("Sample")
        ax.set_ylabel("Channel index")
        ax.set_xticks(range(len(sample_ids)))
        ax.set_xticklabels(sample_ids, rotation=45, ha="right")
        ax.set_yticks(range(len(channels)))
        ax.set_yticklabels([f"{ch} ({channel_names[ch]})" for ch in channels])
        fig.colorbar(image, ax=ax, label="log(scale)")
        fig.tight_layout()

        out_path = qc_dir / "pixel_scale_factor_heatmap.png"
        fig.savefig(out_path)
        plt.close(fig)
        saved.append(str(out_path))

    return saved


def plot_pixel_histograms_qc(
    qc_dir,
    paths,
    sample_ids,
    channel_names,
    diagnostics,
    ranked_channels,
    min_value,
    sample_size,
    top_n_channels,
):
    try:
        import matplotlib.pyplot as plt
    except Exception as error:
        print(f"Plotting skipped: matplotlib unavailable ({error})")
        return []

    saved = []
    selected_channels = ranked_channels[: max(0, top_n_channels)]
    for channel in selected_channels:
        diag = diagnostics[channel]
        global_min = diag["global_min"]
        global_max = diag["global_max"]

        fig, axes = plt.subplots(1, 2, figsize=(14, 5), dpi=140)
        for sample_idx, (sid, path) in enumerate(zip(sample_ids, paths)):
            image = load_image_as_cyx(path)
            sampled = sample_pixels(image[channel], sample_size, min_value=min_value)

            before_vals = log_transform(sampled, min_value=min_value)
            before_vals = before_vals[np.isfinite(before_vals)]

            scale = float(diag["scales"][sample_idx])
            after_vals = log_transform(sampled.astype(np.float64) * scale, min_value=min_value)
            after_vals = after_vals[np.isfinite(after_vals)]

            if before_vals.size:
                axes[0].hist(before_vals, bins=80, range=(global_min, global_max), histtype="step", linewidth=1.2, label=sid)
            if after_vals.size:
                axes[1].hist(after_vals, bins=80, range=(global_min, global_max), histtype="step", linewidth=1.2, label=sid)

        channel_label = f"{channel} ({channel_names[channel]})"
        axes[0].set_title(f"Before normalization\nchannel {channel_label}")
        axes[1].set_title(f"After normalization\nchannel {channel_label}")
        axes[0].set_xlabel("log(value)")
        axes[1].set_xlabel("log(value)")
        axes[0].set_ylabel("Count")
        axes[1].set_ylabel("Count")

        handles, labels = axes[1].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc="center right", frameon=False)
        fig.tight_layout(rect=[0, 0, 0.88, 1])

        safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", channel_names[channel])
        out_path = qc_dir / f"pixel_hist_before_after_channel_{channel}_{safe_name}.png"
        fig.savefig(out_path)
        plt.close(fig)
        saved.append(str(out_path))

    return saved


def plot_adata_qc(
    qc_dir,
    paths,
    group_records,
    feature_name_map,
    diagnostics,
    ranked_features,
    min_value,
    adata_sample_size,
    top_n_features,
    max_heatmap_features,
):
    ad = get_anndata_module()
    try:
        import matplotlib.pyplot as plt
    except Exception as error:
        print(f"Plotting skipped: matplotlib unavailable ({error})")
        return []

    adatas = [ad.read_h5ad(path) for path in paths]
    saved = []
    sample_ids = [record["label"] for record in group_records]
    rng = np.random.default_rng(0)

    selected_features = ranked_features[: max(0, top_n_features)]
    for feature_idx in selected_features:
        diag = diagnostics[feature_idx]
        global_min = diag["global_min"]
        global_max = diag["global_max"]

        fig, axes = plt.subplots(1, 2, figsize=(14, 5), dpi=140)
        for sample_idx, (sid, record) in enumerate(zip(sample_ids, group_records)):
            adata_obj = adatas[record["adata_idx"]]
            values = sample_adata_group_feature_values(
                adata_obj,
                record["obs_idx"],
                feature_idx,
                sample_size=adata_sample_size,
                min_value=min_value,
                rng=rng,
            )
            scale = float(diag["scales"][sample_idx])

            before_vals = log_transform(values, min_value=min_value)
            before_vals = before_vals[np.isfinite(before_vals)]

            after_vals = log_transform(values.astype(np.float64) * scale, min_value=min_value)
            after_vals = after_vals[np.isfinite(after_vals)]

            if before_vals.size:
                axes[0].hist(before_vals, bins=80, range=(global_min, global_max), histtype="step", linewidth=1.2, label=sid)
            if after_vals.size:
                axes[1].hist(after_vals, bins=80, range=(global_min, global_max), histtype="step", linewidth=1.2, label=sid)

        feature_name = feature_name_map.get(feature_idx, str(feature_idx))
        axes[0].set_title(f"Before normalization\nfeature {feature_name}")
        axes[1].set_title(f"After normalization\nfeature {feature_name}")
        axes[0].set_xlabel("log(value)")
        axes[1].set_xlabel("log(value)")
        axes[0].set_ylabel("Count")
        axes[1].set_ylabel("Count")

        handles, labels = axes[1].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc="center right", frameon=False)
        fig.tight_layout(rect=[0, 0, 0.88, 1])

        safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", feature_name)
        out_path = qc_dir / f"adata_hist_before_after_feature_{feature_idx}_{safe_name}.png"
        fig.savefig(out_path)
        plt.close(fig)
        saved.append(str(out_path))

    if ranked_features:
        selected = ranked_features[: max(1, max_heatmap_features)]
        scale_matrix = np.array([diagnostics[idx]["scales"] for idx in selected], dtype=float)
        log_scale = np.log(scale_matrix + 1e-12)

        ylabels = [feature_name_map.get(idx, str(idx)) for idx in selected]

        fig, ax = plt.subplots(figsize=(max(7, len(sample_ids) * 1.2), max(6, len(selected) * 0.25)), dpi=140)
        image = ax.imshow(log_scale, aspect="auto", cmap="coolwarm")
        ax.set_title("Log scale factors by AnnData feature and sample")
        ax.set_xlabel("Sample")
        ax.set_ylabel("Feature")
        ax.set_xticks(range(len(sample_ids)))
        ax.set_xticklabels(sample_ids, rotation=45, ha="right")
        ax.set_yticks(range(len(selected)))
        ax.set_yticklabels(ylabels)
        fig.colorbar(image, ax=ax, label="log(scale)")
        fig.tight_layout()

        out_path = qc_dir / "adata_scale_factor_heatmap.png"
        fig.savefig(out_path)
        plt.close(fig)
        saved.append(str(out_path))

    return saved


def main():
    args = parse_args()
    mode = args.mode

    generate_plots = parse_bool(args.generate_plots)
    qc_dir = Path(args.qc_dir)
    qc_dir.mkdir(parents=True, exist_ok=True)

    outputs = []

    if mode == "geojson":
        geojson_inputs = [path for path in args.inputs if path.lower().endswith(".geojson")]
        if geojson_inputs:
            print(f"[geojson] Starting normalization for {len(geojson_inputs)} input file(s)")
            records = load_geojsons(geojson_inputs)
            records_before, records_after, sample_ids, diagnostics, ranked_keys, geo_outputs = normalize_geojson(
                records,
                num_bins=args.num_bins,
                min_value=args.min_value,
                exclude_pattern=args.exclude_pattern,
                output_suffix=args.output_suffix,
            )
            outputs.extend(geo_outputs)

            summary_path, run_summary_path = write_qc_summary(qc_dir, diagnostics, ranked_keys, item_prefix="key")
            print(f"Wrote {summary_path}")
            print(f"Wrote {run_summary_path}")

            if generate_plots:
                for plot_path in plot_geojson_qc(
                    qc_dir,
                    records_before,
                    records_after,
                    sample_ids,
                    diagnostics,
                    ranked_keys,
                    args.min_value,
                    args.qc_top_n_keys,
                    args.qc_max_heatmap_keys,
                ):
                    print(f"Wrote {plot_path}")
        else:
            print("[geojson] No geojson inputs found for selected mode")

    if mode == "ome_tiff":
        tiff_inputs = [
            path
            for path in args.inputs
            if path.lower().endswith(".tiff") or path.lower().endswith(".tif")
        ]
        if tiff_inputs:
            print(f"[pixel] Starting normalization for {len(tiff_inputs)} input image(s)")
            sample_ids, channel_names, diagnostics, ranked_channels, tiff_outputs = normalize_pixel(
                tiff_inputs,
                num_bins=args.num_bins,
                min_value=args.min_value,
                pixel_output_suffix=args.pixel_output_suffix,
                pixel_sample_size=args.pixel_sample_size,
            )
            outputs.extend(tiff_outputs)

            channel_name_map = {idx: name for idx, name in enumerate(channel_names)}
            summary_path, run_summary_path = write_qc_summary(
                qc_dir,
                diagnostics,
                ranked_channels,
                item_prefix="channel",
                item_name_map=channel_name_map,
            )
            print(f"Wrote {summary_path}")
            print(f"Wrote {run_summary_path}")

            if generate_plots:
                for plot_path in plot_pixel_qc(
                    qc_dir,
                    sample_ids,
                    channel_names,
                    diagnostics,
                    ranked_channels,
                    args.qc_max_heatmap_keys,
                ):
                    print(f"Wrote {plot_path}")

                for plot_path in plot_pixel_histograms_qc(
                    qc_dir,
                    tiff_inputs,
                    sample_ids,
                    channel_names,
                    diagnostics,
                    ranked_channels,
                    args.min_value,
                    args.pixel_sample_size,
                    args.qc_top_n_keys,
                ):
                    print(f"Wrote {plot_path}")
        else:
            print("[pixel] No TIFF inputs found for selected mode")

    if mode == "adata":
        adata_inputs = [path for path in args.inputs if path.lower().endswith(".h5ad")]
        if adata_inputs:
            print(f"[adata] Starting normalization for {len(adata_inputs)} input file(s)")
            sample_ids, group_records, feature_name_map, diagnostics, ranked_features, adata_outputs = normalize_adata(
                adata_inputs,
                num_bins=args.num_bins,
                min_value=args.min_value,
                output_suffix=args.output_suffix,
                group_by=args.adata_group_by,
                adata_sample_size=args.adata_sample_size,
                adata_target=args.adata_target,
                adata_filter_column=args.adata_filter_column,
                adata_filter_regex=args.adata_filter_regex,
            )
            outputs.extend(adata_outputs)

            summary_path, run_summary_path = write_qc_summary(
                qc_dir,
                diagnostics,
                ranked_features,
                item_prefix="feature",
                item_name_map=feature_name_map,
            )
            print(f"Wrote {summary_path}")
            print(f"Wrote {run_summary_path}")

            if generate_plots:
                for plot_path in plot_adata_qc(
                    qc_dir,
                    adata_inputs,
                    group_records,
                    feature_name_map,
                    diagnostics,
                    ranked_features,
                    args.min_value,
                    args.adata_sample_size,
                    args.qc_top_n_keys,
                    args.qc_max_heatmap_keys,
                ):
                    print(f"Wrote {plot_path}")
        else:
            print("[adata] No h5ad inputs found for selected mode")

    if not outputs:
        raise SystemExit("No compatible inputs found for selected mode.")

    for out in outputs:
        print(f"Wrote {out}")


if __name__ == "__main__":
    main()
