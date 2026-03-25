#!/usr/bin/env python3
"""
UniFORM-style normalization for GeoJSON measurements and OME-TIFF/TIFF images.
"""

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import tifffile
from scipy.signal import correlate


def parse_args():
    parser = argparse.ArgumentParser(description="Normalize GeoJSON measurements and/or OME-TIFF intensities")
    parser.add_argument("--inputs", nargs="+", required=True, help="Input files")
    parser.add_argument("--mode", choices=["geojson", "ome_tiff", "both"], default="geojson")
    parser.add_argument("--num-bins", type=int, default=1024)
    parser.add_argument("--min-value", type=float, default=1.0)
    parser.add_argument("--exclude-pattern", default=r"^(kronos_|emb_)")
    parser.add_argument("--output-suffix", default="_uniform")
    parser.add_argument("--pixel-output-suffix", default="_uniform")
    parser.add_argument("--pixel-sample-size", type=int, default=200000)
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

    key_diagnostics = {}
    ranked_keys = []

    n_samples = len(records)
    if keys and n_samples >= 2:
        scales_by_key = {}
        for key in keys:
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
    flat = channel_arr.reshape(-1).astype(float)
    valid = flat[flat >= min_value]
    if valid.size == 0:
        return np.array([0.0], dtype=float)
    if valid.size <= sample_size:
        return valid
    rng = np.random.default_rng(0)
    idx = rng.choice(valid.size, size=sample_size, replace=False)
    return valid[idx]


def normalize_pixel(paths, num_bins, min_value, pixel_output_suffix, pixel_sample_size):
    sample_ids = [Path(path).stem.replace('.ome', '') for path in paths]
    images = [load_image_as_cyx(path) for path in paths]

    n_samples = len(images)
    n_channels = min(img.shape[0] for img in images)

    scales_by_channel = {}
    ch_diagnostics = {}

    if n_samples >= 2:
        for channel in range(n_channels):
            sample_values = []
            global_min = float("inf")
            global_max = float("-inf")

            for image in images:
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
    for sample_idx, (path, image) in enumerate(zip(paths, images)):
        corrected = image.astype(np.float64, copy=True)
        for channel in range(n_channels):
            corrected[channel] *= float(scales_by_channel[channel][sample_idx])

        orig_dtype = image.dtype
        if np.issubdtype(orig_dtype, np.integer):
            info = np.iinfo(orig_dtype)
            corrected = np.clip(corrected, info.min, info.max).astype(orig_dtype)
        else:
            corrected = corrected.astype(orig_dtype)

        out_path = output_path_for_image(path, pixel_output_suffix)
        tifffile.imwrite(out_path, corrected)
        outputs.append(str(out_path))

    ranked_channels = sorted(
        range(n_channels),
        key=lambda ch: float(np.std(np.log(np.asarray(scales_by_channel[ch], dtype=float) + 1e-12))),
        reverse=True,
    )

    return sample_ids, ch_diagnostics, ranked_channels, outputs


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


def write_qc_summary(qc_dir, diagnostics, ranked_items, item_prefix="key"):
    qc_dir.mkdir(parents=True, exist_ok=True)

    summary_path = qc_dir / f"uniform_{item_prefix}_summary.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                item_prefix,
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

        for item in ranked_items:
            diag = diagnostics[item]
            scales = np.asarray(diag["scales"], dtype=float)
            shifts = np.asarray(diag["shifts"], dtype=float)
            writer.writerow(
                [
                    item,
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


def plot_pixel_qc(qc_dir, sample_ids, diagnostics, ranked_channels, max_heatmap_keys):
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
        ax.set_yticklabels(channels)
        fig.colorbar(image, ax=ax, label="log(scale)")
        fig.tight_layout()

        out_path = qc_dir / "pixel_scale_factor_heatmap.png"
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

    if mode in {"geojson", "both"}:
        geojson_inputs = [path for path in args.inputs if path.lower().endswith(".geojson")]
        if geojson_inputs:
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

    if mode in {"ome_tiff", "both"}:
        tiff_inputs = [
            path
            for path in args.inputs
            if path.lower().endswith(".tiff") or path.lower().endswith(".tif")
        ]
        if tiff_inputs:
            sample_ids, diagnostics, ranked_channels, tiff_outputs = normalize_pixel(
                tiff_inputs,
                num_bins=args.num_bins,
                min_value=args.min_value,
                pixel_output_suffix=args.pixel_output_suffix,
                pixel_sample_size=args.pixel_sample_size,
            )
            outputs.extend(tiff_outputs)

            summary_path, run_summary_path = write_qc_summary(qc_dir, diagnostics, ranked_channels, item_prefix="channel")
            print(f"Wrote {summary_path}")
            print(f"Wrote {run_summary_path}")

            if generate_plots:
                for plot_path in plot_pixel_qc(
                    qc_dir,
                    sample_ids,
                    diagnostics,
                    ranked_channels,
                    args.qc_max_heatmap_keys,
                ):
                    print(f"Wrote {plot_path}")

    if not outputs:
        raise SystemExit("No compatible inputs found for selected mode.")

    for out in outputs:
        print(f"Wrote {out}")


if __name__ == "__main__":
    main()
