#!/usr/bin/env python3
"""
UniFORM-style measurement normalization for cellmeasurement GeoJSON outputs.

Applies cohort-level histogram/correlation shift alignment to numeric cell
measurement features across samples and writes per-sample normalized GeoJSONs.
KRONOS/embedding keys are excluded by pattern.
"""

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.signal import correlate


def parse_args():
    parser = argparse.ArgumentParser(description="Normalize cellmeasurement GeoJSON measurements across samples")
    parser.add_argument("--inputs", nargs="+", required=True, help="Input GeoJSON files")
    parser.add_argument("--num-bins", type=int, default=1024, help="Histogram bin count")
    parser.add_argument("--min-value", type=float, default=1.0, help="Minimum value before log transform")
    parser.add_argument(
        "--exclude-pattern",
        default=r"^(kronos_|emb_)",
        help="Regex for measurement keys to exclude from normalization",
    )
    parser.add_argument(
        "--output-suffix",
        default="_uniform",
        help="Suffix to append before .geojson for normalized outputs",
    )
    parser.add_argument("--qc-dir", default="qc", help="Directory for QC artifacts")
    parser.add_argument("--generate-plots", default="true", help="Generate QC plots (true/false)")
    parser.add_argument("--qc-top-n-keys", type=int, default=12, help="Number of keys to plot")
    parser.add_argument("--qc-max-heatmap-keys", type=int, default=40, help="Max keys shown in scale heatmap")
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


def normalize_records(records, keys, num_bins, min_value, sample_ids):
    n_samples = len(records)
    if n_samples < 2:
        print("Only one sample found; writing pass-through normalized files.")
        return {}, []

    scales_by_key = {}
    key_diagnostics = {}

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

    return key_diagnostics, ranked_keys


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


def write_qc_summary(qc_dir, key_diagnostics, ranked_keys):
    qc_dir.mkdir(parents=True, exist_ok=True)

    summary_path = qc_dir / "uniform_key_summary.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "key",
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

        for key in ranked_keys:
            diag = key_diagnostics[key]
            scales = np.asarray(diag["scales"], dtype=float)
            shifts = np.asarray(diag["shifts"], dtype=float)
            writer.writerow(
                [
                    key,
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

    run_summary_path = qc_dir / "uniform_run_summary.json"
    with run_summary_path.open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "n_keys": len(ranked_keys),
                "top_changed_keys": ranked_keys[:20],
            },
            handle,
            indent=2,
        )

    return summary_path, run_summary_path


def plot_qc(
    qc_dir,
    records_before,
    records_after,
    sample_ids,
    key_diagnostics,
    ranked_keys,
    min_value,
    top_n_keys,
    max_heatmap_keys,
):
    try:
        import matplotlib.pyplot as plt
    except Exception as error:
        print(f"Plotting skipped: matplotlib unavailable ({error})")
        return []

    qc_dir.mkdir(parents=True, exist_ok=True)
    saved = []

    selected_keys = ranked_keys[: max(0, top_n_keys)]
    for key in selected_keys:
        before_vals = collect_per_sample_values(records_before, key)
        after_vals = collect_per_sample_values(records_after, key)

        diag = key_diagnostics[key]
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
        scale_matrix = np.array([key_diagnostics[key]["scales"] for key in heatmap_keys], dtype=float)
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


def write_outputs(records, output_suffix):
    output_paths = []
    for rec in records:
        in_path = Path(rec["path"])
        out_path = in_path.with_name(f"{in_path.stem}{output_suffix}.geojson")
        with open(out_path, "w", encoding="utf-8") as handle:
            json.dump(rec["data"], handle)
        output_paths.append(str(out_path))
    return output_paths


def main():
    args = parse_args()
    records = load_geojsons(args.inputs)
    sample_ids = [rec["sample_id"] for rec in records]
    records_before = load_geojsons(args.inputs)

    if not records:
        raise SystemExit("No input GeoJSON files provided.")

    keys = collect_measurement_keys(records, args.exclude_pattern)
    print(f"Loaded {len(records)} GeoJSON files")
    print(f"Normalizing {len(keys)} shared numeric measurement keys")

    key_diagnostics = {}
    ranked_keys = []
    if keys:
        key_diagnostics, ranked_keys = normalize_records(
            records,
            keys,
            num_bins=args.num_bins,
            min_value=args.min_value,
            sample_ids=sample_ids,
        )

    outputs = write_outputs(records, args.output_suffix)
    for out in outputs:
        print(f"Wrote {out}")

    qc_dir = Path(args.qc_dir)
    summary_path, run_summary_path = write_qc_summary(qc_dir, key_diagnostics, ranked_keys)
    print(f"Wrote {summary_path}")
    print(f"Wrote {run_summary_path}")

    if parse_bool(args.generate_plots):
        plot_paths = plot_qc(
            qc_dir=qc_dir,
            records_before=records_before,
            records_after=records,
            sample_ids=sample_ids,
            key_diagnostics=key_diagnostics,
            ranked_keys=ranked_keys,
            min_value=args.min_value,
            top_n_keys=args.qc_top_n_keys,
            max_heatmap_keys=args.qc_max_heatmap_keys,
        )
        for plot_path in plot_paths:
            print(f"Wrote {plot_path}")


if __name__ == "__main__":
    main()
