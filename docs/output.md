# Output

## Main output

- `uniformnormalize/*_uniform.geojson`
  - Input GeoJSON files with normalized numeric cellmeasurement keys.
  - Keys matching `uniform_exclude_pattern` are left unchanged.

- `uniformnormalize/*_uniform.ome.tiff` and/or `uniformnormalize/*_uniform.tiff`
  - Input OME-TIFF/TIFF files with per-channel cohort-level scaling applied.
  - Output format follows the input extension where possible.

## QC output

- `uniformnormalize/qc/uniform_key_summary.csv`
  - Per-key normalization diagnostics (scale ranges, shift ranges, chosen reference sample).
- `uniformnormalize/qc/uniform_run_summary.json`
  - Run-level summary with the most changed keys.
- `uniformnormalize/qc/hist_before_after_<key>.png`
  - Before/after log-distribution overlays per sample for the top affected keys.
- `uniformnormalize/qc/scale_factor_heatmap.png`
  - Heatmap of log scale factors (keys × samples) to quickly spot unstable settings.

## Versions

- `uniformnormalize/versions.yml`
