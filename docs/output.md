# Output

## Main output

- `uniformnormalize/*_uniform.geojson`
  - Input GeoJSON files with normalized numeric cellmeasurement keys.
  - Keys matching `uniform_exclude_pattern` are left unchanged.

- `uniformnormalize/*_uniform.ome.tiff` and/or `uniformnormalize/*_uniform.tiff`
  - Input OME-TIFF/TIFF files with per-channel cohort-level scaling applied.
  - Output format follows the input extension where possible.
- `uniformnormalize/*_unifrom.ome.tiff` and/or `uniformnormalize/*_unifrom.tiff`
  - Same as above when using the current default `uniform_pixel_output_suffix`.

## QC output

- `uniformnormalize/qc/uniform_key_summary.csv`
  - Per-key normalization diagnostics (scale ranges, shift ranges, chosen reference sample).
- `uniformnormalize/qc/uniform_run_summary.json`
  - Run-level summary with the most changed keys.
- `uniformnormalize/qc/hist_before_after_<key>.png`
  - Before/after log-distribution overlays per sample for the top affected keys.
- `uniformnormalize/qc/scale_factor_heatmap.png`
  - Heatmap of log scale factors (keys × samples) to quickly spot unstable settings.
- `uniformnormalize/qc/uniform_channel_summary.csv`
  - Per-channel normalization diagnostics, including channel index and channel name.
- `uniformnormalize/qc/uniform_channel_run_summary.json`
  - Run-level summary with the most changed channels.
- `uniformnormalize/qc/pixel_scale_factor_heatmap.png`
  - Heatmap of log scale factors (channels × samples) for pixel mode.
- `uniformnormalize/qc/pixel_hist_before_after_channel_<idx>_<name>.png`
  - Before/after log-intensity overlays per sample for top changed channels in pixel mode.

## Versions

- `uniformnormalize/versions.yml`
