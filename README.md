# uniform_batch_correction

Standalone Nextflow pipeline for UniFORM-style batch normalization of cellmeasurement GeoJSON outputs.

## Input

Provide a CSV/YAML samplesheet with:

- `sample`: sample ID
- `geojson`: path to input cellmeasurement GeoJSON (for feature-level mode)
- `ome_tiff`: path to input OME-TIFF/TIFF image (for pixel-level mode)

Example CSV is provided in [assets/samplesheet.csv](assets/samplesheet.csv).

## Usage

```bash
nextflow run . \
  -profile conda,large \
  --input assets/samplesheet.csv \
  --outdir results
```

The `large` profile targets SLURM large-memory nodes (default queue `lrg`) and increases process resources for large images.

## Key parameters

- `--run_uniform` (default: `true`)
- `--uniform_apply_to` (default: `geojson`; options: `geojson`, `ome_tiff`, `both`)
- `--uniform_num_bins` (default: `1024`)
- `--uniform_min_value` (default: `1.0`)
- `--uniform_exclude_pattern` (default: `^(kronos_|emb_)`)
- `--uniform_output_suffix` (default: `_uniform`)
- `--uniform_pixel_output_suffix` (default: `_uniform`)
- `--uniform_pixel_sample_size` (default: `200000`)
- `--uniform_generate_plots` (default: `true`)
- `--uniform_qc_top_n_keys` (default: `12`)
- `--uniform_qc_max_heatmap_keys` (default: `40`)

Normalized GeoJSON files are published to `results/uniformnormalize/`.
QC files are published to `results/uniformnormalize/qc/`.
