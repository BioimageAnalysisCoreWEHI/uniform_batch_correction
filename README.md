# uniform_batch_correction

Standalone Nextflow pipeline for UniFORM-style batch normalization of cellmeasurement GeoJSON outputs.

## Input

Provide a CSV/YAML samplesheet with:

- `sample`: sample ID
- `geojson`: path to input cellmeasurement GeoJSON (for feature-level mode)
- `ome_tiff`: path to input OME-TIFF/TIFF image (for pixel-level mode)
- `adata`: path to input AnnData `.h5ad` file (for feature-level matrix mode)

For `adata` mode, normalization is performed per group inside each `.h5ad` using an `obs` column (default: `image`).
Each unique group value is treated as one sample for cohort alignment.

Example CSV is provided in [assets/samplesheet.csv](assets/samplesheet.csv).

## Usage

```bash
nextflow run . \
  -profile conda,large \
  --input assets/samplesheet.csv \
  --outdir results
```

The `large` profile targets SLURM large-memory nodes (default queue `regular`) and increases process resources for large images.

## Key parameters

- `--run_uniform` (default: `true`)
- `--uniform_apply_to` (default: `geojson`; options: `geojson`, `ome_tiff`, `adata`)
- `--uniform_num_bins` (default: `1024`)
- `--uniform_min_value` (default: `1.0`)
- `--uniform_exclude_pattern` (default: `^(kronos_|emb_)`)
- `--uniform_output_suffix` (default: `_uniform`)
- `--uniform_pixel_output_suffix` (default: `_unifrom`)
- `--uniform_pixel_sample_size` (default: `200000`)
- `--uniform_adata_group_by` (default: `image`)
- `--uniform_adata_sample_size` (default: `200000`)
- `--uniform_adata_filter_column` (default: empty; e.g. `statistic`)
- `--uniform_adata_filter_regex` (default: empty; only matching features are normalized)
- `--uniform_generate_plots` (default: `true`)
- `--uniform_qc_top_n_keys` (default: `12`)
- `--uniform_qc_max_heatmap_keys` (default: `40`)

Normalized GeoJSON files are published to `results/uniformnormalize/`.
Normalized AnnData files are published to `results/uniformnormalize/`.
QC files are published to `results/uniformnormalize/qc/`.
