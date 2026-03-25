# uniform_batch_correction

Standalone Nextflow pipeline for UniFORM-style batch normalization of cellmeasurement GeoJSON outputs.

## Input

Provide a CSV/YAML samplesheet with:

- `sample`: sample ID (must match GeoJSON filename stem for normalized mapping)
- `geojson`: path to input cellmeasurement GeoJSON

Example CSV is provided in [assets/samplesheet.csv](assets/samplesheet.csv).

## Usage

```bash
nextflow run . \
  -profile conda \
  --input assets/samplesheet.csv \
  --outdir results
```

## Key parameters

- `--run_uniform` (default: `true`)
- `--uniform_num_bins` (default: `1024`)
- `--uniform_min_value` (default: `1.0`)
- `--uniform_exclude_pattern` (default: `^(kronos_|emb_)`)
- `--uniform_output_suffix` (default: `_uniform`)
- `--uniform_generate_plots` (default: `true`)
- `--uniform_qc_top_n_keys` (default: `12`)
- `--uniform_qc_max_heatmap_keys` (default: `40`)

Normalized GeoJSON files are published to `results/uniformnormalize/`.
QC files are published to `results/uniformnormalize/qc/`.
