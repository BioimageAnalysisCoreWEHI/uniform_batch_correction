# Usage

Run the pipeline with a samplesheet containing `sample` and `geojson` columns.

```bash
nextflow run . \
  -profile conda \
  --input /path/to/samplesheet.csv \
  --outdir /path/to/results
```

Optional tuning:

```bash
--uniform_num_bins 1024 \
--uniform_min_value 1.0 \
--uniform_exclude_pattern '^(kronos_|emb_)' \
--uniform_output_suffix _uniform \\
--uniform_pixel_output_suffix _uniform \\
--uniform_pixel_sample_size 200000 \\
--uniform_generate_plots true \\
--uniform_qc_top_n_keys 12 \\
--uniform_qc_max_heatmap_keys 40
```

Pixel-level mode (OME-TIFF/TIFF):

```bash
nextflow run . \\
  -profile conda \\
  --input /path/to/samplesheet.csv \\
  --outdir /path/to/results \\
  --uniform_apply_to ome_tiff
```
