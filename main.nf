#!/usr/bin/env nextflow

include { UNIFORM_BATCH_CORRECTION } from './workflows/uniform_batch_correction'

workflow {
    def mode = (params.uniform_apply_to ?: 'geojson').toString()

    if (!params.input) {
        error "Please provide --input <samplesheet.csv>"
    }

    def ch_samplesheet = channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample) {
                error "Samplesheet must contain at least a sample column"
            }

            def geojson = row.geojson ? file(row.geojson.toString()) : null
            def omeTiff = row.ome_tiff ? file(row.ome_tiff.toString()) : null
            def adata = row.adata ? file(row.adata.toString()) : null

            if (mode == 'geojson' && !geojson) {
                error "Mode '${mode}' requires geojson column/path for sample ${row.sample}"
            }
            if (mode == 'ome_tiff' && !omeTiff) {
                error "Mode '${mode}' requires ome_tiff column/path for sample ${row.sample}"
            }
            if (mode == 'adata' && !adata) {
                error "Mode '${mode}' requires adata column/path for sample ${row.sample}"
            }

            if (!geojson && !omeTiff && !adata) {
                error "Each samplesheet row must include geojson and/or ome_tiff and/or adata"
            }

            [ [id: row.sample.toString()], geojson, omeTiff, adata ]
        }

    UNIFORM_BATCH_CORRECTION(
        ch_samplesheet
    )
}
