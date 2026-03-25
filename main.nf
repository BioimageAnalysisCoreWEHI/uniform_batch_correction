#!/usr/bin/env nextflow

include { UNIFORM_BATCH_CORRECTION } from './workflows/uniform_batch_correction'

workflow {
    def mode = (params.uniform_apply_to ?: 'geojson').toString()

    if (!params.input) {
        error "Please provide --input <samplesheet.csv>"
    }

    def ch_samplesheet = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample) {
                error "Samplesheet must contain at least a sample column"
            }

            def geojson = row.geojson ? file(row.geojson.toString()) : null
            def omeTiff = row.ome_tiff ? file(row.ome_tiff.toString()) : null

            if ((mode == 'geojson' || mode == 'both') && !geojson) {
                error "Mode '${mode}' requires geojson column/path for sample ${row.sample}"
            }
            if ((mode == 'ome_tiff' || mode == 'both') && !omeTiff) {
                error "Mode '${mode}' requires ome_tiff column/path for sample ${row.sample}"
            }

            if (!geojson && !omeTiff) {
                error "Each samplesheet row must include geojson and/or ome_tiff"
            }

            [ [id: row.sample.toString()], geojson, omeTiff ]
        }

    UNIFORM_BATCH_CORRECTION(
        ch_samplesheet
    )
}
