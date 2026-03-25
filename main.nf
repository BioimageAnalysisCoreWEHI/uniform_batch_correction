#!/usr/bin/env nextflow

include { UNIFORM_BATCH_CORRECTION } from './workflows/uniform_batch_correction'

workflow {
    if (!params.input) {
        error "Please provide --input <samplesheet.csv>"
    }

    def ch_samplesheet = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample || !row.geojson) {
                error "Samplesheet must contain columns: sample,geojson"
            }
            [ [id: row.sample.toString()], file(row.geojson.toString()) ]
        }

    UNIFORM_BATCH_CORRECTION(
        ch_samplesheet
    )
}
