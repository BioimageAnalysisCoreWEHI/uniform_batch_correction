include { UNIFORMNORMALIZE } from '../modules/local/uniformnormalize/main.nf'

workflow UNIFORM_BATCH_CORRECTION {

    take:
    ch_samplesheet // channel: [ val(meta), path(geojson) ]

    main:

    ch_versions = Channel.empty()
    ch_annotations = ch_samplesheet.map { meta, geojson -> [meta, geojson] }
    ch_qc = Channel.empty()

    if (params.run_uniform) {
        def uniformSuffix = (params.uniform_output_suffix ?: '_uniform') as String

        UNIFORMNORMALIZE(
            ch_samplesheet.map { meta, geojson -> geojson }.collect()
        )
        ch_versions = ch_versions.mix(UNIFORMNORMALIZE.out.versions.first())
        ch_qc = UNIFORMNORMALIZE.out.qc

        ch_samplesheet
            .map { meta, geojson -> [meta.id.toString(), meta] }
            .join(
                UNIFORMNORMALIZE.out.normalized_annotations.map { geojson ->
                    def base = geojson.getBaseName()
                    def sample = base.endsWith(uniformSuffix)
                        ? base.substring(0, base.length() - uniformSuffix.length())
                        : base
                    [sample, geojson]
                },
                by: 0
            )
            .map { sample_id, meta, normalized_geojson ->
                [meta, normalized_geojson]
            }
            .set { ch_annotations }
    }

    emit:
    annotations = ch_annotations   // channel: [ val(meta), *.geojson ]
    qc          = ch_qc            // channel: [ qc files ]
    versions    = ch_versions      // channel: [ versions.yml ]
}
