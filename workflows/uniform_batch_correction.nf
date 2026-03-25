include { UNIFORMNORMALIZE as UNIFORMNORMALIZE_GEOJSON } from '../modules/local/uniformnormalize/main.nf'
include { UNIFORMNORMALIZE as UNIFORMNORMALIZE_PIXEL   } from '../modules/local/uniformnormalize/main.nf'

workflow UNIFORM_BATCH_CORRECTION {

    take:
    ch_samplesheet // channel: [ val(meta), path(geojson?) , path(ome_tiff?) ]

    main:

    ch_versions = Channel.empty()
    def mode = (params.uniform_apply_to ?: 'geojson').toString()

    ch_geojson_input = ch_samplesheet
        .filter { meta, geojson, ome_tiff -> geojson != null }
        .map { meta, geojson, ome_tiff -> [meta, geojson] }

    ch_ome_tiff_input = ch_samplesheet
        .filter { meta, geojson, ome_tiff -> ome_tiff != null }
        .map { meta, geojson, ome_tiff -> [meta, ome_tiff] }

    ch_annotations = ch_geojson_input
    ch_images = ch_ome_tiff_input
    ch_qc = Channel.empty()

    if (params.run_uniform && (mode == 'geojson' || mode == 'both')) {
        def uniformSuffix = (params.uniform_output_suffix ?: '_uniform') as String

        UNIFORMNORMALIZE_GEOJSON(
            ch_geojson_input.map { meta, geojson -> geojson }.collect()
        )
        ch_versions = ch_versions.mix(UNIFORMNORMALIZE_GEOJSON.out.versions.first())
        ch_qc = UNIFORMNORMALIZE_GEOJSON.out.qc

        ch_geojson_input
            .map { meta, geojson -> [meta.id.toString(), meta] }
            .join(
                UNIFORMNORMALIZE_GEOJSON.out.normalized_annotations.map { geojson ->
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

    if (params.run_uniform && (mode == 'ome_tiff' || mode == 'both')) {
        def pixelSuffix = (params.uniform_pixel_output_suffix ?: '_uniform') as String

        UNIFORMNORMALIZE_PIXEL(
            ch_ome_tiff_input.map { meta, ome_tiff -> ome_tiff }.collect()
        )
        ch_versions = ch_versions.mix(UNIFORMNORMALIZE_PIXEL.out.versions.first())
        ch_qc = ch_qc.mix(UNIFORMNORMALIZE_PIXEL.out.qc)

        ch_normalized_images = UNIFORMNORMALIZE_PIXEL.out.normalized_images.mix(UNIFORMNORMALIZE_PIXEL.out.normalized_images_ome)

        ch_ome_tiff_input
            .map { meta, ome_tiff -> [meta.id.toString(), meta] }
            .join(
                ch_normalized_images.map { image ->
                    def base = image.getBaseName()
                    def sample = base.endsWith(pixelSuffix)
                        ? base.substring(0, base.length() - pixelSuffix.length())
                        : base
                    [sample, image]
                },
                by: 0
            )
            .map { sample_id, meta, normalized_image ->
                [meta, normalized_image]
            }
            .set { ch_images }
    }

    emit:
    annotations = ch_annotations   // channel: [ val(meta), *.geojson ]
    images      = ch_images        // channel: [ val(meta), *.tiff/*.ome.tiff ]
    qc          = ch_qc            // channel: [ qc files ]
    versions    = ch_versions      // channel: [ versions.yml ]
}
