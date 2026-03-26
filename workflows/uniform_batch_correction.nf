include { UNIFORMNORMALIZE as UNIFORMNORMALIZE_GEOJSON } from '../modules/local/uniformnormalize/main.nf'
include { UNIFORMNORMALIZE as UNIFORMNORMALIZE_PIXEL   } from '../modules/local/uniformnormalize/main.nf'
include { UNIFORMNORMALIZE as UNIFORMNORMALIZE_ADATA   } from '../modules/local/uniformnormalize/main.nf'

workflow UNIFORM_BATCH_CORRECTION {

    take:
    ch_samplesheet // channel: [ val(meta), path(geojson?) , path(ome_tiff?) , path(adata?) ]

    main:

    ch_versions = channel.empty()
    def mode = (params.uniform_apply_to ?: 'geojson').toString()

    ch_geojson_input = ch_samplesheet
        .filter { _meta, geojson, _ome_tiff, _adata -> geojson != null }
        .map { meta, geojson, _ome_tiff, _adata -> [meta, geojson] }

    ch_ome_tiff_input = ch_samplesheet
        .filter { _meta, _geojson, ome_tiff, _adata -> ome_tiff != null }
        .map { meta, _geojson, ome_tiff, _adata -> [meta, ome_tiff] }

    ch_adata_input = ch_samplesheet
        .filter { _meta, _geojson, _ome_tiff, adata -> adata != null }
        .map { meta, _geojson, _ome_tiff, adata -> [meta, adata] }

    ch_annotations = ch_geojson_input
    ch_images = ch_ome_tiff_input
    ch_adatas = ch_adata_input
    ch_qc = channel.empty()

    if (params.run_uniform && mode == 'geojson') {
        def uniformSuffix = (params.uniform_output_suffix ?: '_uniform') as String

        UNIFORMNORMALIZE_GEOJSON(
            ch_geojson_input.map { _meta, geojson -> geojson }.collect()
        )
        ch_versions = ch_versions.mix(UNIFORMNORMALIZE_GEOJSON.out.versions.first())
        ch_qc = UNIFORMNORMALIZE_GEOJSON.out.qc

        ch_geojson_input
            .map { meta, _geojson -> [meta.id.toString(), meta] }
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
            .map { _sample_id, meta, normalized_geojson ->
                [meta, normalized_geojson]
            }
            .set { ch_annotations }
    }

    if (params.run_uniform && mode == 'adata') {
        def uniformSuffix = (params.uniform_output_suffix ?: '_uniform') as String

        UNIFORMNORMALIZE_ADATA(
            ch_adata_input.map { _meta, adata -> adata }.collect()
        )
        ch_versions = ch_versions.mix(UNIFORMNORMALIZE_ADATA.out.versions.first())
        ch_qc = ch_qc.mix(UNIFORMNORMALIZE_ADATA.out.qc)

        ch_adata_input
            .map { meta, _adata -> [meta.id.toString(), meta] }
            .join(
                UNIFORMNORMALIZE_ADATA.out.normalized_adata.map { adata ->
                    def base = adata.getBaseName()
                    def sample = base.endsWith(uniformSuffix)
                        ? base.substring(0, base.length() - uniformSuffix.length())
                        : base
                    [sample, adata]
                },
                by: 0
            )
            .map { _sample_id, meta, normalized_adata ->
                [meta, normalized_adata]
            }
            .set { ch_adatas }
    }

    if (params.run_uniform && mode == 'ome_tiff') {
        def pixelSuffix = (params.uniform_pixel_output_suffix ?: '_unifrom') as String

        UNIFORMNORMALIZE_PIXEL(
            ch_ome_tiff_input.map { _meta, ome_tiff -> ome_tiff }.collect()
        )
        ch_versions = ch_versions.mix(UNIFORMNORMALIZE_PIXEL.out.versions.first())
        ch_qc = ch_qc.mix(UNIFORMNORMALIZE_PIXEL.out.qc)

        ch_normalized_images = UNIFORMNORMALIZE_PIXEL.out.normalized_images_tif_unifrom
            .mix(UNIFORMNORMALIZE_PIXEL.out.normalized_images_tif_uniform)
            .mix(UNIFORMNORMALIZE_PIXEL.out.normalized_images_tiff_unifrom)
            .mix(UNIFORMNORMALIZE_PIXEL.out.normalized_images_tiff_uniform)
            .mix(UNIFORMNORMALIZE_PIXEL.out.normalized_images_ome_tif_unifrom)
            .mix(UNIFORMNORMALIZE_PIXEL.out.normalized_images_ome_tif_uniform)
            .mix(UNIFORMNORMALIZE_PIXEL.out.normalized_images_ome_tiff_unifrom)
            .mix(UNIFORMNORMALIZE_PIXEL.out.normalized_images_ome_tiff_uniform)

        ch_ome_tiff_input
            .map { meta, _ome_tiff -> [meta.id.toString(), meta] }
            .join(
                ch_normalized_images.map { image ->
                    def name = image.getName()
                    def sample = name
                    def endings = [
                        "${pixelSuffix}.ome.tiff",
                        "${pixelSuffix}.ome.tif",
                        "${pixelSuffix}.tiff",
                        "${pixelSuffix}.tif"
                    ]

                    def matchedEnding = endings.find { ending -> sample.endsWith(ending) }
                    if (matchedEnding) {
                        sample = sample.substring(0, sample.length() - matchedEnding.length())
                    }

                    [sample, image]
                },
                by: 0
            )
            .map { _sample_id, meta, normalized_image ->
                [meta, normalized_image]
            }
            .set { ch_images }
    }

    emit:
    annotations = ch_annotations   // channel: [ val(meta), *.geojson ]
    images      = ch_images        // channel: [ val(meta), *.tiff/*.ome.tiff ]
    adata       = ch_adatas        // channel: [ val(meta), *.h5ad ]
    qc          = ch_qc            // channel: [ qc files ]
    versions    = ch_versions      // channel: [ versions.yml ]
}
