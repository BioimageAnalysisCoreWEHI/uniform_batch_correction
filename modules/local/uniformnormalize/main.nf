process UNIFORMNORMALIZE {
    tag "uniform"
    label 'process_multi'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/python_tifffile_scikit-image_scikit-learn_pruned:593e00ba324c12b3'

    input:
    path(geojsons)

    output:
    path "*_uniform.geojson"            , emit: normalized_annotations, optional: true
    path "*_uniform.h5ad"               , emit: normalized_adata, optional: true
    path "*_unifrom.tif"                , emit: normalized_images, optional: true
    path "*_uniform.tif"                , emit: normalized_images, optional: true
    path "*_unifrom.tiff"               , emit: normalized_images, optional: true
    path "*_uniform.tiff"               , emit: normalized_images, optional: true
    path "*_unifrom.ome.tif"            , emit: normalized_images_ome, optional: true
    path "*_uniform.ome.tif"            , emit: normalized_images_ome, optional: true
    path "*_unifrom.ome.tiff"           , emit: normalized_images_ome, optional: true
    path "*_uniform.ome.tiff"           , emit: normalized_images_ome, optional: true
    path "qc/*"                         , emit: qc
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def geojson_list = geojsons.join(' ')
    def uniform_script = "${projectDir}/bin/uniform_normalize_geojson.py"
    """
    python ${uniform_script} \
        --inputs ${geojson_list} \\
        --qc-dir qc \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    for f in ${geojsons.join(' ')}; do
        base=\$(basename "\$f" .geojson)
        touch "\${base}_uniform.geojson"
    done
    mkdir -p qc
    touch qc/uniform_key_summary.csv
    touch qc/uniform_run_summary.json
    touch qc/scale_factor_heatmap.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
