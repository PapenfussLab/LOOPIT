process SPLIT_BED {
    tag "${meta.id}_${chrom}"
    label 'process_low'

    input:
    tuple val(meta), val(chrom), path(tumour_bed), path(normal_bed, stageAs: "normal/*")

    output:
    tuple val(meta), val(chrom), path("tumour_chrom/${chrom}.bed.gz"), path("normal_chrom/${chrom}.bed.gz"), emit: chrom_beds
    path "versions.yml", emit: versions

    script:
    def has_normal = normal_bed && normal_bed.name != 'null'
    """
    if [ ! -f ${tumour_bed}.tbi ]; then
        tabix -0 -p bed ${tumour_bed}
    fi

    mkdir -p tumour_chrom normal_chrom

    tabix -0 ${tumour_bed} ${chrom} | bgzip > tumour_chrom/${chrom}.bed.gz
    tabix -0 -p bed tumour_chrom/${chrom}.bed.gz

    ${has_normal ? """
    if [ ! -f ${normal_bed}.tbi ]; then
        tabix -0 -p bed ${normal_bed}
    fi

    tabix -0 ${normal_bed} ${chrom} | bgzip > normal_chrom/${chrom}.bed.gz
    tabix -0 -p bed normal_chrom/${chrom}.bed.gz
    """ : """
    echo "" | bgzip > normal_chrom/${chrom}.bed.gz
    """}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(tabix --version 2>&1 | head -1 | sed 's/tabix (htslib) //')
    END_VERSIONS
    """
}