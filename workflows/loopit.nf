/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SETUP_REFERENCES } from '../modules/local/setup_references/main'
include { GET_REPEAT_REGIONS     } from '../modules/local/get_repeat_regions/main'
include { GET_REFERENCE_REGIONS  } from '../modules/local/get_reference_regions/main'
include { GET_BREAKPOINTS        } from '../modules/local/get_breakpoints/main'
include { GET_SEGMENTS           } from '../modules/local/get_segments/main'
include { SPLIT_BED        } from '../modules/local/split_bed/main'
include { COMPUTE_COVERAGE } from '../modules/local/compute_coverage/main'
include { MERGE_COVERAGE         } from '../modules/local/merge_coverage/main'
include { DETECT_CYCLES          } from '../modules/local/detect_cycles/main'
include { DRAW_CYCLES            } from '../modules/local/draw_cycles/main'

workflow LOOPIT {

    take:
    ch_samplesheet

    main:
    ch_versions = Channel.empty()

    //setup references
    if (!params.reference_dir || !file(params.reference_dir).exists()) {
        log.info "Reference directory not found — downloading reference data..."
        SETUP_REFERENCES()
        ch_reference_dir = SETUP_REFERENCES.out.reference_dir
        ch_versions      = ch_versions.mix(SETUP_REFERENCES.out.versions.first())
    } else {
        log.info "Using existing reference directory: ${params.reference_dir}"
        ch_reference_dir = Channel.value(file(params.reference_dir))
    }

    //process repeat mask
    GET_REPEAT_REGIONS(ch_reference_dir)
    ch_repeat_bed = GET_REPEAT_REGIONS.out.repeat_bed
    ch_versions   = ch_versions.mix(GET_REPEAT_REGIONS.out.versions.first())

    //get reference genome
    GET_REFERENCE_REGIONS(ch_reference_dir)
    ch_reference_regions = GET_REFERENCE_REGIONS.out.reference_regions
    ch_versions          = ch_versions.mix(GET_REFERENCE_REGIONS.out.versions.first())

    //get breakpoints
    ch_vcf = ch_samplesheet.map { meta, t_bed, n_bed, vcf, cnv -> [ meta, vcf ] }
    GET_BREAKPOINTS(ch_vcf.combine(ch_repeat_bed))
    ch_versions = ch_versions.mix(GET_BREAKPOINTS.out.versions.first())

    //get segments (breakpoints + reference regions)
    GET_SEGMENTS(GET_BREAKPOINTS.out.breakpoints.combine(ch_reference_regions))
    ch_versions = ch_versions.mix(GET_SEGMENTS.out.versions.first())

    //split BED files by chromosomes
    def chroms = (1..22).collect { "chr${it}" } + ["chrX", "chrY"]

    ch_beds_per_chrom = ch_samplesheet
        .map { meta, t_bed, n_bed, vcf, cnv -> [ meta, t_bed, n_bed ] }
        .combine(Channel.from(chroms))
        .map { meta, t_bed, n_bed, chrom -> [ meta, chrom, t_bed, n_bed ] }

    SPLIT_BED(ch_beds_per_chrom)
    ch_versions = ch_versions.mix(SPLIT_BED.out.versions.first())

    //compute coverage per chromosome 
    ch_blacklist = ch_reference_dir.map { d ->
        params.blacklist ? file(params.blacklist) : file("${d}/encode_blacklist.bed")
    }

    ch_chrom_input = SPLIT_BED.out.chrom_beds
    .combine(GET_SEGMENTS.out.segments, by: 0)
    .combine(ch_blacklist)
    .combine(ch_repeat_bed)
    .map { meta, chrom, t_bed, n_bed, segments, blacklist, repeat_bed ->
        [ meta + [chrom: chrom], t_bed, n_bed, segments, blacklist, repeat_bed ]
    }

    COMPUTE_COVERAGE(ch_chrom_input)
    ch_versions = ch_versions.mix(COMPUTE_COVERAGE.out.versions.first())

    //merge coverage across all chromosomes
    ch_merged_input = COMPUTE_COVERAGE.out.coverage
        .map { meta, file -> [ meta.subMap(['id','mode','purity','ploidy','sex']), file ] }
        .groupTuple(by: 0, size: 24)

    MERGE_COVERAGE(ch_merged_input)
    ch_versions = ch_versions.mix(MERGE_COVERAGE.out.versions.first())
    //detect ecDNA cycles
    ch_cnv = ch_samplesheet.map { meta, t_bed, n_bed, vcf, cnv ->
        [ meta.subMap(['id','mode','purity','ploidy','sex']), cnv ]
    }

    ch_detect_input = MERGE_COVERAGE.out.coverage
        .join(GET_BREAKPOINTS.out.breakpoints, by: 0)
        .join(ch_cnv,                          by: 0)

    DETECT_CYCLES(ch_detect_input)
    ch_versions = ch_versions.mix(DETECT_CYCLES.out.versions.first())

    //draw cycles
     ch_annotation = ch_reference_dir.map { d ->
        params.annotation
            ? file(params.annotation)
            : file("${d}/gencode.v42.primary_assembly.basic.annotation.gtf")
    }

    ch_draw_input = DETECT_CYCLES.out.results
        .map { meta, cycles_txt, decomposed_csv, rds, debug_log -> tuple(meta.subMap(['id','mode','purity','ploidy','sex']), meta, rds) }
        .join(ch_cnv, by: 0)
        .map { key, meta, rds, cnv -> tuple(meta, rds, cnv) }
        .combine(ch_annotation)
        .map { meta, rds, cnv, annotation -> tuple(meta, rds, annotation, cnv) }

    DRAW_CYCLES(ch_draw_input)
    ch_versions = ch_versions.mix(DRAW_CYCLES.out.versions.first())


emit:
    cycles   = DETECT_CYCLES.out.results
    plots    = DRAW_CYCLES.out.plots
    versions = ch_versions

}