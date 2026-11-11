/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { fromSamplesheet           } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: display version and exit
    validate_params   // boolean: validate parameters against schema at runtime
    monochrome_logs   // boolean: do not use coloured log outputs
    nextflow_cli_args //   array: list of positional nextflow CLI args
    outdir            //  string: output directory
    input             //  string: path to input samplesheet

    main:
    ch_versions = Channel.empty()

    // print version and exit if required
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile == 'conda' || workflow.profile == 'mamba'
    )

    // validate params
    UTILS_NFVALIDATION_PLUGIN(
    false,
    workflow.commandLine,
    null,
    null,
    validate_params,
    "${projectDir}/nextflow_schema.json"
    )

    // create channel from samplesheet
    Channel
        .fromSamplesheet("input")
        .map { meta, tumour_bed, normal_bed, vcf, cnv ->
            def new_meta = meta + [
                purity: meta.purity ?: 1.0,
                ploidy: meta.ploidy ?: 2.0
            ]
            // cnv is optional
            def cnv_file = cnv ? file(cnv) : []
            def normal_file = normal_bed  ? file(normal_bed)  : [] 
            [ new_meta, tumour_bed, normal_file, vcf, cnv_file ] 
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION TASKS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: send plain-text email instead of HTML
    outdir          //    path: path to output directory
    monochrome_logs // boolean: disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //     any: placeholder (not used — no MultiQC in LOOPIT)

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                []
            )
        }
        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}
