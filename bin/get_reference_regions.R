#!/usr/bin/env Rscript
# get_reference_regions.R
# Computes reference genome windows from FASTA
# Called by GET_REFERENCE_REGIONS Nextflow process

suppressPackageStartupMessages(library(optparse))

option_list = list(
    make_option("--reference",     type="character", help="Reference FASTA file"),
    make_option("--reference-dir", type="character", help="Reference directory"),
    make_option("--outfile",       type="character", help="Output RDS file")
)

opt = parse_args(OptionParser(option_list = option_list))

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))
source(file.path(script_dir, "ecdna_method.R"))

message("Computing reference regions from: ", opt$reference)

outdir = dirname(opt$outfile)
config = list(
    output_dir    = outdir,
    manifest_path = file.path(outdir, "manifest.csv")
)

result = get_reference_regions(
    reference_fasta  = opt$reference,
    output_directory = outdir
)

saveRDS(result, opt$outfile)
message("Reference regions saved to ", opt$outfile)
