#!/usr/bin/env Rscript
# get_coverage_chrom.R
# Computes mean coverage for one chromosome for tumour AND normal


suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list = list(
    make_option("--sample",     type="character", help="Sample name"),
    make_option("--chrom",      type="character", help="Chromosome e.g. chr7"),
    make_option("--tumour-bed", type="character", help="Tumour per-chrom BED.gz"),
    make_option("--normal-bed", type="character", help="Normal per-chrom BED.gz"),
    make_option("--segments",   type="character", help="Whole-genome segment regions BED file"),
    make_option("--blacklist",  type="character", help="Blacklist BED file"),
    make_option("--repeat-bed", type="character", help="Repeat mask BED file"),
    make_option("--outfile",    type="character", help="Output RDS file path")
)

opt = parse_args(OptionParser(option_list = option_list))

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))
source(file.path(script_dir, "ecdna_method.R"))

message("Computing coverage for ", opt$sample, " ", opt$chrom)


bed_dir = "."


repeat_dest = file.path(bed_dir, "repeat_mask.bed")
if (!file.exists(repeat_dest)) {
    file.copy(opt[["repeat-bed"]], repeat_dest)
}


message("Splitting segments for ", opt$chrom)
all_segs = data.table::fread(opt$segments, header = FALSE, col.names = c("chr", "start", "end"))
chrom_segs = all_segs[chr == opt$chrom]

for (px in c(paste0(opt$sample, "_t"), paste0(opt$sample, "_n"))) {
    chrom_seg_file = file.path(bed_dir, paste0(px, "_", opt$chrom, "_segments.bed"))
    if (nrow(chrom_segs) > 0) {
        data.table::fwrite(chrom_segs, chrom_seg_file, sep = "\t", col.names = FALSE)
        message(opt$chrom, ": ", nrow(chrom_segs), " segments written to ", chrom_seg_file)
    } else {
        message(opt$chrom, ": no segments — tumour/normal coverage for this chromosome will be skipped")
    }
}
rm(all_segs)

message("Computing tumour coverage for ", opt$chrom)
tumour_path = get_coverage_single_chrom(
    prefix          = paste0(opt$sample, "_t"),
    chrom_bed       = opt[["tumour-bed"]],
    bed_dir         = bed_dir,
    segment_regions = opt$segments,
    blacklist_file  = opt$blacklist,
    chrom           = opt$chrom
)
tumour_result = if (!is.null(tumour_path) && file.exists(tumour_path)) {
    data.table::fread(tumour_path)
} else {
    NULL
}


has_normal = FALSE
if (!is.null(opt[["normal-bed"]]) && opt[["normal-bed"]] != "" && 
    opt[["normal-bed"]] != "NA" && file.exists(opt[["normal-bed"]])) {
    first_line = tryCatch(
        system(paste("zcat", shQuote(opt[["normal-bed"]]), "| head -1"), intern = TRUE),
        error = function(e) ""
    )
    n_fields = if (length(first_line) > 0 && nzchar(first_line)) length(strsplit(first_line, "\t")[[1]]) else 0
    has_normal = n_fields >= 4
}

writeLines(
  c(
    paste("Sample:", opt$sample),
    paste("Chrom:", opt$chrom),
    paste("normal-bed arg:", if (is.null(opt[["normal-bed"]])) "NULL" else opt[["normal-bed"]]),
    paste("file.exists check:", if (!is.null(opt[["normal-bed"]])) file.exists(opt[["normal-bed"]]) else NA),
    paste("file size:", if (!is.null(opt[["normal-bed"]]) && file.exists(opt[["normal-bed"]])) file.info(opt[["normal-bed"]])$size else NA),
    paste("first line content:", if (exists("first_line")) paste0("[", paste(first_line, collapse="|"), "]") else "not reached"),
    paste("n_fields:", if (exists("n_fields")) n_fields else NA),
    paste("has_normal final:", has_normal)
  ),
  file.path(".", paste0(opt$sample, "_", opt$chrom, "_has_normal_debug.txt"))
)

normal_result = NULL
if (has_normal) {
    message("Computing normal coverage for ", opt$chrom)
    normal_path = get_coverage_single_chrom(
        prefix          = paste0(opt$sample, "_n"),
        chrom_bed       = opt[["normal-bed"]],
        bed_dir         = bed_dir,
        segment_regions = opt$segments,
        blacklist_file  = opt$blacklist,
        chrom           = opt$chrom
    )
    normal_result = if (!is.null(normal_path) && file.exists(normal_path)) {
        data.table::fread(normal_path)
    } else {
        NULL
    }
} else {
    message(opt$chrom, ": no normal BED supplied — skipping normal coverage computation")
}
saveRDS(list(
    chrom  = opt$chrom,
    tumour = tumour_result,  # path to tumour regions_mean.bed for this chrom
    normal = normal_result   # path to normal regions_mean.bed for this chrom
), opt$outfile)

message("Saved coverage paths for ", opt$chrom, " to ", opt$outfile)
