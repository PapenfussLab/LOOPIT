#!/usr/bin/env Rscript
# merge_coverage.R
# Merges per-chromosome coverage files into tumour and normal combined files
# Called by MERGE_COVERAGE Nextflow process

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list = list(
    make_option("--sample",  type="character", help="Sample name"),
    make_option("--infiles", type="character", help="Space-separated per-chrom RDS files"),
    make_option("--outfile", type="character", help="Output merged RDS file")
)

opt = parse_args(OptionParser(option_list = option_list))

infiles = strsplit(opt$infiles, " ")[[1]]
infiles = infiles[file.exists(infiles)]

message("Merging ", length(infiles), " chromosome coverage files for ", opt$sample)

chrom_order = paste0("chr", c(1:22, "X", "Y"))

#load all per-chrom RDS files
all_data = lapply(infiles, function(f) {
    message("Loading: ", f)
    readRDS(f)
})
#combine tumour across all chromosomes
tumour_files = lapply(all_data, function(x) x$tumour)
tumour_files = tumour_files[!sapply(tumour_files, is.null)]

if (length(tumour_files) == 0) stop("No tumour coverage files found")

tumour_combined = data.table::rbindlist(tumour_files, use.names = TRUE, fill = TRUE)
if (nrow(tumour_combined) == 0 || !"chrom" %in% names(tumour_combined)) {
    stop("Merged tumour coverage is empty or missing 'chrom' column")
}
tumour_combined = tumour_combined[chrom %in% chrom_order]
tumour_combined[, chrom := factor(chrom, levels = chrom_order)]
tumour_combined = tumour_combined[order(chrom, start)]
tumour_combined[, chrom := as.character(chrom)]

#combine normal across all chromosomes
normal_files = lapply(all_data, function(x) x$normal)
normal_files = normal_files[!sapply(normal_files, is.null)]

normal_combined = if (length(normal_files) > 0) {
    dt = data.table::rbindlist(normal_files, use.names = TRUE, fill = TRUE)
    if (nrow(dt) > 0 && "chrom" %in% names(dt)) {
        dt = dt[chrom %in% chrom_order]
        dt[, chrom := factor(chrom, levels = chrom_order)]
        dt = dt[order(chrom, start)]
        dt[, chrom := as.character(chrom)]
        dt
    } else {
        NULL
    }
} else {
    NULL
}

message("Tumour segments: ", nrow(tumour_combined))
message("Normal segments: ", if (!is.null(normal_combined)) nrow(normal_combined) else "NA")


writeLines(
  c(
    paste("Sample:", opt$sample),
    paste("N input files:", length(infiles)),
    paste("N files with non-null $tumour:", length(tumour_files)),
    paste("N files with non-null $normal:", length(normal_files)),
    paste("Tumour segments:", nrow(tumour_combined)),
    paste("Normal segments:", if (!is.null(normal_combined)) nrow(normal_combined) else "NA")
  ),
  file.path(dirname(opt$outfile), paste0(opt$sample, "_merge_coverage_debug.txt"))
)

message("Saved merged coverage to ", opt$outfile)
saveRDS(list(
    tumour = tumour_combined,
    normal = normal_combined
), opt$outfile)


