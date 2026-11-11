library(data.table)
library(gGnome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(gTrack)
library(stringr)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(tidyverse)
library(StructuralVariantAnnotation)
library(dplyr)
library(fastmap)
library(rtracklayer)
library(optparse)
library(this.path)
library(igraph)

gr_bp_to_bed = function(gr){
  bedpe = data.frame(chrom = GenomeInfoDb::seqnames(gr),
                     start = floor((start(gr)+end(gr))/2) - 1,
                     end = floor((start(gr)+end(gr))/2),
                     strand = strand(gr))
  return(bedpe)
}

bp_to_bedpe = function(gr,mode,vaf_field){
  if(mode=="SR"){
    bedpe = data.frame(  chrom1=seqnames(gr),
                         start1=floor((start(gr)+end(gr))/2) - 1,
                         end1=floor((start(gr)+end(gr))/2),
                         chrom2=seqnames(partner(gr)),
                         start2=floor((start(partner(gr))+end(partner(gr)))/2) - 1,
                         end2=floor((start(partner(gr))+end(partner(gr)))/2),
                         name=names(gr),
                         partner.name = names(partner(gr)),
                         score=gr$QUAL,
                         strand1=strand(gr),
                         strand2=strand(partner(gr)),
                         TAF=mcols(gr)[[vaf_field]]) 
  }else{
    bedpe = data.frame(  chrom1 = seqnames(gr),
                         start1 = floor((start(gr)+end(gr))/2) - 1,
                         end1 = floor((start(gr)+end(gr))/2),
                         chrom2 = seqnames(partner(gr)),
                         start2 = floor((start(partner(gr))+end(partner(gr)))/2) - 1,
                         end2 = floor((start(partner(gr))+end(partner(gr)))/2),
                         name = names(gr),
                         partner.name = names(partner(gr)),
                         score = gr$QUAL,
                         strand1 = strand(gr),
                         strand2 = strand(partner(gr)),
                         VAF = mcols(gr)[[vaf_field]]
    )
  }
  
  return(bedpe)
}


get_coverage_single_chrom = function(prefix, chrom_bed, bed_dir,
                                     segment_regions, blacklist_file, chrom) {
  
    repeat_file              = file.path(bed_dir, "repeat_mask.bed")
    bed_bl_repeat_removed_gz = file.path(bed_dir, paste0(prefix, "_", chrom, "_bl_repeat_removed_regions.bed.gz"))  
 

    regions_to_remove        = file.path(bed_dir, paste0(prefix, "_", chrom, "_removed_regions.bed"))
    
    blacklist_df = data.table::fread(blacklist_file)
    repeat_df    = data.table::fread(repeat_file)
    
    subtract_df = data.table::rbindlist(list(blacklist_df, repeat_df))[
      V1 == chrom  #only this chromosome
    ] |> unique() |> _[order(V1, V2)]
    
    data.table::fwrite(subtract_df, regions_to_remove, sep = "\t", col.names = FALSE)
    
    message("Removing blacklist and repeat regions for ", chrom)
    exit_code = system(paste(
      "bedtools subtract -a", chrom_bed, "-b", regions_to_remove,
      "| bgzip >", bed_bl_repeat_removed_gz,
      "&& 
      tabix -p bed", bed_bl_repeat_removed_gz
    ))
    
    if (exit_code != 0 || !file.exists(bed_bl_repeat_removed_gz)) {
      stop("bedtools subtract failed for ", chrom)
    }  
  
  chrom_segments  = file.path(bed_dir, paste0(prefix, "_", chrom, "_segments.bed"))
  if (!file.exists(chrom_segments) || file.size(chrom_segments) == 0) {
    message("No segments for ", chrom, " — skipping")
    return(NULL)
  }

  bed_filtered = file.path(bed_dir, paste0(prefix, "_", chrom, "_filtered_regions.bed.gz"))

  exit_code = system(paste(
      "bedtools intersect -sorted -wa -wb",
      "-a", chrom_segments,
      "-b", bed_bl_repeat_removed_gz,
      "| bgzip >", bed_filtered,
      "&& tabix -p bed", bed_filtered
    ))
    
    if (exit_code != 0 || !file.exists(bed_filtered)) {
      stop("bedtools intersect failed for ", chrom)
    }

    regions_mean = file.path(bed_dir, paste0(prefix, "_", chrom, "_regions_mean.bed"))

    message("Calculating mean coverage for ", chrom)
    covered_tmp = file.path(bed_dir, paste0(prefix, "_", chrom, "_covered_means.tsv"))
    
    exit_code = system(paste(
      "zcat", bed_filtered,
      "| awk 'BEGIN{OFS=\"\\t\"} {w=$6-$5; print $1,$2,$3, $7*w, w}'",
      "| bedtools groupby -g 1,2,3 -c 4,5 -o sum,sum",
      "| awk 'BEGIN{OFS=\"\\t\"} {print $1,$2,$3, $4/$5}'",
      ">", covered_tmp
    ))
    
    if (exit_code != 0) {
      stop("Coverage mean calculation failed for ", chrom)
    }
    
    if (!file.exists(covered_tmp) || file.size(covered_tmp) == 0) {
      message("WARNING: empty coverage for ", chrom, " — all regions masked")
      cov_per_region = data.table::data.table(
        chrom  = character(),
        start  = integer(),
        end    = integer(),
        tcn.em = numeric()
      )
    } else {
      cov_per_region = data.table::fread(covered_tmp,
                                         col.names = c("chrom", "start", "end", "tcn.em")
      )
    }
    
    #restore segments lost entirely to masking
    segments = data.table::fread(chrom_segments, select = c(1, 2, 3))
    data.table::setnames(segments, c("V1", "V2", "V3"), c("chrom", "start", "end"))
    
    missing = segments[
      !cov_per_region, on = c("chrom", "start", "end")
    ][, tcn.em := -1]
    
    cov_per_region = data.table::rbindlist(list(cov_per_region, missing))
    data.table::fwrite(cov_per_region, regions_mean, sep = ",", col.names = TRUE)
   

  
  return(regions_mean)
}


get_repeat_regions = function(reference_directory,output_directory){
  REPEAT_FILE     = file.path(reference_directory, "repeat_mask.tsv")
  REPEAT_BED_FILE = file.path(output_directory,    "repeat_mask.bed")
  print("Getting the repeat regions of the reference")
  repeat_df = data.table::fread(
    input   = REPEAT_FILE,
    select  = c("genoName", "genoStart", "genoEnd", "repClass")  
  )[repClass %in% c("Simple_repeat", "Satellite", "Low_complexity") 
  ][, repClass := NULL]                                          
  data.table::fwrite(x = repeat_df, file = REPEAT_BED_FILE, sep = "\t", col.names = FALSE)
  return (REPEAT_BED_FILE)
}



get_reference_regions = function(reference_fasta, output_directory){
  CHROM_NAMES = paste0("chr", c(1:22, "X", "Y"))
  ref_filename = tools::file_path_sans_ext(basename(reference_fasta)) 
  ref_bed = file.path(output_directory, paste0(ref_filename, ".bed"))
    #get chromosome regions of the reference
    chr_reg = data.table::fread(
      input  = paste0(reference_fasta, ".fai"),
      select = c(1, 2),                          
      col.names = c("chrom", "end")
    )[chrom %in% CHROM_NAMES                    
    ][, start := 0                            
    ][, .(chrom, start, end)]                    
    
    data.table::fwrite(chr_reg, file = ref_bed, sep = "\t", col.names = FALSE)
    

  return(list(
    region_ct = nrow(chr_reg),
    chr_regions =  chr_reg))
}



get_breakpoints = function(sample_name, vcf_file, mode, repeat_file, output_directory,chromosomes = NA) {
  CHROM_NAMES = paste0("chr", c(1:22, "X", "Y"))
  message("get_breakpoints starting for sample: ", sample_name)
  message("mode: ", mode)
  
  
  message("Reading VCF file: ", vcf_file)
  
  vcf = readVcf(vcf_file)
  
  VAF_FIELD = {
    vcf_header = VariantAnnotation::header(vcf)
    info_fields = rownames(VariantAnnotation::info(vcf_header))
    if (mode == "SR" && "TAF" %in% info_fields) "TAF"
    else if ("VAF" %in% info_fields) "VAF"
    else if ("AF" %in% info_fields) "AF"
    else "VAF"
  }
  
  bp_for_bed = breakpointRanges(vcf)
  bp         = data.table::as.data.table(gr_bp_to_bed(bp_for_bed))
  
  bp[, start := as.integer(floor(start))]
  bp[, end   := as.integer(floor(end))]
  
  bp = bp[chrom %in% CHROM_NAMES & start >= 0L & end >= 0L]
 
 
  bp_no_infer = breakpointRanges(vcf, info_columns = c(VAF_FIELD, "EVENTTYPE"))
  bp_full = breakpointRanges(vcf,
                             info_columns          = c(VAF_FIELD, "EVENTTYPE"),
                             inferMissingBreakends = TRUE
  )
  
  
  
  #breakpoints        = data.table::as.data.table(bp_to_bedpe(bp_full, mode, VAF_FIELD))
  #breakpoints_df_raw = as.data.frame(bp_full)
  breakpoints = data.table::as.data.table(bp_to_bedpe(bp_full, mode, VAF_FIELD))
  breakpoints_df_raw = as.data.frame(bp_no_infer)
  source_ids         = rownames(breakpoints_df_raw)   
  
  breakpoints_df     = data.table::as.data.table(breakpoints_df_raw)
  breakpoints_df[, start := floor((start + end) / 2L)][, end := start]
  
  if (mode == "LR") {
    breakpoints_df[, sourceId := source_ids]
  }
  
  
  
  if (VAF_FIELD %in% colnames(breakpoints_df)) {
    VAF_COL = VAF_FIELD
  } else {
    #try to find column matching VAF_FIELD prefix with .value suffix
    matches = grep(paste0("^", VAF_FIELD, ".*\\.value$"),
                   colnames(breakpoints_df), value = TRUE)
    if (length(matches) == 0) {
    #fall back to any column starting with VAF_FIELD
      matches = grep(paste0("^", VAF_FIELD),
                     colnames(breakpoints_df), value = TRUE)
    }
    if (length(matches) == 0) {
      stop("Cannot find VAF column matching '", VAF_FIELD,
           "' in breakpoints_df. Available columns: ",
           paste(colnames(breakpoints_df), collapse=", "))
    }
    VAF_COL = matches[1]
    
  }
  
  
  breakpoints_df[, VAF_num := sapply(
    get(VAF_COL),
    function(x) as.numeric(as.character(unlist(x)[[1]]))
  )]
  
  
  chrom_suffix = if (identical(chromosomes, NA)) "wg" else paste(chromosomes, collapse="_")
  
  sample_bed  = file.path(output_directory, paste0(sample_name, "_", chrom_suffix, ".bed"))
  data.table::fwrite(
      x         = bp[, .(chrom, start, end, strand)],
      file      = sample_bed,
      sep       = "\t",
      col.names = FALSE
    )
 
 bp_no_repeat_bed_file = file.path(output_directory, paste0(sample_name, "_", chrom_suffix, "_bp_no_repeat.bed"))
 if (nrow(bp) == 0L) {
  message("No breakpoints for sample: ", sample_name, " — skipping repeat-region subtraction")
  bp_no_repeat = data.table::data.table(chrom=character(), start=integer(), end=integer(), strand=character())
  } else { 
    message("Removing breakpoints overlapping repeat regions")
    exit_code = system(paste(
      "bedtools subtract -a", sample_bed, "-b", repeat_file,
      ">", bp_no_repeat_bed_file
    ))
    if (exit_code != 0 || !file.exists(bp_no_repeat_bed_file) || file.size(bp_no_repeat_bed_file) == 0) {
      stop("bedtools subtract failed or produced empty output for sample: ", sample_name,
          " (sample_bed=", sample_bed, ", repeat_file=", repeat_file, ")")
    }
      
      bp_no_repeat = data.table::fread(
        input     = bp_no_repeat_bed_file,
        col.names = c("chrom", "start", "end", "strand")
      )
  }
  #print(head(bp_no_repeat))
  #print(head(breakpoints_df[, .(seqnames, start, end, strand)]))
  
  breakpoints_df = breakpoints_df[
    bp_no_repeat,
    on      = .(seqnames = chrom, start = end, strand),
    nomatch = NULL
  ]
  breakpoints_df[, svLen := abs(svLen)]
  
  breakpoints = breakpoints[name %in% breakpoints_df$partner]
  
  valid_ids   = breakpoints_df$sourceId
  partner_ids = breakpoints_df$partner
  
  breakpoints_df = breakpoints_df[sourceId %in% partner_ids & partner %in% valid_ids]
  breakpoints    = breakpoints[name %in% breakpoints_df$partner]
  
  valid_ids_final   = breakpoints_df$sourceId
  partner_ids_final = breakpoints_df$partner
  bed = bp_no_repeat[!(end %in% breakpoints_df[!sourceId %in% partner_ids_final, start])]
  
  
  
  bed = unique(bed, by = c("chrom", "start", "end"))
  data.table::setorder(bed, chrom, start, end)
  
  return(list(
    breakpoint_bed = bed,
    breakpoint_df  = breakpoints_df
  ))
}

get_segments = function(region_ct, chr_regions, breakpoints, sample_name, output_directory,chromosomes = NA){
  
  
    chrom_suffix = if (identical(chromosomes, NA)) "wg" else paste(chromosomes, collapse="_")
    bed_regions = file.path(output_directory, paste0(sample_name, "_", chrom_suffix, "_regions.bed"))
    chr_regions = data.table::as.data.table(chr_regions)
    breakpoints = data.table::as.data.table(breakpoints)
    
    
    bp_by_chrom = split(breakpoints, by = "chrom")
    
    segments_list = vector("list", region_ct)
    
    for (i in seq_len(region_ct)) {
      chrom     = chr_regions$chrom[i]
      chrom_end = chr_regions$end[i]
      bp_df     = bp_by_chrom[[chrom]]
      
      if (!is.null(bp_df) && nrow(bp_df) > 0) {
        seg_starts = c(0L,          bp_df$end)
        seg_ends   = c(bp_df$start, chrom_end)
      } else {
        seg_starts = 0L
        seg_ends   = chrom_end
      }
      
      segments_list[[i]] = data.table::data.table(
        chrom = chrom,
        start = seg_starts,
        end   = seg_ends
      )
    }
    
    segments_df = data.table::rbindlist(segments_list)
    segments_df = unique(segments_df)
    
    data.table::fwrite(segments_df, file = bed_regions, sep = "\t", col.names = FALSE)
 
  
    return(bed_regions)
  
}


normalize_coverage = function(cov) {
  is_sex = cov$chr %in% c("X", "Y", "chrX", "chrY")
  tn     = cov[!is_sex]
  tnxy   = cov[is_sex]
  
  tn[, chr := as.numeric(gsub("chr", "", chr, fixed = TRUE))]
  data.table::setorder(tn, chr, start)
  
  data.table::rbindlist(list(tn, tnxy), use.names = TRUE)
}


calculate_cn = function(normal_cov, tumor_cov){
  has_normal = !is.null(normal_cov)   
  
  tum_cov = tumor_cov[, c(1,2,3,4), with = FALSE]
  data.table::setnames(tum_cov, c("chr", "start", "end", "doc"))
  if (has_normal) {
    normal_cov = normal_cov[, c(1,2,3,4), with = FALSE]
    data.table::setnames(normal_cov, c("chr", "start", "end", "doc"))
  }
  
  
  contigs    = c(1:22, "X", "Y", paste0("chr", 1:22), "chrX", "chrY")
  tum_cov    = tum_cov[chr %in% contigs]
  if (has_normal) {
    normal_cov = normal_cov[chr %in% contigs]
  }
  
  
  tum_cov = normalize_coverage(tum_cov)
  if (has_normal) {
    normal_cov = normalize_coverage(normal_cov)
  }
  
  #med_cov = median(tum_cov[!chr %in% c("X", "Y", "chrX", "chrY"), doc], na.rm = TRUE)
  tum_autosome = tum_cov[!chr %in% c("X", "Y", "chrX", "chrY")]
  med_cov = weighted_median(tum_autosome$doc, tum_autosome$end - tum_autosome$start)
  
  data.table::setnames(tum_cov,
                       old = c("chr", "start", "end"),
                       new = c("Chromosome", "Start_Position", "End_Position")
  )
  
  
  chr_lens_dt = tum_cov[, .(len = max(End_Position, na.rm = TRUE)), by = Chromosome]
  chr_lens    = setNames(chr_lens_dt$len, chr_lens_dt$Chromosome)
  
  
  if (has_normal) {
    data.table::setnames(normal_cov,
                         old = c("chr", "start", "end"),
                         new = c("Chromosome", "Start_Position", "End_Position")
    )
    
    
    data.table::setkey(tum_cov,    Chromosome, Start_Position, End_Position)
    data.table::setkey(normal_cov, Chromosome, Start_Position, End_Position)
    
    all_depth = normal_cov[tum_cov,
                           on = c("Chromosome", "Start_Position", "End_Position")
    ]
    
    data.table::setnames(all_depth, "i.doc", "doc_t")
    
    autosome_mask = !all_depth$Chromosome %in% c("X", "Y", "chrX", "chrY")
    ratios    = all_depth$doc_t[autosome_mask] / all_depth$doc[autosome_mask]
    widths    = all_depth$End_Position[autosome_mask] - all_depth$Start_Position[autosome_mask]
    map_ratio = weighted_median(ratios, widths)

    #map_ratio = median(ratios[is.finite(ratios)], na.rm = TRUE)
    message("Coverage ratio T/N: ", round(map_ratio, digits = 3))
    
    all_depth[, doc   := doc * map_ratio        
    ][, logR  := log2(doc_t + 1) - log2(doc + 1) 
    ][, doc   := NULL]                            
    
  } else {
    all_depth = tum_cov
    all_depth[, logR := log2(doc) - log2(med_cov)]
  }
  
  
  
  all_depth[, Chromosome := data.table::fcase(
    Chromosome %in% c("X", "chrX"), "chrX",
    Chromosome %in% c("Y", "chrY"), "chrY",
    !startsWith(as.character(Chromosome), "chr"), 
    paste0("chr", Chromosome),
    rep(TRUE, .N), as.character(Chromosome)
  )]
  
  return(all_depth)
}


filter_chr_x = function(cn_df, x_label, reference){
  #gets all chr x regions without PAR1/2
  #what if the segment includes PAR region (larger than PAR)  
PAR_REGIONS = data.table::data.table(
  ref    = c("grch38", "grch38", "grch37", "grch37"),
  region = c("PAR1X",  "PAR2X",  "PAR1X",  "PAR2X"),
  start  = c(10000L,   155701382L, 60000L,  154931043L),
  end    = c(2781479L, 156030895L, 2699520L, 155260560L)
)


  par = PAR_REGIONS[ref == reference]
  par1 = par[region == "PAR1X"]
  par2 = par[region == "PAR2X"]
  
  
  cn_df[Chromosome == x_label &
          !(Start_Position >= par1$start & End_Position <= par1$end) &
          !(Start_Position >= par2$start & End_Position <= par2$end)]
  
}


filter_chr_y = function(cov, y_label, reference) {
PAR_REGIONS_Y = data.table::data.table(
  ref    = c("grch38", "grch38", "grch37", "grch37"),
  region = c("PAR1Y",  "PAR2Y",  "PAR1Y",  "PAR2Y"),
  start  = c(10000L,   56887902L, 10000L,  59034049L),
  end    = c(2781479L, 57217415L, 2649520L, 59363566L)
) 

  par = PAR_REGIONS_Y[ref == reference]
  par1 = par[region == "PAR1Y"]
  par2 = par[region == "PAR2Y"]
  
  cov[Chromosome == y_label &
        !(Start_Position >= par1$start & End_Position <= par1$end) &
        !(Start_Position >= par2$start & End_Position <= par2$end)]
}



log2_ratio_to_int=function(row, purity=NA,...)
{ 
  #Purity and Ploidy code is adapted from https://github.com/etal/cnvkit/blob/5c4fd715dc9c103566d32a543de7ef8a3f0d3011/cnvlib/call.py
  #converts log 2 ratio to integer CN
  # Math::
  # 
  #     log2_ratio = log2(ncopies / ploidy)
  #     2^log2_ratio = ncopies / ploidy
  #     ncopies = ploidy * 2^log2_ratio
  # 
  # With rescaling for purity:
  # 
  #     let v = log2 ratio value, p = tumor purity,
  #         r = reference ploidy, x = expected ploidy,
  #         n = tumor ploidy ("ncopies" above);
  # 
  #     v = log_2(p*n/r + (1-p)*x/r)
  #     2^v = p*n/r + (1-p)*x/r
  #     n*p/r = 2^v - (1-p)*x/r
  #     n = (r*2^v - x*(1-p)) / p
  # if purity adjustment is skipped (p=1, ex. germline)
  #     n = r&2^v
  
  logR      = as.numeric(row["logR"])
  reference = as.numeric(row["reference"])
  expected  = as.numeric(row["expected"])
  
  if (!is.na(purity)) {
    cn = (reference * 2^logR - expected * (1 - purity)) / purity
  } else {
    cn = reference * 2^logR
  }
  return(as.numeric(cn))
}



infer_sample_sex = function(cov, y_label = "chrY", autosome_labels = paste0("chr", 1:22)) {
  cov = cov[, c(1,2,3,4), with = FALSE]
  data.table::setnames(cov, c("chr", "start", "end", "doc"))

  y_rows    = cov[chr == y_label]
  auto_rows = cov[chr %in% autosome_labels]

  y_depth    = weighted_median(y_rows$doc, y_rows$end - y_rows$start)
  auto_depth = weighted_median(auto_rows$doc, auto_rows$end - auto_rows$start)

  if (is.na(y_depth) || is.na(auto_depth) || auto_depth == 0) {
    message("Could not infer sex from coverage — defaulting to sex = 'F'")
    return("F")
  }

  y_ratio = y_depth / auto_depth
  sex = if (y_ratio < 0.15) "F" else "M"
  message("Inferred sex: ", sex, " (chrY/autosome depth ratio = ", round(y_ratio, 3), ")")
  sex
}

get_ref_exp_copies = function(cn_df, ploidy, is_haploid_x_ref, sex, labelx, labely, ref){
  #determines the number of copies expected and in the reference
  #reference count is the chr ploidy in the reference
  #expected count is the chr neutral ploidy of the given sample based on the specified sex
  #ex.female sample but male reference 
  # chr X expected = 2 reference = 1
  # chr Y reference is always 1 to avoid division by 0
  # default reference is XXY 
  #--------
  
  half_ploidy = ploidy %/% 2L
  
  
  #cn_df[, `:=`(reference = as.integer(ploidy),
  #             expected  = as.integer(ploidy))]
  cn_df[, `:=`(reference = 2L,    
               expected  = 2L)]   
  

  chrx_filtered = filter_chr_x(cn_df, labelx, ref)
  cn_df[, is_chrx := FALSE]
  cn_df[chrx_filtered,
        on = .(Chromosome, Start_Position, End_Position),
        is_chrx := TRUE
  ]
  
  chry_filtered = filter_chr_y(cn_df, labely, ref)
  cn_df[, is_chry := FALSE]
  cn_df[chry_filtered,
        on = .(Chromosome, Start_Position, End_Position),
        is_chry := TRUE
  ]

  if (sex == "M") {
    cn_df[is_chrx == TRUE, expected := half_ploidy]
    cn_df[is_chry == TRUE, expected := half_ploidy]
  } else {
    cn_df[is_chry == TRUE, expected := 0L]
  }
  

  cn_df[is_chrx == TRUE, reference := 2L]
  if (is_haploid_x_ref) {
    cn_df[is_chrx == TRUE, reference := 1L]
  }
  cn_df[is_chry == TRUE, reference := 1L]
  
  cn_df[, `:=`(is_chrx = NULL, is_chry = NULL)]
  cn_df
  
}


weighted_median = function(x, w) {
  valid = is.finite(x) & is.finite(w) & w > 0
  x = x[valid]; w = as.numeric(w[valid])
  if (length(x) == 0) return(NA_real_)
  ord = order(x)
  x = x[ord]; w = w[ord]
  cum_w = cumsum(w)
  half  = sum(w) / 2
  x[which(cum_w >= half)[1]]
}

rescale_cn = function(cn,purity, ploidy,reference,  sex){
  
  ref_config = list(
    GRCh38 = list(labelx = "chrX", labely = "chrY", ref = "grch38",
                   is_haploid_x_ref = FALSE)
  )
  
  cfg = ref_config[[reference]]
  if (is.null(cfg)) stop("Unknown reference: ", reference)
  
  message("Rescaling sample with purity ", purity, " and ploidy ", ploidy)
  
  cn_df = get_absolute_cn(cn, purity, ploidy,
                          cfg$is_haploid_x_ref, sex, cfg$labelx, cfg$labely, cfg$ref
  )
  
  cn_df[, `:=`(
    cn     = as.integer(cn),
    cn_int = as.integer(cn)
  )][is.na(cn), `:=`(
    cn     = -1L,
    cn_int = -1L
  )]
  
  cn_df
}


get_absolute_cn = function(cn_df, purity, ploidy, is_haploid_x_ref, sex, labelx, labely, ref){
  copies_df = get_ref_exp_copies(cn_df, ploidy, is_haploid_x_ref, sex, labelx, labely, ref)
  
  cn_abs = apply(
    as.data.frame(copies_df), 
    1, 
    log2_ratio_to_int, 
    purity = purity
  )
  
  cn_df[, cn := cn_abs]
  cn_df
  
}

get_hlfa = function(cn_df, min_cn, min_width, max_width, chromosomes){

  hlfa_df = data.table::copy(cn_df)
  hlfa_df[, width := End_Position - Start_Position + 1L]  
  hlfa_df = hlfa_df[(cn >= min_cn & width >= min_width & width <= max_width) |    (cn == -1 & width >= 100 & width < min_width)]
  hlfa_df[, Chromosome := data.table::fcase(
    Chromosome %in% c("chrX", "X", "23", 23L), "chrX",
    Chromosome %in% c("chrY", "Y", "24", 24L), "chrY",
    !startsWith(Chromosome, "chr"),              paste0("chr", Chromosome),
    rep(TRUE, .N),                               Chromosome   # already has chr prefix
  )]
  if (!identical(chromosomes, NA)) {
    hlfa_df = hlfa_df[Chromosome %in% chromosomes]
  }
  
  hlfa_df
}


build_breakpoint_hashmaps = function(breakpoints_df){
  n = as.integer(nrow(breakpoints_df))
  
  return(list(
    bp_start = with(breakpoints_df,
                    setNames(seq_len(n), paste0(seqnames, "|", start))),
    bp_end   = with(breakpoints_df,
                    setNames(seq_len(n), paste0(seqnames, "|", end))),
    bp_partner = with(breakpoints_df,
                      setNames(seq_len(n), partner)),
    breakpoints_df = breakpoints_df
  ))
}

get_breakend = function(chr, position, strand, idx) {
  be_type  = 0
  qual     = -1
  breakend = NA
  partner  = NA
  
  if (strand == "+") {
    key = paste0(chr, "|", position + 1)
    ii  = as.integer(unname(idx$bp_start[key])) 
    if (!is.na(ii)) {
      breakend = idx$breakpoints_df[ii, ]
      be_type  = 2
      qual     = breakend[["QUAL"]]
      jj       = as.integer(unname(idx$bp_partner[breakend[["sourceId"]]]))
      partner  = if (!is.na(jj)) idx$breakpoints_df[jj, ] else NA
    }
  } else {
    key = paste0(chr, "|", position)
    ii  = as.integer(unname(idx$bp_end[key]))    
    if (!is.na(ii)) {
      breakend = idx$breakpoints_df[ii, ]
      be_type  = -2
      qual     = breakend[["QUAL"]]
      jj       = as.integer(unname(idx$bp_partner[breakend[["sourceId"]]]))
      partner  = if (!is.na(jj)) idx$breakpoints_df[jj, ] else NA
    }
  }
  
  list(breakend_type = be_type,
       breakend      = breakend,
       partner       = partner,
       qual          = qual,
       strand        = strand)
}


get_node_pos = function(chr, pos, hlfa_df, full_segments_df, tolerance = 200000) {
  if (length(chr) != 1L || length(pos) != 1L) {
    message("get_node_pos called with non-scalar input: chr=[", paste(chr, collapse=","),
            "] (length ", length(chr), "), pos=[", paste(pos, collapse=","),
            "] (length ", length(pos), ")")
    chr = chr[1]
    pos = pos[1]
  }
  exact = hlfa_df[Chromosome == chr &
                    Start_Position <= pos &
                    End_Position   >= pos]
  if (nrow(exact) > 0) {
    if (abs(exact$Start_Position[1] - pos) <= abs(exact$End_Position[1] - pos)) {
      return(exact$Start_Position[1])
    } else {
      return(exact$End_Position[1])
    }
  }
  
  nearby_debug = tryCatch({
  message("hlfa_df dims: ", nrow(hlfa_df), " x ", ncol(hlfa_df))
  message("class(hlfa_df$Chromosome): ", class(hlfa_df$Chromosome), " length: ", length(hlfa_df$Chromosome))
  message("class(hlfa_df$cn_int): ", class(hlfa_df$cn_int), " length: ", length(hlfa_df$cn_int))
  message("class(hlfa_df$Start_Position): ", class(hlfa_df$Start_Position))
  message("tolerance: ", tolerance, " length: ", length(tolerance))
  NULL
}, error = function(e) message("debug check itself failed: ", conditionMessage(e)))

  nearby = hlfa_df[Chromosome == chr & cn_int >= 0 &
                     (abs(Start_Position - pos) <= tolerance |
                        abs(End_Position   - pos) <= tolerance)]
  if (nrow(nearby) == 0) return(NULL)
  
  starts = abs(nearby$Start_Position - pos)
  ends   = abs(nearby$End_Position   - pos)
  dists  = pmin(starts, ends)
  i      = which.min(dists)
  best   = nearby[i]
  
  true_segment = full_segments_df[Chromosome == chr &
                                    Start_Position <= pos &
                                    End_Position   >= pos]
  if (nrow(true_segment) == 0) return(NULL)
  
  adjacent = true_segment$Start_Position[1] == best$End_Position[1] ||
             true_segment$End_Position[1]   == best$Start_Position[1]
  if (!adjacent) {
    message("Rejected snap: ", chr, ":", pos, " -> nearest candidate ",
            best$Start_Position, "-", best$End_Position,
            " (", round(dists[i]/1000, 1), "kb, not contiguous)")
    return(NULL)
  }
  
  if (starts[i] <= ends[i]) {
    snapped = best$Start_Position
  } else {
    snapped = best$End_Position
  }
  
  return(snapped)
}


create_graph_components = function(hlfa_df, breakpoints_df, full_segments_df, min_cn, min_width) {
  message("full_segments_df columns: ", paste(names(full_segments_df), collapse=", "))
  message("full_segments_df dims: ", nrow(full_segments_df), " x ", ncol(full_segments_df))


  get_node_pos = function(chr, pos, hlfa_df, full_segments_df, tolerance = 200000) {
    
    message("hlfa_df dims: ", nrow(hlfa_df), " x ", ncol(hlfa_df))
    message("class(hlfa_df$Chromosome): ", class(hlfa_df$Chromosome), " length: ", length(hlfa_df$Chromosome))
    message("class(hlfa_df$cn_int): ", class(hlfa_df$cn_int), " length: ", length(hlfa_df$cn_int))
    message("class(hlfa_df$Start_Position): ", class(hlfa_df$Start_Position))
    message("tolerance: ", tolerance, " length: ", length(tolerance))
    message("chr=[", paste(chr, collapse=","), "] length=", length(chr))
    message("pos=[", paste(pos, collapse=","), "] length=", length(pos))
    
    
    exact = hlfa_df[Chromosome == chr &
                      Start_Position <= pos &
                      End_Position   >= pos]
    if (nrow(exact) > 0) {
      if (abs(exact$Start_Position[1] - pos) <= abs(exact$End_Position[1] - pos)) {
        return(exact$Start_Position[1])
      } else {
        return(exact$End_Position[1])
      }
    }
    
    nearby = hlfa_df[Chromosome == chr & cn_int >= 0 &
                       (abs(Start_Position - pos) <= tolerance |
                          abs(End_Position   - pos) <= tolerance)]
    if (nrow(nearby) == 0) return(NULL)
    
    starts = abs(nearby$Start_Position - pos)
    ends   = abs(nearby$End_Position   - pos)
    dists  = pmin(starts, ends)
    i      = which.min(dists)
    best   = nearby[i]
    
    if (starts[i] <= ends[i]) {
      snapped = best$Start_Position
    } else {
      snapped = best$End_Position
    }
    
    return(snapped)
  }
  
  hm = build_breakpoint_hashmaps(breakpoints_df)
  
  idx = list(
    bp_start       = hm$bp_start,
    bp_end         = hm$bp_end,
    bp_partner     = hm$bp_partner,
    breakpoints_df = breakpoints_df
  )
  
  bp_start = hm$bp_start
  bp_end   = hm$bp_end
  
  n = nrow(hlfa_df)
  
  
  nodes_list = vector("list", n * 2L)
  edges_list = vector("list", n * 4L)  
  edge_idx   = 0L
  
  
  
  for (i in seq_len(n)) {
    
    seg       = hlfa_df[i]
    chr       = seg$Chromosome
    seg_start = seg$Start_Position
    seg_end   = seg$End_Position
    cn_val    = seg$cn_int
    width     = seg$width
    seg_type  = if (cn_val >= min_cn && width >= min_width) 1L else -1L
    
    if ((chr == "chr14" && seg_start == 33545797) || (chr == "chr11" && seg_start == 78469228)) {
    message("FOUND TARGET SEGMENT: chr=", chr, " start=", seg_start, " end=", seg_end,
            " cn_int=", cn_val, " width=", width, " seg_type=", seg_type,
            " min_cn=", min_cn, " min_width=", min_width)
    }
    
    sv_start_plus  = get_breakend(chr, seg_start, "+", idx)
    sv_start_minus = get_breakend(chr, seg_start, "-", idx)
    sv_start = if (sv_start_plus$breakend_type != 0) sv_start_plus else sv_start_minus
    
    
    node_start = data.table::data.table(
      ID           = paste0(chr, ":", seg_start),
      chr          = chr,
      pos          = seg_start,
      parent       = as.integer(i - 1L),
      width        = as.integer(width),
      start_flag   = TRUE,
      text         = paste0(chr, ":", seg_start),
      segment_type = seg_type,
      be_type      = as.numeric(sv_start$breakend_type),
      be_qual      = as.numeric(sv_start$qual),
      cn           = as.numeric(cn_val)
    )
    
    node_end = data.table::data.table(
      ID           = paste0(chr, ":", seg_end),
      chr          = chr,
      pos          = seg_end,
      parent       = as.integer(i - 1L),
      width        = as.integer(width),
      start_flag   = FALSE,
      text         = paste0(chr, ":", seg_end),
      segment_type = seg_type,
      be_type      = as.numeric(sv_start$breakend_type),
      be_qual      = as.numeric(sv_start$qual),
      cn           = as.numeric(cn_val)
    )
    
    nodes_list[[2L * i - 1L]] = node_start
    nodes_list[[2L * i]]      = node_end
    
    be_plus  = get_breakend(chr, seg_end, "+", idx)
    be_minus = get_breakend(chr, seg_end, "-", idx)
    be = if (be_plus$breakend_type != 0) be_plus else be_minus
    
    
    
    if (be$breakend_type != 0 && !identical(be$partner, NA)) {
      partner_row = be$partner
      
      
      if (!identical(partner_row, NA) && nrow(partner_row) > 0) {
        p_chr = as.character(partner_row$seqnames[1])
        p_pos = as.integer(partner_row$start[1])
        
        p_pos_snapped = get_node_pos(p_chr, p_pos, hlfa_df, full_segments_df)
        
        if (!is.null(p_pos_snapped)) {
          edge_idx = edge_idx + 1L
          edges_list[[edge_idx]] = data.table::data.table(
            ID           = paste0(chr, ":", seg_end, "_", p_chr, ":", p_pos_snapped),
            source       = paste0(chr, ":", seg_end),
            target       = paste0(p_chr, ":", p_pos_snapped),
            source_chr   = chr,
            source_start = seg_start,
            source_end   = seg_end,
            target_chr   = p_chr,
            target_start = p_pos_snapped,
            target_end   = p_pos_snapped,
            qual         = as.numeric(be$qual),
            VAF          = as.numeric(partner_row$VAF_num[1]),
            svLen        = as.numeric(partner_row$svLen[1]),
            svtype       = as.character(partner_row$svtype[1]),
            strand       = as.character(be$strand),
            filter       = NA_character_
          )
        }
      }
    }
    
    if (seg_type == 1L) {
      seg_end_num = as.numeric(seg_end)
      
      nearby_bp = breakpoints_df[
        as.character(seqnames) == chr &
          as.numeric(start) > seg_end_num &
          as.numeric(start) - seg_end_num <= 200000
      ]
      
      if (nrow(nearby_bp) > 0) {
        nearby_bp[, partner_chr := sapply(partner, function(pid) {
          pr = breakpoints_df[sourceId == pid]
          if (nrow(pr) > 0) as.character(pr$seqnames[1]) else NA_character_
        })]
        
        nearby_bp = nearby_bp[!is.na(partner_chr) & partner_chr != chr]
      }
      
      if (nrow(nearby_bp) > 0) {
        nearby_bp = nearby_bp[which.min(abs(as.numeric(nearby_bp$start) - seg_end_num))]
        partner_id  = nearby_bp$partner[1]
        partner_row = breakpoints_df[sourceId == partner_id]
        
        
        
        if (nrow(partner_row) > 0) {
          p_chr = as.character(partner_row$seqnames[1])
          p_pos = as.integer(partner_row$start[1])
          p_pos_snapped = get_node_pos(p_chr, p_pos, hlfa_df, full_segments_df)
          
          if (!is.null(p_pos_snapped)) {
            edge_idx = edge_idx + 1L
            edges_list[[edge_idx]] = data.table::data.table(
              ID           = paste0(chr, ":", seg_end, "_nearby_", p_chr, ":", p_pos_snapped),
              source       = paste0(chr, ":", seg_end),
              target       = paste0(p_chr, ":", p_pos_snapped),
              source_chr   = chr,
              source_start = seg_start,
              source_end   = seg_end,
              target_chr   = p_chr,
              target_start = p_pos_snapped,
              target_end   = p_pos_snapped,
              qual         = as.numeric(nearby_bp$QUAL[1]),
              VAF          = as.numeric(partner_row$VAF_num[1]),
              svLen        = as.numeric(partner_row$svLen[1]),
              svtype       = as.character(partner_row$svtype[1]),
              strand       = as.character(nearby_bp$strand[1]),
              filter       = NA_character_
            )
          }
        }
      }
    }
  }
    
  nodes_dt = data.table::rbindlist(nodes_list[!sapply(nodes_list, is.null)])
  edges_dt = data.table::rbindlist(edges_list[seq_len(edge_idx)])
  edges_dt = unique(edges_dt, by = c("source", "target"))
  
  return(list("nodes" = nodes_dt, "edges" = edges_dt))
}




apply_filters = function(nodes_df, edges_df, breakpoints_df, filters, min_cn, min_width){
  keys = names(filters)
  
  
  low_qual_bp = data.frame()
  #get filtered breakends
  high_qual_bp = data.table::as.data.table(breakpoints_df)
  for(i in 1:length(keys)){
    if(keys[i] == 'MIN_QUAL_BREAKPOINT'){
      MIN_QUAL_BREAKPOINT  =  filters[['MIN_QUAL_BREAKPOINT']]
      #to get breakpoints that pass the filter but have the same location as breakpoints that don't pass the filter
      high_qual_bp = dplyr::filter(breakpoints_df,QUAL >= MIN_QUAL_BREAKPOINT)
      #get breakpoints 
      low_qual_bp = dplyr::filter(breakpoints_df,QUAL < MIN_QUAL_BREAKPOINT)
    }
    
    if(keys[i] == 'MIN_QUAL_VAF'){
      MIN_QUAL_VAF  =  filters[['MIN_QUAL_VAF']]
      #to get breakpoints that pass the filter but have the same location as breakpoints that don't pass the filter
      #get breakpoints
      low_VAF_bp = dplyr::filter(breakpoints_df,VAF_num < MIN_QUAL_VAF)
      low_qual_bp = rbind(low_qual_bp,low_VAF_bp)
      high_qual_bp = dplyr::filter(high_qual_bp,VAF_num >= MIN_QUAL_VAF & !(sourceId %in% low_VAF_bp$partner))
    }
    
    if(keys[i] == 'MIN_SV_LEN'){
      MIN_SV_LEN =  filters[['MIN_SV_LEN']]
      low_svlen_bp = dplyr::filter(breakpoints_df,svLen < MIN_SV_LEN & (svtype == "INS" | svtype == "DEL"))
      low_qual_bp = rbind(low_qual_bp,low_svlen_bp)
      high_qual_bp = dplyr::filter(high_qual_bp, !(svLen < MIN_SV_LEN & (svtype == "INS" | svtype == "DEL"))  & !(sourceId %in%  low_svlen_bp$partner))
    }
    
    
  }
  
  
  #get partners of low quality breakpoints
  
  
  low_qual_bp = unique(
    data.table::as.data.table(low_qual_bp),
    by = c("seqnames", "start", "end", "strand", "sourceId")
  )
  
  merged_nodes_df = nodes_df
  
  #join segments for filtered breakpoints
  
  if(nrow(low_qual_bp) > 0){
    for(i in 1:nrow(low_qual_bp)){
      
      if(nrow(dplyr::filter(high_qual_bp,seqnames==low_qual_bp[[i,"seqnames"]] & start==low_qual_bp[[i,"start"]] ))==0){ # &&
        source = NA
        target = NA
        source = which(merged_nodes_df$pos+1==low_qual_bp[[i,"start"]] & merged_nodes_df$chr==low_qual_bp[[i,"seqnames"]],arr.ind = TRUE )[1]
        target = which(merged_nodes_df$pos==low_qual_bp[[i,"end"]] & merged_nodes_df$chr==low_qual_bp[[i,"seqnames"]],arr.ind = TRUE )[1]
        
        if(!is.na(source)&!is.na(target)) #if there is a segment between the breakpoint then there will be a target
        { 
          
          #merge by removing right end of source node and left end of target node
          nodes_to_remove =  c();
          nodes_to_remove = c(nodes_to_remove,source)
          nodes_to_remove = c(nodes_to_remove,target)
          merged_nodes_df[[target+1,"parent"]] = merged_nodes_df[[source-1,"parent"]]
          #recompute copy number
          if(merged_nodes_df[[source,"cn"]]!=-1 & merged_nodes_df[[target,"cn"]]!=-1){
            new_cn= as.integer((merged_nodes_df[[source,"cn"]] * merged_nodes_df[[source,"width"]] + merged_nodes_df[[target,"cn"]] * merged_nodes_df[[target,"width"]]) / (merged_nodes_df[[source,"width"]] + merged_nodes_df[[target,"width"]]))
          }else if(merged_nodes_df[[source,"cn"]]!=-1) {
            new_cn= as.integer(merged_nodes_df[[source,"cn"]])
          }else if(merged_nodes_df[[target,"cn"]]!=-1) {
            new_cn= as.integer(merged_nodes_df[[target,"cn"]])
          }else{
            new_cn = -1
          }
          merged_nodes_df[[source-1,"cn"]] = new_cn
          merged_nodes_df[[target+1,"cn"]] = new_cn
          #recompute width
          new_width = merged_nodes_df[[target+1,"pos"]] -  merged_nodes_df[[source-1,"pos"]] + 1
          merged_nodes_df[[source-1,"width"]] = new_width
          merged_nodes_df[[target+1,"width"]] = new_width
          
          if(new_cn < min_cn | new_width < min_width){
            merged_nodes_df[[source-1,"segment_type"]] = -1
            merged_nodes_df[[target+1,"segment_type"]] = -1
          }
          else{
            merged_nodes_df[[source-1,"segment_type"]] = 1
            merged_nodes_df[[target+1,"segment_type"]] = 1
          }
          
          #update text and ID
          ID = paste( merged_nodes_df[[source-1,"chr"]],  merged_nodes_df[[source-1,"pos"]],sep=":")
          merged_nodes_df[[source-1,"text"]]= ID
          merged_nodes_df[[source+1,"text"]]= ID
          merged_nodes_df[[source-1,"ID"]]= ID
          merged_nodes_df[[source+1,"ID"]]= ID
          
          merged_nodes_df =  merged_nodes_df[-nodes_to_remove,]
        }
        else if(!is.na(source)) #if there is no segment between the breakpoints (ends of the breakends of the breakpoints are 1bp apart)
        {
          #get the partner
          partner = dplyr::filter(low_qual_bp,low_qual_bp[[i,"partner"]]==sourceId)[1]
        }
        
        
      }
    }
  }
  
  edges_filtered = edges_df
  for(i in 1:length(keys)){
    if(keys[i] == 'MIN_QUAL_BREAKPOINT'){
      MIN_QUAL_BREAKPOINT  =  filters[['MIN_QUAL_BREAKPOINT']]
      edges_filtered =  dplyr::filter(edges_filtered, qual>= MIN_QUAL_BREAKPOINT)
      
    }
    else if(keys[i] == 'MIN_QUAL_VAF'){
      MIN_QUAL_VAF  =  filters[['MIN_QUAL_VAF']]
      
      edges_filtered =  dplyr::filter(edges_filtered, VAF>= MIN_QUAL_VAF)
      
      
    }
    else if(keys[i] == 'MIN_SV_LEN'){
      MIN_SV_LEN  =  filters[['MIN_SV_LEN']]
      edges_filtered = dplyr::filter(edges_filtered, 
                                     is.na(svLen) | !(svLen < MIN_SV_LEN & (svtype == "INS" | svtype == "DEL")))
    }
  }  
  
  
  #remove edges with no matching nodes because they were removed ex. CN below minimum CN
  edges_filtered = dplyr::filter(edges_filtered,source %in% merged_nodes_df$ID & target %in% merged_nodes_df$ID )
  return(
    list("filtered_edges" = edges_filtered,
         "merged_nodes" = merged_nodes_df)
  )

  
}



add_adjacency_edges = function(nodes_df, edges_df, min_width) {
  
  start_nodes = nodes_df[start_flag == TRUE]
  end_nodes   = nodes_df[start_flag == FALSE]
  
  #sort both by chromosome then position
  data.table::setorder(start_nodes, chr, pos)
  data.table::setorder(end_nodes,   chr, pos)
  
  nodes_ct = nrow(start_nodes)
  if (nodes_ct < 2L) return(edges_df)
  
  #get end position for each segment via parent lookup
  end_pos = setNames(end_nodes$pos, as.character(end_nodes$parent))
  
  chr_vec       = start_nodes$chr
  start_pos_vec = start_nodes$pos
  end_pos_vec   = end_pos[as.character(start_nodes$parent)]
  
  #consecutive segments on same chromosome
  same_chr = chr_vec[-nodes_ct] == chr_vec[-1L]
  idx      = which(same_chr)
  
  
  if (length(idx) == 0L) return(edges_df)
  
  #end of segment N to start of segment N+1
  src   = paste0(chr_vec[idx],      ":", end_pos_vec[idx])
  tgt   = paste0(chr_vec[idx + 1L], ":", start_pos_vec[idx + 1L])
  
  #remove self-loops 
  valid = src != tgt
  idx   = idx[valid]; src = src[valid]; tgt = tgt[valid]
  
  if (length(idx) == 0L) return(edges_df)
  
  adj_edges = data.table::data.table(
    ID           = seq(nrow(edges_df), nrow(edges_df) + length(idx) - 1L),
    source       = src,
    target       = tgt,
    source_chr   = chr_vec[idx],
    source_start = end_pos_vec[idx],
    source_end   = end_pos_vec[idx],
    target_chr   = chr_vec[idx + 1L],
    target_start = start_pos_vec[idx + 1L],
    target_end   = start_pos_vec[idx + 1L],
    filter       = FALSE,
    qual         = 100,
    svtype       = "ADJ",
    svLen        = min_width + 10L,
    VAF          = 1.0
  )
  
  data.table::rbindlist(list(edges_df, adj_edges), use.names = TRUE, fill = TRUE)
}


get_cycles_from_df = function(edges_df, nodes_df) {
  
  if (is.null(edges_df) || nrow(edges_df) == 0) return(list())

  endpoint_to_segment = nodes_df[nodes_df$segment_type != -1L, c('ID', 'parent')]
  endpoint_to_segment = endpoint_to_segment[!duplicated(endpoint_to_segment$ID), ]
  
  valid_edges = edges_df[edges_df$source %in% endpoint_to_segment$ID &
                           edges_df$target %in% endpoint_to_segment$ID, ]
  
  if (nrow(valid_edges) == 0) return(list())
  
  edge_pairs = valid_edges |>
    dplyr::inner_join(endpoint_to_segment, by = c('source' = 'ID')) |>
    dplyr::rename(from = parent) |>
    dplyr::inner_join(endpoint_to_segment, by = c('target' = 'ID')) |>
    dplyr::rename(to = parent) |>
    dplyr::mutate(from = as.character(from), to = as.character(to)) |>
    dplyr::select(from, to) |>
    dplyr::distinct()
  
  self_loops = edge_pairs[edge_pairs$from == edge_pairs$to, ]
  edge_pairs  = edge_pairs[edge_pairs$from != edge_pairs$to, ]

  #self_loop_cycles = lapply(seq_len(nrow(self_loops)), function(i) {
  #  list(as.numeric(self_loops$from[i]))
  #})
  

  g = igraph::graph_from_data_frame(edge_pairs, directed = TRUE)

  #SCC
  scc = igraph::components(g, mode = "strong")
  scc_sizes = table(scc$membership)
  
  scc_multi_members = unlist(lapply(names(scc_sizes[scc_sizes >= 2L]), function(comp_id) {
  names(scc$membership[scc$membership == as.integer(comp_id)])
  }))

  self_loops = self_loops[!(self_loops$from %in% scc_multi_members), ]

  self_loop_cycles = lapply(seq_len(nrow(self_loops)), function(i) {
    list(as.numeric(self_loops$from[i]))
  })

  #if (nrow(edge_pairs) == 0) return(lapply(seq_len(nrow(self_loops)), function(i) list(as.numeric(self_loops$from[i]))))
  

  scc_cycles = lapply(names(scc_sizes[scc_sizes >= 2L]), function(comp_id) {
    members = names(scc$membership[scc$membership == as.integer(comp_id)])
    as.list(as.numeric(members))
  })
  

  c(self_loop_cycles, scc_cycles)
}

find_closed_walk_cn = function(scc_members, edges_filtered, merged_nodes_df) {
  
  segs = data.table::copy(
    data.table::as.data.table(merged_nodes_df)[
      parent %in% as.numeric(scc_members) & start_flag == TRUE,
      .(parent, chr, pos, cn)
    ]
  )
  
  if (nrow(segs) == 0 || all(segs$cn <= 0)) {
    fallback_cn = min(merged_nodes_df[parent %in% as.numeric(scc_members) & cn > 0, cn], na.rm=TRUE)
    return(list(list(cn = fallback_cn, segments = as.numeric(scc_members))))
  }
  
  node_ids  = merged_nodes_df[parent %in% as.numeric(scc_members), ID]
  scc_edges = data.table::as.data.table(
    edges_filtered[edges_filtered$source %in% node_ids &
                     edges_filtered$target %in% node_ids, ]
  )
  
  if (nrow(scc_edges) == 0) {
    fallback_cn = min(segs[cn > 0, cn], na.rm=TRUE)
    return(list(list(cn = fallback_cn, segments = as.numeric(scc_members))))
  }
  
  ep_to_seg = setNames(merged_nodes_df$parent, merged_nodes_df$ID)
  scc_edges[, from_parent := ep_to_seg[source]]
  scc_edges[, to_parent   := ep_to_seg[target]]
  scc_edges = scc_edges[from_parent != to_parent]
  
  if (nrow(scc_edges) == 0) {
    fallback_cn = min(segs[cn > 0, cn], na.rm=TRUE)
    return(list(list(cn = fallback_cn, segments = as.numeric(scc_members))))
  }
  
  find_walk = function(start_seg, segs_cn) {

    max_stack_iter = 10000L  
    stack_iter     = 0L
    stack = list(list(current = start_seg, walk = integer(0), tried = integer(0)))
    
    while (length(stack) > 0) {
      frame   = stack[[length(stack)]]
      current = frame$current
      walk    = c(frame$walk, current)
      tried   = frame$tried
      stack_iter = stack_iter + 1L
      if (stack_iter > max_stack_iter) return(integer(0))
      if (length(walk) > 1 && current == start_seg) return(walk)
      
      outgoing     = scc_edges[from_parent == current]
      next_parents = unique(outgoing$to_parent)
      next_parents = next_parents[next_parents %in% segs_cn[cn > 0, parent] & !is.na(next_parents)]
      
      if (start_seg %in% next_parents && length(walk) > 1) return(c(walk, start_seg))
      
      candidates = next_parents[!next_parents %in% walk & !next_parents %in% tried]
      
      if (length(candidates) == 0) {
        stack = stack[-length(stack)]
        if (length(stack) > 0) {
          parent_frame = stack[[length(stack)]]
          parent_frame$tried = c(parent_frame$tried, current)
          stack[[length(stack)]] = parent_frame
        }
        next
      }
      
      cn_vals  = segs_cn[parent %in% candidates, .(parent, cn)]
      next_seg = cn_vals[which.max(cn), parent]
      
      frame$tried = c(tried, next_seg)
      stack[[length(stack)]] = frame
      stack = c(stack, list(list(current = next_seg, walk = walk, tried = integer(0))))
    }
    return(integer(0))
  }
  
  walks_found = list()  
  max_iter    = 500L
  iter        = 0L
  
  repeat {
    iter = iter + 1L
    if (iter > max_iter) break
    if (all(segs$cn <= 0)) break
    
    start_seg = segs[cn > 0][which.max(cn), parent]
    walk      = find_walk(start_seg, segs)
    
    if (length(walk) < 2) break
    
    walk_unique = unique(walk)
    bottleneck  = min(segs[parent %in% walk_unique, cn], na.rm=TRUE)
    
    if (bottleneck <= 0) break
    
    segs[parent %in% walk_unique, cn := cn - bottleneck]
    
    walks_found[[length(walks_found) + 1L]] = list(  
      cn       = as.integer(round(bottleneck)),
      segments = walk_unique
    )
    
    message("Walk ", length(walks_found), ": ", paste(walk_unique, collapse="→"),
            " bottleneck=", bottleneck)
  }
  
  if (length(walks_found) == 0) {
    fallback_cn = min(merged_nodes_df[parent %in% as.numeric(scc_members) & cn > 0, cn], na.rm=TRUE)
    return(list(list(cn = fallback_cn, segments = as.numeric(scc_members))))
  }
  
  return(walks_found)   
}

decompose_cycles = function(cycles_df, merged_nodes, threshold = 4.5){
  
  if (is.null(cycles_df) || nrow(cycles_df) == 0L) {
    message("No cycles to decompose")
    return(data.table::data.table())
  }
  
  cycles_dt  = data.table::as.data.table(cycles_df)
  segment_cn = unique(cycles_dt[, .(segment_id, cn)], by = "segment_id")
  data.table::setkey(segment_cn, segment_id)
  
  seg_chr = unique(merged_nodes[, .(segment_id = parent, chr = chr)])
  data.table::setkey(seg_chr, segment_id)
  
  cycle_nchroms = seg_chr[cycles_dt[, .(segment_id, cycle_id)], on = "segment_id"]
  cycle_nchroms = cycle_nchroms[, .(
    n_chroms = data.table::uniqueN(chr[!is.na(chr)])
  ), by = cycle_id]
  
  
  message("cycle_nchroms top 5:")
  print(head(cycle_nchroms[order(-n_chroms)]))
  
  remaining  = data.table::copy(cycles_dt)
  results    = vector("list", nrow(cycles_dt))
  result_idx = 0L
  
  message("Starting cycle decomposition with ", nrow(cycles_dt), " rows")
  message("cycle_nchroms sample:")
  print(head(cycle_nchroms[order(-n_chroms)]))
  
  selected_cycles = character(0)
  
  repeat {
    
    bottlenecks = segment_cn[remaining, on = "segment_id"
    ][cn > 0, .(bottleneck = min(cn)), by = cycle_id]
    
    best_bn = max(bottlenecks$bottleneck)
    message("Best bottleneck: ", best_bn)
    
    if (best_bn <= 0) break
    
    #among cycles with highest bottleneck, prefer most chromosomes
    high_bn = bottlenecks[bottleneck == best_bn & !cycle_id %in% selected_cycles]
    
    
    if (nrow(high_bn) == 0) break
    
    high_bn = cycle_nchroms[high_bn, on = "cycle_id"]
    high_bn = high_bn[order(-n_chroms)][1L]
    
    selected_cycles = c(selected_cycles, high_bn$cycle_id)
    
    message("Selected cycle: ", high_bn$cycle_id,
            " with bottleneck: ", high_bn$bottleneck,
            " n_chroms: ", high_bn$n_chroms)
    
    segments = remaining[cycle_id == high_bn$cycle_id, segment_id]
    message("Segments in cycle: ", paste(segments, collapse = ", "))
    
    segment_cn[segment_id %in% segments, cn := cn - high_bn$bottleneck]
    
    result_idx = result_idx + 1L
    results[[result_idx]] = data.table::data.table(
      cycle_id    = high_bn$cycle_id,
      assigned_cn = high_bn$bottleneck,
      ecdna       = high_bn$bottleneck > threshold
    )
  }
  
  if (result_idx == 0L) return(data.table::data.table())
  dc = data.table::rbindlist(results[seq_len(result_idx)])

  MAX_ECDNA_TOTAL_LENGTH = 50000000  

  merged_dt = data.table::as.data.table(merged_nodes)
  seg_widths = unique(merged_dt[, .(segment_id = parent, width)], by = "segment_id")

  cycle_total_length = cycles_dt[seg_widths, on = "segment_id", nomatch = NULL][
    , .(total_cycle_length = sum(width, na.rm = TRUE)), by = cycle_id]

  dc = cycle_total_length[dc, on = "cycle_id"]
  dc[is.na(total_cycle_length), total_cycle_length := 0L]
  dc[total_cycle_length > MAX_ECDNA_TOTAL_LENGTH, ecdna := FALSE]

  dc
}


detect_cycles = function(outdir, all_cn, breakpoints_df, min_width, chromosomes, filter_grid,mode="SR",ploidy=2.0) {


  
  breakpoints_df = dplyr::mutate(breakpoints_df, seqnames_int = as.numeric(gsub("chr", "", seqnames)))
  breakpoints_df = dplyr::mutate(breakpoints_df, seqnames_int = if_else(breakpoints_df$seqnames == "chrX", 22, breakpoints_df$seqnames_int))
  breakpoints_df = dplyr::mutate(breakpoints_df, seqnames_int = if_else(breakpoints_df$seqnames == "chrY", 23, breakpoints_df$seqnames_int))
  breakpoints_df = arrange(breakpoints_df, seqnames_int, start)
  
  filter_grid = filter_grid[order(
    -filter_grid$MIN_CN,
    filter_grid$MIN_QUAL_BREAKPOINT,
    -filter_grid$MIN_QUAL_VAF,
    filter_grid$MIN_SV_LEN,
    -filter_grid$ECDNA_THRESHOLD
  ), ]
  
  prev_min_cn = -1L
  nodes_iter  = NULL
  edges_iter  = NULL
  all_results = vector("list", nrow(filter_grid))
  
  for (i in 1:nrow(filter_grid)) {
    filters              = as.list(filter_grid[i, ])
    current_min_cn       = as.integer(filters$MIN_CN)
    ecdna_threshold_used = ploidy + as.numeric(filters$ECDNA_THRESHOLD)
    
    message("Testing row ", i,
            " MIN_CN=", current_min_cn,
            " VAF=", filters$MIN_QUAL_VAF,
            " ECDNA_THRESHOLD=", ecdna_threshold_used)
    
    if (current_min_cn != prev_min_cn) {
      message("Rebuilding graph for MIN_CN=", current_min_cn)
      hlfa_iter = get_hlfa(all_cn, current_min_cn, min_width, 1e7, chromosomes)
      
      if (nrow(hlfa_iter) == 0 || !any(hlfa_iter$cn_int >= current_min_cn)) {
        message("No amplified segments at MIN_CN=", current_min_cn, " — skipping")
        prev_min_cn = current_min_cn
        next
      }
      
      graph_iter = create_graph_components(hlfa_iter, breakpoints_df, all_cn, current_min_cn, min_width)
      nodes_iter = graph_iter$nodes
      edges_iter = graph_iter$edges
      nodes_iter = dplyr::mutate(nodes_iter, seqnames_int = as.numeric(gsub("chr", "", chr)))
      nodes_iter = dplyr::mutate(nodes_iter, seqnames_int = if_else(nodes_iter$chr == "chrX", 22, nodes_iter$seqnames_int))
      nodes_iter = dplyr::mutate(nodes_iter, seqnames_int = if_else(nodes_iter$chr == "chrY", 23, nodes_iter$seqnames_int))
      nodes_iter = arrange(nodes_iter, seqnames_int, pos)
      prev_min_cn    = current_min_cn
      message("raw edges: ", nrow(edges_iter))

    }
    
    
    if (is.null(nodes_iter)) {
      message("No graph available — skipping row ", i)
      next
    }
    
    
    if (is.null(nodes_iter) || nrow(nodes_iter) == 0) {
      message("No nodes available — skipping row ", i)
      next
    }

    if (is.null(edges_iter) || nrow(edges_iter) == 0) {
       message("No edges in graph — skipping row ", i)
       next
    }

    
    filtered_graph = apply_filters(nodes_iter, edges_iter,
                                   breakpoints_df, filters,
                                   current_min_cn, min_width)

    
    merged_nodes_df = filtered_graph[["merged_nodes"]]
    edges_filtered  = filtered_graph[["filtered_edges"]]
    edges_filtered  = edges_filtered[!duplicated(edges_filtered$ID), ]
    edges_filtered  = add_adjacency_edges(merged_nodes_df, edges_filtered, filters$MIN_SV_LEN)
    
    if (!identical(chromosomes, NA)) {
      edges_filtered = dplyr::filter(edges_filtered, source_chr %in% chromosomes)
    }
    
    scc_result = get_cycles_from_df(edges_filtered, merged_nodes_df)
    
    if (length(scc_result) == 0) {
      message("No cycles at this filter combination")
      next
    }
    
    min_sv_vaf = if (mode == "LR") 0.2 else 0.5
    real_sv_sources = unique(c(
      edges_filtered[svtype %in% c("BND", "DUP") & VAF >= min_sv_vaf, source],
      edges_filtered[svtype %in% c("BND", "DUP") & VAF >= min_sv_vaf, target]
    ))
    
    real_cycles = Filter(function(cycle) {
      if (length(cycle) > 1) return(TRUE)
      seg_nodes = merged_nodes_df[parent == as.numeric(cycle[[1]]), ID]
      any(seg_nodes %in% real_sv_sources)
    }, scc_result)
    
    if (length(real_cycles) == 0) {
      message("Only spurious cycles at this filter combination")
      next
    }
    
    message("Found ", length(real_cycles), " real cycles")
    
    cycles_df  = data.frame()
    cycles_out = ""
    for (j in seq_along(real_cycles)) {
      cycle = real_cycles[[j]]
      cycles_out = paste(cycles_out, "Cycle", j, "\n")
      
      walk_results = if (length(cycle) > 1) {
  tryCatch({
    find_closed_walk_cn(cycle, edges_filtered, merged_nodes_df)
  }, error = function(e) {
    message("Walk CN failed for cycle ", j, ": ", e$message)
    fallback_cn = min(merged_nodes_df[parent %in% as.numeric(unlist(cycle)) & cn > 0, cn], na.rm=TRUE)
    list(list(cn = fallback_cn, segments = as.numeric(unlist(cycle))))
  })
} else {
  list(list(cn = NULL, segments = as.numeric(unlist(cycle))))
}

      walk_cn = if (length(walk_results) > 0) {
        cns = sapply(walk_results, function(w) w$cn)
        cns = cns[!sapply(cns, is.null)]
        if (length(cns) > 0) max(unlist(cns)) else NULL
      } else {
        NULL
      }

            for (k in seq_along(cycle)) {
              parent_id = as.numeric(cycle[[k]])
              
              node = dplyr::filter(merged_nodes_df, parent == parent_id)
              if (nrow(node) == 0) next
              
              start_node = node[node$start_flag == TRUE,  "ID"][[1]]
              end_node = node[node$start_flag == FALSE, "ID"][[1]]
              
              sv_from_end   = any(edges_filtered$source == end_node   & 
                                    edges_filtered$svtype %in% c("BND","DEL","DUP"))
              sv_from_start = any(edges_filtered$source == start_node & 
                                    edges_filtered$svtype %in% c("BND","DEL","DUP"))
              
              strand = if (sv_from_end) "+" else if (sv_from_start) "-" else "+"
              
              cycles_out = paste0(cycles_out, " Segment ", parent_id, " : ",
                                  node[1, "ID"], " - ", node[2, "ID"], "\t", strand, "\n")
              
              cn_val = if (!is.null(walk_cn)) as.integer(round(walk_cn)) else node[[1, "cn"]]
              
              cycles_df = rbind(cycles_df, data.frame(
                cycle_id   = paste0("cycle", j),
                segment_id = node[[1, "parent"]],
                cn         = cn_val,
                strand     = strand
              ))
            }
    }
    
    if (nrow(cycles_df) == 0) next
    
    dc = decompose_cycles(cycles_df, merged_nodes = merged_nodes_df,
                          threshold = ecdna_threshold_used)
    
    if (is.null(dc) || nrow(dc) == 0) next
    
    ecdna_cycles = dc[dc$ecdna == TRUE, ]
    if (nrow(ecdna_cycles) == 0) {
      message("No ecDNA cycles after decomposition")
      next
    }
    
    cycles_dt = data.table::as.data.table(cycles_df)
    merged_dt = data.table::as.data.table(merged_nodes_df)
    
    seg_lengths = merged_dt[start_flag == TRUE, .(
      segment_id = parent, seg_start = pos
    )][merged_dt[start_flag == FALSE, .(
      segment_id = parent, seg_end = pos
    )], on = "segment_id"][, length := seg_end - seg_start + 1L]
    
    cycle_lengths = cycles_dt[seg_lengths, on = "segment_id", nomatch = NULL][
      , .(total_length = sum(length, na.rm = TRUE)), by = cycle_id]
    
    ecdna_dt = data.table::as.data.table(ecdna_cycles)
    ecdna_dt = cycle_lengths[ecdna_dt, on = "cycle_id"]
    ecdna_dt[is.na(total_length), total_length := 0L]
    
    all_results[[i]] = list(
      filters         = filters,
      cycles_df       = cycles_df,
      cycles_out      = cycles_out,
      merged_nodes_df = merged_nodes_df,
      edges_filtered  = edges_filtered,
      dc              = dc,
      ecdna_dt        = ecdna_dt,
      ecdna_threshold = ecdna_threshold_used,
      min_cn          = current_min_cn,
      score           = 0,
      grid_row        = i
    )
  }
  
  # collect valid results
  valid_results = all_results[!sapply(all_results, is.null)]
  
  if (length(valid_results) == 0) {
    message("No cycles found across all filter combinations")
    return(list(
      filtered_edges  = data.table::data.table(),
      merged_nodes    = data.table::data.table(),
      cycles_df       = data.frame(),
      cycles_out      = "",
      ecdna_threshold = ploidy + 2.5,
      min_cn          = as.integer(filter_grid$MIN_CN[1])
    ))
  }
  

  #global_max_total_len = max(sapply(valid_results, function(r) 
  #  sum(r$ecdna_dt$total_length, na.rm=TRUE)))
  #global_max_n_ecdna   = max(sapply(valid_results, function(r) 
  #  nrow(r$dc[r$dc$ecdna == TRUE, ])))
  
  #message("Global max_total_length: ", global_max_total_len)
  #message("Global max_n_ecdna: ", global_max_n_ecdna)
  global_max_total_len = max(sapply(valid_results, function(r) 
  if (nrow(r$ecdna_dt) > 0) median(r$ecdna_dt$total_length, na.rm=TRUE) else 0))

  message("Global max_median_length: ", global_max_total_len) 

  for (k in seq_along(valid_results)) {
    r = valid_results[[k]]
    
    #total_ecdna_len  = sum(r$ecdna_dt$total_length, na.rm=TRUE)
    #total_ecdna_len = if (nrow(r$ecdna_dt) > 0) median(r$ecdna_dt$total_length, na.rm = TRUE) else 0

    #best_cycle_score = total_ecdna_len / global_max_total_len
    
    total_ecdna_len  = if (nrow(r$ecdna_dt) > 0) median(r$ecdna_dt$total_length, na.rm=TRUE) else 0
    best_cycle_score = if (global_max_total_len > 0) total_ecdna_len / global_max_total_len else 0


    #n_ecdna        = nrow(r$dc[r$dc$ecdna == TRUE, ])
    #coverage_score = n_ecdna / global_max_n_ecdna
    
    vaf_s = as.numeric(r$filters$MIN_QUAL_VAF)
    
    bp_range = max(filter_grid$MIN_QUAL_BREAKPOINT) - min(filter_grid$MIN_QUAL_BREAKPOINT)
    bp_s = if (bp_range == 0) 1.0 else
      (as.numeric(r$filters$MIN_QUAL_BREAKPOINT) - min(filter_grid$MIN_QUAL_BREAKPOINT)) / bp_range
    
    sv_range = max(filter_grid$MIN_SV_LEN) - min(filter_grid$MIN_SV_LEN)
    sv_s = if (sv_range == 0) 1.0 else
      (as.numeric(r$filters$MIN_SV_LEN) - min(filter_grid$MIN_SV_LEN)) / sv_range
    

    strictness_score = vaf_s * 0.625 +  bp_s * 0.25 + sv_s * 0.125
    valid_results[[k]]$score = best_cycle_score * 0.75 +
    strictness_score  * 0.25
    
    message("Row ", r$grid_row,
            " MIN_CN=", r$min_cn,
            " VAF=", r$filters$MIN_QUAL_VAF,
            " total_ecdna_len=", total_ecdna_len,
            #" n_ecdna=", n_ecdna,
            " best_cycle_score=", round(best_cycle_score, 3),
            " strictness=", round(strictness_score, 3),
            #" coverage=", round(coverage_score, 3),
            " total=", round(valid_results[[k]]$score, 3))
  }
  
  scores    = sapply(valid_results, function(r) r$score)
  grid_rows = sapply(valid_results, function(r) r$grid_row)
  best_idx  = order(-scores, grid_rows)[1]
  best      = valid_results[[best_idx]]
  
  message("Best filter combination:")
  message(" MIN_CN=", best$min_cn,
          " MIN_QUAL_VAF=", best$filters$MIN_QUAL_VAF,
          " ECDNA_THRESHOLD=", best$ecdna_threshold,
          " score=", round(best$score, 4))
  
  
  message("Cycles were detected")
  message(best$cycles_df)
  
  message("edges_filtered rows: ", nrow(best$edges_filtered))
  

  
  return(list(
    filtered_edges  = best$edges_filtered,
    merged_nodes    = best$merged_nodes_df,
    cycles_df       = best$cycles_df,
    cycles_out      = best$cycles_out,
    ecdna_threshold = best$ecdna_threshold,
    min_cn          = best$min_cn,
    dc              = best$dc
  ))
}



draw_cycles = function(edges_df,nodes_df, cycles_df, annotation_file, cnv_df, output_directory, min_cn){

  graph_dir = output_directory
  #dir.create(graph_dir, recursive = TRUE, showWarnings = FALSE)
  edges_df   = data.table::as.data.table(edges_df)
  nodes_df   = data.table::as.data.table(nodes_df)
  cycles_df  = data.table::as.data.table(cycles_df)
  
  
  message("Loading annotation file")
  gtf        = rtracklayer::import(annotation_file)
  gencode_gr = gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"]
  gene_track = create_gene_track("hg38", gencode_gr)
  rm(gtf) 
  
  if (!is.null(cnv_df) && nrow(cnv_df) > 0) {
    cnv_external_df = data.table::as.data.table(cnv_df)[,
                                                        .(chr = chrom, start, end, cn)
    ]
    cnv_external_df[, chr := data.table::fcase(
      chr %in% c("chrX", "X", "23", 23L), "chrX",
      chr %in% c("chrY", "Y", "24", 24L), "chrY",
      !startsWith(chr, "chr"), paste0("chr", chr),
      rep(TRUE, .N),chr
    )]
  } else {
    cnv_external_df = NULL
  }
  
  nodes_by_parent = split(nodes_df, by = "parent")
  cycles_by_id    = split(cycles_df, by = "cycle_id")
  
  cycle_ids = unique(cycles_df$cycle_id)
  
  corrected_strands = list()  
  #process each cycle
  for(cycle in cycle_ids){
    #get segments belonging to this cycle
    
    cycle_segments = cycles_by_id[[as.character(cycle)]]
    segment_ids    = cycle_segments$segment_id
    
    #get the corresponding nodes
    all_cycle_nodes = data.table::rbindlist(
      nodes_by_parent[as.character(segment_ids)],
      use.names = TRUE
    )
    
    if (nrow(all_cycle_nodes) == 0L) next
    if (!any(all_cycle_nodes$cn >= min_cn)) next
    
    #only amplified nodes for visualization
    cycle_nodes = all_cycle_nodes[cn >= min_cn]
    if (nrow(cycle_nodes) == 0L) next
    start_nodes = cycle_nodes[start_flag == TRUE,  .(parent, chr, cn, start = pos)]
    end_nodes   = cycle_nodes[start_flag == FALSE, .(parent, end = pos)]
    
    
    nodes_segments_df = end_nodes[start_nodes, on = "parent", nomatch = NULL]

    
    message("nodes_segments_df parents: ", paste(nodes_segments_df$parent, collapse=","))
    message("duplicated parents: ", paste(nodes_segments_df$parent[duplicated(nodes_segments_df$parent)], collapse=","))
    message("start_nodes nrow: ", nrow(start_nodes))
    message("end_nodes nrow: ", nrow(end_nodes))
    
    if (any(duplicated(nodes_segments_df$parent))) {
      nodes_segments_df = nodes_segments_df[!duplicated(nodes_segments_df$parent)]
    }
    
    nodes_segments_df[, `:=`(
      strand      = "+",
      circle_id   = cycle,
      fragment_id = paste0("seg_", parent - 1L),
      idx         = .I,
      start_id    = paste0(chr, ":", start),
      end_id      = paste0(chr, ":", end)
    )]
    
    all_start = all_cycle_nodes[start_flag == TRUE  & cn >= min_cn, .(parent, chr, pos)]
    all_end   = all_cycle_nodes[start_flag == FALSE & cn >= min_cn, .(parent, chr, pos)]
    
    amplified_parents = nodes_segments_df$parent
    parent_to_idx = setNames(nodes_segments_df$idx, nodes_segments_df$parent)


    all_id_to_parent = setNames(all_cycle_nodes$parent, all_cycle_nodes$ID)
    all_id_to_side    = setNames(ifelse(all_cycle_nodes$start_flag, 0L, 1L), all_cycle_nodes$ID)
   
   cycle_edges = edges_df[
      source %in% all_cycle_nodes$ID & target %in% all_cycle_nodes$ID
    ]

    resolve_to_amplified = function(node_id, visited = character(0)) {
    if (node_id %in% visited) return(NA_character_) 
    parent = all_id_to_parent[node_id]
    if (is.na(parent)) return(NA_character_)
    if (parent %in% amplified_parents) return(node_id)  

    visited = c(visited, node_id)

    other_end = all_cycle_nodes[parent == as.numeric(parent) & ID != node_id, ID]
    if (length(other_end) == 0L) return(NA_character_)
    other_end = other_end[1]
 
    next_hop = cycle_edges[source == other_end, target]
    if (length(next_hop) == 0L) next_hop = cycle_edges[target == other_end, source]
    if (length(next_hop) == 0L) return(NA_character_)

    resolve_to_amplified(next_hop[1], visited)
    }
    resolved_source = sapply(edges_df$source, resolve_to_amplified)
    resolved_target = sapply(edges_df$target, resolve_to_amplified)

    src      = parent_to_idx[as.character(all_id_to_parent[resolved_source])]
    tgt      = parent_to_idx[as.character(all_id_to_parent[resolved_target])]
    src_side = all_id_to_side[resolved_source]
    tgt_side = all_id_to_side[resolved_target]
    # snap_to_amplified = function(chr_val, pos_val, nodes_segments_df) {
    #   candidates = nodes_segments_df[chr == chr_val]
    #   if (nrow(candidates) == 0L) return(NA_integer_)
    #   dists = pmin(abs(candidates$start - pos_val), abs(candidates$end - pos_val))
    #   candidates$idx[which.min(dists)]
    # }
    
    # start_ids  = paste0(all_start$chr, ":", all_start$pos)
    # start_idxs = mapply(snap_to_amplified,
    #                     all_start$chr, all_start$pos,
    #                     MoreArgs = list(nodes_segments_df = nodes_segments_df))
    # start_lookup = setNames(as.integer(start_idxs), start_ids)
    
    # end_ids  = paste0(all_end$chr, ":", all_end$pos)
    # end_idxs = mapply(snap_to_amplified,
    #                   all_end$chr, all_end$pos,
    #                   MoreArgs = list(nodes_segments_df = nodes_segments_df))
    # end_lookup = setNames(as.integer(end_idxs), end_ids)
    
    # src_from_start    = start_lookup[edges_df$source]
    # src_from_end      = end_lookup[edges_df$source]
    # target_from_start = start_lookup[edges_df$target]
    # target_from_end   = end_lookup[edges_df$target]
    
    # src      = ifelse(!is.na(src_from_start), src_from_start, src_from_end)
    # src_side = ifelse(!is.na(src_from_start), 0L, 1L)
    # tgt      = ifelse(!is.na(target_from_start), target_from_start, target_from_end)
    # tgt_side = ifelse(!is.na(target_from_start), 0L, 1L)
    
    valid = !is.na(src) & !is.na(tgt)
    
    if (!any(valid)) next
    
    matched_edges = data.table::data.table(
      circle_id = cycle,
      n1        = src[valid],
      n2        = tgt[valid],
      n1.side   = src_side[valid],
      n2.side   = tgt_side[valid]
    )

    
    non_self   = matched_edges[n1 != n2]
    self_loops = unique(matched_edges[n1 == n2], by = "n1")
    
  
    real_selfloop_positions = unique(edges_df[source == target, source])
    
    message("non_self BEFORE norm: ", sum(non_self$n2.side == 1L))
    non_self[n1.side == 1L & n2.side == 1L, n2.side := 0L]
    message("non_self AFTER norm: ", sum(non_self$n2.side == 1L))
    
    
    self_loops[, is_real_sv := (
      nodes_segments_df$start_id[n1] %in% real_selfloop_positions |
        nodes_segments_df$end_id[n1]   %in% real_selfloop_positions
    )]
    
    nodes_with_any_edge = unique(c(non_self$n1, non_self$n2))
    self_loops = self_loops[!n1 %in% nodes_with_any_edge | is_real_sv == TRUE]
    self_loops[, is_real_sv := NULL]
     
    self_loops[n1 == n2, `:=`(n1.side = 1L, n2.side = 0L)]
   
    matched_edges = data.table::rbindlist(list(non_self, self_loops))
    
    message("matched_edges AFTER rebuild n2.side==1: ", sum(matched_edges$n2.side == 1L))
    
    matched_edges[, col := alpha("orange", 0.8)]
    
    nodes_gr = GenomicRanges::makeGRangesFromDataFrame(
      as.data.frame(nodes_segments_df),
      keep.extra.columns = TRUE,
      ignore.strand      = FALSE
    )
    nodes_gr$coverage = nodes_gr$cn
    
    if (!is.null(cnv_external_df)) {
      cnv_gr = GenomicRanges::makeGRangesFromDataFrame(
        as.data.frame(cnv_external_df),
        keep.extra.columns = TRUE,
        ignore.strand      = FALSE
      )
        cnv_gr = IRanges::subsetByOverlaps(cnv_gr, nodes_gr + 500000)
    } else {
      cnv_gr = NULL
    }
    
    
    cn_track = create_copy_number_track (nodes_gr)
    cnv_track = create_external_copy_number_track (cnv_gr)
    
    
    GenomeInfoDb::seqinfo(nodes_gr) = GenomeInfoDb::seqinfo(
      BSgenome.Hsapiens.UCSC.hg38
    )[GenomeInfoDb::seqlevels(nodes_gr)]
    
    graph_filename =  paste(graph_dir,"/",cycle,".png",sep="")
    
    message("columns: ", paste(names(matched_edges), collapse=","))
    for (i in seq_len(nrow(matched_edges))) {
      message("edge ", i, ": n1=", matched_edges$n1[i], 
              " n2=", matched_edges$n2[i],
              " n1.side=", matched_edges$n1.side[i],
              " n2.side=", matched_edges$n2.side[i],
              " col=", matched_edges$col[i])
    }
    print(nodes_segments_df[nodes_segments_df$idx %in% 1:15, 
                            c("idx","parent","chr","start","end","strand")])
    #matched_edges = matched_edges[matched_edges$n1 != matched_edges$n2, ]
    
    result = draw_amplicon( graph_filename, cycle, nodes_gr,  matched_edges, cn_track , cnv_track, gene_track, outdir = output_directory)
    
  }
  
  if (length(corrected_strands) > 0) {
 
    cycles_file = file.path(output_directory, "cycles.txt")
    lines = readLines(cycles_file)
    
    current_cycle = NA
    new_lines = character(length(lines))
    
    for (i in seq_along(lines)) {
      line = lines[i]
 
      cycle_match = regmatches(line, regexpr("Cycle\\s+(\\S+)", line))
      if (length(cycle_match) > 0) {
        current_cycle = trimws(gsub("Cycle\\s+", "", cycle_match))
      }
      
    
      if (!is.na(current_cycle) && 
          current_cycle %in% names(corrected_strands) &&
          grepl("^\\s*Segment", line)) {
        
        
        parent_match = regmatches(line, regexpr("Segment\\s+(\\d+)", line))
        if (length(parent_match) > 0) {
          parent_id = as.integer(gsub("Segment\\s+", "", parent_match))
          if (parent_id %in% corrected_strands[[current_cycle]]) {
          
            line = sub("\\+$", "-", trimws(line))
            line = paste0(" ", line)
          }
        }
      }
      new_lines[i] = line
    }
    
    writeLines(new_lines, cycles_file)
  }
  
}


#VISUALIZATION FUNCTIONS

draw_amplicon = function(output_filename, circle_id, nodes_gr, edges_df, cn_track, cnv_track, gene_track, maxcov = 500,outdir = ".") {
  
  
  filtered_nodes = nodes_gr[elementMetadata(nodes_gr)[,"circle_id"] == circle_id]
  filtered_edges = edges_df[edges_df$circle_id == circle_id,]
  
  idx_map = setNames(seq_len(length(filtered_nodes)), as.character(filtered_nodes$idx))
  filtered_edges$n1 = as.integer(idx_map[as.character(filtered_edges$n1)])
  filtered_edges$n2 = as.integer(idx_map[as.character(filtered_edges$n2)])
  filtered_edges = filtered_edges[!is.na(filtered_edges$n1) & !is.na(filtered_edges$n2), ]
  
  filtered_edges = data.table::as.data.table(filtered_edges)
  filtered_edges[n1 == n2, `:=`(n1.side = 1L, n2.side = 0L)]
  filtered_edges[n1.side == 1L & n2.side == 1L & n1 != n2, n2.side := 0L]
  
  
  filtered_nodes$col = alpha("gray40", 0.7)
  filtered_edges$col = alpha("orange", 0.8)
  
  
  filtered_nodes_gg = filtered_nodes
  GenomicRanges::end(filtered_nodes_gg) = GenomicRanges::end(filtered_nodes_gg) - 1L



  chroms = unique(as.character(GenomeInfoDb::seqnames(filtered_nodes)))
  is_interchromosomal = length(chroms) > 1
  
  if (is_interchromosomal) {
    
    GenomeInfoDb::seqinfo(filtered_nodes) = GenomeInfoDb::seqinfo(
      BSgenome.Hsapiens.UCSC.hg38
    )[GenomeInfoDb::seqlevels(filtered_nodes)]

    graph = gG(nodes = filtered_nodes_gg,
               edges = filtered_edges,
               meta  = list(sep.bg.col = alpha("white", 0),
                            sep.draw   = FALSE,
                            bundle        = FALSE,
                            draw.ref   = FALSE))

    gt_amp = graph$gt
    message("all formatting slots:")
    print(names(formatting(gt_amp)))
    
    gt_amp$name       = circle_id
    gt_amp$sep.bg.col = NA
    gt_amp$xaxis.cex  = 0.2
    gt_amp$sep.draw   = FALSE
    
    formatting(gt_amp)$refcol = "transparent"
    formatting(gt_amp)$ref.col = "transparent" 
    formatting(gt_amp)$junction.col = alpha("orange", 0.8)
    formatting(gt_amp)$col = alpha("orange", 0.8)
    
    formatting(gt_amp)$ref.col     = alpha("white", 0)  
    formatting(gt_amp)$edge.ref.col = alpha("white", 0)
    
    formatting(gt_amp)$height = if (nrow(filtered_edges) <= 2) 5 else 20
    formatting(gt_amp)$y0       = NA
    formatting(gt_amp)$sep.draw = FALSE
    formatting(gt_amp)$ygap     = 20
    formatting(gt_amp)$ywid = max(0.1, min(0.5, 5 / length(filtered_nodes)))
    #formatting(gt_amp)$ywid = max(0.3, min(1.5, length(filtered_nodes) * 0.1))

    windows = do.call(c, lapply(chroms, function(chrom) {
      chrom_nodes = filtered_nodes[
        as.character(GenomeInfoDb::seqnames(filtered_nodes)) == chrom
      ]
      GenomicRanges::reduce(chrom_nodes + 200000)
    }))
    windows = sort(windows)
    

    GenomeInfoDb::seqinfo(windows) = GenomeInfoDb::seqinfo(
      BSgenome.Hsapiens.UCSC.hg38
    )[GenomeInfoDb::seqlevels(windows)]
    

    cn_nodes = cn_track$gr
    if (!is.null(cn_nodes)) {
      common_seqlevels = intersect(
        GenomeInfoDb::seqlevels(cn_nodes),
        GenomeInfoDb::seqlevels(BSgenome.Hsapiens.UCSC.hg38)
      )
      GenomeInfoDb::seqinfo(cn_track$gr) = GenomeInfoDb::seqinfo(
        BSgenome.Hsapiens.UCSC.hg38
      )[common_seqlevels]
    }
    
    png_width = max(5000, 3000 * length(chroms))
    if (png_width > 20000) {
       message("Skipping plot for ", circle_id, ": png_width (", png_width, ") exceeds safe rendering limit — likely reflects an oversized/merged cycle, not a rendering issue")
        file.create(file.path(outdir, paste0(circle_id, "_plot_skipped_too_large.txt")))
   } else {
    png(output_filename, width = png_width, height = 4000, res = 300)
    par(mar = c(8, 2, 0, 20))
    
    formatting(gene_track)$ygap = 10
    formatting(cn_track)$ygap   = 10
    
    if (!is.null(cnv_track)) {
      formatting(cnv_track)$ygap = 20
      gTrack::plot.gTrack(
        c(gene_track, cn_track, cnv_track, gt_amp),
        windows,
        sep.draw  = FALSE,
        y.heights = c(0.1, 0.15, 0.15, 0.4)
      )
    } else {
      gTrack::plot.gTrack(
        c(gene_track, cn_track, gt_amp),
        windows,
        sep.draw  = FALSE,
        y.heights = c(0.1, 0.15, 0.4)
      )
    }
    
    dev.off()
    return(invisible(NULL))
    } 
  } else {
  
    graph = gG(nodes = filtered_nodes_gg,
               edges = filtered_edges,
               meta  = list(sep.bg.col = alpha("white", 0),
                            ref.col    = alpha("white", 0)))
    
    gt_amp = graph$gt
    gt_amp$name       = circle_id
    gt_amp$sep.bg.col = NA
    gt_amp$xaxis.cex  = 0.2
    windows = GenomicRanges::reduce(copy(nodes_gr) + 200000)
    
    formatting(gt_amp)$ref.col     = alpha("white", 0)  
    formatting(gt_amp)$edge.ref.col = alpha("white", 0)
    
    formatting(gt_amp)$height = if (nrow(filtered_edges) <= 2) 5 else 20
    formatting(gt_amp)$y0       = NA
    formatting(gt_amp)$sep.draw = FALSE
    formatting(gt_amp)$ygap     = 20
    formatting(gt_amp)$ywid = max(0.5, nrow(filtered_nodes) * 0.5)
    
    png(output_filename, width = 5000, height = 4000, res = 300)
    formatting(gene_track)$ygap = 10
    formatting(cn_track)$ygap   = 10
    
    if (!is.null(cnv_track)) {
      formatting(cnv_track)$ygap = 20
      gTrack::plot.gTrack(
        c(gene_track, cn_track, cnv_track, gt_amp),
        windows,
        sep.draw  = FALSE,
        y.heights = c(0.1, 0.15, 0.15, 0.4)
      )
    } else {
      gTrack::plot.gTrack(
        c(gene_track, cn_track, gt_amp),
        windows,
        sep.draw  = FALSE,
        y.heights = c(0.1, 0.15, 0.4)
      )
    }
    
    dev.off()
    return(gt_amp)
  }
}




create_copy_number_track = function(nodes_gr, maxcn = NULL){
  if(is.null(maxcn)){
    maxcn <- max(nodes_gr$cn, na.rm = TRUE) * 1.2
  }
  gt_cov= gTrack(nodes_gr,
                 y.field = "cn",
                 line = TRUE,
                 name = "CN",
                 y0 = 0,
                 y1 = maxcn,
                 height = 30)
  
  return(gt_cov)
}

create_external_copy_number_track = function(cnv_gr, maxcn = NULL){
  if(!is.null(cnv_gr))
  {
    if(is.null(maxcn)){
      maxcn <- max(cnv_gr$cn, na.rm = TRUE) * 1.2
    }
  
    gt_cnv= gTrack(cnv_gr,
                   y.field = "cn",
                   line = TRUE,
                   name = "FACETS",
                   y0 = 0,
                   y1 = maxcn,
                   height = 30)
    
    return(gt_cnv)
  }else{
    return (NULL)
  }
}


create_gene_track <- function(build = "hg38", gencode_gr) {
  gencode_gr$col= alpha("#0B6623",0.3)
  gencode = gTrack(gencode_gr,
                    name = "Genes",
                    gr.labelfield = "gene_name",
                    gr.cex.label = 0.5)
  gencode$height = 40 

  return(gencode)
}