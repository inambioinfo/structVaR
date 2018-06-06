
#' Function that maps the translocation mates, and assigns translocation name to variable ABERRATION (using standard nomenclature, e.g. t(14;18)(p21;q33))
#'
#' @param manta_df Manta data frame with translocation data
#'
map_translocation_mate <- function(manta_df){

  mate_mapped_translocations <- data.frame()
  if('ID' %in% colnames(manta_df) && 'MATEID' %in% colnames(manta_df) && 'CYTOBAND' %in% colnames(manta_df) && 'VCF_SAMPLE_ID' %in% colnames(manta_df) && 'ALT' %in% colnames(manta_df)){
    manta_mate_data <- dplyr::select(manta_df, MATEID, ALT, VCF_SAMPLE_ID, CYTOBAND) %>%
      dplyr::rename(CYTOBAND_PARTNER = CYTOBAND, ALT_PARTNER = ALT)

    translocations <- dplyr::left_join(manta_df, manta_mate_data, by=c("ID" = "MATEID", "VCF_SAMPLE_ID" = "VCF_SAMPLE_ID"))
    translocations <- tidyr::separate(translocations, CYTOBAND, c('cyto_chrom','cyto_band'),sep=":")
    translocations <- tidyr::separate(translocations, CYTOBAND_PARTNER, c('cyto_chrom_partner','cyto_band_partner'),sep=":")

    translocations$ABERRATION <- paste0("t(",translocations$cyto_chrom,";",translocations$cyto_chrom_partner,")(",translocations$cyto_band,";",translocations$cyto_band_partner,")")
    translocations <- dplyr::filter(translocations, !(cyto_chrom_partner == 'X' | cyto_chrom_partner == 'Y' | cyto_chrom_partner == 'M'))
    translocations_sex_mito <- dplyr::filter(translocations, cyto_chrom == 'X' | cyto_chrom == 'Y' | cyto_chrom == 'M') %>%
      dplyr::select(-c(chromosome, cyto_chrom, cyto_chrom_partner, cyto_band, cyto_band_partner))
    translocations_auto <- dplyr::filter(translocations, !(cyto_chrom == 'X' | cyto_chrom == 'Y' | cyto_chrom == 'M'))
    translocations_auto$cyto_chrom <- as.integer(translocations_auto$cyto_chrom)
    translocations_auto$cyto_chrom_partner <- as.integer(translocations_auto$cyto_chrom_partner)
    translocations_auto <- dplyr::filter(translocations_auto, cyto_chrom <= cyto_chrom_partner) %>%
      dplyr::select(-c(chromosome, cyto_chrom, cyto_chrom_partner, cyto_band, cyto_band_partner))

    mate_mapped_translocations <- rbind(translocations_sex_mito, translocations_auto)
  }
  return(mate_mapped_translocations)
}


valid_manta_data <- function(manta_df_raw){
  for(col in c("CHROM","POS","ALT","VCF_SAMPLE_ID","SVTYPE","ID","FILTER","IMPRECISE","INV3","INV5","BND_DEPTH","MATE_BND_DEPTH","MATEID","SOMATICSCORE","PR","SR")){
    if(!(col %in% colnames(manta_df_raw))){
      rlogging::stop(paste0('Manta input VCF data is missing a required column - ', col, ' - EXITING'))
    }
  }
}


#' Function that maps the chromosome bands of chromosomal breakpoints (SVs)
#'
#' @param break_gr genomic ranges object with genomic breakpoints
#' @param build genomic assembly (grch37/grch38)
#' @param sub_bands logical indicating if subcytobands should be included
#'
map_cytoband <- function(break_gr, build = 'grch37', sub_bands = F){

  cytoband_gr <- structVaR::cytoband_gr_hg38
  if(build == 'grch37'){
    cytoband_gr <- structVaR::cytoband_gr_hg19
  }
  cyto_hits <- GenomicRanges::findOverlaps(break_gr, cytoband_gr, type="any", select="all")
  ranges <- cytoband_gr[subjectHits(cyto_hits)]
  mcols(ranges) <- c(mcols(ranges),mcols(break_gr[queryHits(cyto_hits)]))
  cyto_df <- as.data.frame(mcols(ranges))
  cyto_df$segment_start <- start(ranges(break_gr[queryHits(cyto_hits)]))
  cyto_df$segment_end <- end(ranges(break_gr[queryHits(cyto_hits)]))
  cyto_df$segment_length <- width(ranges(break_gr[queryHits(cyto_hits)]))

  cytoband_map <- as.data.frame(dplyr::group_by(cyto_df,segmentID,segment_length) %>% dplyr::summarise(CYTOBAND = paste(name2, collapse=", "), chromosome_arm = paste(unique(arm), collapse=",")))
  if(sub_bands == T){
    cytoband_map <- as.data.frame(dplyr::group_by(cyto_df,segmentID,segment_length) %>% dplyr::summarise(CYTOBAND = paste(name, collapse=", "), chromosome_arm = paste(unique(arm), collapse=",")))
  }
  cytoband_map$CYTOBAND <- stringr::str_replace(cytoband_map$CYTOBAND, pattern = ", (\\S+, ){0,}", replacement = " - ")
  cytoband_map <- tidyr::separate(cytoband_map,segmentID,sep=":",into = c('chrom','start','stop'),remove=F)
  cytoband_map$CYTOBAND <- paste0(cytoband_map$chrom,":",cytoband_map$CYTOBAND)
  cytoband_map$CYTOBAND <- stringr::str_replace_all(cytoband_map$CYTOBAND,"chr","")

  cytoband_map <- dplyr::select(cytoband_map, segmentID, CYTOBAND)

  return(cytoband_map)
}

#' Function that maps the chromosome bands of chromosomal breakpoints (SVs)
#'
#' @param manta_vcf2tsv_fname Manta somatic SV candidates, VCF to TSV converted with https://github.com/sigven/vcf2tsv
#' @param build genomic build (grch37/grch38)
#' @param translocations_only consider translocations only (SVTYPE = BND)
#' @param control_sample_pattern regular expression that matches control sample in VCF

parse_manta_sv <- function(manta_vcf2tsv_fname, build = 'grch37', translocations_only = T, control_sample_pattern = '-N01'){

  rlogging::message('Reading somatic Manta data (i.e. \'somaticSV.vcf.gz\'), VCF converted to tab-separated output with https://github.com/sigven/vcf2tsv')
  rlogging::message('Filename: ', manta_vcf2tsv_fname)
  manta_tsv_raw <- read.table(file=manta_vcf2tsv_fname,skip = 1,header=T,stringsAsFactors = F,quote = "")
  rlogging::message('Validating Manta data')
  structVaR::valid_manta_data(manta_tsv_raw)

  manta_df <- dplyr::mutate(manta_tsv_raw, chromosome = CHROM) %>%
    dplyr::select(-c(CHROM)) %>%
    dplyr::distinct()

  rlogging::message('Excluding non-valid chromosomes')
  manta_df <- dplyr::filter(manta_df, chromosome != 'hs37d5')
  if(build == 'grch37'){
    manta_df <- pcgrr::get_valid_chromosomes(manta_df, chromosome_column = 'chromosome', bsg = BSgenome.Hsapiens.UCSC.hg19)
  }
  if(build == 'grch38'){
    manta_df <- pcgrr::get_valid_chromosomes(manta_df, chromosome_column = 'chromosome', bsg = BSgenome.Hsapiens.UCSC.hg38)
  }

  manta_df$segmentID <- paste0(manta_df$chromosome,":",manta_df$POS,":",manta_df$POS)
  for(v in c('SOMATIC','IMPRECISE','INV3','INV5')){
    if(v %in% colnames(manta_df)){
      manta_df[,v] <- as.logical(dplyr::recode(manta_df[,v], True = TRUE, False = FALSE))
    }
  }

  manta_df <- manta_df %>%
    dplyr::select(chromosome,POS,ALT,segmentID,VCF_SAMPLE_ID,SVTYPE,ID,FILTER,BND_DEPTH,MATE_BND_DEPTH,MATEID,SOMATICSCORE,PR,SR) %>%
    #dplyr::filter(FILTER == 'PASS') %>%
    dplyr::filter(!stringr::str_detect(VCF_SAMPLE_ID,control_sample_pattern)) %>%
    dplyr::rename(PAIRED_READ_SUPPORT = PR, SPLIT_READ_SUPPORT = SR) %>%
    dplyr::select(-FILTER)


  if(translocations_only){
    manta_df <- dplyr::filter(manta_df, SVTYPE == 'BND')
  }

  manta_gr <- NULL
  if(nrow(manta_df) > 0){
    if(build == 'grch37'){
      manta_gr <- GenomicRanges::makeGRangesFromDataFrame(manta_df, keep.extra.columns = T, seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19), seqnames.field = 'chromosome',start.field = 'POS', end.field = 'POS', ignore.strand = T, starts.in.df.are.0based = T)
    }
    if(build == 'grch38'){
      manta_gr <- GenomicRanges::makeGRangesFromDataFrame(manta_df, keep.extra.columns = T, seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38), seqnames.field = 'chromosome',start.field = 'segment_start', end.field = 'segment_end', ignore.strand = T, starts.in.df.are.0based = T)
    }

    cytoband_map <- structVaR::map_cytoband(manta_gr, build = build)
    manta_df <- dplyr::left_join(manta_df, cytoband_map, by=c("segmentID"))
    manta_df <- structVaR::map_translocation_mate(manta_df)
    manta_df$ID_MATEID <- paste0(manta_df$ID, '_', manta_df$MATEID)
    manta_df$BREAKENDS <- paste0(manta_df$ALT, ' - ', manta_df$ALT_PARTNER)
    manta_df <- dplyr::select(manta_df, -c(ALT, ALT_PARTNER, ID, MATEID, SVTYPE))

    manta_df <- dplyr::left_join(manta_df, structVaR::mitelman_aberrations, by=c("ABERRATION"))
  }

  return(list('df' = manta_df, 'gr' = manta_gr))

}
