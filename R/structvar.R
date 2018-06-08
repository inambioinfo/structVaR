
#' Function that maps the translocation mates, and assigns translocation name to variable ABERRATION (using standard nomenclature, e.g. t(14;18)(p21;q33))
#'
#' @param manta_df Manta data frame with translocation data
#'
map_translocation_mate <- function(manta_df){

  mate_mapped_translocations <- data.frame()
  if('ID' %in% colnames(manta_df) && 'MATEID' %in% colnames(manta_df) && 'CYTOBAND' %in% colnames(manta_df) && 'VCF_SAMPLE_ID' %in% colnames(manta_df) && 'ALT' %in% colnames(manta_df)){
    manta_mate_data <- dplyr::select(manta_df, MATEID, ALT, VCF_SAMPLE_ID, CYTOBAND, GENE_LOCUS) %>%
      dplyr::rename(CYTOBAND_PARTNER = CYTOBAND, ALT_PARTNER = ALT, GENE_LOCUS_PARTNER = GENE_LOCUS)

    translocations <- dplyr::left_join(manta_df, manta_mate_data, by=c("ID" = "MATEID", "VCF_SAMPLE_ID" = "VCF_SAMPLE_ID"))
    translocations <- tidyr::separate(translocations, CYTOBAND, c('cyto_chrom','cyto_band'),sep=":")
    translocations <- tidyr::separate(translocations, CYTOBAND_PARTNER, c('cyto_chrom_partner','cyto_band_partner'),sep=":")

    translocations$ABERRATION <- paste0("t(",translocations$cyto_chrom,";",translocations$cyto_chrom_partner,")(",translocations$cyto_band,";",translocations$cyto_band_partner,")")
    translocations <- dplyr::filter(translocations, !(cyto_chrom_partner == 'X' | cyto_chrom_partner == 'Y' | cyto_chrom_partner == 'M' | stringr::str_detect(cyto_chrom_partner,"Un_|_random|_hap")))
    translocations_sex_mito <- dplyr::filter(translocations, cyto_chrom == 'X' | cyto_chrom == 'Y' | cyto_chrom == 'M') %>%
      dplyr::select(-c(chromosome, cyto_chrom, cyto_chrom_partner, cyto_band, cyto_band_partner))
    translocations_auto <- dplyr::filter(translocations, !(cyto_chrom == 'X' | cyto_chrom == 'Y' | cyto_chrom == 'M' | stringr::str_detect(cyto_chrom,"Un_|_random|_hap")))
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
      rlogging::stop(paste0('Manta input data (vcf2tsv output) is missing a required column - ', col, ' - EXITING'))
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

  if(!file.exists(manta_vcf2tsv_fname)){
    rlogging::stop(paste0('Input file ',manta_vcf2tsv_fname, ' does not exist - EXITING'))
  }
  if(!(build == 'grch37' || build == 'grch38')){
    rlogging::stop(paste0('Genome build must be one of grch37/grch38 - current input is ',build))
  }
  rlogging::message('Reading somatic Manta data (i.e. \'somaticSV.vcf.gz\'), VCF converted to tab-separated output with https://github.com/sigven/vcf2tsv')
  rlogging::message('Filename: ', manta_vcf2tsv_fname)
  manta_tsv_raw <- read.table(file=manta_vcf2tsv_fname,skip = 1,header=T,stringsAsFactors = F,quote = "")
  rlogging::message('Validating Manta data - check that all required columns are present')
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

  ## check the existence of genotypes for control sample
  control_genotypes_found <- TRUE
  manta_df_control <- manta_df %>% dplyr::filter(stringr::str_detect(VCF_SAMPLE_ID,control_sample_pattern))
  if(nrow(manta_df_control) == 0){
    rlogging::warning(paste0('Not able to retrieve genotype data for the control sample - regular expression (',control_sample_pattern,') not correctly specified?'))
    control_genotypes_found <- FALSE
  }else{
    manta_df_control <- manta_df_control %>%
      dplyr::rename(PAIRED_READ_SUPPORT_CONTROL = PR, SPLIT_READ_SUPPORT_CONTROL = SR) %>%
      dplyr::select(ID,PAIRED_READ_SUPPORT_CONTROL,SPLIT_READ_SUPPORT_CONTROL)
  }

  manta_df <- manta_df %>%
    dplyr::select(chromosome,POS,ALT,segmentID,VCF_SAMPLE_ID,SVTYPE,ID,FILTER,BND_DEPTH,MATE_BND_DEPTH,MATEID,SOMATICSCORE,PR,SR) %>%
    dplyr::filter(!stringr::str_detect(VCF_SAMPLE_ID,control_sample_pattern)) %>%
    dplyr::rename(PAIRED_READ_SUPPORT_TUMOR = PR, SPLIT_READ_SUPPORT_TUMOR = SR) %>%
    dplyr::select(-FILTER)

  if(control_genotypes_found){
    manta_df <- dplyr::left_join(manta_df, manta_df_control, by=c("ID"))
  }else{
    manta_df$PAIRED_READ_SUPPORT_CONTROL <- NA
    manta_df$SPLIT_READ_SUPPORT_CONTROL <- NA
  }

  if(translocations_only){
    manta_df <- dplyr::filter(manta_df, SVTYPE == 'BND')
  }

  manta_gr <- NULL
  if(nrow(manta_df) > 0){
    manta_df$SOMATIC_CALL_CONFIDENCE <- 'High'
    if(nrow(manta_df[as.integer(manta_df$SOMATICSCORE) < 30,]) > 0){
      manta_df[as.integer(manta_df$SOMATICSCORE) < 30,]$SOMATIC_CALL_CONFIDENCE <- 'Low'
    }

    gencode_gr <- NULL
    if(build == 'grch37'){
      manta_gr <- GenomicRanges::makeGRangesFromDataFrame(manta_df, keep.extra.columns = T, seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19), seqnames.field = 'chromosome',start.field = 'POS', end.field = 'POS', ignore.strand = T, starts.in.df.are.0based = T)
      gencode_gr <- structVaR::gencode_grch37
    }
    if(build == 'grch38'){
      manta_gr <- GenomicRanges::makeGRangesFromDataFrame(manta_df, keep.extra.columns = T, seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38), seqnames.field = 'chromosome',start.field = 'POS', end.field = 'POS', ignore.strand = T, starts.in.df.are.0based = T)
      gencode_gr <- structVaR::gencode_grch38
    }

    #nearest_gene_locus_index <- GenomicRanges::nearest(manta_gr, gencode_gr)
    #nearest_gene_locus <- as.data.frame(mcols(gencode_gr[nearest_gene_locus_index,]))

    nearest_gene_locus_hits <- GenomicRanges::nearest(manta_gr, gencode_gr, select="all")
    ranges <- gencode_gr[subjectHits(nearest_gene_locus_hits)]
    mcols(ranges) <- c(mcols(ranges),mcols(manta_gr[queryHits(nearest_gene_locus_hits)]))
    gene_df <- as.data.frame(mcols(ranges))
    gene_df <- dplyr::select(gene_df, symbol, segmentID,tsgene,ts_oncogene) %>% dplyr::group_by(segmentID) %>% dplyr::summarise(GENE_LOCUS = paste(sort(unique(symbol)),collapse="|"))
    manta_df <- dplyr::left_join(manta_df, gene_df, by=c("segmentID"))
    cytoband_map <- structVaR::map_cytoband(manta_gr, build = build)
    manta_df <- dplyr::left_join(manta_df, cytoband_map, by=c("segmentID"))
    manta_df <- structVaR::map_translocation_mate(manta_df)
    manta_df$ID_MATEID <- paste0(manta_df$ID, '_', manta_df$MATEID)
    manta_df$BREAKENDS <- paste0(manta_df$ALT, ' - ', manta_df$ALT_PARTNER)
    manta_df$GENE_LOCI_BREAKPOINTS <- paste0(manta_df$GENE_LOCUS, ';', manta_df$GENE_LOCUS_PARTNER)
    manta_df <- dplyr::select(manta_df, -c(ALT, ALT_PARTNER, ID, MATEID, GENE_LOCUS, GENE_LOCUS_PARTNER, SVTYPE))
    manta_df <- dplyr::left_join(manta_df, structVaR::mitelman_aberrations, by=c("ABERRATION"))
  }else{
    rlogging::message('Zero translocation events detected')
  }

  return(list('df' = manta_df, 'gr' = manta_gr))

}


