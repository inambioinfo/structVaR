library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)


## prepare cytoband into a data frame (hg19)
cytobandURL <- 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBandIdeo.txt.gz'
cytobandFile <- paste0('data-raw/',stringr::str_replace(basename(cytobandURL),pattern = 'txt\\.gz',replacement = paste('grch37.txt.gz')))
if(!file.exists(cytobandFile)){
  download.file(cytobandURL, cytobandFile,quiet = T)
}
cytobandData <- read.table(cytobandFile,sep="\t",stringsAsFactors = F,quote="",comment.char="")
colnames(cytobandData) <- c('chrom','start','end','name','gieStain')
cytobandData$name2 <- stringr::str_replace(cytobandData$name,"\\.[0-9]{1,}$","")
cytobandData$arm <- paste0(cytobandData$chrom,stringr::str_replace(cytobandData$name,'[0-9]{1,}(\\.[0-9]{1,}){0,}',''))
cytobandData[stringr::str_detect(cytobandData$arm,'random|hap|Un'),]$arm <- NA

# make genomic ranges
cytoband_gr_hg19 <- GenomicRanges::makeGRangesFromDataFrame(cytobandData, keep.extra.columns = T, seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19), seqnames.field = 'chrom',start.field = 'start', end.field = 'end', ignore.strand = T, starts.in.df.are.0based = T)


## prepare cytoband into a data frame (hg38)
cytobandURL <- 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBandIdeo.txt.gz'
cytobandFile <- paste0('data-raw/',stringr::str_replace(basename(cytobandURL),pattern = 'txt\\.gz',replacement = paste('grch38.txt.gz')))
if(!file.exists(cytobandFile)){
  download.file(cytobandURL, cytobandFile,quiet = T)
}
cytobandData <- read.table(cytobandFile,sep="\t",stringsAsFactors = F,quote="",comment.char="")
colnames(cytobandData) <- c('chrom','start','end','name','gieStain')
cytobandData$name2 <- stringr::str_replace(cytobandData$name,"\\.[0-9]{1,}$","")
cytobandData$arm <- paste0(cytobandData$chrom,stringr::str_replace(cytobandData$name,'[0-9]{1,}(\\.[0-9]{1,}){0,}',''))
cytobandData[stringr::str_detect(cytobandData$arm,'random|hap|Un'),]$arm <- NA

# make genomic ranges
cytoband_gr_hg38 <- GenomicRanges::makeGRangesFromDataFrame(cytobandData, keep.extra.columns = T, seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38), seqnames.field = 'chrom',start.field = 'start', end.field = 'end', ignore.strand = T, starts.in.df.are.0based = T)


devtools::use_data(cytoband_gr_hg19, overwrite = T)
devtools::use_data(cytoband_gr_hg38, overwrite = T)



## MITELMAN DATABASE

library(magrittr)

##ftp://ftp1.nci.nih.gov/pub/CGAP/mitelman.tar.gz

mitelman_db_release <- '20180214'

pmids <- read.table(file="data-raw/mitelman/reference.dat",header = T,stringsAsFactors = F,sep="\t", quote="",comment.char="")
pmids <- dplyr::select(pmids, REFNO, PUBMED) %>% dplyr::distinct()

disease_codes <- read.table(file="data-raw/mitelman/koder.dat",header = T,stringsAsFactors = F,sep="\t", quote="",comment.char="")
disease_codes <- disease_codes %>%
  dplyr::filter(KODTYP == 'MORPH') %>%
  dplyr::select(KOD, BENAMNING, KORTNAMN) %>%
  dplyr::rename(CANCERCODE = KOD, CANCERTYPE = BENAMNING, CANCERTYPE_ABBR = KORTNAMN)
disease_codes$CANCERTYPE_ABBR <- stringr::str_replace(disease_codes$CANCERTYPE_ABBR,"SA ","SARC ")
disease_codes$CANCERTYPE_ABBR <- stringr::str_replace(disease_codes$CANCERTYPE_ABBR,"CA ","CANCER ")
disease_codes[disease_codes$CANCERTYPE_ABBR == 'DLBL',]$CANCERTYPE_ABBR <- 'DLBCL'
disease_codes[disease_codes$CANCERTYPE_ABBR == 'Mesothel',]$CANCERTYPE_ABBR <- 'MESO'
disease_codes[disease_codes$CANCERTYPE_ABBR == 'MBC NOS',]$CANCERTYPE_ABBR <- 'MBCN NOS'
disease_codes$CANCERTYPE_ABBR <- toupper(disease_codes$CANCERTYPE_ABBR)
mitelman_disease_codes <- disease_codes

molclinabnorm <- read.table(file="data-raw/mitelman/molclinabnorm.dat",header = T,stringsAsFactors = F,sep="\t", quote="",comment.char="")
molclinabnorm <- molclinabnorm %>% dplyr::select(REFNO, INVNO, ORDERNO, ABNORMALITY) %>% dplyr::rename(ABERRATION = ABNORMALITY) %>% dplyr::distinct()
molclinabnorm <- dplyr::filter(molclinabnorm, !stringr::str_detect(ABERRATION,"\\?|\\+|\\-[0-9]{1,}"))

molclingene <- read.table(file="data-raw/mitelman/molclingene.dat",header = T,stringsAsFactors = F,sep="\t", quote="",comment.char="")
molclingene <- molclingene %>% dplyr::select(REFNO, INVNO, ORDERNO, GENE) %>% dplyr::distinct()
molclingene <- dplyr::inner_join(molclinabnorm, molclingene) %>% dplyr::distinct()
molclingene <- dplyr::left_join(molclingene, pmids, by=c("REFNO"))
molclingene <- as.data.frame(dplyr::group_by(molclingene, ABERRATION, GENE) %>% dplyr::summarise(PMID = paste(unique(PUBMED), collapse=",")))

aberration_events <- read.table(file="data-raw/mitelman/recurrent_data.dat",header = T,stringsAsFactors = F,sep="\t", quote="",comment.char="")
aberration_events <- aberration_events %>% dplyr::select(ABERRATION, CODE, TOTAL_CASES) %>% dplyr::rename(CANCERCODE = CODE) %>% dplyr::distinct()
aberration_events <- as.data.frame(dplyr::group_by(aberration_events, ABERRATION, CANCERCODE) %>% dplyr::summarise(N_CASES = sum(TOTAL_CASES)))
aberration_events$CANCERCODE <- as.character(aberration_events$CANCERCODE)
aberration_events <- dplyr::left_join(aberration_events, disease_codes, by=c("CANCERCODE"))
aberration_events <- dplyr::left_join(aberration_events, molclingene, by=c("ABERRATION"))
aberration_events <- dplyr::arrange(aberration_events, ABERRATION, desc(N_CASES))

aberrations_slim <- as.data.frame(aberration_events %>% dplyr::group_by(ABERRATION, GENE) %>% dplyr::summarise(MITELMAN_DB_CANCERTYPE = paste(unique(CANCERTYPE_ABBR), collapse=","), MITELMAN_PMID = paste(unique(PMID), collapse=",")))
mitelman_aberrations <- as.data.frame(aberrations_slim %>% dplyr::group_by(ABERRATION, MITELMAN_DB_CANCERTYPE) %>% dplyr::summarise(MITELMAN_DB_PMID = paste(unique(MITELMAN_PMID), collapse=","), MITELMAN_DB_GENE = paste(GENE, collapse=",")))

write.table(mitelman_aberrations, file=paste0("data-raw/mitelman/mitelman_aberrations.",mitelman_db_release,".tsv"),col.names = T, quote = F, row.names = F,sep="\t")
write.table(disease_codes, file=paste0("data-raw/mitelman/mitelman_disease_codes.",mitelman_db_release,".tsv"),col.names = T, quote = F, row.names = F,sep="\t")


#system(paste0('ln -s -f output/mitelman_aberrations.',mitelman_db_release,'.tsv mitelman_aberrations.tsv'))
#system(paste0('ln -s -f output/mitelman_disease_codes.',mitelman_db_release,'.tsv mitelman_disease_codes.tsv'))

devtools::use_data(mitelman_aberrations, overwrite = T)
devtools::use_data(mitelman_disease_codes, overwrite = T)




get_gencode_ensembl_transcripts <- function(basedir = "/Users/sigven/research/structVaR/data-raw/", build = "grch37", append_regulatory_region = FALSE, gencode_version = "28", gene_info = NULL){
  if(build == 'grch37'){
    gencode_version <- '19'
  }
  rlogging::message(paste0("Retrieving GENCODE transcripts - version ", gencode_version,", build ",build))
  if(build == 'grch37'){
    if(!file.exists(paste0(basedir,"/gencode/",build,"/gencode.",build,".annotation.gtf.gz"))){
      rlogging::message(paste0("Downloading ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"))
      download.file("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz", destfile = paste0(basedir,"/gencode/",build,"/gencode.",build,".annotation.gtf.gz"),quiet = T)
    }
  }
  if(build == 'grch38'){
    if(!file.exists(paste0(basedir,"/gencode/",build,"/gencode.",build,".annotation.gtf.gz"))){
      rlogging::message(paste0("Downloading ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz"))
      download.file("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz", destfile = paste0(basedir,"/gencode/",build,"/gencode.",build,".annotation.gtf.gz"),quiet = T)
    }
  }
  gencode_gtf <- read.table(gzfile(paste0(basedir,"/gencode/",build,"/gencode.",build,".annotation.gtf.gz")), sep="\t", stringsAsFactors = F,na.strings=".", skip=5, comment.char = "#", quote = "")
  gencode_gtf <- dplyr::filter(gencode_gtf, V3 == 'transcript')

  gencode_transcript_start <- gencode_gtf$V4
  gencode_transcript_end <- gencode_gtf$V5
  gencode_transcript_strand <- gencode_gtf$V7
  gencode_chrom <- gencode_gtf$V1
  gencode_gtf_annotation <- stringr::str_replace_all(gencode_gtf$V9,'"','')
  gencode_ensembl_gene_id <- stringr::str_match(gencode_gtf_annotation,"gene_id ENSG[0-9]{1,}")[,1]
  gencode_ensembl_gene_id <- stringr::str_replace(gencode_ensembl_gene_id,"gene_id ","")
  gencode_ensembl_trans_id <- stringr::str_match(gencode_gtf_annotation,"transcript_id ENST[0-9]{1,}\\.[0-9]{1,}")[,1]
  gencode_ensembl_trans_id <- stringr::str_replace(gencode_ensembl_trans_id,"transcript_id ","")
  gencode_ensembl_trans_id_original <- gencode_ensembl_trans_id
  gencode_ensembl_trans_id <- stringr::str_replace(gencode_ensembl_trans_id,"\\.[0-9]{1,}$","")
  gencode_gene_type <- stringr::str_match(gencode_gtf_annotation,"gene_type \\S+")[,1]
  gencode_gene_type <- stringr::str_replace_all(gencode_gene_type,"gene_type |;$","")
  gencode_trans_type <- stringr::str_match(gencode_gtf_annotation,"transcript_type \\S+")[,1]
  gencode_trans_type <- stringr::str_replace_all(gencode_trans_type,"transcript_type |;$","")
  gencode_tag <- stringr::str_match(gencode_gtf_annotation,"tag \\S+")[,1]
  gencode_tag <- stringr::str_replace_all(gencode_tag,"tag |;$","")
  gencode_symbol <- stringr::str_match(gencode_gtf_annotation,"gene_name \\S+")[,1]
  gencode_symbol <- stringr::str_replace_all(gencode_symbol,"gene_name |;$","")
  gencode <- data.frame('chrom' = gencode_chrom, 'start' = gencode_transcript_start, 'end' = gencode_transcript_end, 'strand' = gencode_transcript_strand, 'ensembl_gene_id' = gencode_ensembl_gene_id, 'gencode_gene_biotype' = gencode_gene_type, 'gencode_tag' = gencode_tag, 'gencode_transcript_type' = gencode_trans_type, 'ensembl_transcript_id' = gencode_ensembl_trans_id, 'ensembl_transcript_id_full' = gencode_ensembl_trans_id_original, stringsAsFactors = F)
  gencode <- gencode %>% dplyr::distinct()
  rlogging::message(paste0("A total of ",nrow(gencode)," transcripts parsed"))

  gencode_symbols <- data.frame('ensembl_gene_id' = gencode_ensembl_gene_id, 'symbol' = gencode_symbol, stringsAsFactors = F) %>% dplyr::distinct()

  ## include regulatory region (for VEP annotation)
  if(append_regulatory_region == T){
    rlogging::message("Parameter 'append_regulatory_region' is TRUE: expanding transcript start/end with 5kb (for VEP consequence compliance)")
    chromosome_lengths <- data.frame('chrom' = head(names(seqlengths(BSgenome.Hsapiens.UCSC.hg19)),24), 'chrom_length' = head(seqlengths(BSgenome.Hsapiens.UCSC.hg19),24), stringsAsFactors = F, row.names = NULL)
    if(build == 'grch38'){
      chromosome_lengths <- data.frame('chrom' = head(names(seqlengths(BSgenome.Hsapiens.UCSC.hg38)),24), 'chrom_length' = head(seqlengths(BSgenome.Hsapiens.UCSC.hg38),24), stringsAsFactors = F, row.names = NULL)
    }
    gencode <- dplyr::left_join(gencode, chromosome_lengths,by=c("chrom"))
    gencode <- as.data.frame(gencode %>% dplyr::rowwise() %>% dplyr::mutate(start = max(1,start - 5000)))
    gencode <- as.data.frame(gencode %>% dplyr::rowwise() %>% dplyr::mutate(end = min(chrom_length, end + 5000)))
    gencode <- dplyr::select(gencode, -chrom_length)
  }
  gencode$gencode_release <- gencode_version

  ensembl_mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL")
  if(build == 'grch37'){
    ensembl_mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "http://grch37.ensembl.org")
  }
  ensembl_genes <- biomaRt::useDataset("hsapiens_gene_ensembl", mart=ensembl_mart)
  queryAttributes <- c('ensembl_gene_id','gene_biotype','ensembl_transcript_id','refseq_mrna')
  ensembl_genes_xref1 <- biomaRt::getBM(attributes = queryAttributes, mart=ensembl_genes)
  ensembl_genes_xref1[ensembl_genes_xref1$refseq_mrna == '',]$refseq_mrna <- NA

  queryAttributes <- c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol','description')
  ensembl_genes_xref2 <- biomaRt::getBM(attributes = queryAttributes, mart=ensembl_genes)
  ensembl_genes_xref2[ensembl_genes_xref2$hgnc_symbol == '',]$hgnc_symbol <- NA
  ensembl_genes_xref2[ensembl_genes_xref2$description == '',]$description <- NA
  ensembl_genes_xref2$description <- stringr::str_replace(ensembl_genes_xref2$description," \\[.+$","")
  ensembl_genes_xref2 <- dplyr::select(ensembl_genes_xref2, -ensembl_transcript_id) %>% dplyr::rename(symbol = hgnc_symbol, name = description) %>% dplyr::distinct()

  ensembl_genes_xref <- dplyr::left_join(ensembl_genes_xref1, ensembl_genes_xref2, by=c("ensembl_gene_id"))
  gencode <- dplyr::left_join(gencode, dplyr::select(ensembl_genes_xref, ensembl_gene_id, ensembl_transcript_id, name, refseq_mrna), by=c("ensembl_gene_id","ensembl_transcript_id")) %>% dplyr::distinct()
  gencode <- dplyr::left_join(gencode, dplyr::select(gene_info, entrezgene, symbol, ensembl_gene_id), by=c("ensembl_gene_id"))

  gencode_not_missing_name <- gencode %>% dplyr::filter(!is.na(symbol) & !is.na(entrezgene) & !is.na(name))
  gencode_missing_symbol_entrez <- gencode %>% dplyr::filter(is.na(symbol) & is.na(entrezgene) & !is.na(name))
  gencode_missing_symbol_entrez <- dplyr::left_join(dplyr::select(gencode_missing_symbol_entrez,-symbol), gencode_symbols,by=c("ensembl_gene_id"))
  gencode_missing_name <- gencode %>% dplyr::filter(!is.na(symbol) & !is.na(entrezgene) & is.na(name))
  gencode_missing_symbol_entrez <- dplyr::left_join(dplyr::select(gencode_missing_symbol_entrez,-entrezgene),dplyr::select(gene_info,symbol,entrezgene),by=c("symbol"))

  gencode_missing_name <- dplyr::left_join(dplyr::select(gencode_missing_name, -name), dplyr::select(gene_info, -gene_biotype),by=c("ensembl_gene_id","entrezgene","symbol"))
  gencode <- rbind(gencode_not_missing_name, gencode_missing_symbol_entrez, gencode_missing_name)

  gencode <- gencode %>% dplyr::filter(!is.na(ensembl_transcript_id)) %>% dplyr::distinct()
  rlogging::message(paste0("A total of ",nrow(gencode)," valid transcripts remaining"))

  return(gencode)

}

get_gene_info_ncbi <- function(basedir = "/Users/sigven/research/structVaR/data-raw/", build = "grch37"){
  rlogging::message("Retrieving gene_info from NCBI/Entrez")
  if(!file.exists(paste0(basedir,'/ncbi_gene/',build,'/Homo_sapiens.gene_info.gz'))){
    download.file("ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz", destfile=paste0(basedir,'ncbi_gene/',build,'/Homo_sapiens.gene_info.gz'), quiet=T)
  }
  gene_info <- read.table(gzfile(paste0(basedir,'/ncbi_gene/',build,'/Homo_sapiens.gene_info.gz')), sep="\t", stringsAsFactors = F,na.strings="-", skip=1, comment.char = "#", quote = "", fill=T)
  gene_info <- gene_info %>% dplyr::filter(V1 == 9606)
  gene_info <- dplyr::select(gene_info, c(V2,V3,V9,V6,V10))
  tmp <- as.data.frame(stringr::str_match(gene_info$V6,"Ensembl:ENSG[0-9]{1,}"))
  gene_info$ensembl_gene_id <- stringr::str_replace(tmp$V1,"Ensembl:","")
  gene_info <- dplyr::select(gene_info, -V6)
  gene_info <- dplyr::rename(gene_info, entrezgene = V2, symbol = V3, name = V9, gene_biotype = V10)
  if(nrow(gene_info[!is.na(gene_info$gene_biotype) & gene_info$gene_biotype == 'protein-coding',])>0){
    gene_info[!is.na(gene_info$gene_biotype) & gene_info$gene_biotype == 'protein-coding',]$gene_biotype <- 'protein_coding'
  }
  return(gene_info)

}


get_tsgene_data <- function(basedir = '/Users/sigven/research/structVaR/data-raw/'){
  rlogging::message("Retrieving known proto-oncogenes/tumor suppressor genes from TSGene2.0")
  download.file("https://bioinfo.uth.edu/TSGene/Human_TSGs.txt", destfile=paste0(basedir,'/tsgene/Human_TSGs.txt'), quiet=T)
  tsgene <- read.table(paste0(basedir,'/tsgene/Human_TSGs.txt'), sep="\t", header = T, stringsAsFactors = F,na.strings="-", comment.char = "#", quote = "", fill=T)
  tsgene <- dplyr::select(tsgene, GeneID) %>% dplyr::rename(entrezgene = GeneID)
  tsgene$tsgene <- 'TRUE'
  rlogging::message("A total of ",nrow(tsgene)," tumor suppressor genes were parsed")
  download.file("https://bioinfo.uth.edu/TSGene/oncogene.txt", destfile=paste0(basedir,'/tsgene/Human_ONCGs.txt'), quiet=T)
  ts_oncogene <- read.table(paste0(basedir,'/tsgene/Human_ONCGs.txt'), sep="\t", header = T, stringsAsFactors = F,na.strings="-", comment.char = "#", quote = "", fill=T)
  ts_oncogene <- dplyr::select(ts_oncogene, GeneID) %>% dplyr::rename(entrezgene = GeneID)
  ts_oncogene$ts_oncogene <- 'TRUE'
  rlogging::message("A total of ",nrow(ts_oncogene)," tumor suppressor genes were parsed")

  tsgene_full <- dplyr::full_join(tsgene, ts_oncogene, by="entrezgene")
  n_onc_ts <- dplyr::filter(tsgene_full, ts_oncogene == T & tsgene == T)
  rlogging::message("A total of ",nrow(n_onc_ts)," genes were annotated with dual roles as tumor suppressor genes and oncogenes")

  return(tsgene_full)

}

for(build in c('grch37','grch38')){
  ## get NCBI gene info (gene symbols, names, and entrez gene identifiers)
  gene_info <- get_gene_info_ncbi(basedir = '/Users/sigven/research/structVaR/data-raw/', build = build)

  ## get GENCODE transcripts (chromosomal start/ends, biotypes etc.)
  gencode <- get_gencode_ensembl_transcripts(basedir = '/Users/sigven/research/structVaR/data-raw/', build = build, append_regulatory_region = T, gene_info = gene_info)

  ## Known oncogenes/tumor suppressor genes (TSgene 2.0)
  tsgene <- get_tsgene_data(basedir = '/Users/sigven/research/structVaR/data-raw/')
  gencode <- dplyr::left_join(gencode,tsgene,by=c("entrezgene"))
  gencode <- dplyr::rename(gencode, transcript_start = start, transcript_end = end)
  gencode <- dplyr::filter(gencode, chrom != 'chrM')

  seqinfo <- GenomeInfoDb::Seqinfo(seqnames = GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg38)), seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg38)), genome = 'hg38')
  if(build == 'grch37'){
    seqinfo <- GenomeInfoDb::Seqinfo(seqnames = GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19)), seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19)), genome = 'hg19')
  }
  gencode <- dplyr::select(gencode, -strand)
  gencode <- dplyr::filter(gencode, is.na(gencode_tag) | gencode_tag == 'basic') %>% dplyr::filter(gencode_transcript_type != 'misc_RNA')

  gencode_gr <- GenomicRanges::makeGRangesFromDataFrame(gencode, keep.extra.columns = T, seqinfo = seqinfo, seqnames.field = 'chrom', start.field = 'transcript_start', end.field = 'transcript_end', ignore.strand = T, starts.in.df.are.0based = T)

  if(build == 'grch37'){
    gencode_grch37 <- gencode_gr
    devtools::use_data(gencode_grch37, overwrite = T)
  }
  if(build == 'grch38'){
    gencode_grch38 <- gencode_gr
    devtools::use_data(gencode_grch38, overwrite = T)
  }
}
