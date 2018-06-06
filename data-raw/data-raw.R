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


