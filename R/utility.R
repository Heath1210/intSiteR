#' @title Read bed file
#'
#' @description Read bed file
#'
#' @param bed a .bed file
#'
#' @importFrom data.table fread
#'
#' @export

readBed <- function(bed){
  col.name = c('chr','start','end','name','quality','strand')
  col.class= c('character','numeric','numeric','character','numeric','character')
  fread(bed,col.names = col.name,colClasses = col.class,header = FALSE)
}


#' @title Find duplicates based on edit distance
#'
#' @description Find duplicates' index based on edit distance
#'
#' @param str_vec a character vector
#' @param len UMI length
#'
#' @export

closeSeqIdx <- function(str_vec,len){
  dup_idx = NULL
  for (i in 2:(len-1)) {
    for (j in (i+1):len) {
      dup_idx = c(dup_idx,which(duplicated(
        paste0(substr(str_vec,1,nchar(str_vec)-j),
               substr(str_vec,nchar(str_vec)-j+2,nchar(str_vec)-i),
               substr(str_vec,nchar(str_vec)-i+2,nchar(str_vec)))
      )))
    }
  }
  unique(dup_idx)
}


#' @title Translate barcode sequences to number
#'
#' @description Translate barcode sequences to number
#'
#' @param barc_vec a vector containing short sequences need to be parsed
#' @param barcode a vector containing the reference barcodes.
#'     Default is c('ATCACG','CGATGT','TTAGGC','TGACCA','ACAGTG',
#'     'GCCAAT','CAGATC','ACTTGA','GATCAG')
#'
#' @importFrom utils adist
#'
#' @export


barcTrans <- function(barc_vec,barcode = def_barcode){
  def_barcode = c('ATCACG','CGATGT','TTAGGC','TGACCA','ACAGTG',
               'GCCAAT','CAGATC','ACTTGA','GATCAG')
  code = rep(0,length(barc_vec))
  for(i in 1:length(barcode)) code[which(adist(barc_vec,barcode[i])<2)] = i
  code
}


#' @title Using bowtie2 to align reads
#'
#' @description Using bowtie2 to align reads
#'
#' @param fa1 a fasta file
#' @param fa2 a fasta file
#' @param outdir the folder of output
#' @param bowtie2 bowtie2 path
#' @param ref reference genome
#' @param threads number of threads used in alignment
#'
#' @export

alignBowtie2 <- function(fa1,
                         fa2 = NULL,
                         outdir,
                         bowtie2 = 'bowtie2',
                         ref,
                         threads = 8){
  # outdir check
  if(is.null(outdir)){
    outdir = dirname(fa1)
  }else stopifnot('Invalid outdir' = dir.exists(outdir))

  if(is.null(fa2)) {
      system(paste0(bowtie2,' -f --end-to-end -t --no-unal -p',threads,' -x ',ref,' -U ',fa1,' -S ',
                    paste0(outdir,'/',strsplit(basename(fa1),"\\.")[[1]][1],'.sam')))
  } else {
      system(paste0(bowtie2,' -f --end-to-end -t --no-unal --no-mixed --no-discordant --dovetail --no-contain --no-overlap  -p',
                    threads,' -x ',ref,' -1 ',fa1,' -2 ',fa2,' -S ',
                    paste0(outdir,'/',strsplit(basename(fa1),"\\.")[[1]][1],'.sam')))
  }
}


#' @title Transform sam to bam
#'
#' @description Using samtools to transform sam to bam
#'
#' @param sam a .sam file
#' @param outdir the folder of output
#' @param samtools samtools path
#'
#' @export

sam2bam <- function(sam,
                    outdir,
                    samtools = 'samtools'){
  # outdir check
  if(is.null(outdir)){
    outdir = dirname(sam)
  }else stopifnot('Invalid outdir' = dir.exists(outdir))

  system(paste0(samtools,' view -Sb ',sam,'> ',
                paste0(outdir,'/',strsplit(basename(sam),"\\.")[[1]][1],'.bam')))
}


#' @title Transform bam to bed
#'
#' @description Using bedtools to transform bam to bed
#'
#' @param bam a .bam file
#' @param outdir the folder of output
#' @param bedtools bedtools path
#'
#' @export

bam2bed <- function(bam,
                    outdir,
                    bedtools = 'bedtools'){
  # outdir check
  if(is.null(outdir)){
    outdir = dirname(bam)
  }else stopifnot('Invalid outdir' = dir.exists(outdir))

  system(paste0(bedtools,' bamtobed -i ',bam,'> ',
                paste0(outdir,'/',strsplit(basename(bam),"\\.")[[1]][1],'.bed')))
}


#' @title Merge paired-end Reads to single reads
#'
#' @description Using fastp to merge paired-end Reads
#'
#' @param fq1 a fastq file
#' @param fq2 a fastq file
#' @param fastp fastp path
#' @param fq1_out filtered fq1 output file
#' @param fq2_out filtered fq2 output file
#' @param mg_out merged reads output file
#' @param threads number of threads used in alignment
#'
#' @export

mergeReads <- function(fq1,
                       fq2,
                       fastp = 'fastp',
                       fq1_out = NULL,
                       fq2_out = NULL,
                       mg_out = NULL,
                       threads = 8){
  if(is.null(fq1_out)){
    fq1_out = paste0(dirname(fq1),'/',strsplit(basename(fq1),"\\.")[[1]][1],'.fq')
  }
  if(is.null(fq2_out)){
    fq2_out = paste0(dirname(fq2),'/',strsplit(basename(fq2),"\\.")[[1]][1],'.fq')
  }
  if(is.null(mg_out)){
    mg_out = paste0(dirname(fq1),'/',strsplit(basename(fq1),"\\.")[[1]][1],'_merge.fq')
  }

  system(paste0(fastp," -g -w ",threads," -i ",fq1," -o ",fq1_out,
                " -I ",fq2," -O ",fq2_out," -m --merged_out ",mg_out))

}
