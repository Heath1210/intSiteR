#' @rdname intSite
#' @title intSite
#'
#' @author Cai Haodong
#'
#' @description Fetch the integration sites from raw fastq files.
#'    This package integrated mergeReads, insTrim, alignBowtie2,
#'    sam2bam, bam2bed and intBed, is used for a complete processing
#'    from raw fastq to integration site coordinates.
#'    Alternatively, you can also analysis your data step by step,
#'    in such mode you can tune lots of parameters to achieve better
#'    performance.
#'
#' @param input a folder containing raw fastq files. Not
#'     supporting file path beginning with dot, as './'!!
#' @param ref genome reference for bowtie2
#' @param bowtie2 bowtie2 path
#' @param samtools samtools path
#' @param bedtools bedtools path
#'
#' @examples
#' if(FALSE){
#'   setwd('~/')
#'
#'   # Not supporting file path beginning with dot, like './'
#'   if(!dir.exists('~/intTest')) dir.create('~/intTest')
#'
#'   file.copy(from = system.file('extdata','raw',package = 'intSiteR'),
#'             to = '~/intTest',
#'             recursive = TRUE)
#'
#'   intSite(input = '~/intTest/raw',
#'           ref = '~/hg38/Homo_sapiens.GRCh38')
#'           # change to your reference folder
#'  }
#' @export


intSite <- function(input,
                    ref,
                    bowtie2 = 'bowtie2',
                    samtools = 'samtools',
                    bedtools = 'bedtools'){
  setwd(dirname(input))
  num_threads = parallel::detectCores()

  rawfiles = dir(paste0('./',basename(input)) ,'q.gz$',full.names = TRUE)
  fileNameCore = unique(unlist(strsplit(basename(rawfiles),'\\_R[1-2]'))[rep(c(TRUE,FALSE),length(rawfiles)/2)])

  # raw data QC and merge reads
  if(!dir.exists('./1fastp')) dir.create('./1fastp')

  fq1_out = paste0('./1fastp/',fileNameCore,'_R1.fq')
  fq2_out = paste0('./1fastp/',fileNameCore,'_R2.fq')
  mg_out = paste0('./1fastp/',fileNameCore,'_merge.fq')

  for (i in 1:length(fileNameCore)) {
    mergeReads(fq1 = rawfiles[2*i-1],
               fq2 = rawfiles[2*i],
               fq1_out = fq1_out[i],
               fq2_out = fq2_out[i],
               mg_out = mg_out[i],
               threads = num_threads
    )
  }


  # trim reads
  if(!dir.exists('./2fasta')) dir.create('./2fasta')

  for (i in 1:length(fileNameCore)) {
    intTrim(fq1_out[i],fq2_out[i],mg_out[i],outdir = './2fasta')
  }


  # align
  if(!dir.exists('./3sam')) dir.create('./3sam')

  fa_merge <- dir(paste0('./2fasta') ,'_merge.fa$',full.names = TRUE)
  fa_r1 <-  dir(paste0('./2fasta') ,'_R1.fa$',full.names = TRUE)
  fa_r1_corename <- unlist(strsplit(basename(fa_r1),'\\_R1'))[rep(c(TRUE,FALSE),length(fa_r1))]

  for (i in 1:length(fa_r1_corename)) {
    alignBowtie2(fa1 = paste0('./2fasta/',fa_r1_corename[i],'_R1.fa'),
                 fa2 = paste0('./2fasta/',fa_r1_corename[i],'_R2.fa'),
                 outdir = './3sam',
                 bowtie2 = 'bowtie2',
                 ref = ref,
                 threads = num_threads)}

  for (i in 1:length(fa_merge)) {
    alignBowtie2(fa1 = fa_merge[i],
                 outdir = './3sam',
                 bowtie2 = 'bowtie2',
                 ref = ref,
                 threads = num_threads)
  }


  # sam to bam
  if(!dir.exists('./4bam')) dir.create('./4bam')

  samfiles <- dir(paste0('./3sam') ,'sam$',full.names = TRUE)

  for (i in 1:length(samfiles)) {
    sam2bam(sam = samfiles[i],
            outdir = './4bam',
            samtools = 'samtools')
  }

  # bam to bed

  if(!dir.exists('./5bed')) dir.create('./5bed')

  bamfiles <- dir(paste0('./4bam') ,'bam$',full.names = TRUE)

  for (i in 1:length(bamfiles)) {
    bam2bed(bam = bamfiles[i],
            outdir = './5bed',
            bedtools = 'bedtools')
  }

  # process bed file

  if(!dir.exists('./6insite')) dir.create('./6insite')

  bed_r1 <- dir(paste0('./5bed') ,'R1.bed$',full.names = TRUE)
  bed_merge <- dir(paste0('./5bed') ,'merge.bed$',full.names = TRUE)
  bed_corename <- unlist(strsplit(basename(bed_r1),'\\_R1'))[rep(c(TRUE,FALSE),length(bed_r1))]


  for (i in 1:length(bed_r1)) {
    intBed(mgBed = bed_merge,
           fqBed = bed_r1,
           outdir = './6insite')
  }

  file.copy(from = dir('./2fasta','.log',full.names = TRUE),
            to   = "./6insite")
}
