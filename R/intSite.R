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

  myfiles = dir(paste0('./',basename(input)) ,'q.gz$',full.names = TRUE)

  fileNameCore = unique(unlist(strsplit(basename(myfiles),'\\_R[1-2]'))[rep(c(TRUE,FALSE),length(myfiles)/2)])

  num_threads = parallel::detectCores()



  # raw data QC and merge reads
  if(!dir.exists('./1fastp')) dir.create('./1fastp')

  fq1_out = paste0('./1fastp/',fileNameCore,'_R1.fq')
  fq2_out = paste0('./1fastp/',fileNameCore,'_R2.fq')
  mg_out = paste0('./1fastp/',fileNameCore,'_merge.fq')

  for (i in 1:length(fileNameCore)) {
    mergeReads(fq1 = myfiles[2*i-1],
               fq2 = myfiles[2*i],
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

  for (i in 1:length(fileNameCore)) {
    alignBowtie2(fa1 = paste0('./2fasta/',strsplit(basename(fq1_out[i]), "\\.")[[1]][1],'.fa'),
                 fa2 = paste0('./2fasta/',strsplit(basename(fq2_out[i]), "\\.")[[1]][1],'.fa'),
                 outdir = './3sam',
                 bowtie2 = 'bowtie2',
                 ref = ref,
                 threads = num_threads)
    alignBowtie2(fa1 = paste0('./2fasta/',strsplit(basename(mg_out[i]), "\\.")[[1]][1],'.fa'),
                 outdir = './3sam',
                 bowtie2 = 'bowtie2',
                 ref = ref,
                 threads = num_threads)
  }


  # sam to bam
  if(!dir.exists('./4bam')) dir.create('./4bam')

  for (i in 1:length(fileNameCore)) {
    sam2bam(sam = paste0('./3sam/',strsplit(basename(fq1_out[i]), "\\.")[[1]][1],'.sam'),
            outdir = './4bam',
            samtools = 'samtools')

    sam2bam(sam = paste0('./3sam/',strsplit(basename(mg_out[i]), "\\.")[[1]][1],'.sam'),
            outdir = './4bam',
            samtools = 'samtools')
  }

  # bam to bed

  if(!dir.exists('./5bed')) dir.create('./5bed')

  for (i in 1:length(fileNameCore)) {
    bam2bed(bam = paste0('./4bam/',strsplit(basename(fq1_out[i]), "\\.")[[1]][1],'.bam'),
            outdir = './5bed',
            bedtools = 'bedtools')

    bam2bed(bam = paste0('./4bam/',strsplit(basename(mg_out[i]), "\\.")[[1]][1],'.bam'),
            outdir = './5bed',
            bedtools = 'bedtools')
  }

  # process bed file

  if(!dir.exists('./6insite')) dir.create('./6insite')

  for (i in 1:length(fileNameCore)) {
    intBed(mgBed = paste0('./5bed/',strsplit(basename(mg_out[i]), "\\.")[[1]][1],'.bed'),
           fqBed = paste0('./5bed/',strsplit(basename(fq1_out[i]), "\\.")[[1]][1],'.bed'),
           outdir = './6insite')
  }

  file.copy(from = dir('./2fasta','.log',full.names = TRUE),
            to   = "./6insite")
}
