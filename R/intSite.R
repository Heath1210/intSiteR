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
#' @param mode output mode, whether by file type or by sample,
#'     value can be 'filetype' or 'sample'
#' @param theta a parameter controlling the looseness of filter by duplicate degree
#' @param fThreshold the minimum duplicate degree allowed for further analysis
#' @param monoSite a threshold that insiteCode with frequency above it
#'     will not be filtered by breaking point
#' @param collapse whether neighboring sites collapsed to one.
#' @param mixBarcode whether samples with different barcodes mixed in one library.
#'
#' @importFrom stringr str_match
#' @importFrom stringr str_replace
#' @importFrom fs path
#' @importFrom fs path_dir
#' @importFrom fs path_file
#' @importFrom fs dir_ls
#' @importFrom fs file_move
#' @importFrom fs dir_create
#' @importFrom parallel detectCores
#' @importFrom parallel mcmapply
#' @importFrom dplyr `%>%`
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
                    bedtools = 'bedtools',
                    mode = 'filetype',
                    theta = 100,
                    fThreshold = 5,
                    monoSite = 50,
                    collapse = 7,
                    mixBarcode = FALSE){

  num_threads = parallel::detectCores()

  raw_fq1_files = dir_ls(input,regexp='.+R1.+gz')
  raw_fq2_files = dir_ls(input,regexp='.+R2.+gz')
  raw_r1 = str_match(raw_fq1_files,'(.+)_R1')[,2] %>% path_file()
  raw_r2 = str_match(raw_fq1_files,'(.+)_R1')[,2] %>% path_file()

  stopifnot('Please check your fastq pairs!'= raw_r1 == raw_r2)

  if(mode != 'sample'){
    if(mode != 'filetype') {
      print('You did not choose a right output mode, it will run as filetype mode.')}

    # raw data QC and merge reads
    fp_out <- path(dirname(input),'1fastp')
    if(!dir.exists(fp_out)) dir.create(fp_out)

    mapply(mergeReads,
           fq1 = raw_fq1_files,
           fq2 = raw_fq2_files,
           fq1_out = paste0(fp_out,'/',raw_r1,'_R1.fq'),
           fq2_out = paste0(fp_out,'/',raw_r1,'_R2.fq'),
           mg_out = paste0(fp_out,'/',raw_r1,'_merge.fq'),
           threads = num_threads)

    # trim reads
    fa_out <- path(dirname(input),'2fasta')
    if(!dir.exists(fa_out)) dir.create(fa_out)

    mcmapply(intTrim,
             path2fq1 = paste0(fp_out,'/',raw_r1,'_R1.fq'),
             path2fq2 = paste0(fp_out,'/',raw_r1,'_R2.fq'),
             path2mg = paste0(fp_out,'/',raw_r1,'_merge.fq'),
             outdir = fa_out,
             mc.cores = num_threads)

    # align
    sam_out <- path(dirname(input),'3sam')
    if(!dir.exists(sam_out)) dir.create(sam_out)

    mapply(alignBowtie2,
           fa1 = dir(fa_out,'_R1.fa$',full.names = TRUE),
           fa2 = dir(fa_out,'_R2.fa$',full.names = TRUE),
           outdir = sam_out,
           bowtie2 = bowtie2,
           ref = ref,
           threads = num_threads)

    mapply(alignBowtie2,
           fa1 = dir(fa_out,'_merge.fa$',full.names = TRUE),
           outdir = sam_out,
           bowtie2 = bowtie2,
           ref = ref,
           threads = num_threads)


    # sam to bam
    bam_out <- path(dirname(input),'4bam')
    if(!dir.exists(bam_out)) dir.create(bam_out)

    mcmapply(sam2bam,
             sam = dir(sam_out,'sam$',full.names = TRUE),
             outdir = bam_out,
             samtools = samtools,
             mc.cores = num_threads)

    # bam to bed
    bed_out <- path(dirname(input),'5bed')
    if(!dir.exists(bed_out)) dir.create(bed_out)

    mcmapply(bam2bed,
             bam = dir(bam_out ,'bam$',full.names = TRUE),
             outdir = bed_out,
             bedtools = bedtools,
             mc.cores = num_threads)

    # process bed file
    intsite_out <- path(dirname(input),'6intsite')
    if(!dir.exists(intsite_out)) dir.create(intsite_out)

    mcmapply(intBed,
             mgBed = dir(bed_out,'merge.bed$',full.names = TRUE),
             fqBed = dir(bed_out,'R1.bed$',full.names = TRUE),
             outdir = intsite_out,
             theta = theta,
             fThreshold = fThreshold,
             monoSite = monoSite,
             collapse = collapse,
             mixBarcode = mixBarcode,
             mc.cores = num_threads)

    file.copy(from = dir(fa_out,'.log',full.names = TRUE),
              to = intsite_out)

    file.rename(from = dir(intsite_out,full.names = TRUE),
                to = str_replace(dir(intsite_out,full.names = TRUE),'_R1_','_'))
  }

  if(mode == 'sample'){
    dir_create(path(input,raw_r1))
    dir_create(path(input,raw_r1,'output'))
    dir_create(path(input,raw_r1,'result'))
    file_move(raw_fq1_files,path(input,raw_r1))
    file_move(raw_fq2_files,path(input,raw_r1))

    # merge fastq reads
    my_pwd <- path(input,raw_r1)

    mapply(mergeReads,
           fq1 = dir_ls(my_pwd, regexp = "_R1"),
           fq2 = dir_ls(my_pwd, regexp = "_R2"),
           fq1_out = paste0(my_pwd,'/output/',raw_r1,'_R1.fq'),
           fq2_out = paste0(my_pwd,'/output/',raw_r1,'_R2.fq'),
           mg_out = paste0(my_pwd,'/output/',raw_r1,'_merge.fq'),
           threads = num_threads)

    # trim the LTR and linker
    mcmapply(intTrim,
             path2fq1 = paste0(my_pwd, "/output/", raw_r1, "_R1.fq"),
             path2fq2 = paste0(my_pwd, "/output/", raw_r1, "_R2.fq"),
             path2mg = paste0(my_pwd, "/output/", raw_r1, "_merge.fq"),
             outdir = paste0(my_pwd,'/output'),
             mc.cores = num_threads)

    # align
    mapply(alignBowtie2,
           fa1 = dir_ls(input, regexp = "R1.fa$",recurse = TRUE),
           fa2 = dir_ls(input, regexp = "R2.fa$",recurse = TRUE),
           outdir = path_dir(dir_ls(input, regexp = "R1.fa$",recurse = TRUE)),
           bowtie2 = 'bowtie2',
           ref = ref,
           threads = num_threads)

    mapply(alignBowtie2,
           fa1 = dir_ls(input, regexp = "merge.fa$",recurse = TRUE),
           outdir = path_dir(dir_ls(input, regexp = "merge.fa$",recurse = TRUE)),
           bowtie2 = 'bowtie2',
           ref = ref,
           threads = num_threads)

    mcmapply(sam2bam,
             sam = dir_ls(input, regexp = ".sam",recurse = TRUE),
             outdir = path_dir(dir_ls(input, regexp =  ".sam",recurse = TRUE)),
             samtools = 'samtools',
             mc.cores = num_threads)

    mcmapply(bam2bed,
             bam = dir_ls(input, regexp = ".bam",recurse = TRUE),
             outdir = path_dir(dir_ls(input, regexp =  ".bam",recurse = TRUE)),
             bedtools = 'bedtools',
             mc.cores = num_threads)

    # parse bed files to extract integration sites
    mcmapply(intBed,
             mgBed = dir_ls(input, regexp = "merge.bed",recurse = TRUE),
             fqBed = dir_ls(input, regexp = "R1.bed",recurse = TRUE),
             outdir = fs::path(dir_ls(input, regexp = "R1.bed",recurse = TRUE) %>%
                                 path_dir %>%
                                 path_dir,
                               'result'),
             theta = theta,
             fThreshold = fThreshold,
             monoSite = monoSite,
             collapse = collapse,
             mixBarcode = mixBarcode,
             mc.cores = num_threads)

    file_move(path = dir_ls(input, regexp = "P1.log",recurse = TRUE),
              new_path = path(dir_ls(input, regexp = "P1.log",recurse = TRUE) %>%
                                    path_dir %>%
                                    path_dir,
                                  'result'))

    file.rename(from = dir(path(input,raw_r1,'result'),full.names = TRUE),
                to = str_replace(dir(path(input,raw_r1,'result'),full.names = TRUE),
                                 '_R1_','_'))
  }

 }


