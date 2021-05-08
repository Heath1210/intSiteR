#' @rdname intTrim
#' @title Reads filtration and trimming
#'
#' @description Filter reads matching primer and linker, dump others.
#'    Adjust reads to same direction that linker is on the left side.
#'    Trim off linker and primer sequence followed by adding UMI, barcode
#'    to reads ID which has been simplied but informative enough to
#'    distinguish each reads. Save fasta files for further analysis.
#'
#' @param path2fq1 the file path to fq1 file.
#' @param path2fq2 the file path to fq2 file.
#' @param path2mg the file path to merged fastq file.
#' @param outdir the file fold where output files put in
#' @param LTR the sequence of vector closest to integrated host genome.
#'     Recommended 18~28bp. Default value is for lentiviral vector.
#'     If you use other vector, change to corresponding sequence.
#' @param linker the sequence of adaptor linker at the end, near the genome part.
#'     Default linker is from INSPIIRED pipeline.
#'     If you use other linker, change to corresponding sequence.
#' @param avoidseq1 Reads coontaining this sequence 1 will bring false positive
#'     integration site results.
#'     Default value is vector sequence close to 5'LTR tail.
#'     Change to NULL if you don't need it.
#' @param avoidseq2 Reads coontaining this sequence 2 will bring false positive
#'     integration site results
#'     Default value is vector sequence close to 3'LTR tail.
#'     Change to NULL if you don't need it.
#'
#' @importFrom Biostrings vcountPattern
#' @importFrom Biostrings vmatchPattern
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings DNAString
#' @importFrom ShortRead readFastq
#' @importFrom ShortRead writeFasta
#' @importFrom Biostrings width
#' @importFrom dplyr `%>%`
#' @importFrom data.table fwrite
#'
#' @export
#'

intTrim <- function(path2fq1,
                    path2fq2,
                    path2mg,
                    outdir = NULL,
                    LTR = 'AGTCAGTGTGGAAAATCTCTAGCA',
                    linker = 'CTCCGCTTAAGGGACT',
                    avoidseq1 = 'GTGGCGCCCGAACAGGGACTTGAAAGCGAAAGGGAAACCAGAGGAGCTCT',
                    avoidseq2 = 'GTAGTAGTTCATGTCATCTTATTATTCAGTATTTATAACT'
                    ){
  # check validity of input arguments
  stopifnot('Error: input file not exists, please check your file path or file name!' =
              file.exists(path2fq1) & file.exists(path2fq2) & file.exists(path2mg))

  # output
  if(is.null(outdir)){
    outdir = paste0(dirname(path2fq1),'/fasta')
    if(!dir.exists(outdir)) dir.create(outdir)
  }

  # convert charactor to DNAString
  LTR = DNAString(LTR)
  linker = DNAString(linker)

  # loading fastq
  mg = readFastq(path2mg)
  fq1 = readFastq(path2fq1)
  fq2 = readFastq(path2fq2)
  # QC total reads
  total_reads = length(fq1)+length(mg)
  merge_percentage = length(mg)/total_reads

  # find out the merged reads that can be matched with reverse and complement LTR sequence
  mg_rc_idx = which(vcountPattern(reverseComplement(LTR),
                                       substr(mg@sread,width(mg)-40,width(mg)))==1)
  # replace reverse and complement reads so that all reads have same orientation
  mg@sread[mg_rc_idx] = reverseComplement(mg@sread[mg_rc_idx])

  # positively filter reads that can be matched with LTR and linker
  mg_match_LTR = which(vcountPattern(LTR,substr(mg@sread,1,40))==1)

  mg_match_linker = which(vcountPattern(reverseComplement(linker),
                                  substr(mg@sread,width(mg)-60,width(mg)-20))==1)

  mg_match_idx = intersect(mg_match_LTR,mg_match_linker)

  # negatively filter reads that can be matched with vector sequences that close to LTR
  mg_match_avs1 = ifelse(is.null(avoidseq1),rep(0,length(mg)),
                         vcountPattern(avoidseq1,mg@sread,max.mismatch=3))

  mg_match_avs2 = ifelse(is.null(avoidseq2),rep(0,length(mg)),
                         vcountPattern(avoidseq2,mg@sread,max.mismatch=5))

  mg_match_avs = which(mg_match_avs1 + mg_match_avs2 > 0)

  # drop unwanted reads
  mg = mg[setdiff(mg_match_idx,mg_match_avs)]

  if(length(mg)!=0) {
    # get the genome sequence start location in each reads
    mg_start = unlist(vmatchPattern(LTR,substr(mg@sread,1,40))@ends)+1

    # get the genome sequence end location in each reads
    mg_end = unlist(vmatchPattern(reverseComplement(linker),
                                  substr(mg@sread,width(mg)-60,width(mg)-20))@ends)+width(mg)-length(linker)-61

    # the start and end locations are calculated in a vectorization way, normally have same length
    # if not, maybe LTR or linker matched twice in one reads, so the match boundaries should be narrowed
    stopifnot('Please adjust your strictness of matching!'=
                length(mg_start)==length(mg_end))

    # trim off LTR and linker sequences in reads, keep only genome sequences
    mg_trim = substr(mg@sread,mg_start,mg_end)

    # name trimed reads the shortest unique strings of reads IDs appending with UMI sequencs
    # normally reads IDs have two patterns: with and without barcode
    # the length of reads ID with barcode is larger than 60
    mg_id = data.frame(matrix(unlist(strsplit(as.character(mg@id),' ')),ncol = 3,byrow = T),
                       stringsAsFactors = FALSE)[,1]

    names(mg_trim) <- paste0(substr(mg_id,24,nchar(mg_id)),':',
                             substr(mg@sread,mg_end+length(linker)+1,
                                    mg_end+length(linker)+18) %>%
                               DNAStringSet() %>%
                               reverseComplement())

    mg_trim = mg_trim[which(nchar(mg_trim)>30)]
  }


  # find out the fq1 reads that can be matched with reverse and complement LTR sequence
  fq_rc_idx = which(vcountPattern(LTR,substr(fq2@sread,1,40))==1)

  # exchange reads so that all reads are in same orientation
  fq1_temp = fq1
  fq1@sread[fq_rc_idx] = fq2@sread[fq_rc_idx]
  fq1@quality@quality[fq_rc_idx] = fq2@quality@quality[fq_rc_idx]

  fq2@sread[fq_rc_idx] = fq1_temp@sread[fq_rc_idx]
  fq2@quality@quality[fq_rc_idx] = fq1_temp@quality@quality[fq_rc_idx]

  # the index of PE reads that can be matched with LTR and linker
  fq_match_LTR = which(vcountPattern(LTR,substr(fq1@sread,1,40))==1)

  fq_match_linker = which(vcountPattern(linker,substr(fq2@sread,20,60))==1)

  fq_match_idx = intersect(fq_match_LTR,fq_match_linker)

  reads_matching_primers = length(mg_match_idx)+length(fq_match_idx)

  # the index of PE reads that can be matched with vector sequences that close to LTR
  fq_match_avs1 = ifelse(is.null(avoidseq1),rep(0,length(fq1)),
                         vcountPattern(avoidseq1,fq1@sread,max.mismatch=3))

  fq_match_avs2 = ifelse(is.null(avoidseq2),rep(0,length(fq1)),
                         vcountPattern(avoidseq2,fq1@sread,max.mismatch=5))

  fq_match_avs = which(fq_match_avs1 + fq_match_avs2 > 0)

  # positive filter and negative filter
  fq1 = fq1[setdiff(fq_match_idx,fq_match_avs)]
  fq2 = fq2[setdiff(fq_match_idx,fq_match_avs)]

  if(length(fq1)!=0){
    # get the genome sequence start and end location in fq1, fq2 reads
    fq_start = unlist(vmatchPattern(LTR,substr(fq1@sread,1,40))@ends)+1
    fq_end = unlist(vmatchPattern(linker,substr(fq2@sread,20,60))@ends)+20

    stopifnot('Please adjust your strictness of matching!' =
                length(fq_start) == length(fq_end))

    # trim off LTR and linker sequences ,get genome sequences and UMI
    fq1_trim = substr(fq1@sread,fq_start,width(fq1))
    fq2_trim = substr(fq2@sread,fq_end,width(fq2))

    # name trimed reads the shortest unique strings of reads IDs appending with UMI sequencs
    names(fq1_trim) = names(fq2_trim) =
      paste0(substr(fq1@id,24,width(fq1@id)-ifelse(width(fq1@id)>60,24,15)),':',
             substr(fq2@sread,fq_end-length(linker)-18,fq_end-length(linker)-1))

    fq_short_idx = which(nchar(fq1_trim) < 20 | nchar(fq2_trim) < 20)

    if(length(fq_short_idx) > 0){
      fq1_trim = fq1_trim[-fq_short_idx]
      fq2_trim = fq2_trim[-fq_short_idx]
    }
  }

  reads_usable_for_sites_detection =
    ifelse(length(fq1)==0,0,length(fq1_trim))+
    ifelse(length(mg)==0,0,length(mg_trim))

  if(length(mg) !=0 ){
    writeFasta(DNAStringSet(mg_trim),
               paste0(outdir,'/',strsplit(basename(path2mg), "\\.f")[[1]][1],'.fa'))
  }

  if(length(fq1) != 0){
    writeFasta(DNAStringSet(fq1_trim),
               paste0(outdir,'/',strsplit(basename(path2fq1), "\\.f")[[1]][1],'.fa'))
    writeFasta(DNAStringSet(fq2_trim),
               paste0(outdir,'/',strsplit(basename(path2fq2), "\\.f")[[1]][1],'.fa'))
  }

  log_info = list(
    total_reads = total_reads,
    merge_percentage = merge_percentage,
    reads_matching_primers = reads_matching_primers,
    reads_usable_for_sites_detection = reads_usable_for_sites_detection)

  fwrite(log_info,
         paste0(outdir,'/',strsplit(basename(path2fq1), "\\.f")[[1]][1],'_P1.log'),
         sep = '\t' )

  log_info

}
