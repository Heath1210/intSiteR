#' @rdname intBed
#' @title Fetch integration sites from bed files
#'
#' @description Combine bed files from alignment of merged and unmerged reads.
#'     Deduplication is based on UMI and breaking point as well as the relative
#'     duplicate degree of reads.
#'
#' @param mgBed the file path of bed from alignment of merged reads
#' @param fqBed the file path of bed from alignment of unmerged reads
#' @param outdir the folder of output
#' @param qThreshold the threshold of alignment quality filter.
#'     Default value 20 means keep only > Q20.
#' @param theta a parameter controlling the looseness of filter by duplicate degree
#' @param fThreshold the minimum duplicate degree allowed for further analysis.
#' @param monoSite a threshold that insiteCode with frequency above it
#'     will not be filtered by breaking point
#' @param collapse whether neighboring sites collapsed to one.
#' @param mixBarcode whether different samples with different barcode mixed in
#'     one library.
#'
#' @importFrom dplyr `%>%`
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr summarise
#' @importFrom dplyr desc
#' @importFrom dplyr group_by
#' @importFrom data.table fwrite
#'
#'
#' @export


intBed <- function(mgBed,
                   fqBed,
                   outdir = NULL,
                   qThreshold = 20,
                   theta = 100,
                   fThreshold = 5,
                   monoSite = 50,
                   collapse = 7,
                   mixBarcode = FALSE){
  # prohibit scientific notation
  options(scipen = 200)
  # outdir check
  if(is.null(outdir)){
    outdir = paste0(dirname(mgBed),'/result')
    if(!dir.exists(outdir)) dir.create(outdir)
  }else stopifnot('Invalid outdir!' = dir.exists(outdir))

  mgbed = readBed(mgBed)
  fqbed = readBed(fqBed)

  # mgbed fqbed check
  stopifnot('Bed for fq1&2 seems wrong.' =
              substring(fqbed$name[1],nchar(fqbed$name[1])) == '1')
  stopifnot('Bed for merged fq seems wrong.' =
              substring(mgbed$name[1],nchar(mgbed$name[1])) != '1')

  fqbed_temp = fqbed[rep(c(TRUE,FALSE),nrow(fqbed)/2)]
  fqbed_temp$start = ifelse(fqbed_temp$strand=='+',fqbed_temp$start,fqbed$start[rep(c(FALSE,TRUE),nrow(fqbed)/2)])
  fqbed_temp$end = ifelse(fqbed_temp$strand=='+',fqbed$end[rep(c(FALSE,TRUE),nrow(fqbed)/2)],fqbed_temp$end)
  fqbed_temp$name = substr(fqbed_temp$name,1,nchar(fqbed_temp$name)-2)

  bed = rbind(mgbed,fqbed_temp)
  bed = bed %>% filter(quality> qThreshold)
  bed = bed %>% mutate(
    barcode = barcTrans(substr(name,nchar(name)-5,nchar(name))),
    insiteCode = paste0('B',barcode,':',chr,':',start,':',end,':',strand,':',
                        substr(name,nchar(name)-17,nchar(name)-6)),
    insite = paste0('chr',chr,':',ifelse(strand=='+',start,end),':',strand)
  ) %>% filter(barcode!=0)

  if(nrow(bed)==0){
    fwrite(list(putative_barc = 0,
                unique_insiteCode = 0,
                removed_insiteCode_count_by_UMI_distance = 0,
                filter_threshold_automatically_calculated = 0,
                max_insiteCode = 0),
           paste0(outdir,'/',strsplit(basename(fqBed),"\\.")[[1]][1],'_P2.log'),
           sep='\t')
  }else{
    barcode_table = data.frame(table(bed$barcode))
    bed_serial = list()

    putative_barc = ifelse(mixBarcode,
                           filter(barcode_table,Freq/sum(Freq) > 0.05)$Var1,
                           arrange(barcode_table,desc(Freq))$Var1[1])

    log_info = list()

    for(i in putative_barc){
      bed_serial[[i]] = filter(bed,barcode==i)
      insiteCode_table = data.frame(table(bed_serial[[i]]$insiteCode))
      rownames(insiteCode_table) = insiteCode_table$Var1

      bed_serial[[i]] = bed_serial[[i]] %>% mutate(count = insiteCode_table[insiteCode,2]) %>% arrange(desc(count))

      bed_serial[[i]] = bed_serial[[i]] %>% filter(!duplicated(insiteCode))

      bed_serial[[i]] = bed_serial[[i]][-closeSeqIdx(bed_serial[[i]]$insiteCode,12),]

      unique_insiteCode = length(unique(bed_serial[[i]]$insiteCode))

      removed_insiteCode_count_by_UMI_distance = nrow(insiteCode_table)-unique_insiteCode

      filter_threshold_automatically_calculated = max(fThreshold,
                                                      mean(bed_serial[[i]]$count[1:10])/theta)

      log_info[[paste0(i)]] = list(
        paste0('putative_barc = ',i,'\n',
               'unique_insiteCode = ',unique_insiteCode,'\n',
               'removed_insiteCode_count_by_UMI_distance = ',removed_insiteCode_count_by_UMI_distance,'\n',
               'filter_threshold_automatically_calculated = ',filter_threshold_automatically_calculated,'\n',
               'max_insiteCode = ',bed_serial[[i]]$count[1])
      )

      bed_serial[[i]] = bed_serial[[i]] %>% filter(count > max(fThreshold,
                                                               mean(sort(count,decreasing = TRUE)[1:10])/theta))

      if(nrow(bed_serial[[i]])==0){
        fwrite(log_info[[paste0(i)]],
               paste0(outdir,'/',strsplit(basename(fqBed),"\\.")[[1]][1],'_B',i,'_P2.log'),
               quote = FALSE)
      }else{
        bed_serial[[i]] = bed_serial[[i]] %>% mutate(breakSite = paste0('B',barcode,':',chr,':',end))

        sum_intBreak = summarise(group_by(bed_serial[[i]],insite,breakSite))

        sum_intBreak_table = data.frame(table(sum_intBreak$insite))

        multi_break = filter(sum_intBreak_table,Freq >50)$Var1 %>% as.character()

        bed_serial[[i]] = bed_serial[[i]] %>%
          mutate(multiBreak = ifelse(insite %in% multi_break,1,0),
                 dupBreak = ifelse(duplicated(breakSite),1,0))

        dup_idx = with(bed_serial[[i]],which(multiBreak==0 & dupBreak==1))

        if(length(dup_idx)>0){
          bed_serial[[i]] = bed_serial[[i]][-dup_idx,]
        }

        bed_merge = as.data.frame(table(bed_serial[[i]]$insite,dnn=list('coordinate')),responseName='frequency')

        bed_merge = bed_merge %>% arrange(desc(frequency))
        bed_merge$coordinate = as.character(bed_merge$coordinate)

        bed_merge_loc = t(data.frame(strsplit(bed_merge$coordinate,':')))

        if (nrow(bed_merge_loc)==0){
          fwrite(list('Warning: no effective integration site had been detected!'),
                 paste0(outdir,'/warning.txt'))
        }else{
          bed_merge = bed_merge %>% mutate(
            chr = bed_merge_loc[,1],
            loc = bed_merge_loc[,2],
            strand = bed_merge_loc[,3]
          )

          bed_merge$loc = as.numeric(bed_merge$loc)
          bed_merge = bed_merge %>% arrange(chr,loc)
          bed_merge$neighbour = 0
          bed_merge$neighbour_head = 0

          j=1
          while(j < nrow(bed_merge)) {
            if(bed_merge$chr[j]==bed_merge$chr[j+1] & bed_merge$strand[j]==bed_merge$strand[j+1]){
              if(bed_merge$loc[j+1]-bed_merge$loc[j] < collapse){
                recorder = 0
                k = 1
                while (j+k<=nrow(bed_merge)) {
                  if(bed_merge$loc[k+j]-bed_merge$loc[k+j-1] < collapse){
                    recorder = c(recorder,k)
                    k=k+1
                  }else break
                }
                bed_merge$neighbour[j+recorder] = 1
                bed_merge$neighbour_head[j+which.max(bed_merge$frequency[j+recorder])-1] = 1
                bed_merge$frequency[j+which.max(bed_merge$frequency[j+recorder])-1] = sum(bed_merge$frequency[j+recorder])
                j=j+max(recorder)
              }else(j=j+1)

            }else(j=j+1)
          }

          rmRows = which(bed_merge$neighbour+bed_merge$neighbour_head==1)
          if(length(rmRows)>0){
            bed_merge = bed_merge[-which(bed_merge$neighbour+bed_merge$neighbour_head==1),]
          }

          fwrite(bed_serial[[i]],
                 paste0(outdir,'/',strsplit(basename(fqBed),"\\.")[[1]][1],'_insite_',i,'.csv'))
          fwrite(bed_merge[,1:2],
                 paste0(outdir,'/',strsplit(basename(fqBed),"\\.")[[1]][1],'_loc_',i,'.csv'))
          fwrite(log_info[[paste0(i)]],
                 paste0(outdir,'/',strsplit(basename(fqBed),"\\.")[[1]][1],'_B',i,'_P2.log'),
                 quote = FALSE)
      }
      }
    }

    log_info
  }
}
