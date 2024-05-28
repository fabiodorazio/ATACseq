#' @include dk
#' @author Fabio M D'Orazio
#' @description Split read count by Genomic Position

#####################################################
### split read counts by Genomic Position ####

library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(dplyr)
library(ChIPseeker)

## bin the genome in 100 bp fragments
bin100 <- tileGenome(seqinfo(BSgenome.Hsapiens.UCSC.hg19), tilewidth = 100)
bin100 <- unlist(bin100)
## add one column for ensembl ids
bin100$BIN <- paste0('BIN', c(1:length(bin100)))

## read csv file
csv.coord <- readRDS('../CAGEpromHs19colon.rds')
csv.coord.range <- GRanges(csv.coord)

subset.overlap.sample <- function(x,y){

  ov.first <- findOverlaps(x, y) 
  x[queryHits(ov.first),]
  ov.first <- y[subjectHits(ov.first),]
  ov.first <- data.frame(ov.first)

  ov.second <- findOverlaps(y, x) 
  y[queryHits(ov.second),]
  ov.second <- x[subjectHits(ov.second),]
  ov.second <- data.frame(ov.second)
  colnames(ov.second) <- paste0(colnames(ov.second), '_gene')
  #merge the two dataframes and convert to Granges
  merged.sample <- cbind(ov.second, ov.first)
  merged.sample.frame <- GRanges(merged.sample)
  
  ## to group peaks belonging to the same bin and sum their tpm
  merged.sample.frame.grouped <- merged.sample.frame %>%
    data.frame() %>% group_by(BIN) %>% summarise(sum_tpm = sum(tpm))
  
  ## remove duplicated BIN from merged.sample
  temp_merged <- merged.sample[!duplicated(merged.sample$BIN),] 
  ## merge the two data frames to retrive coordinates
  final_merged <- merge(as.data.frame(merged.sample.frame.grouped), temp_merged, by = 'BIN')
  rownames(final_merged) <- final_merged$BIN
  final_merged$BIN <- NULL
  ## annotate bins
  
  return(final_merged)
}

# run function
binned_tpm1 <- subset.overlap.sample(csv.coord.range, bin100)
