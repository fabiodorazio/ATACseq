#' @include dk
#' @author Fabio M D'Orazio
#' @description ATAC Coverage Analysis

#####################################################
### Adjust cut sites and calculate Coverage for ATAC sequencing ####

library(plyr)
library(dplyr)
library(tidyr)
library(magrittr)
library(GenomicRanges)
library(GenomicAlignments)

## read bam files, adjust for Tn5 cut sites and calculate the coverage
readAlignBAM <- function(x){
  genome.size <- 1.42e9
  name <- strsplit(x, "\\.Fabio")
  name <- name[[1]][1]
  ## import bam as galigment pair
  alignment = readGAlignmentPairs(file = x)
  alignment_range = granges(alignment, on.discordant.seqnames = 'drop', use.names=T, use.mcols=F)
  
  expected.cov <- sum(width(alignment_range))/genome.size
  
  ## adjust for Tn5 cut sites
  tn1 = GRanges(seqnames(alignment_range), IRanges(start(alignment_range)+5, start(alignment_range)+5), strand(alignment_range))
  tn2 = GRanges(seqnames(alignment_range), IRanges(end(alignment_range)-4, end(alignment_range)-4), strand(alignment_range))
  alignment_tn = GRanges(seqnames(alignment_range), IRanges(start(tn1), end(tn2)), strand(alignment_range)) #sum of cut0 and cut1
  
  seqlevels(alignment_tn) = seqlevels(BSgenome.Drerio.UCSC.danRer7)
  seqinfo(alignment_tn)   = seqinfo(BSgenome.Drerio.UCSC.danRer7)
  alignment_tn = trim(alignment_tn)
  
  ## select open chromatin
  open_chromatin <- alignment_tn[width(alignment_tn) < 120]
  Cov <- coverage(open_chromatin) / expected.cov
  
  saveRDS(Cov, paste0('ATAC_FoldChange_Coverage120', name, '.rds'))
}

## resize to 1bp resolution for Tn5 cut sites
cov.function <- function(a, resize = 1){
  #a <- subset(a, width(a) < 120)
  seqlevels(a) = seqlevels(BSgenome.Drerio.UCSC.danRer7)
  seqinfo(a)   = seqinfo(BSgenome.Drerio.UCSC.danRer7)
  a = trim(a)
  ## remove blacklisted regions
  a <- subsetByOverlaps(a, blacklist, invert = T)
  ## Tn5 cut sites: resize = 1
  a <- resize(a, width = resize, fix = "start", ignore.strand = FALSE)

  a_plus <- a[strand(a) == '+']
  a_minus <- a[strand(a) == '-']

  genome.size <- 1.42e9
  expected.cov <- sum(width(a))/genome.size
  ## select only nucleosomes
  FoldChange.minus <- GRanges(coverage(a_minus) / expected.cov, strand = '-')
  FoldChange.plus <- GRanges(coverage(a_plus) / expected.cov, strand = '+')
  FoldChange.strand <- c(FoldChange.minus, FoldChange.plus)
  
  #export.bed(FoldChange, paste0(x, '.bed'))
  return(FoldChange.strand)
}
