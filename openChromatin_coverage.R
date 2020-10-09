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
  
  ## select only nucleosomes
  open_chromatin <- alignment_tn[width(alignment_tn) < 120]
  Cov <- coverage(open_chromatin) / expected.cov
  
  saveRDS(Cov, paste0('~/ATAC_FoldChange_Coverage120', name, '.rds'))
}
