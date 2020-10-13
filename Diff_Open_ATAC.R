### differential openess 
#atacUp <- read.csv('../All_ATAC_Upreg_PGC_Prim5.txt', sep = '\t')
#atacDown <- read.csv('../All_ATAC_Upreg_Soma_Prim5.txt', sep = '\t')

## create chr and range column from merged rownames to generate a Grange
convert.to.grange <- function(x){
  x$chr <- as.character(lapply(strsplit(rownames(x), "\\:"), function(u) u[[1]][1]))
  x$range <- as.character(lapply(strsplit(rownames(x), "\\:"), function(u) u[[2]][1]))
  x$start <- as.character(lapply(strsplit(x$range, "\\-"), function(u) u[[1]][1]))
  x$end <- as.character(lapply(strsplit(x$range, "\\-"), function(u) u[[2]][1]))
  x$range <- NULL
  return(x)
}
atacUpRange <- convert.to.grange(atacUp)
atacDownRange <- convert.to.grange(atacDown)

library(AnnotationDbi)
txdb <- loadDb('../annotation/txdb_DanRer7.sqlite')
ensembl.ids.atac <- function(x){
  grange.obj <- toGRanges(x)
  anno.ensembl.id <- annotatePeakInBatch(grange.obj, AnnotationData=toGRanges(txdb), 
                                         output="nearestBiDirectionalPromoters",
                                         bindingRegion=c(-10000, 10000))
  anno.ensembl.id.order <- anno.ensembl.id[order(abs(anno.ensembl.id$log2FoldChange)),]
  anno.ensembl.id.200 <- tail(anno.ensembl.id.order, 50)
  ## ENSEMBL IDs of genes with Upreg ATAC peaks
  ensembl.ids <- anno.ensembl.id.200$feature
  return(ensembl.ids)
  
}
