### Bisulfite seq analysis ###
## data downloaded from :https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41923 ###

library(ChIPseeker)
library(ChIPpeakAnno)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(org.Mm.eg.db)
library(BSgenome.Mmusculus.UCSC.mm9)
library(biomaRt)

args <- commandArgs(trailingOnly = TRUE)

# location
# '/mnt/biggley/home/fabio/Fabio_Projects/CAGE_Fabio/CAGE_Mouse_ESC/'
path.to.files <- args[1]

coverage <- read.csv(paste0(path.to.files, 'GSM1027571_DNA_CpG_coverage_E14_serum_LIF.bedGraph'),
                     sep = '\t', header = F, skip = 1)
meth.level <- read.csv(paste0(path.to.files, 'GSM1027571_DNA_CpG_methcounts_E14_serum_LIF.bedGraph'), 
                       sep = '\t', header = F, skip = 1)

## merge the two dataframes if the first three column match

merge.meth.data <- function(x, y){
  ifelse(x[,c(1:3)] == y[,c(1:3)],
         test <- data.frame(x,y[,4]), stop('x and y are not equal'))
  return(test)
}
bis.seq <- merge.meth.data(coverage, meth.level)
colnames(bis.seq) <- c('chr', 'start', 'end', 'coverage', 'meth')

## set a threshold to the coverage and convert it to GRanges
bis.seq.threshold <- subset(bis.seq, bis.seq$coverage > 10)
bis.seq.threshold.range <- GRanges(seqnames = bis.seq.threshold$chr,
                                   ranges = IRanges(start = bis.seq.threshold$start, end = bis.seq.threshold$end),
                                   coverage = bis.seq.threshold$coverage, meth = bis.seq.threshold$meth,
                                   seqlengths = seqlengths(BSgenome.Mmusculus.UCSC.mm9))
seqinfo(bis.seq.threshold.range) <- seqinfo(BSgenome.Mmusculus.UCSC.mm9)[seqlevels(bis.seq.threshold.range)]

## take all CpG with a score higher than half and viceversa
bis.seq.high <- subset(bis.seq, bis.seq$meth > 0.5)
bis.seq.low <- subset(bis.seq, bis.seq$meth <= 0.5)


### aim 1:
## plot DNA meth levels for sharp and broad promoters based on mouse CAGE

#### START ####

## retrieve mouse promoters
mouse67_mart <- useMart(host='may2012.archive.ensembl.org',
                        biomart='ENSEMBL_MART_ENSEMBL',
                        dataset = "mmusculus_gene_ensembl")

promoters_from_biomart <- function (mart, upstream=500, downstream=500) {
  # retrieves biomart attributes
  mouse_genes <- getBM(attributes = c("ensembl_gene_id",
                                      "ensembl_transcript_id",
                                      "external_gene_id",
                                      "description",
                                      "chromosome_name",
                                      "transcript_start",
                                      "transcript_end",
                                      "strand",
                                      "gene_biotype"),
                       mart = mart)
  mouse_genes <- mouse_genes[mouse_genes$gene_biotype == "protein_coding",] # subset protein coding transcripts only
  mouse_genes$strand <- ifelse(mouse_genes$strand==1, "+",
                               ifelse(mouse_genes$strand==-1, "-", "*")) # replaces +-1 with +-
  
  mouse_genes_gr <- makeGRangesFromDataFrame(mouse_genes,
                                             start.field = "transcript_start",
                                             end.field = "transcript_end",
                                             keep.extra.columns = TRUE)
  seqlevelsStyle(mouse_genes_gr) <- "UCSC"
  
  #make promoters with total width (-upstream, downstream) to use for CAGE clustering
  mouse_promoters <- promoters(mouse_genes_gr,
                               upstream = upstream,
                               downstream = downstream)
  promoters_df <- GenomicRanges::as.data.frame(mouse_promoters)
  names(promoters_df)[1] <- "chr"
  
  # select only unique promoters for each gene - there can be multiple promoters
  # from different transcripts, but they should not have the same start and end
  promoters_df <- promoters_df[order(promoters_df$chr,
                                     promoters_df$start,
                                     promoters_df$strand,
                                     promoters_df$ensembl_gene_id
  ),
  ]
  promoters_df <- promoters_df[!duplicated(promoters_df[c("chr", "start", "strand")]),] # discards lines duplicated in chr, start and end
  return(promoters_df)
}
mouse.promoters <- promoters_from_biomart(mouse67_mart, 1000, 1000)

#### END ####

### aim 2:
## define sharp and broad promoters 

#### START ####

cageESC <- readRDS('CAGEsets/CAGEsetMouseESCfantom5.rds')
CAGEset <- tagClusters(cageESC, sample = 'CAGE_mouseMm9_ESC', returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)

## alternatively annotate the cage set
CAGEset <- readRDS("CAGEsetMouseESCfantom5.rds")
return.cage.anno <- function(x){
  cage.set.ESC <- tagClusters(x,sample = "CAGE_mouseMm9_ESC", returnInterquantileWidth = TRUE,
                              qLow = 0.1,qUp = 0.9)
  ## annotate peaks
  txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
  cage.range <- toGRanges(cage.set.ESC)
  anno.cage <- annotatePeak(cage.range, TxDb = txdb,  annoDb = "org.Mm.eg.db", sameStrand = TRUE, verbose = FALSE)
  cage.esc.anno <- data.frame(anno.cage@anno)
  cage.esc.anno <- subset(cage.esc.anno, cage.esc.anno$ENSEMBL != 'NA')
  return(cage.esc.anno)
}
cage.esc.anno <- return.cage.anno(CAGEset)

## define sharp and broad promoters in ESC
sIdx <- subset(cage.esc.anno, CAGEset$interquantile_width <= 9)
bIdx <- subset(cage.esc.anno, CAGEset$interquantile_width > 9)

## function for verifying the methylation levels of promoters with CAGE peaks
## overlaps cage and dna methylation, returns a list of mouse promoters with a CAGE peak
## and promoters with CpGs
define.sharp.broad <- function(CAGE_width, dna_meth){
  ## overlap with dna methylation
  CAGEset <- makeGRangesFromDataFrame(CAGE_width,
                                      start.field = "start",
                                      end.field = "end",
                                      keep.extra.columns = TRUE)
  ## sites where transcription starts
  dna_meth_width <- subsetByOverlaps(dna_meth, CAGEset)
  
  ## promoter regions
  mouse.promoters <- GRanges(mouse.promoters)
  promoters <- subsetByOverlaps(mouse.promoters, CAGEset)
  dna_meth_prom <- subsetByOverlaps(dna_meth, promoters)
  
  dna_meth_width <- data.frame(dna_meth_width)
  dna_meth_prom <- data.frame(dna_meth_prom)
  
  list.width <- list(dna_meth_width, dna_meth_prom)
  return(dna_meth_prom)
  # mean DNA methylation ## (sum(dna_meth_prom$meth))/(nrow(dna_meth_prom))
  mean(dna_meth_prom)
}

broad.meth <- define.sharp.broad(bIdx, bis.seq.threshold.range)
sharp.meth <- define.sharp.broad(sIdx, bis.seq.threshold.range)
# mean DNA methylation ## (sum(broad.meth$meth))/(nrow(broad.meth))
mean(broad.meth)
mean(sharp.meth)

#### START ####
## define meth level for up and downregulated promoters by DNMTKO
up.down.promoters <- function(x,y,z){
  x <- x[rownames(y),]
  x <- na.omit(x)
  promotersTSS<-GRanges(seqnames = x$chr,
                        ranges=IRanges(start = x$dominant_ctss, end = x$dominant_ctss),
                        strand = x$strand,
                        IQwidth = x$interquantile_width,
                        tpm = x$tpm,
                        seqlengths = seqlengths(BSgenome.Mmusculus.UCSC.mm9))
  
  seqinfo(promotersTSS)<-seqinfo(BSgenome.Mmusculus.UCSC.mm9) ## adds chromosome lengths
  promflank <- promoters(promotersTSS,upstream=1000,downstream=500) ## extends promoter range 1000 and 500 bp
  
  promflank = promflank[width(trim(promflank)) == 1500] ## ensures all widths are of the same lenght
  ## to order by coordinate
  promflank.ordered <- promflank
  ## to order by IQwidth
  promflank.ordered<-promflank[order(promflank$IQwidth),] ## orders by IQ width
  promflank.ordered<-keepStandardChromosomes(promflank.ordered)
  promflank.ordered<-dropSeqlevels(promflank.ordered,"chrM",pruning.mode = "coarse")
  ## overlap with DNA methylated CpGs
  dna_meth <- subsetByOverlaps(z, promflank.ordered)

  ## plot heatmap
  weight <- promflank.ordered$tpm
  
  hm.l = CoverageHeatmap(
    promflank.ordered,
    dna_meth,
    weight=weight,
    coords = c(-1000,500),
    label = "tpm.dominant_ctss")
  
  #smooth heatmap and scale 
  hm_smoothed.l <- smoothHeatmap(hm.l, sigma = c(2, 3), output.size=c(5000, 500))
  scale(hm_smoothed.l) <- quantile(weight, c(0.1, 0.9))
  
  # plot heatmaps
  samplename <- deparse(substitute(y))
  pdf(paste('SmoothHeatmapCpGmeth', samplename, '.pdf'), height = 8, width = 6)
  plotHeatmapList(hm_smoothed.l,
                  cex.label = 1.5,
                  legend = TRUE,
                  legend.pos = "l",
                  legend.width = 0.15,
                  color = c("white", "blue"))
  dev.off()
  ## plot density heatmap
  PromotersTSSflank <- promoters(promotersTSS, upstream = 600, downstream = 600)
  PromotersTSSflank <- PromotersTSSflank[order(PromotersTSSflank$IQwidth),]
  PromotersTSSflankSeq <- getSeq(Mmusculus, PromotersTSSflank)
  
  pdf(paste('DensityHeatmap', samplename, '.pdf'), height = 8, width = 6)
  plotPatternDensityMap(PromotersTSSflankSeq, c("TA"))
  dev.off()
  return(PromotersTSSflank)
}

up.promoters <- up.down.promoters(cage.esc.anno, dl1, bis.seq.threshold.range)
down.promoters <- up.down.promoters(cage.esc.anno, dl2, bis.seq.threshold.range)

#### START ####
## define promoters
define.promoters <- function(CAGE_width, dna_meth){
  ## overlap with dna methylation
  CAGEset <- makeGRangesFromDataFrame(CAGE_width,
                                      start.field = "start",
                                      end.field = "end",
                                      keep.extra.columns = TRUE)
  ## sites where transcription starts
  dna_meth_width <- subsetByOverlaps(dna_meth, CAGEset)
  
  ## promoter regions
  mouse.promoters <- GRanges(mouse.promoters)
  promoters <- subsetByOverlaps(mouse.promoters, CAGEset)
  dna_meth_prom <- subsetByOverlaps(dna_meth, promoters)
  
  dna_meth_width <- data.frame(dna_meth_width)
  dna_meth_prom <- data.frame(dna_meth_prom)
  
  list.width <- list(dna_meth_width, dna_meth_prom)
  return(dna_meth_prom)
  # mean DNA methylation ## (sum(dna_meth_prom$meth))/(nrow(dna_meth_prom))
  mean(dna_meth_prom)
}

## upregulated and downregulated by DNMTKO

## plot heatmap
plot.heatmap <- function(x,y){
  promotersTSS<-GRanges(seqnames = x$chr,
                        ranges=IRanges(start = x$dominant_ctss, end = x$dominant_ctss),
                        strand = x$strand,
                        IQwidth = x$interquantile_width,
                        tpm = x$tpm,
                        seqlengths = seqlengths(BSgenome.Mmusculus.UCSC.mm9))
  seqinfo(promotersTSS)<-seqinfo(BSgenome.Mmusculus.UCSC.mm9) ## adds chromosome lengths
  promflank <- promoters(promotersTSS,upstream=1000,downstream=500) ## extends promoter range 1000 and 500 bp
  
  promflank = promflank[width(trim(promflank)) == 1500] ## ensures all widths are of the same lenght
  ## to order by coordinate
  promflank.ordered <- promflank
  ## to order by IQwidth
  promflank.ordered<-promflank[order(promflank$IQwidth),] ## orders by IQ width
  promflank.ordered<-keepStandardChromosomes(promflank.ordered)
  promflank.ordered<-dropSeqlevels(promflank.ordered,"chrM",pruning.mode = "coarse")
 
  dna_meth <- subsetByOverlaps(y, promflank.ordered)
  ## plot heatmap
  weight <- promflank.ordered$tpm
  
  hm.l = CoverageHeatmap(
    promflank.ordered,
    dna_meth,
    weight=weight,
    coords = c(-1000,500),
    label = "tpm.dominant_ctss")
  
  #smooth heatmap and scale 
  hm_smoothed.l <- smoothHeatmap(hm.l, sigma = c(2, 3), output.size=c(5000, 500))
  scale(hm_smoothed.l) <- quantile(weight, c(0.1, 0.9))
  
  # plot heatmaps
  samplename <- deparse(substitute(x))
  pdf(paste('SmoothHeatmap', samplename, '.pdf'), height = 8, width = 6)
  plotHeatmapList(hm_smoothed.l,
                  cex.label = 1.5,
                  legend = TRUE,
                  legend.pos = "l",
                  legend.width = 0.15,
                  color = c("white", "blue"))
  dev.off()
  return(promflank.ordered)
}

plot.heatmap(CAGEset, bis.seq.threshold.range)

#### END ####
annotate.cgs <- function(x, txdb, genome){
  
  grange <- subsetByOverlaps(promoter, grange, ignore.strand = T)
  
  grange.anno <- annotatePeak(grange, tssRegion = c(-1000, 1000), TxDb = txdb,
                              annoDb = genome, sameStrand = T, verbose = F)
  grange.anno2 <- data.frame(grange.anno@anno)
  grange.anno3 <- subset(grange.anno2, grange.anno2$annotation == 'Promoter')
  return(grange.anno3)
}
anno.bis.seq <- annotate.cgs(bis.seq.threshold, txdb9, org.Mm.eg.db)

