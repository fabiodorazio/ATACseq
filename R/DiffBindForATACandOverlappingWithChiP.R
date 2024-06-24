#generate a sample sheet inidicating the PATH of bed file outputs of macs2

source("https://bioconductor.org/biocLite.R")
biocLite('RDAVIDWebService')

library(DiffBind)
library(ChIPpeakAnno)
library(ChIPseeker)
library(GenomicFeatures)

setwd(system.file("extra", package="DiffBind"))

Tdrd7 = dba(sampleSheet="ATACTdrd7PGCTotalPeaks.csv")

Tdrd7 = dba(Tdrd7, mask = Tdrd7$mask$PGC) #subsetting dba object

Tdrd7 = dba.count(Tdrd7) #counts reads from bam files
Tdrd7 = dba.contrast(Tdrd7, categories = DBA_CONDITION, minMembers = 2) # defines diff condition for differential analyses
Tdrd7 = dba.contrast(Tdrd7, group1 = Tdrd7$masks$PGC, group2 = Tdrd7$masks$Somatic)
Tdrd7 = dba.analyze(Tdrd7)

Tdrd7.DB <- dba.report(Tdrd7) #subset the differential bound samples in a grange object

totalSig <- subset(Tdrd7.DB, Tdrd7.DB$FDR < 0.01)
Downreg <- subset(totalSig, totalSig$Fold < 0) # downreg in PGCs vs Soma (or group1 vs group2)
Upreg <- subset(totalSig, totalSig$Fold > 0) # upreg in PGCs vs Soma (or group1 vs group2)

#save txt
DownregEx = data.frame(Downreg)
write.table(DownregEx, 'ATAC_Peaks_UpregInSoma_prim5.txt', row.names = F, sep = '\t')

#find promoters vs distal elements
txdb <- makeTxDbFromUCSC('danRer7', 'ensGene')
promoter = getPromoters(TxDb=txdb, upstream=200, downstream=200)

distalelementsSoma = subsetByOverlaps(Downreg, promoter, invert = T, ignore.strand = T)

#plot distribution on genomic elements

aCRUp<-assignChromosomeRegion(Upreg, nucleotideLevel=FALSE, 
                            precedence=c("Promoters", "immediateDownstream", 
                                         "fiveUTRs", "threeUTRs", 
                                         "Exons", "Introns"), 
                            TxDb=txdb)
barplot(aCRUp$percentage, las=3)

aCRDown<-assignChromosomeRegion(Downreg, nucleotideLevel=FALSE, 
                              precedence=c("Promoters", "immediateDownstream", 
                                           "fiveUTRs", "threeUTRs", 
                                           "Exons", "Introns"), 
                              TxDb=txdb)
barplot(aCRDown$percentage, las=3)

# overlap with enhancers
active_enhancers = read.csv('allstage_activeEnhancer_merged_sorted.bed', header = F, sep = '\t')
colnames(active_enhancers) <- c('chromosome', 'start', 'end')
active_enhancers <- toGRanges(active_enhancers)

active_enhancers = read.csv('Enhancer Yixuan/allstage_activeEnhancer_merged_sorted_Zv9.bed', header = F, sep = '\t')
colnames(active_enhancers) = c('chromosome', 'start', 'end')
active_enhancersGR = toGRanges(active_enhancers)

Top_Downreg = Downreg[order(Downreg$Fold),]
Top_Downreg = head(Top_Downreg, 2000)

overlapSomaEnhancers <- findOverlapsOfPeaks(active_enhancers, distalelementsSoma)
makeVennDiagram(overlapSomaEnhancers)
SomaEnhancers = subsetByOverlaps(distalelementsSoma, active_enhancers, ignore.strand = T)

# correlation with RNA-Seq
annoData <- toGRanges(txdb, feature="gene")

SomaEnhancersID <- annotatePeakInBatch(SomaEnhancers, 
                                   AnnotationData=annoData, 
                                   output="nearestBiDirectionalPromoters",
                                   bindingRegion=c(-30000, 30000))
distalelementsSoma.frame <- annotatePeakInBatch(distalelementsSoma,
AnnotationData=annoData,
output="nearestBiDirectionalPromoters",
bindingRegion=c(-30000, 30000))
distalelementsSoma.frame <- data.frame(distalelementsSoma.frame)

enhancers <- read.csv('allstage_activeEnhancer_merged_sorted_Zv9.bed', sep = '\t', header = F)
colnames(enhancers) = c('chromosome', 'start', 'end')
active_enhancersGR = toGRanges(enhancers)
overlapSomaEnhancers <- subsetByOverlaps(active_enhancersGR, distalelementsSoma)
overlapSomaEnhancersID <- annotatePeakInBatch(overlapSomaEnhancers, 
                                              AnnotationData=annoData, 
                                              output="nearestBiDirectionalPromoters",
                                              bindingRegion=c(-50000, 50000))
overlapSomaEnhancersID.frame <- data.frame(overlapSomaEnhancersID)
overlap.genes <- merge(overlapSomaEnhancersID.frame, genesSoma, by = 'feature')

enhancersIds <- as.vector(overlapSomaEnhancersID$feature)

tpm$feature <- rownames(tpm)
tpm.enhancer <- merge(tpm, distalelementsSoma.frame, by = 'feature')
tpm.enhancer <- tpm.enhancer[order(tpm.enhancer$Fold),]

tpm.enhancer$meanSoma <- rowMeans(tpm.enhancer[c('sPrim5Soma1', 'sPrim5Soma2')],)
tpm.enhancer$meanPGC <- rowMeans(tpm.enhancer[c('sPrim5PGC1', 'sPrim5PGC2')],)
pheatmap(log(newframe+1), cluster_rows = F, cluster_cols = F)



setwd("//adf/Storage/F/D/FMD523/Desktop/FabioRNAseq (2)/2108/RNAseqtables")
tpm <- read.csv('tpmStagesFiltered.txt', header = T, sep = '\t')

SomaGenesIDs <- as.vector(SomaEnhancersID$feature)
tpmSomaGenes <- tpm[SomaGenesIDs,]
boxplot(log(tpmSomaGenes[,c(17:20)]))

# venn of overlapping upregulated peaks and genes
ids <- as.vector(rownames(SomaEnhancersID.frame))
overlap.genes.peaks.soma <- na.omit(overlap.genes.peaks.soma)
overlap.genes.peaks.soma <- overlap.genes.peaks.soma[,c(17:20)]

overlap.genes.peaks.soma <- tpm[ids,]
overlap.genes.peaks.soma$mean <- rowMeans(overlap.genes.peaks.soma[,c(3:4)])

overlap.genes.peaks.soma <- overlap.genes.peaks.soma[order(overlap.genes.peaks.soma$mean),]

allpeaks.frame <- annotatePeakInBatch(Tdrd7.DB, 
                                       AnnotationData=annoData, 
                                       output="nearestBiDirectionalPromoters",
                                       bindingRegion=c(-30000, 30000))

Tdrd7.DB.frame <- data.frame(allpeaks.frame)
library(plyr)
ddply(Tdrd7.DB.frame, 'feature', numcolwise(sum))
rnaseqDiffExpressed <- read.csv('log2FoldChange_padj0.1_s24PGCVSs24Soma (1).txt', sep = '\t')
# GO germ cell GO_REF:0000107

#filter for developmental genes and distance from tss < 10000
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset='drerio_gene_ensembl')
mart <- useDataset("drerio_gene_ensembl", useMart("ensembl"), verbose = FALSE)
genes <- df$genes
df<-df[,-4]
s4 <- rownames(s3)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "go_id", "external_gene_name"),values='',mart= mart)
grep("go", listAttributes(mart)[,1], value = TRUE) #check term for go in attributes
row.names(G_list) <- G_list$ensembl_gene_id
G_list$ensembl_gene_id <- NULL
s5 <- merge(s3,G_list, by = 'row.names')
row.names(s5) <- s5$external_gene_name
s5$Row.names <- NULL
s5$external_gene_name <- NULL
s6 <- order(rowMeans(s5), decreasing=TRUE)

#density map for enhancers and promoters

#tagMatrix <- getTagMatrix(PGCMO1, windows=promoter)

#ovlp = subsetByOverlaps(PGC5mm1, promoter, invert = T, ignore.strand = T)

#tagHeatmap(tagMatrix, xlim=c(-1000, 1000), color="red")

#finds peaks that overlap in Tdrd7 MO and mismatch and somatic with enhancer marks
Over_lap <- findOverlapsOfPeaks(Tdrd7.DB, H3K4me1)
Over_lap <- addMetadata(Over_lap, colNames="score", FUN=mean) 
Over_lap$peaklist[["///gr2"]][1:2] #shows the peak list
makeVennDiagram(Over_lap)

path.to.chip <- "~../ATAC analysis/ChIPseqDataForEnhancers/GSM915199_H3K27ac_24hpf_Zv9_reads.bed"
h3k27 <- read.csv(paste0(path.to.chip, '/GSM915199_H3K27ac_24hpf_Zv9_reads.bed'))

colnames(h3k27) <- c('chr', 'start', 'end', 'V4', 'V5', 'strand')

path.to.atac <- "~/Desktop/PhD-March-2019-backup/ATAC analysis"
up.pgc <- read.csv(paste0(path.to.atac, "/UpregInPGCsATACpeaksENSEMBLprim5.txt"), sep = '\t')


