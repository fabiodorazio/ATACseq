# ATACseq

Small R sandbox for **ATAC-seq** analyses, with an emphasis on integrating **DNA methylation (bisulfite-seq)** with promoter classes defined from **CAGE tag clusters**

## Requirements

### R / Bioconductor packages

The script imports (at least) the following libraries:

- `ChIPseeker`
- `ChIPpeakAnno`
- `GenomicRanges`
- `TxDb.Mmusculus.UCSC.mm9.knownGene`
- `org.Mm.eg.db`
- `BSgenome.Mmusculus.UCSC.mm9`
- `biomaRt`

  
# Scripts in repo

## 1) Bisulfite Integration

Logic:
1. **Loads bisulfite bedGraphs** (coverage + methylation counts), merges them by genomic coordinates, and filters CpGs by coverage (default `> 10`)
2. **Builds promoter windows from BioMart** (mm9; using the May 2012 Ensembl archive host) and keeps unique promoters per gene. 
3. **Loads a CAGEset** (RDS) and computes tag clusters with **interquantile width** (`qLow=0.1`, `qUp=0.9`) then defines:
   - **sharp promoters:** `interquantile_width <= 9`
   - **broad promoters:** `interquantile_width > 9`
4. **Overlaps methylation with promoter regions / CAGE clusters** and summarizes methylation.
5. Generates **coverage heatmaps** (smoothed) and **pattern density plots** (example motif `"TA"`), writing PDFs
---
## Usage

```bash
Rscript Bisulfite_Integration.R /path/to/bisulfite_bedgraphs/
```

## 2) Coverage

Logic
1. **Read BAM**
2. **wResize Tn5 cut sizes**
3. **Calculte coverage** as peak width / genome size
4. **Fold Change** vs expected coverage

## Split reads by GP
Workflows for promoter-centered regions

Logic:
1. **Read BAM**
2. **Bin reads** by genomic position

