# Rsamstool installeren voor indexeren
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsamtools")
library(Rsamtools)

# BAI bestanden maken
bam_files <- list.files(pattern = "\\.sorted\\.bam$", full.names = TRUE)
sapply(bam_files, indexBam)

# FAI bestand maken
indexFa('GCF_000001405.40_GRCh38.p14_genomic.fna')
