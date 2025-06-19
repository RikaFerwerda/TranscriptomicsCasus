setwd("C:/Users/rikaf/Downloads/TranscriptomicsCasus")
getwd()
unzip("Data_RA_Raw.zip", exdir = "RA_data")
install.packages('BiocManager')
BiocManager::install('Rsubread')
library(Rsubread)
buildindex(
  basename = 'ref_hg38',
  reference = 'GCF_000001405.40_GRCh38.p14_genomic.fna',
  memory = 4000,
  indexSplit = TRUE)

# RA monsters
align.RA1 <- align(index = "ref_hg38", readfile1 = "SRR4785979_1_subset40k.fastq", readfile2 = "SRR4785979_2_subset40k.fastq", output_file = "RA1.BAM")
align.RA2 <- align(index = "ref_hg38", readfile1 = "SRR4785980_1_subset40k.fastq", readfile2 = "SRR4785980_2_subset40k.fastq", output_file = "RA2.BAM")
align.RA3 <- align(index = "ref_hg38", readfile1 = "SRR4785986_1_subset40k.fastq", readfile2 = "SRR4785986_2_subset40k.fastq", output_file = "RA3.BAM")
align.RA4 <- align(index = "ref_hg38", readfile1 = "SRR4785988_1_subset40k.fastq", readfile2 = "SRR4785988_2_subset40k.fastq", output_file = "RA4.BAM")

# Controle monsters
align.Controle1 <- align(index = "ref_hg38", readfile1 = "SRR4785819_1_subset40k.fastq", readfile2 = "SRR4785819_2_subset40k.fastq", output_file = "Controle1.BAM")
align.Controle2 <- align(index = "ref_hg38", readfile1 = "SRR4785820_1_subset40k.fastq", readfile2 = "SRR4785820_2_subset40k.fastq", output_file = "Controle2.BAM")
align.Controle3 <- align(index = "ref_hg38", readfile1 = "SRR4785828_1_subset40k.fastq", readfile2 = "SRR4785828_2_subset40k.fastq", output_file = "Controle3.BAM")
align.Controle4 <- align(index = "ref_hg38", readfile1 = "SRR4785831_1_subset40k.fastq", readfile2 = "SRR4785831_2_subset40k.fastq", output_file = "Controle4.BAM")

# Sorteren en indexeren
library(Rsamtools)
  # Bestandsnamen van de monsters
samples <- c('Controle1', 'Controle2', 'Controle3', 'Controle4', 'RA1', 'RA2', 'RA3', 'RA4')
lapply(samples, function(s) {sortBam(file = paste0(s, '.BAM'), destination = paste0(s, '.sorted'))
})

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

# Count Matrix maken
# Je definieert een vector met namen van BAM-bestanden. Elke BAM bevat reads van een RNA-seq-experiment (bijv. behandeld vs. controle).

allsamples <- c('Controle1.BAM', 'Controle2.BAM', 'Controle3.BAM', 'Controle4.BAM', 'RA1.BAM', 'RA2.BAM', 'RA3.BAM', 'RA4.BAM')

count_matrix <- featureCounts(
  files = allsamples,
  annot.ext = "genomic.gtf",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE,
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE
)

# Volledige countmatrix 

counts <- count_matrix$counts
colnames(count_matrix_volledig) <- c('Gene ID', 'Controle1', 'Controle2', 'Controle3', 'Controle4', 'RA1', 'RA2', 'RA3', 'RA4')
rownames(count_matrix_volledig) <- count_matrix_volledig[,1]
count_matrix_volledig <- count_matrix_volledig[,-1]
write.csv(counts, "bewerkt_countmatrix_DEF.csv")

count_matrix_volledig <- round(count_matrix_volledig)

# Statistiek
treatment <- c('Controle', 'Controle', 'Controle', 'Controle', 'RA', 'RA', 'RA', 'RA')
treatment_table <- data.frame(treatment)
rownames(treatment_table) <- c('Controle1', 'Controle2', 'Controle3', 'Controle4', 'RA1', 'RA2', 'RA3', 'RA4')

# Packages installeren en inladen
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("KEGGREST")
library(DESeq2)
library(KEGGREST)

# DESeqDataSet aanmaken
dds <- DESeqDataSetFromMatrix(countData = count_matrix_volledig, colData = treatment_table, design = ~ treatment)
dds <- DESeq(dds)
resultaten <- results(dds)
write.table(resultaten, file = 'ResultatenCasusDEF.csv', row.names = TRUE, col.names = TRUE)

# Resultaten tellen die statistisch significant zijn
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange > 1, na.rm = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange < -1, na.rm = TRUE)

# Genen die eruit springen
hoogste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde <- resultaten[order(resultaten$padj, decreasing = FALSE), ]

# Visualisatie --> Volcano Plot
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  BiocManager::install("EnhancedVolcano")
}
library(EnhancedVolcano)

EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj')
# Alternatieve plot zonder p-waarde cutoff (alle genen zichtbaar)
EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0)

dev.copy(png, 'VolcanoplotCasus.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()

# GO-analyse
BiocManager::install('goseq')
BiocManager::install('geneLenDataBase')
BiocManager::install("org.Dm.eg.db")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("AnnotationDbi")
library("goseq")
library("geneLenDataBase")
library("org.Dm.eg.db")
library('org.Hs.eg.db')
library('TxDb.Hsapiens.UCSC.hg38.knownGene')
library('AnnotationDbi')
library(tidyverse)
library(dplyr)

# ALL bestand maken met alle genen
res_all <- results(dds, alpha = 0.05)
res_all_df <- as.data.frame(res_all)
write.csv(res_all_df, "ALL_genes_DESeq2.csv")

# Read in van de twee files
DEG <- read.table("ResultatenCasusDEF.csv", header = TRUE, sep = "\t", comment.char = "#", check.names = FALSE)
ALL <- read.table("ALL_genes_DESeq2.csv", header = TRUE, sep = "\t", comment.char = "#", check.names = FALSE)
# Vectoren maken van DEG en ALL object
res_df <- as.data.frame(resultaten)
filter_res <- filter(res_df, res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1 )

class(DEG)
DEG.vector <- rownames(filter_res)
ALL.vector <- rownames(resultaten)

gene.vector=as.integer(ALL.vector%in%DEG.vector)
names(gene.vector)=ALL.vector


pwf=nullp(gene.vector, 'hg19', 'geneSymbol')
GO.wall=goseq(pwf,"hg19", 'geneSymbol')

goResults <- goseq(pwf, "hg19", 'geneSymbol', test.cats = c("GO:BP"))

# Visualiseren
goResults %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")

library(GO.db)
GOTERM[[goResults$category[1]]]

# Significante genen in biologische context
if (!requireNamespace("pathview", quietly = TRUE)) {
  BiocManager::install("pathview")
}
library(pathview)

# KEGG-pathway visualiseren
resultaten[1] <- NULL
resultaten[2:5] <- NULL
gene_vector <- resultaten$log2FoldChange
names(gene_vector) <- rownames(resultaten)
pathview(
  gene.data = gene_vector,
  pathway.id = "hsa04062",  # KEGG ID voor cell cyclus – Hompo Sapiens
  species = "hsa",          # 'hsa' = Humaan in KEGG
  gene.idtype = "SYMBOL",     # Geef aan dat het KEGG-ID's zijn
  limit = list(gene = 5)    # Kleurbereik voor log2FC van -5 tot +5
)

BiocManager::install('clusterProfiler')
BiocManager::install('enrichplot')
BiocManager::install('pheatmap')
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(pheatmap)

# Stel: je hebt je lijst met DEG's als ENTREZ ID's
entrez_df <- bitr(DEG.vector,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)

entrez_gene_list <- entrez_df$ENTREZID

kk <- enrichKEGG(gene = entrez_gene_list,
                 organism = 'hsa')

dotplot(kk) + ggtitle("KEGG Pathway Enrichment")

# Selecteer 3 pathways
my_pathways <- c("hsa04062", "hsa04060", "hsa04660")

# Haal genen uit alleen deze pathways
genes_in_my_pathways <- subset(pathway2gene, ID %in% my_pathways)
genes_in_my_pathways_list <- unique(unlist(strsplit(genes_in_my_pathways$geneID, "/")))

# Map naar SYMBOLS
genes_symbols <- bitr(genes_in_my_pathways_list,
                      fromType = "ENTREZID",
                      toType = "SYMBOL",
                      OrgDb = org.Hs.eg.db)

# Filter matrix
heatmap_genes <- count_matrix_volledig[rownames(count_matrix_volledig) %in% genes_symbols$SYMBOL, ]

# Plot nettere heatmap
pheatmap(heatmap_genes, 
         scale = "row", 
         show_rownames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 8)

# Variance Stabilizing Transformation
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

# Zelfde filtering:
heatmap_genes <- vsd_mat[rownames(vsd_mat) %in% genes_symbols$SYMBOL, ]

# Nettere heatmap:
pheatmap(heatmap_genes, 
         scale = "row", 
         show_rownames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 8,
         fontsize_col = 10,
         main = "Pathway genes expression heatmap")
dotplot(kk) + ggtitle("KEGG pathway enrichment in RA vs Control")
