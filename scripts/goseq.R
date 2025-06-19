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

# Visualiseren van GO-analyse
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
