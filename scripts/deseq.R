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
