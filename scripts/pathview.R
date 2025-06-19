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