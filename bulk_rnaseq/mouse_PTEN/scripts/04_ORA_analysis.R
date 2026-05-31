# Load Libraries
library(here)
library(enrichplot)
library(org.Mm.eg.db)
library(clusterProfiler)

# Load DEG results
DEG_raw <- read.csv(
  here(
    "bulk_rnaseq", "mouse_PTEN", "result",
    "differential_gene_expression", "DGE_results_raw_annotated.csv"
  ),
  row.names = 1
)

sig_genes <- subset(
  DEG_raw,
  DEG_raw$padj < 0.05 & abs(DEG_raw$log2FoldChange) > 1
)

# Upregulated significant genes
up_genes <- sig_genes[sig_genes$log2FoldChange > 1,]

# Downregulated significant genes
down_genes <-  sig_genes[sig_genes$log2FoldChange < -1, ]
  
# Get Ensembl ids of upregulated, downregulated and background genes
up_genes_id <- up_genes$GENEID
down_genes_id <- down_genes$GENEID
bg_genes_id <- DEG_raw$GENEID

# Convert gene IDs (Ensembl -> Entrez)
up_entrez <- bitr(
  up_genes_id,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

down_entrez <- bitr(
  down_genes_id,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

bg_entrez <- bitr(
  bg_genes_id,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

# Unregulated go
ego_up <- enrichGO(
  gene          = up_entrez$ENTREZID,
  universe      = bg_entrez$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable = TRUE
)

# Downregulated go
ego_down <- enrichGO(
  gene          = down_entrez$ENTREZID,
  universe      = bg_entrez$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable = TRUE
)

# KEGG pathway enrichment
# Upregulated KEGG

ekegg_up <- setReadable(
  enrichKEGG(
    gene          = up_entrez$ENTREZID,
    universe      = bg_entrez$ENTREZID,
    organism      = "mmu",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH"
  ),
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID"
)


# Downregulated KEGG
ekegg_down <- setReadable(
  enrichKEGG(
    gene          = down_entrez$ENTREZID,
    universe      = bg_entrez$ENTREZID,
    organism      = "mmu",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH"
  ),
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID"
)


# Visualisation

# Dotplot upregulated (GO)
dotplot(ego_up, showCategory = 10, title="GO Enrichment Analysis of Upregulated Genes")


# Dotplot Downregulated (GO)
dotplot(ego_down, showCategory = 10, title="GO Enrichment Analysis of Downregulated Genes")

# Bar plot Downregulated (GO)
barplot(ego_down, showCategory = 10, title="GO Enrichment Analysis of Downregulated Genes")

# Dotplot upregulated (KEGG)
dotplot(ekegg_up, showCategory = 10, title="KEGG Enrichment Analysis of Upregulated Genes")

# Dotplot downregulated (KEGG)
dotplot(ekegg_down, showCategory = 10, title="KEGG Enrichment Analysis of Downregulated Genes")

# Bar plot Downregulated (KEGG)
barplot(ekegg_down, showCategory = 10, title="KEGG Enrichment Analysis of Downregulated Genes")

# Convert Enrichment Results to Data Frames
ego_up_df <- as.data.frame(ego_up)

ego_down_df <- as.data.frame(ego_down)

ekegg_up_df <- as.data.frame(ekegg_up)

ekegg_down_df <- as.data.frame(ekegg_down)

# Sort by Adjusted p-value
ego_up_df <- ego_up_df[order(ego_up_df$p.adjust), ]
ego_down_df <- ego_down_df[order(ego_down_df$p.adjust), ]

ekegg_up_df <- ekegg_up_df[order(ekegg_up_df$p.adjust), ]
ekegg_down_df <- ekegg_down_df[order(ekegg_down_df$p.adjust), ]

# Save results
write.csv(
  ego_down_df,
  here("bulk_rnaseq", "mouse_PTEN", "result", "enrichment", "GO_down.csv")
)

write.csv(
  ekegg_down_df,
  here("bulk_rnaseq", "mouse_PTEN", "result", "enrichment", "KEGG_down.csv")
)



