# Load library
library(here)
library(DESeq2)
library(apeglm)

# Load count data
count <- read.csv(
  here("data", "bulk_rnaseq", "GSE225688_PTEN_Raw_count_cleaned.csv"),
  row.names = 1
)

# Load metadata
metadata <- read.csv(
  here("data", "metadata", "GSE225688_PTEN_metadata.csv"),
  stringsAsFactors = T
)

rownames(metadata) <- metadata$sample

all(colnames(count) == rownames(metadata))

# Set reference level
metadata$treatment <-  relevel(metadata$treatment, ref = "CTRL_AB")

# Create deseq2 object
dds <- DESeqDataSetFromMatrix(
  countData = count,
  colData = metadata,
  design = ~treatment
)

# Run deseq2
dds <- DESeq(dds)

# Check coefficient names
resultsNames(dds)

# Extract results
res <- results(
  dds,
  contrast = c("treatment", "NETRIN1_AB", "CTRL_AB")
)

# Get summary 
summary(res)

# Shrink log2 fold changes
res_shrink <- lfcShrink(
  dds,
  coef = "treatment_NETRIN1_AB_vs_CTRL_AB",
  type = "apeglm"
)

# Convert to data frame
res_df <- as.data.frame(res_shrink)

# Order by adjusted p-value
res_df <- res_df[order(res_df$padj), ]

sig_genes <- subset(
  res_df,
  padj < 0.05 & abs(log2FoldChange) > 1
)

# Save raw DGE analysis result
write.csv(
  as.data.frame(res),
  here("bulk_rnaseq", "mouse_PTEN", "result", "differential_gene_expression", "DGE_results_raw.csv")
)

# Save all genes with shrunken LFC
write.csv(
  res_df,
  here("bulk_rnaseq", "mouse_PTEN", "result", "differential_gene_expression", "DGE_results_all_shrunken.csv")
)

# Save significant genes only
write.csv(
  sig_genes,
  here("bulk_rnaseq", "mouse_PTEN", "result", "differential_gene_expression", "DGE_lfcs_significant.csv")
)

