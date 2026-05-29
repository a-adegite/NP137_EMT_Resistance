# Load library
library(here)
library(DESeq2)
library(ggplot2)
library(pheatmap)

# Load count data
count <- read.csv(
  here("data", "bulk_rnaseq", "GSE225688_PTEN_Raw_count_cleaned.csv"),
  row.names = 1
)

# Load metadata
metadata <- read.csv(
  here("data", "metadata", "GSE225688_PTEN_metadata.csv"),
  stringsAsFactors = T,
)

rownames(metadata) <- metadata$sample

# Checks if col names in count matrix matches row names in metadata
all(rownames(metadata) == colnames(count))

# Create deseq object
dds <- DESeqDataSetFromMatrix(
  countData = count,
  colData = metadata,
  design = ~treatment
)

# PCA
# Perform regularized log transform
rld <- rlog(dds, blind = T)

pca <- plotPCA(rld, intgroup = "treatment")
pca

# Sample-to-sample correlation heatmap
# Calculate sample correlation
sample_corr = cor(assay(rld))

# Create an annotation dataframe
annotation_df <- data.frame(
  treatment = metadata$treatment
)

# Match row names to samples
rownames(annotation_df) <- colnames(sample_corr)


# Plot Sample-to-sample heatmap
pheatmap(
  sample_corr,
  annotation_col = annotation_df
)


