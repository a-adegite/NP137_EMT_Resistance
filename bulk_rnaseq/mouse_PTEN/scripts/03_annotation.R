# Load libraries
library(here)
library(biomaRt)

# Load count data
count <- read.csv(
  here("data", "bulk_rnaseq", "GSE225688_PTEN_Raw_count_cleaned.csv"),
  row.names = 1
)

# Load DGE data

# Raw DGE analysis result
DGE_raw <- read.csv(
  here(
    "bulk_rnaseq", "mouse_PTEN", "result",
    "differential_gene_expression", "DGE_results_raw.csv"
  ),
  row.names = 1
)

#Raw DGE analysis result (Log2FoldChange Shrunken)
DGE_raw_lfcs <- read.csv(
  here(
    "bulk_rnaseq", "mouse_PTEN", "result",
    "differential_gene_expression", "DGE_results_all_shrunken.csv"
  )
)

# Significant genes (Log2FoldChange Shrunken)
sig_gene_lfcs <- read.csv(
  here(
    "bulk_rnaseq", "mouse_PTEN", "result", "differential_gene_expression",
    "DGE_lfcs_significant.csv"
  ),
  row.names = 1
)

# Clean ensembl id
ensemblID <-  rownames(count)
ensemblID <- gsub("\\..*$", "", ensemblID)

# Select database and create connection to ensembl
ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl"
)

gene_info_all <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "gene_biotype",
    "description"
  ),
  filters = "ensembl_gene_id",
  values = unique(ensemblID),
  mart = ensembl
)

