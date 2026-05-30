# Load libraries
library(here)
library(AnnotationDbi)
library(EnsDb.Mmusculus.v79)

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
  ),
  row.names = 1
)

# Significant genes (Log2FoldChange Shrunken)
sig_gene_lfcs <- read.csv(
  here(
    "bulk_rnaseq", "mouse_PTEN", "result", "differential_gene_expression",
    "DGE_lfcs_significant.csv"
  ),
  row.names = 1
)

# Create EnsemblID field
DGE_raw$GENEID <- gsub("\\..*$", "", rownames(DGE_raw))
DGE_raw_lfcs$GENEID <- gsub("\\..*$", "", rownames(DGE_raw_lfcs))
sig_gene_lfcs$GENEID <- gsub("\\..*$", "", rownames(sig_gene_lfcs))

# Clean ensembl id
ensemblID <-  rownames(count)
ensemblID <- gsub("\\..*$", "", ensemblID)

# Create EnsDb.Mmusculus.v79 object
edb <- EnsDb.Mmusculus.v79

#keytypes(edb)

gene_info_all <- AnnotationDbi::select(
  edb,
  keys = unique(ensemblID),
  keytype = "GENEID",
  columns = c(
    "GENEID",
    "SYMBOL",
    "GENEBIOTYPE"
  )
)

# Merge annotations with DGE Results
DGE_raw_annot <- merge(DGE_raw, gene_info_all, by = "GENEID")
DGE_raw_lfcs_annot <- merge(DGE_raw_lfcs, gene_info_all, by = "GENEID")
sig_gene_lfcs_annot <-merge(sig_gene_lfcs, gene_info_all, by = "GENEID")

# Save annotated DGE result
write.csv(
  DGE_raw_annot,
  here("bulk_rnaseq", "mouse_PTEN", "result", "differential_gene_expression", "DGE_results_raw_annotated.csv")
)

write.csv(
  DGE_raw_lfcs_annot,
  here("bulk_rnaseq", "mouse_PTEN", "result", "differential_gene_expression", "DGE_results_all_shrunken_annotated.csv")
)

write.csv(
  sig_gene_lfcs_annot,
  here("bulk_rnaseq", "mouse_PTEN", "result", "differential_gene_expression", "DGE_lfcs_significant_annotated.csv")
)
