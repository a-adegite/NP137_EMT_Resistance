# Load libraries
library(here)
library(GEOquery)

# Load count data
count <- read.delim(
  here("data", "bulk_rnaseq", "GSE225688_PTEN_Raw_count.csv"),
  sep = ";",
  row.names = 1
)

# Reformat count fields
colnames(count) <-  c("CTL_A", "CTL_B", "CTL_C", "NP137_D", "NP137_E", "NP137_F")

# Retrieve metadata
gse <- getGEO("GSE225688",GSEMatrix = T)
gse <- gse[[1]]
metadata <- pData(gse)

# Clean metadata
metadata <- metadata[, c(1,2,9,39)]

# Rename field names
colnames(metadata)[c(1,3,4)] <- c("sample", "organism", "treatment")

# Clean metadata values
metadata$treatment <- gsub("control antibody", "CTRL_AB", metadata$treatment)
metadata$treatment <- gsub("anti-Netrin-1 antibody", "NETRIN1_AB", metadata$treatment)
metadata$organism <- gsub("Mus musculus", "mus_musculus", metadata$organism)
metadata$sample <- c("CTL_A", "CTL_B", "CTL_C", "NP137_D", "NP137_E", "NP137_F")

# Save cleaned metadata
write.csv(
  metadata,
  here("data", "metadata", "GSE225688_PTEN_metadata.csv"),
  row.names = F
)
# Save renamed count data
write.csv(
  count,
  here("data", "bulk_rnaseq", "GSE225688_PTEN_Raw_count_cleaned.csv"),
)

