library(dplyr)
library(tidyr)
library(readr)
library(tidyverse)

# Create output directory
outdir <- "results/correlation"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Load the data
df_meth <- read.csv("data/wood_smoke_percent_methylation.csv", sep = ",")

# Fix column names (remove leading/trailing spaces)
colnames(df_meth) <- trimws(colnames(df_meth))

# Define sample columns (these are all between the 5th and 11th columns based on the preview)
sample_cols <- colnames(df_meth)[5:12]

# Remove rows where gene Symbol is NA
df_meth_filtered <- df_meth %>% filter(!is.na(gene.Symbol))

# Convert to long format: one row per sample-site-gene
df_meth_long <- df_meth_filtered %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "Sample", values_to = "Percent_Methylation")

# Group by sample and gene, average methylation across DMRs
df_meth_avg <- df_meth_long %>%
  group_by(Sample, gene.Symbol) %>%
  filter_all(all_vars(. != "")) %>%
  dplyr::summarise(Average_Methylation = mean(Percent_Methylation, na.rm = TRUE), .groups = "drop")

# Pivot to wide format: samples as rows, genes as columns
df_meth_wide <- df_meth_avg %>% 
  pivot_wider(names_from = gene.Symbol, values_from = Average_Methylation)

# Read in the expression data CSV file
df_exp <- read.csv("data/normalized_wood_smoke_counts.csv", sep = ",")

# Process the expression data:
df_exp_filtered <- df_exp %>%
  mutate(X = sub(".*?_", "", X)) %>%                 # Remove any prefix before underscore in 'X' column (to extract sample ID)
  filter(X %in% colnames(df_meth_wide)) %>%          # Keep only rows where sample IDs are present in methylation data
  group_by(X) %>%                                    # Group by sample ID
  dplyr::summarise(across(where(is.numeric),                # Compute the mean of all numeric expression columns
                          mean, na.rm = TRUE), 
                   .groups = "drop")                        # Do not retain grouping structure in the output

# Transpose the expression matrix: genes as columns, samples as rows
df_exp_wide <- t(df_exp_filtered)

# Set the first row as column names (sample IDs)
colnames(df_exp_wide) <- df_exp_wide[1,]

# Remove the first row now that it has been used as headers
df_exp_wide <- as.data.frame(df_exp_wide[-1,])

# Keep only samples whose names start with 'FA5D' or 'WS5D' (e.g., treatment groups)
df_exp_wide <- df_exp_wide[grepl("^(FA5D|WS5D)", rownames(df_exp_wide)), ]

# Convert all values in the data frame to numeric
df_exp_wide <- as.data.frame(lapply(df_exp_wide, as.numeric))

# Ensure methylation matrix is a data frame
df_meth_wide <- as.data.frame(df_meth_wide)

# Set sample IDs as row names (assumes first column contains sample IDs)
rownames(df_meth_wide) <- df_meth_wide[, 1]

# Remove the first column now that it’s used as row names
df_meth_wide <- df_meth_wide[, -1]

# Keep only the genes that are also in the expression dataset
df_meth_wide <- df_meth_wide[, colnames(df_meth_wide) %in% colnames(df_exp_wide)]

# Keep only the genes in the expression data that are also in the methylation data
df_exp_wide <- df_exp_wide[, colnames(df_exp_wide) %in% colnames(df_meth_wide)]

# Sanity check: Ensure column order and names are identical between the datasets
stopifnot(all(colnames(df_meth_wide) == colnames(df_exp_wide)))

# Combine methylation and expression columns side-by-side for each gene
meth_exp_correlation <- do.call(cbind, lapply(seq_along(df_meth_wide), function(i) {
  colname <- colnames(df_meth_wide)[i]
  data.frame(
    setNames(df_meth_wide[i], paste0(colname, "_M")),  # Methylation column
    setNames(df_exp_wide[i], paste0(colname, "_E"))    # Expression column
  )
}))

meth_exp_correlation <- do.call(cbind, lapply(seq_along(df_meth_wide), function(i) {
  colname <- colnames(df_meth_wide)[i]
  
  data.frame(
    setNames(df_meth_wide[i], paste0(colname, "_M")),
    setNames(df_exp_wide[i], paste0(colname, "_E"))
  )
}))

# Reorder rows to pair FA and WS samples
n_rows <- nrow(meth_exp_correlation)
half <- n_rows / 2

row_order <- as.vector(rbind(
  seq_len(half),
  seq_len(half) + half
))

meth_exp_correlation <- meth_exp_correlation[row_order, , drop = FALSE]

# Remove the first two rows (likely headers or invalid samples)
meth_exp_correlation <- meth_exp_correlation[-c(1, 2),]

# Initialize storage structures for correlation results
datalist <- list()
gene_rows <- c()

# Compute Spearman correlation for each gene’s methylation and expression
results <- for (i in seq(from = 1, to = ncol(meth_exp_correlation), by = 2)) {
  dat <- cor.test(meth_exp_correlation[, i], meth_exp_correlation[, i + 1], method = "spearman")
  dat$i <- i
  datalist[[i]] <- dat
  gene_rows <- c(gene_rows, sub("_.*", "", colnames(meth_exp_correlation[i])))  # Extract gene name
}

# Combine results into a single data frame
big_data <- do.call(rbind, datalist)
big_data_df <- data.frame(big_data)

# Extract p-values and correlation coefficients
p_val_spearman <- data.frame(cbind(big_data_df$p.value, big_data_df$estimate))
colnames(p_val_spearman) <- c("p_value", "spearman")

# Adjust p-values using FDR correction
fdr <- p.adjust(p_val_spearman$p_value, method = "fdr")

# Create final correlation result table
fdr_spearman <- data.frame(cbind(p_val_spearman$p_value, fdr, as.numeric(p_val_spearman$spearman)))
colnames(fdr_spearman) <- c("p_val", "fdr", "spearman")

# Convert all values to character (possibly for downstream formatting)
fdr_spearman <- apply(fdr_spearman, 2, as.character)
fdr_spearman <- as.data.frame(fdr_spearman)
rownames(fdr_spearman) <- gene_rows

# Optional: write to CSV
write.csv(
  fdr_spearman,
  file.path(outdir, "final_fdr_spearman_meth_exp_correlation_test.csv")
)

