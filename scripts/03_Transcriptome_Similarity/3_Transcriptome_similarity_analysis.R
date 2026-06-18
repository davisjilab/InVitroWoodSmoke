required_packages <- c(
  "dplyr", "tidyr", "readr", "stringr",
  "tibble", "ggplot2", "pheatmap", "cluster"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
  library(cluster)
})

input_dir <- "data/differential_expression"
outdir <- "results/transcriptome_similarity"
file_pattern <- "\\.csv$"
gene_column <- "X"
padj_column <- "padj"
lfc_column <- "log2FoldChange"
padj_threshold <- 0.05
absolute_lfc_threshold <- 0
gene_clusters_requested <- 2
heatmap_width <- 9
heatmap_height <- 8

# Contrast names to omit. Leave empty to keep all contrasts.
exclude_contrasts <- character(0)

# Optional regex containing one capture group used to derive contrast names.
# Leave NULL to use each filename without its extension.
contrast_regex <- NULL

# Optional interactive directory selection is not provided by base R.
# Enter the folder path above, or use RStudio's Files pane to copy its path.

if (!dir.exists(input_dir)) stop("Input directory not found: ", input_dir, call. = FALSE)
if (!is.numeric(padj_threshold) || length(padj_threshold) != 1L || is.na(padj_threshold)) stop("padj_threshold must be one number.", call. = FALSE)
if (!is.numeric(absolute_lfc_threshold) || absolute_lfc_threshold < 0) stop("absolute_lfc_threshold must be zero or greater.", call. = FALSE)


dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
table_dir <- file.path(outdir, "tables")
figure_dir <- file.path(outdir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(
  input_dir,
  pattern = file_pattern,
  full.names = TRUE
)

if (length(files) < 2) {
  stop("At least two CSV files are required.")
}

read_de_file <- function(path) {
  
  data <- readr::read_csv(
    path,
    show_col_types = FALSE,
    name_repair = "minimal"
  )
  
  # Rename a blank first column
  if (is.na(names(data)[1]) || trimws(names(data)[1]) == "") {
    names(data)[1] <- gene_column
  }
  
  if (!gene_column %in% names(data)) {
    stop(
      basename(path),
      " is missing the gene column: ",
      gene_column
    )
  }
  
  if (!padj_column %in% names(data)) {
    stop(
      basename(path),
      " is missing the adjusted P-value column: ",
      padj_column
    )
  }
  
  if (!lfc_column %in% names(data)) {
    data[[lfc_column]] <- NA_real_
  }
  
  contrast_name <- tools::file_path_sans_ext(
    basename(path)
  )
  
  data.frame(
    gene = as.character(data[[gene_column]]),
    padj = suppressWarnings(as.numeric(data[[padj_column]])),
    log2FoldChange = suppressWarnings(
      as.numeric(data[[lfc_column]])
    ),
    contrast = contrast_name,
    stringsAsFactors = FALSE
  )
}

de_list <- lapply(files, read_de_file)

de_results <- dplyr::bind_rows(de_list)

de_results <- de_results %>%
  dplyr::filter(
    !is.na(gene),
    gene != ""
  ) %>%
  dplyr::group_by(contrast, gene) %>%
  dplyr::summarise(
    padj = if (all(is.na(padj))) {
      NA_real_
    } else {
      min(padj, na.rm = TRUE)
    },
    log2FoldChange = if (all(is.na(log2FoldChange))) {
      NA_real_
    } else {
      log2FoldChange[
        which.min(replace(padj, is.na(padj), Inf))[1]
      ]
    },
    .groups = "drop"
  )

print(basename(files))
print(unique(de_results$contrast))
print(table(de_results$contrast))

excluded <- unique(trimws(as.character(exclude_contrasts)))
excluded <- excluded[!is.na(excluded) & nzchar(excluded)]
if (length(excluded) > 0L) {
  de_results <- de_results %>% filter(!contrast %in% excluded)
}

if (n_distinct(de_results$contrast) < 2L) {
  stop("Fewer than two contrasts remain after exclusion.", call. = FALSE)
}

de_results <- de_results %>%
  mutate(
    significant = !is.na(padj) &
      padj < padj_threshold &
      (
        absolute_lfc_threshold <= 0 |
          (!is.na(log2FoldChange) & abs(log2FoldChange) >= absolute_lfc_threshold)
      )
  )

all_significant_genes <- de_results %>%
  filter(significant) %>%
  distinct(gene) %>%
  pull(gene)

if (length(all_significant_genes) == 0L) {
  stop("No genes pass the requested DEG thresholds.", call. = FALSE)
}

binary_matrix <- de_results %>%
  filter(gene %in% all_significant_genes) %>%
  select(contrast, gene, significant) %>%
  complete(contrast, gene = all_significant_genes, fill = list(significant = FALSE)) %>%
  mutate(significant = as.integer(significant)) %>%
  pivot_wider(names_from = gene, values_from = significant) %>%
  column_to_rownames("contrast") %>%
  as.matrix()

# Jaccard(A,B) = |A intersection B| / |A union B|.
jaccard_pair <- function(x, y) {
  union_size <- sum(x == 1L | y == 1L)
  if (union_size == 0L) return(NA_real_)
  sum(x == 1L & y == 1L) / union_size
}

overlap_pair <- function(x, y) sum(x == 1L & y == 1L)

contrasts <- rownames(binary_matrix)
jaccard_matrix <- outer(
  seq_along(contrasts),
  seq_along(contrasts),
  Vectorize(function(i, j) jaccard_pair(binary_matrix[i, ], binary_matrix[j, ]))
)
overlap_matrix <- outer(
  seq_along(contrasts),
  seq_along(contrasts),
  Vectorize(function(i, j) overlap_pair(binary_matrix[i, ], binary_matrix[j, ]))
)
dimnames(jaccard_matrix) <- list(contrasts, contrasts)
dimnames(overlap_matrix) <- list(contrasts, contrasts)

pairwise_results <- expand_grid(contrast_1 = contrasts, contrast_2 = contrasts) %>%
  rowwise() %>%
  mutate(
    jaccard_similarity = jaccard_matrix[contrast_1, contrast_2],
    overlap_genes = overlap_matrix[contrast_1, contrast_2],
    union_genes = sum(
      binary_matrix[contrast_1, ] == 1L |
        binary_matrix[contrast_2, ] == 1L
    )
  ) %>%
  ungroup()

deg_counts <- tibble(
  contrast = contrasts,
  significant_genes = rowSums(binary_matrix)
)

write_csv(de_results, file.path(table_dir, "combined_differential_expression_results.csv"))
write_csv(as.data.frame(binary_matrix) %>% rownames_to_column("contrast"),
          file.path(table_dir, "binary_deg_matrix.csv"))
write_csv(as.data.frame(jaccard_matrix) %>% rownames_to_column("contrast"),
          file.path(table_dir, "jaccard_similarity_matrix.csv"))
write_csv(as.data.frame(overlap_matrix) %>% rownames_to_column("contrast"),
          file.path(table_dir, "overlap_count_matrix.csv"))
write_csv(pairwise_results, file.path(table_dir, "pairwise_similarity_long.csv"))
write_csv(deg_counts, file.path(table_dir, "deg_counts.csv"))

labels <- matrix(
  paste0(
    overlap_matrix,
    "\n(",
    sprintf("%.2f", jaccard_matrix),
    ")"
  ),
  nrow = nrow(jaccard_matrix),
  dimnames = dimnames(jaccard_matrix)
)

png(
  file.path(figure_dir, "jaccard_similarity_heatmap.png"),
  width = heatmap_width,
  height = heatmap_height,
  units = "in",
  res = 300
)

pheatmap::pheatmap(
  jaccard_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = labels,
  number_format = "%s",
  border_color = NA,
  main = "Transcriptomic similarity\nOverlap count (Jaccard similarity)"
)

dev.off()

# Cluster genes by their binary significance pattern across contrasts.
gene_pattern_matrix <- t(binary_matrix)
informative <- rowSums(gene_pattern_matrix) > 0L &
  apply(gene_pattern_matrix, 1, function(x) length(unique(x)) > 1L)

gene_clusters <- tibble(gene = rownames(gene_pattern_matrix), cluster = NA_integer_)
ordered_genes <- rownames(gene_pattern_matrix)

if (sum(informative) >= 2L) {
  gene_distance <- dist(gene_pattern_matrix[informative, , drop = FALSE], method = "manhattan")
  gene_hclust <- hclust(gene_distance, method = "complete")
  requested_k <- max(1L, min(gene_clusters_requested, sum(informative)))
  assignments <- cutree(gene_hclust, k = requested_k)
  gene_clusters$cluster[match(names(assignments), gene_clusters$gene)] <- assignments

  ordered_informative <- names(assignments)[order(assignments, gene_hclust$order)]
  ordered_genes <- c(
    ordered_informative,
    setdiff(rownames(gene_pattern_matrix), ordered_informative)
  )

  capture.output(gene_hclust, file = file.path(outdir, "gene_hclust_summary.txt"))
}

write_csv(gene_clusters, file.path(table_dir, "gene_pattern_clusters.csv"))

ordered_binary <- binary_matrix[, ordered_genes, drop = FALSE]
row_labels <- paste0(rownames(ordered_binary), " (", rowSums(ordered_binary), " DEGs)")

png(
  file.path(figure_dir, "binary_deg_pattern_heatmap.png"),
  width = heatmap_width,
  height = heatmap_height,
  units = "in",
  res = 300
)
pheatmap(
  ordered_binary,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  labels_row = row_labels,
  legend = FALSE,
  border_color = NA,
  main = sprintf(
    "Differentially expressed genes (adjusted P < %s%s)",
    padj_threshold,
    if (absolute_lfc_threshold > 0)
      paste0("; |log2FC| ≥ ", absolute_lfc_threshold)
    else ""
  )
)
dev.off()