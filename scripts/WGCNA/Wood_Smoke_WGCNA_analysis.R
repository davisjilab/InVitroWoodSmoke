options(stringsAsFactors = FALSE)
set.seed(20240109)

counts_file <- "data/wood_smoke_exposure_rnaseq_counts.csv"
output_dir <- "results/wgcna"

exclude_samples <- c("FA5D_1")
selected_modules <- c("darkorange", "darkorange2", "lightyellow")

soft_threshold_power <- "auto"  # Use "auto" or a positive number, such as 13
variance_percentile <- 0.25
minimum_expression <- 1
network_type <- "signed hybrid"
correlation_method <- "pearson"
number_of_threads <- 2
run_enrichment <- FALSE

# Validate settings and create output folders.
if (!file.exists(counts_file)) {
  stop(
    "Count file not found: ", counts_file,
    "\nUpdate counts_file in the User settings section.",
    call. = FALSE
  )
}

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

fig_dir <- file.path(output_dir, "figures")
tab_dir <- file.path(output_dir, "tables")
obj_dir <- file.path(output_dir, "objects")
for (d in c(fig_dir, tab_dir, obj_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

selected_modules <- trimws(selected_modules)
selected_modules <- selected_modules[nzchar(selected_modules)]
exclude_samples <- trimws(exclude_samples)
exclude_samples <- exclude_samples[nzchar(exclude_samples)]

# ------------------------------- Dependencies ---------------------------------
required_packages <- c(
  "BioNERO", "DESeq2", "SummarizedExperiment", "S4Vectors", "WGCNA",
  "limma", "ggplot2", "dplyr", "tibble", "tidyr", "patchwork",
  "igraph", "ggraph", "scales"
)
optional_enrichment_packages <- c(
  "biomaRt", "clusterProfiler", "org.Hs.eg.db", "ReactomePA", "enrichplot"
)

missing_required <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required)) {
  stop(
    "Missing required packages: ", paste(missing_required, collapse = ", "),
    "\nInstall them before running the analysis.", call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(BioNERO)
  library(DESeq2)
  library(SummarizedExperiment)
  library(WGCNA)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(igraph)
  library(ggraph)
})

WGCNA::allowWGCNAThreads(nThreads = number_of_threads)

# --------------------------------- Helpers -------------------------------------
message_time <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
}

save_plot <- function(plot, filename, width, height, dpi = 300) {
  output_file <- file.path(fig_dir, filename)
  extension <- tolower(tools::file_ext(output_file))
  
  if (inherits(plot, c("Heatmap", "HeatmapList", "AdditiveUnit"))) {
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      stop(
        "ComplexHeatmap is required to save BioNERO heatmaps.",
        call. = FALSE
      )
    }
    
    if (extension == "pdf") {
      grDevices::pdf(
        output_file,
        width = width,
        height = height,
        onefile = FALSE
      )
    } else if (extension %in% c("jpg", "jpeg")) {
      grDevices::jpeg(
        output_file,
        width = width,
        height = height,
        units = "in",
        res = dpi,
        quality = 100,
        bg = "white"
      )
    } else {
      grDevices::png(
        output_file,
        width = width,
        height = height,
        units = "in",
        res = dpi,
        bg = "white"
      )
    }
    
    on.exit(grDevices::dev.off(), add = TRUE)
    ComplexHeatmap::draw(plot)
    
    return(invisible(output_file))
  }
  
  ggplot2::ggsave(
    filename = output_file,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  
  invisible(output_file)
}

write_table <- function(x, filename) {
  utils::write.csv(x, file.path(tab_dir, filename), row.names = FALSE, quote = TRUE)
}

parse_sample_metadata <- function(sample_names) {
  pattern <- "^(FA|WS)(1D|5D|R)_([[:alnum:]-]+)$"
  matched <- regexec(pattern, sample_names)
  pieces <- regmatches(sample_names, matched)
  valid <- lengths(pieces) == 4L

  if (!all(valid)) {
    stop(
      "The following sample names do not match ^(FA|WS)(1D|5D|R)_ID$: ",
      paste(sample_names[!valid], collapse = ", "), call. = FALSE
    )
  }

  condition_code <- vapply(pieces, `[[`, character(1), 2L)
  time_code <- vapply(pieces, `[[`, character(1), 3L)
  animal_id <- vapply(pieces, `[[`, character(1), 4L)

  metadata <- data.frame(
    sample = sample_names,
    condition = factor(condition_code, levels = c("FA", "WS"), labels = c("Filtered_air", "Wood_smoke")),
    time = factor(time_code, levels = c("1D", "5D", "R"), labels = c("Day_1", "Day_5", "Recovery")),
    animal_id = factor(animal_id),
    stringsAsFactors = FALSE,
    row.names = sample_names
  )
  metadata$condition_time <- interaction(metadata$condition, metadata$time, sep = "_", drop = TRUE)
  metadata
}

format_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, format(p, scientific = TRUE, digits = 2), sprintf("%.3f", p)))
}

plot_trait_heatmap <- function(me_trait, group_name, base_size = 14, text_size = 3.5) {
  dat <- dplyr::filter(me_trait, .data$group == group_name)
  if (!nrow(dat)) stop("No module-trait results found for group: ", group_name, call. = FALSE)

  dat$label <- paste0(sprintf("%.2f", dat$cor), "\n(", format_p(dat$pvalue), ")")
  ggplot(dat, aes(x = trait, y = ME, fill = cor)) +
    geom_tile(color = "grey85", linewidth = 0.2) +
    geom_text(aes(label = label), size = text_size) +
    scale_fill_gradient2(
      low = "blue3", mid = "white", high = "red3",
      midpoint = 0, limits = c(-1, 1), name = "Correlation"
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black")
    )
}

scale01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (!all(is.finite(r)) || diff(r) == 0) return(rep(0, length(x)))
  (x - r[1]) / diff(r)
}

build_module_network <- function(
  dat_expr, module_colors, module, power,
  n_inner = 10, n_outer = 15, n_edges = 300,
  network_type = "signed hybrid", tom_type = "signed"
) {
  overlap <- intersect(colnames(dat_expr), names(module_colors))
  dat_expr <- dat_expr[, overlap, drop = FALSE]
  module_colors <- module_colors[overlap]
  genes <- names(module_colors)[module_colors == module]

  if (length(genes) < 5L) {
    warning("Skipping network for ", module, ": fewer than five genes.")
    return(NULL)
  }

  dat_mod <- dat_expr[, genes, drop = FALSE]
  me <- WGCNA::moduleEigengenes(dat_mod, colors = rep(module, ncol(dat_mod)))$eigengenes[, 1]
  kme_matrix <- suppressWarnings(
    stats::cor(
      dat_mod,
      me,
      use = "pairwise.complete.obs",
      method = "pearson"
    )
  )
  
  kme <- abs(kme_matrix[, 1])
  names(kme) <- rownames(kme_matrix)
  
  if (is.null(names(kme)) || any(!nzchar(names(kme)))) {
    names(kme) <- colnames(dat_mod)
  }
  
  kme <- kme[
    is.finite(kme) &
      !is.na(names(kme)) &
      nzchar(names(kme))
  ]
  
  ranked_genes <- names(
    sort(kme, decreasing = TRUE)
  )
  
  hub_genes <- head(
    ranked_genes,
    min(length(ranked_genes), n_inner + n_outer)
  )
  
  message(
    "Module: ", module,
    " | genes before ranking: ", ncol(dat_mod),
    " | genes with finite kME: ", length(kme),
    " | selected hub genes: ", length(hub_genes)
  )

  dat_hub <- dat_mod[
    ,
    hub_genes,
    drop = FALSE
  ]
  
  tom <- WGCNA::TOMsimilarityFromExpr(
    dat_hub, power = power, networkType = network_type, TOMType = tom_type,
    verbose = 0
  )
  diag(tom) <- 0
  dimnames(tom) <- list(hub_genes, hub_genes)

  upper <- which(upper.tri(tom), arr.ind = TRUE)
  weights <- tom[upper]
  keep <- order(weights, decreasing = TRUE)[seq_len(min(n_edges, length(weights)))]
  upper <- upper[keep, , drop = FALSE]

  edge_df <- data.frame(
    from = rownames(tom)[upper[, 1]],
    to = colnames(tom)[upper[, 2]],
    weight = tom[upper],
    stringsAsFactors = FALSE
  )
  edge_df <- edge_df[is.finite(edge_df$weight) & edge_df$weight > 0, , drop = FALSE]
  if (!nrow(edge_df)) return(NULL)

  graph <- igraph::graph_from_data_frame(
    edge_df,
    vertices = data.frame(name = hub_genes, stringsAsFactors = FALSE),
    directed = FALSE
  )

  n_inner_eff <- min(n_inner, length(hub_genes))
  n_outer_eff <- length(hub_genes) - n_inner_eff
  theta_inner <- seq(0, 2 * pi, length.out = n_inner_eff + 1)[-(n_inner_eff + 1)]
  theta_outer <- if (n_outer_eff > 0) seq(0, 2 * pi, length.out = n_outer_eff + 1)[-(n_outer_eff + 1)] else numeric()
  coords <- rbind(
    cbind(0.55 * cos(theta_inner), 0.55 * sin(theta_inner)),
    if (n_outer_eff > 0) cbind(cos(theta_outer), sin(theta_outer)) else NULL
  )
  rownames(coords) <- hub_genes
  coords <- coords[igraph::V(graph)$name, , drop = FALSE]

  list(graph = graph, coordinates = coords, edges = edge_df, genes = hub_genes)
}

plot_module_network <- function(network, module) {
  coords <- as.data.frame(network$coordinates)
  colnames(coords) <- c("x", "y")
  coords$name <- rownames(coords)

  valid_color <- tryCatch({ grDevices::col2rgb(module); TRUE }, error = function(e) FALSE)
  module_color <- if (valid_color) module else "grey70"

  ggraph(network$graph, layout = "manual", x = coords$x, y = coords$y) +
    geom_edge_link(
      aes(edge_alpha = weight), edge_colour = "grey55", show.legend = FALSE
    ) +
    geom_node_point(shape = 21, fill = module_color, color = "black", size = 5) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3.5) +
    scale_edge_alpha_continuous(range = c(0.05, 0.75)) +
    labs(title = paste(module, "module")) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

map_macaque_to_human <- function(ensembl_ids) {
  if (!all(vapply(optional_enrichment_packages, requireNamespace, logical(1), quietly = TRUE))) {
    warning("Optional enrichment packages are unavailable; enrichment will be skipped.")
    return(NULL)
  }

  clean_ids <- unique(sub("\\..*$", "", sub("_.*$", "", ensembl_ids)))
  macaque <- biomaRt::useEnsembl(biomart = "genes", dataset = "mmulatta_gene_ensembl")
  mapping <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_associated_gene_name"),
    filters = "ensembl_gene_id",
    values = clean_ids,
    mart = macaque
  )
  mapping <- mapping[nzchar(mapping$hsapiens_homolog_associated_gene_name), , drop = FALSE]
  unique(mapping)
}

run_module_enrichment <- function(module_gene_ids, background_gene_ids, module_name) {
  mapping <- map_macaque_to_human(c(module_gene_ids, background_gene_ids))
  if (is.null(mapping) || !nrow(mapping)) return(NULL)

  module_clean <- unique(sub("\\..*$", "", sub("_.*$", "", module_gene_ids)))
  background_clean <- unique(sub("\\..*$", "", sub("_.*$", "", background_gene_ids)))

  module_symbols <- unique(mapping$hsapiens_homolog_associated_gene_name[
    mapping$ensembl_gene_id %in% module_clean
  ])
  background_symbols <- unique(mapping$hsapiens_homolog_associated_gene_name[
    mapping$ensembl_gene_id %in% background_clean
  ])

  if (length(module_symbols) < 5L) {
    warning("Too few mapped genes for enrichment of module ", module_name)
    return(NULL)
  }

  go <- clusterProfiler::enrichGO(
    gene = module_symbols,
    universe = background_symbols,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )

  entrez <- clusterProfiler::bitr(
    module_symbols, fromType = "SYMBOL", toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )
  reactome <- if (nrow(entrez)) {
    ReactomePA::enrichPathway(
      gene = unique(entrez$ENTREZID), organism = "human",
      pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
      readable = TRUE
    )
  } else NULL

  list(mapping = mapping, go = go, reactome = reactome)
}

# ------------------------------ 1. Read counts ---------------------------------
count_df <- utils::read.csv(counts_file, check.names = FALSE)
if (ncol(count_df) < 3L) stop("Count file must contain a gene column and at least two samples.")

rownames(count_df) <- make.unique(as.character(count_df[[1]]))
counts <- as.matrix(count_df[, -1, drop = FALSE])
storage.mode(counts) <- "numeric"

if (any(!is.finite(counts))) stop("Count matrix contains non-finite values.")
if (any(counts < 0)) stop("Count matrix contains negative values.")
if (anyDuplicated(colnames(counts))) stop("Sample names are duplicated.")

exclude_samples <- intersect(exclude_samples, colnames(counts))
if (length(exclude_samples)) {
  counts <- counts[, setdiff(colnames(counts), exclude_samples), drop = FALSE]
}

metadata <- parse_sample_metadata(colnames(counts))

metadata$sample_number <- paste0(
  "S",
  sub(".*_", "", rownames(metadata))
)

write_table(tibble::rownames_to_column(metadata, "sample_id"), "sample_metadata.csv")

# --------------------- 2. VST and paired-sample correction ---------------------
count_matrix <- round(as.matrix(counts))
storage.mode(count_matrix) <- "integer"

expr <- DESeq2::vst(
  count_matrix,
  blind = TRUE
)

expr <- as.matrix(expr)
storage.mode(expr) <- "double"

expr_corrected <- limma::removeBatchEffect(
  expr,
  batch = metadata$sample_number
)

bionero_metadata <- metadata[
  ,
  c("condition", "time", "condition_time"),
  drop = FALSE
]

se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(vst_batch_corrected = expr_corrected),
  colData = S4Vectors::DataFrame(bionero_metadata)
)

# --------------------------- 3. Expression filtering ---------------------------
exp_filt <- BioNERO::replace_na(se)
exp_filt <- BioNERO::remove_nonexp(
  exp_filt, method = "median", min_exp = minimum_expression
)
exp_filt <- BioNERO::filter_by_variance(
  exp_filt, percentile = variance_percentile
)

filtered_expr <- SummarizedExperiment::assay(exp_filt)
utils::write.csv(
  data.frame(gene_id = rownames(filtered_expr), filtered_expr, check.names = FALSE),
  file.path(tab_dir, "filtered_vst_expression.csv"), row.names = FALSE
)
saveRDS(exp_filt, file.path(obj_dir, "filtered_expression_se.rds"))

variance_df <- data.frame(variance = apply(filtered_expr, 1, stats::var))
p_variance <- ggplot(variance_df, aes(x = variance)) +
  geom_histogram(bins = 50, color = "white") +
  labs(x = "Gene-expression variance", y = "Number of genes") +
  theme_classic(base_size = 12)
save_plot(p_variance, "gene_variance_distribution.png", 7, 5)

p_sample_cor <- BioNERO::plot_heatmap(exp_filt, type = "samplecor", show_rownames = FALSE)
save_plot(p_sample_cor, "sample_correlation_heatmap.png", 8, 7)

p_pca <- BioNERO::plot_PCA(
  SummarizedExperiment::assay(exp_filt),
  metadata[
    colnames(SummarizedExperiment::assay(exp_filt)),
    c("animal_id", "condition_time"),
    drop = FALSE
  ]
)
save_plot(p_pca, "sample_PCA.png", 8, 6)

# ----------------------- 4. Soft-threshold selection ---------------------------
sft <- BioNERO::SFT_fit(
  exp_filt,
  net_type = network_type,
  cor_method = correlation_method
)

power <- if (tolower(soft_threshold_power) == "auto") {
  sft$power
} else {
  as.numeric(soft_threshold_power)
}
if (!is.finite(power) || power <= 0) stop("Invalid soft-threshold power: ", soft_threshold_power)

save_plot(sft$plot, "scale_free_topology_fit.png", 8, 6)
write_table(data.frame(selected_power = power), "selected_soft_threshold.csv")

# -------------------------- 5. Network construction ----------------------------
net <- BioNERO::exp2gcn(
  exp_filt,
  net_type = network_type,
  SFTpower = power,
  cor_method = correlation_method
)
saveRDS(net, file.path(obj_dir, "wgcna_network.rds"))

module_assignments <- net$genes_and_modules
colnames(module_assignments)[1:2] <- c("gene_id", "module")
write_table(module_assignments, "module_assignments.csv")

p_dendro <- BioNERO::plot_dendro_and_colors(net)
p_module_size <- BioNERO::plot_ngenes_per_module(net)
save_plot(p_dendro, "module_dendrogram.png", 12, 7)
save_plot(p_module_size, "genes_per_module.png", 8, 6)

# ---------------------- 6. Module-trait associations ---------------------------
me_trait <- BioNERO::module_trait_cor(exp = exp_filt, MEs = net$MEs)
write_table(me_trait, "module_trait_correlations.csv")

p_condition <- plot_trait_heatmap(me_trait, "condition")
p_time <- plot_trait_heatmap(me_trait, "time")
p_condition_time <- plot_trait_heatmap(me_trait, "condition_time")

save_plot(p_condition, "module_trait_condition.png", 10, 8)
save_plot(p_time, "module_trait_time.png", 10, 8)
save_plot(p_condition_time, "module_trait_condition_time.png", 12, 8)

main_summary <- (sft$plot | p_module_size) / (p_condition | p_condition_time)
save_plot(main_summary, "wgcna_manuscript_summary.png", 18, 12)
save_plot(main_summary, "wgcna_manuscript_summary.pdf", 18, 12)

# ----------------------------- 7. Hub genes -------------------------------------
hubs <- BioNERO::get_hubs_gcn(exp_filt, net)
write_table(hubs, "hub_genes_all_modules.csv")

for (module in selected_modules) {
  module_hubs <- dplyr::filter(hubs, .data$Module == module)
  if (nrow(module_hubs)) write_table(module_hubs, paste0(module, "_hub_genes.csv"))
}

# ---------------------- 8. Selected module profiles ----------------------------
condition_time_results <- dplyr::filter(me_trait, .data$group == "condition_time")
for (module in selected_modules) {
  me_name <- paste0("ME", module)
  module_results <- dplyr::filter(condition_time_results, .data$ME == me_name)
  if (!nrow(module_results)) next

  module_results$trait <- factor(
    module_results$trait,
    levels = c(
      "Filtered_air_Day_1", "Filtered_air_Day_5", "Filtered_air_Recovery",
      "Wood_smoke_Day_1", "Wood_smoke_Day_5", "Wood_smoke_Recovery"
    )
  )
  module_results$exposure <- ifelse(grepl("^Filtered_air", module_results$trait), "Filtered air", "Wood smoke")

  valid_color <- tryCatch({ grDevices::col2rgb(module); TRUE }, error = function(e) FALSE)
  module_color <- if (valid_color) module else "black"

  p_profile <- ggplot(module_results, aes(x = trait, y = cor, group = exposure)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
    geom_line(linewidth = 1.1, color = module_color) +
    geom_point(aes(shape = pvalue <= 0.05), size = 3.5, color = module_color) +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 19), name = "P ≤ 0.05") +
    scale_y_continuous(limits = c(-1, 1)) +
    labs(x = NULL, y = "Module-trait correlation", title = paste(module, "module")) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_plot(p_profile, paste0(module, "_trait_profile.png"), 9, 6)
}

# ------------------------ 9. Module network figures ----------------------------
dat_expr <- t(SummarizedExperiment::assay(exp_filt))
module_colors <- setNames(module_assignments$module, module_assignments$gene_id)
network_plots <- list()

for (module in selected_modules) {
  network <- build_module_network(
    dat_expr = dat_expr,
    module_colors = module_colors,
    module = module,
    power = power,
    network_type = network_type
  )
  if (is.null(network)) next

  write_table(network$edges, paste0(module, "_network_edges.csv"))
  network_plots[[module]] <- plot_module_network(network, module)
  save_plot(network_plots[[module]], paste0(module, "_hub_network.png"), 7, 7)
  save_plot(network_plots[[module]], paste0(module, "_hub_network.pdf"), 7, 7)
}

# ---------------------- 10. Optional functional enrichment ---------------------
if (isTRUE(run_enrichment)) {
  background_ids <- rownames(exp_filt)

  for (module in selected_modules) {
    module_ids <- module_assignments$gene_id[module_assignments$module == module]
    if (!length(module_ids)) next

    enrichment <- tryCatch(
      run_module_enrichment(module_ids, background_ids, module),
      error = function(e) {
        warning("Enrichment failed for ", module, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(enrichment)) next

    write_table(enrichment$mapping, paste0(module, "_macaque_to_human_mapping.csv"))

    if (!is.null(enrichment$go)) {
      go_df <- as.data.frame(enrichment$go)
      write_table(go_df, paste0(module, "_GO_BP_enrichment.csv"))
      if (nrow(go_df)) {
        p_go <- enrichplot::dotplot(enrichment$go, showCategory = min(15, nrow(go_df))) +
          ggtitle(paste("GO biological process:", module))
        save_plot(p_go, paste0(module, "_GO_BP_dotplot.png"), 9, 7)
      }
    }

    if (!is.null(enrichment$reactome)) {
      reactome_df <- as.data.frame(enrichment$reactome)
      write_table(reactome_df, paste0(module, "_Reactome_enrichment.csv"))
      if (nrow(reactome_df)) {
        p_reactome <- enrichplot::dotplot(enrichment$reactome, showCategory = min(15, nrow(reactome_df))) +
          ggtitle(paste("Reactome pathways:", module))
        save_plot(p_reactome, paste0(module, "_Reactome_dotplot.png"), 9, 7)
      }
    }
  }
}