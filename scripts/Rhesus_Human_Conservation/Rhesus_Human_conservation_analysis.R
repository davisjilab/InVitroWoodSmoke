required_packages <- c(
  "dplyr", "readr", "purrr", "stringr",
  "tibble", "ggplot2", "biomaRt", "GenomicRanges",
  "IRanges", "rtracklayer"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install required CRAN/Bioconductor packages before running:\n  ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(biomaRt)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
})

mode <- "both"  # "orthology", "liftover", or "both"
outdir <- "results/conservation"

# Orthology settings
orthology_input <- "data/orthology_inputs"  # one CSV file or a directory of CSV files
gene_column <- "gene"
id_type <- "ensembl_gene_id"  # or "external_gene_name"
strip_gene_suffix_regex <- "_.*$"  # use "" to disable suffix removal
minimum_percent_identity <- 0
ensembl_host <- "https://www.ensembl.org"

# liftOver settings
interval_input <- "data/rhesus_intervals.csv"
chain_file <- "data/rheMac10ToHg38.over.chain.gz"
chrom_column <- "chr"
start_column <- "start"
end_column <- "end"
coordinate_system <- "one_based_closed"  # or "zero_based_half_open"
mapping_selection <- "longest"  # or "all"

# Optional interactive file selection examples:
# orthology_input <- file.choose()
# interval_input <- file.choose()
# chain_file <- file.choose()

valid_modes <- c("orthology", "liftover", "both")
if (!mode %in% valid_modes) stop("mode must be one of: ", paste(valid_modes, collapse = ", "), call. = FALSE)
if (!id_type %in% c("ensembl_gene_id", "external_gene_name")) stop("Invalid id_type.", call. = FALSE)
if (!coordinate_system %in% c("one_based_closed", "zero_based_half_open")) stop("Invalid coordinate_system.", call. = FALSE)
if (!mapping_selection %in% c("longest", "all")) stop("Invalid mapping_selection.", call. = FALSE)


dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
table_dir <- file.path(outdir, "tables")
figure_dir <- file.path(outdir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

run_orthology <- mode %in% c("orthology", "both")
run_liftover <- mode %in% c("liftover", "both")

if (run_orthology && (!nzchar(orthology_input) || (!file.exists(orthology_input) && !dir.exists(orthology_input)))) {
  stop("Set orthology_input to an existing CSV file or directory.", call. = FALSE)
}
if (run_liftover && (!file.exists(interval_input) || !file.exists(chain_file))) {
  stop("Set interval_input and chain_file to existing files.", call. = FALSE)
}

orthology_rhesus_to_human <- function(
    genes,
    id_type,
    minimum_percent_identity = 0,
    host = "https://www.ensembl.org") {

  genes <- unique(genes[!is.na(genes) & genes != ""])
  if (length(genes) == 0L) stop("No valid gene identifiers were supplied.", call. = FALSE)

  ensembl_mart <- useEnsembl(
    biomart = "genes",
    host = host
  )
  
  available_datasets <- listDatasets(ensembl_mart)
  
  rhesus_matches <- available_datasets[
    grepl(
      "Macaca mulatta|macaque|mulatta",
      paste(
        available_datasets$dataset,
        available_datasets$description
      ),
      ignore.case = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  if (nrow(rhesus_matches) == 0L) {
    stop(
      "No rhesus macaque dataset was found on the selected Ensembl host: ",
      host,
      call. = FALSE
    )
  }
  
  if ("mmulatta_gene_ensembl" %in% rhesus_matches$dataset) {
    rhesus_dataset <- "mmulatta_gene_ensembl"
  } else {
    rhesus_dataset <- rhesus_matches$dataset[1]
  }
  
  rhesus_mart <- useDataset(
    dataset = rhesus_dataset,
    mart = ensembl_mart
  )
  attributes <- c(
    "ensembl_gene_id",
    "external_gene_name",
    "hsapiens_homolog_ensembl_gene",
    "hsapiens_homolog_associated_gene_name",
    "hsapiens_homolog_orthology_type",
    "hsapiens_homolog_perc_id",
    "hsapiens_homolog_orthology_confidence",
    "hsapiens_homolog_wga_coverage"
  )

  mapping <- getBM(
    attributes = attributes,
    filters = id_type,
    values = genes,
    mart = rhesus_mart
  ) %>%
    as_tibble() %>%
    rename(
      rhesus_ensembl = ensembl_gene_id,
      rhesus_symbol = external_gene_name,
      human_ensembl = hsapiens_homolog_ensembl_gene,
      human_symbol = hsapiens_homolog_associated_gene_name,
      orthology_type = hsapiens_homolog_orthology_type,
      protein_percent_identity = hsapiens_homolog_perc_id,
      orthology_confidence = hsapiens_homolog_orthology_confidence,
      whole_genome_alignment_coverage = hsapiens_homolog_wga_coverage
    ) %>%
    mutate(
      protein_percent_identity = as.numeric(protein_percent_identity),
      orthology_confidence = as.integer(orthology_confidence),
      input_id = if (id_type == "ensembl_gene_id") rhesus_ensembl else rhesus_symbol
    ) %>%
    filter(
      is.na(protein_percent_identity) |
        protein_percent_identity >= minimum_percent_identity
    )

  mapped_inputs <- unique(mapping$input_id[mapping$input_id != ""])
  unmapped <- setdiff(genes, mapped_inputs)

  # Explicitly retain unmapped inputs in the exported table.
  unmapped_rows <- tibble(
    rhesus_ensembl = if (id_type == "ensembl_gene_id") unmapped else NA_character_,
    rhesus_symbol = if (id_type == "external_gene_name") unmapped else NA_character_,
    human_ensembl = NA_character_,
    human_symbol = NA_character_,
    orthology_type = NA_character_,
    protein_percent_identity = NA_real_,
    orthology_confidence = NA_integer_,
    whole_genome_alignment_coverage = NA_real_,
    input_id = unmapped
  )

  bind_rows(mapping, unmapped_rows)
}

summarize_orthology <- function(mapping, comparison) {
  input_ids <- unique(mapping$input_id)
  mapped_ids <- mapping %>%
    filter(!is.na(human_ensembl), human_ensembl != "") %>%
    distinct(input_id) %>%
    pull(input_id)
  one_to_one_ids <- mapping %>%
    filter(!is.na(orthology_type), orthology_type == "ortholog_one2one") %>%
    distinct(input_id) %>%
    pull(input_id)
  high_confidence_ids <- mapping %>%
    filter(
      !is.na(orthology_type),
      orthology_type == "ortholog_one2one",
      !is.na(orthology_confidence),
      orthology_confidence == 1L
    ) %>%
    distinct(input_id) %>%
    pull(input_id)

  tibble(
    comparison = comparison,
    total_input_genes = length(input_ids),
    mapped_any = length(mapped_ids),
    one_to_one = length(one_to_one_ids),
    high_confidence_one_to_one = length(high_confidence_ids)
  ) %>%
    mutate(
      percent_mapped_any = 100 * mapped_any / total_input_genes,
      percent_one_to_one = 100 * one_to_one / total_input_genes,
      percent_high_confidence_one_to_one =
        100 * high_confidence_one_to_one / total_input_genes
    )
}

if (run_orthology) {
  orthology_files <- if (dir.exists(orthology_input)) {
    list.files(orthology_input, pattern = "\\.csv$", full.names = TRUE)
  } else {
    orthology_input
  }

  if (length(orthology_files) == 0L) {
    stop("No orthology input CSV files were found.", call. = FALSE)
  }

  mapped_tables <- map(orthology_files, function(path) {
    data <- read_csv(
      path,
      show_col_types = FALSE,
      name_repair = "minimal"
    )
    
    if (is.na(names(data)[1]) || trimws(names(data)[1]) == "") {
      names(data)[1] <- gene_column
    }
    
    if (!gene_column %in% names(data)) {
      stop(
        basename(path),
        " is missing gene column: ",
        gene_column,
        call. = FALSE
      )
    }
    
    genes <- as.character(data[[gene_column]])

    genes <- as.character(data[[gene_column]])
    if (nzchar(strip_gene_suffix_regex)) {
      genes <- sub(strip_gene_suffix_regex, "", genes)
    }

    comparison <- tools::file_path_sans_ext(basename(path))
    mapping <- orthology_rhesus_to_human(
      genes,
      id_type = id_type,
      minimum_percent_identity = minimum_percent_identity,
      host = ensembl_host
    ) %>%
      mutate(comparison = comparison)

    write_csv(
      mapping,
      file.path(table_dir, paste0(comparison, "_orthology_mapping.csv"))
    )
    mapping
  })

  all_mappings <- bind_rows(mapped_tables)
  orthology_summary <- map_dfr(
    split(all_mappings, all_mappings$comparison),
    ~ summarize_orthology(.x, unique(.x$comparison))
  )

  write_csv(all_mappings, file.path(table_dir, "all_orthology_mappings.csv"))
  write_csv(orthology_summary, file.path(table_dir, "orthology_summary.csv"))

  mapping_plot <- ggplot(
    orthology_summary,
    aes(x = reorder(comparison, percent_high_confidence_one_to_one),
        y = percent_high_confidence_one_to_one)
  ) +
    geom_col(width = 0.65) +
    geom_text(
      aes(label = sprintf(
        "%.1f%% (%d/%d)",
        percent_high_confidence_one_to_one,
        high_confidence_one_to_one,
        total_input_genes
      )),
      hjust = -0.05,
      size = 3.5
    ) +
    coord_flip(clip = "off") +
    scale_y_continuous(limits = c(0, 105), expand = expansion(mult = c(0, 0.02))) +
    labs(
      x = NULL,
      y = "High-confidence one-to-one ortholog mapping rate (%)"
    ) +
    theme_classic(base_size = 11)

  ggsave(
    file.path(figure_dir, "orthology_mapping_rate.png"),
    mapping_plot,
    width = 7,
    height = max(4, 0.55 * nrow(orthology_summary) + 1.5),
    dpi = 300
  )

  identity_data <- all_mappings %>%
    filter(
      orthology_type == "ortholog_one2one",
      orthology_confidence == 1L,
      !is.na(protein_percent_identity)
    )

  if (nrow(identity_data) > 0L) {
    identity_plot <- ggplot(identity_data, aes(x = protein_percent_identity)) +
      geom_histogram(bins = 20, color = "black", fill = "grey75") +
      facet_wrap(~ comparison, scales = "free_y") +
      labs(
        x = "Protein sequence identity (%)",
        y = "Number of genes"
      ) +
      theme_classic(base_size = 11)

    ggsave(
      file.path(figure_dir, "ortholog_protein_identity.png"),
      identity_plot,
      width = max(7, 3 * n_distinct(identity_data$comparison)),
      height = 4.5,
      dpi = 300
    )
  }
}

select_lifted_regions <- function(lifted, selection = "longest") {
  all_ranges <- unlist(lifted, use.names = FALSE)
  
  if (length(all_ranges) == 0L) {
    return(GRanges())
  }
  
  all_input_index <- rep(seq_along(lifted), lengths(lifted))
  
  if (selection == "all") {
    selected <- all_ranges
    input_index <- all_input_index
  } else {
    range_number <- seq_along(all_ranges)
    ranges_by_input <- split(range_number, all_input_index)
    
    selected_number <- vapply(
      ranges_by_input,
      function(i) i[which.max(width(all_ranges)[i])],
      integer(1)
    )
    
    selected <- all_ranges[selected_number]
    input_index <- all_input_index[selected_number]
  }
  
  selected$input_index <- input_index
  selected
}

if (run_liftover) {
  intervals <- read_csv(interval_input, show_col_types = FALSE, name_repair = "minimal")
  required_interval_columns <- c(chrom_column, start_column, end_column)
  absent <- setdiff(required_interval_columns, names(intervals))
  if (length(absent) > 0L) {
    stop("Interval input is missing: ", paste(absent, collapse = ", "), call. = FALSE)
  }

  starts <- as.integer(intervals[[start_column]])
  ends <- as.integer(intervals[[end_column]])
  if (coordinate_system == "zero_based_half_open") starts <- starts + 1L

  if (anyNA(starts) || anyNA(ends) || any(ends < starts)) {
    stop("Invalid interval coordinates detected.", call. = FALSE)
  }

  rhesus_ranges <- GRanges(
    seqnames = intervals[[chrom_column]],
    ranges = IRanges(start = starts, end = ends)
  )
  
  # GRanges reserves these names and does not allow them as metadata columns.
  reserved_granges_names <- c(
    "seqnames",
    "ranges",
    "strand",
    "seqlevels",
    "seqlengths",
    "isCircular",
    "start",
    "end",
    "width",
    "element"
  )
  
  metadata_columns <- setdiff(
    names(intervals),
    required_interval_columns
  )
  
  # Rename any reserved metadata columns rather than discarding them.
  renamed_metadata_columns <- metadata_columns
  
  reserved_hits <- renamed_metadata_columns %in% reserved_granges_names
  
  renamed_metadata_columns[reserved_hits] <- paste0(
    "input_",
    renamed_metadata_columns[reserved_hits]
  )
  
  interval_metadata <- intervals[
    ,
    metadata_columns,
    drop = FALSE
  ]
  
  names(interval_metadata) <- renamed_metadata_columns
  
  mcols(rhesus_ranges) <- S4Vectors::DataFrame(
    interval_metadata,
    check.names = FALSE
  )

  chain <- import.chain(chain_file)
  lifted <- rtracklayer::liftOver(rhesus_ranges, chain)
  
  map_counts <- lengths(lifted)
  mapped <- map_counts > 0L
  
  # Calculate total non-overlapping lifted width for every input interval
  reduced_lifted <- GenomicRanges::reduce(lifted)
  
  reduced_aligned_width <- vapply(
    seq_along(reduced_lifted),
    function(i) {
      gr <- reduced_lifted[i]
      
      if (length(gr) == 0L) {
        return(0)
      }
      
      sum(as.numeric(IRanges::width(unlist(gr, use.names = FALSE))))
    },
    numeric(1)
  )
  
  original_width <- as.numeric(IRanges::width(rhesus_ranges))
  
  width_retention <- ifelse(
    original_width > 0,
    pmin(reduced_aligned_width / original_width, 1),
    NA_real_
  )
  
  liftover_qc <- intervals %>%
    mutate(
      input_index = seq_len(nrow(intervals)),
      original_width = original_width,
      number_of_human_fragments = map_counts,
      lifted = mapped,
      aligned_width = reduced_aligned_width,
      width_retention_ratio = width_retention
    )

  selected_human <- select_lifted_regions(lifted, mapping_selection)
  if (length(selected_human) > 0L) {
    original_metadata <- intervals[selected_human$input_index, , drop = FALSE]
    human_table <- as.data.frame(selected_human) %>%
      bind_cols(original_metadata %>% rename_with(~ paste0("rhesus_", .x)))
  } else {
    human_table <- tibble()
  }

  write_csv(liftover_qc, file.path(table_dir, "liftover_qc_per_interval.csv"))
  write_csv(human_table, file.path(table_dir, "lifted_hg38_intervals.csv"))

  liftover_summary <- tibble(
    total_intervals = length(rhesus_ranges),
    lifted_intervals = sum(mapped),
    unlifted_intervals = sum(!mapped),
    percent_lifted = 100 * mean(mapped),
    median_width_retention = median(width_retention[mapped], na.rm = TRUE),
    percent_with_at_least_50_percent_width =
      100 * mean(width_retention >= 0.5),
    percent_with_at_least_80_percent_width =
      100 * mean(width_retention >= 0.8)
  )
  write_csv(liftover_summary, file.path(table_dir, "liftover_summary.csv"))

  status_plot_data <- tibble(
    status = c("Lifted", "Unlifted"),
    count = c(sum(mapped), sum(!mapped))
  ) %>%
    mutate(percent = 100 * count / sum(count))

  status_plot <- ggplot(status_plot_data, aes(x = status, y = count)) +
    geom_col(width = 0.55, color = "black", fill = "grey75") +
    geom_text(
      aes(label = sprintf("%d (%.1f%%)", count, percent)),
      vjust = -0.4
    ) +
    labs(x = NULL, y = "Number of rhesus intervals") +
    theme_classic(base_size = 12)

  ggsave(
    file.path(figure_dir, "liftover_mapping_rate.png"),
    status_plot,
    width = 5.5,
    height = 4.5,
    dpi = 300
  )

  retention_plot <- ggplot(
    filter(liftover_qc, lifted),
    aes(x = width_retention_ratio)
  ) +
    geom_histogram(bins = 40, color = "black", fill = "grey75") +
    labs(
      x = "Width-retention ratio after liftOver",
      y = "Number of rhesus intervals"
    ) +
    theme_classic(base_size = 12)

  ggsave(
    file.path(figure_dir, "liftover_width_retention.png"),
    retention_plot,
    width = 6,
    height = 4.5,
    dpi = 300
  )
}