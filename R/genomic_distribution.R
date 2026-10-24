# REQUIREMENT:
# All annotation libraries provided to this script must contain non-overlapping genomic intervals.


# Combine two libraries with priority
cut_overlap <- function(gr1, gr2, ignore_strand = TRUE) {
  gr2_trim <- GenomicRanges::setdiff(gr2, gr1, ignore.strand = ignore_strand)
  if (length(gr2_trim) > 0L) {
    hits <- GenomicRanges::findOverlaps(
      gr2_trim, gr2,
      type = "any", ignore.strand = ignore_strand
    )
    mcols(gr2_trim) <- mcols(gr2)[rep(1, length(gr2_trim)), ]
    mcols(gr2_trim)[, ] <- NA

    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)
    if (length(qh) > 0) {
      if (anyDuplicated(qh)) {
        warning("gr2 has overlapping regions. Taking first match for each trimmed fragment.")
        unique_qh <- unique(qh)
        for (q in unique_qh) {
          idx <- which(qh == q)[1]
          mcols(gr2_trim)[q, ] <- mcols(gr2)[sh[idx], ]
        }
      } else {
        mcols(gr2_trim)[qh, ] <- mcols(gr2)[sh, ]
      }
    }
  }
  combined <- c(gr1, gr2_trim)
  combined <- GenomeInfoDb::sortSeqlevels(combined)
  sort(combined)
}


# Assigns feature annotations to query regions based on overlaps with a library.
# Returns a data.frame with feature counts and percentages.
annotate_GRanges <- function(queries, library, feature_col, mode, feature_order = NULL) {
  n_queries <- length(queries)

  # Find overlaps between queries and library
  overlaps <- GenomicRanges::findOverlaps(queries, library)
  if (length(overlaps) == 0) {
    result <- data.frame(feature = "other", count = n_queries, percentage = 100, stringsAsFactors = FALSE)
    return(result)
  }

  # Extract overlapping regions
  query_hits <- S4Vectors::queryHits(overlaps)
  lib_hits <- S4Vectors::subjectHits(overlaps)
  query_gr <- queries[query_hits]
  lib_gr <- library[lib_hits]
  hit_features <- as.character(GenomicRanges::mcols(lib_gr)[[feature_col]])

  if (mode == "nearest") {
    # Assign each query to closest feature by distance to center
    query_centers <- GenomicRanges::start(query_gr) +
      (GenomicRanges::end(query_gr) - GenomicRanges::start(query_gr)) %/% 2
    dist_to_start <- abs(query_centers - GenomicRanges::start(lib_gr))
    dist_to_end <- abs(query_centers - GenomicRanges::end(lib_gr))
    min_dists <- pmin(dist_to_start, dist_to_end)

    dt <- data.table::data.table(query_idx = query_hits, feature = hit_features, min_dist = min_dists)
    dt_nearest <- dt[dt[, .I[which.min(min_dist)], by = query_idx]$V1]

    features <- rep("other", n_queries)
    features[dt_nearest$query_idx] <- dt_nearest$feature
    feature_counts <- table(features)

    result <- data.frame(
      feature = names(feature_counts),
      count = as.vector(feature_counts),
      percentage = as.vector(feature_counts) / n_queries * 100,
      stringsAsFactors = FALSE
    )
  } else {
    # Weight features by overlap length proportion
    overlap_ranges <- GenomicRanges::pintersect(query_gr, lib_gr)
    overlap_widths <- GenomicRanges::width(overlap_ranges)
    query_widths <- GenomicRanges::width(query_gr)
    overlap_props <- overlap_widths / query_widths

    dt <- data.table::data.table(feature = hit_features, weight = overlap_props)
    feature_weights <- dt[, .(count = sum(weight)), by = feature]
    other_weight <- n_queries - sum(feature_weights$count)

    result <- data.frame(
      feature = c(feature_weights$feature, "other"),
      count = c(feature_weights$count, other_weight),
      percentage = c(feature_weights$count, other_weight) / n_queries * 100,
      stringsAsFactors = FALSE
    )
  }

  # Sort by feature_order if provided
  if (!is.null(feature_order)) {
    full_order <- feature_order[feature_order %in% result$feature]
    result <- result[match(full_order, result$feature), ]
  } else {
    result <- result[order(result$percentage, decreasing = TRUE), ]
  }

  rownames(result) <- NULL
  return(result)
}

# Wrapper function that processes GRanges or GRangesList inputs and returns
# feature annotation percentages in table format (samples × features).
annotate_by_overlap <- function(queries, library, feature_col, mode = c("nearest", "weighted"), feature_order = NULL, BPPARAM = BiocParallel::SerialParam()) {
  mode <- match.arg(mode)

  mcols_lib <- GenomicRanges::mcols(library)
  if (!feature_col %in% colnames(mcols_lib)) {
    stop("Column '", feature_col, "' not found in library metadata")
  }
  feature <- as.character(mcols_lib[[feature_col]])
  bad <- is.na(feature) | feature == ""
  if (any(bad)) {
    library <- library[!bad]
  }

  # Prepare feature order
  all_features <- unique(as.character(GenomicRanges::mcols(library)[[feature_col]]))
  if (!is.null(feature_order)) {
    feature_order <- feature_order[feature_order %in% all_features]
    if (!"other" %in% feature_order) {
      feature_order <- c(feature_order, "other")
    }
  }

  if (inherits(queries, "GRanges")) {
    result_df <- annotate_GRanges(
      queries = queries,
      library = library,
      feature_col = feature_col,
      mode = mode,
      feature_order = feature_order
    )
    pct_vec <- setNames(result_df$percentage, result_df$feature)

    if (!is.null(feature_order)) {
      full_vec <- setNames(rep(0, length(feature_order)), feature_order)
      full_vec[names(pct_vec)] <- pct_vec
      pct_vec <- full_vec
    }

    result_table <- as.data.frame(t(pct_vec))
    rownames(result_table) <- "query"
    return(result_table)
  }

  if (!inherits(queries, "GRangesList")) {
    stop("queries must be GRanges or GRangesList")
  }

  result_list <- BiocParallel::bplapply(queries, function(query_gr) {
    result <- annotate_GRanges(
      queries = query_gr,
      library = library,
      feature_col = feature_col,
      mode = mode,
      feature_order = feature_order
    )

    if (!is.null(feature_order)) {
      pct_vec <- setNames(rep(0, length(feature_order)), feature_order)
      pct_vec[result$feature] <- result$percentage
    } else {
      pct_vec <- setNames(result$percentage, result$feature)
    }
    pct_vec
  }, BPPARAM = BPPARAM)
  
  n_queries <- length(result_list)
  n_features <- length(result_list[[1]])
  result_matrix <- matrix(
    unlist(result_list, use.names = FALSE),
    nrow = n_queries,
    ncol = n_features,
    byrow = TRUE
  )
  colnames(result_matrix) <- names(result_list[[1]])
  rownames(result_matrix) <- names(queries)
  result_matrix[is.na(result_matrix)] <- 0
  
  result_df <- as.data.frame(result_matrix)
  return(result_df)
}

# Create stacked barplot for genomic distribution annotation
plot_stacked_barplot <- function(df, out_dir, prefix, title, palette = "Set3", ylabel = "Percentage (%)", plot_width = 6, plot_height = 5) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(tidyr)
    library(dplyr)
    library(grid)
  })

  if (nrow(df) == 1) {
    df_long <- df %>% pivot_longer(cols = everything(), names_to = "Feature", values_to = "Frequency")
    df_long$cluster <- 1
  } else {
    df_with_rownames <- tibble::rownames_to_column(df, var = "cluster")
    df_long <- df_with_rownames %>% pivot_longer(cols = -cluster, names_to = "Feature", values_to = "Frequency")
    cluster_levels <- rownames(df)
    df_long$cluster <- factor(df_long$cluster, levels = rev(cluster_levels))
  }
  feature_levels <- colnames(df)
  df_long$Feature <- factor(df_long$Feature, levels = feature_levels)
  df_long <- df_long %>% arrange(cluster, Feature)

  if (length(palette) == 1) {
    n_colors <- length(feature_levels)
    colors <- RColorBrewer::brewer.pal(n_colors, palette)
    color_mapping <- setNames(colors, feature_levels)
  } else {
    color_mapping <- palette[feature_levels]
    missing_colors <- is.na(color_mapping)
    if (any(missing_colors)) {
      warning("Some features missing colors, using gray")
      color_mapping[missing_colors] <- "#CCCCCC"
    }
  }

  p <- ggplot(df_long, aes(x = cluster, y = Frequency, fill = Feature)) +
    geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
    coord_flip() +
    scale_fill_manual(values = color_mapping, breaks = feature_levels) +
    labs(title = title, x = NULL, y = ylabel) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 18, hjust = 0.5),
      axis.text.x = element_text(size = 12),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.5, "cm"),
      axis.title.x = element_text(size = 20)
    ) +
    guides(fill = guide_legend(title = NULL, ncol = 1))

  if (nrow(df) == 1) {
    p <- p + theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  } else {
    p <- p + theme(axis.text.y = element_text(size = 12))
  }

  output_file <- file.path(out_dir, paste0(prefix, "_distribution.pdf"))
  ggsave(output_file, plot = p, width = plot_width, height = plot_height)
  message("Plot saved to: ", output_file)
}


# Repeat Element Annotation Function
#   Annotates genomic regions with repeat elements using RepeatMasker data,
#   computes overlap-based repeat class percentages, and generates a summary
#   table along with an optional visualization.
repeat_distribution <- function(query, out_dir = "./", ref_genome = "hg38", mode = "nearest", plot = TRUE, BPPARAM = BiocParallel::SerialParam()) {
  repeat_order <- c(
    "SINE", "LINE", "LTR", "Retroposon", "RC", "DNA",
    "Satellite", "Simple_repeat", "Low_complexity",
    "rRNA", "tRNA", "snRNA", "scRNA", "srpRNA", "RNA",
    "Unknown"
  )
  repeat_colors <- setNames(
    c(
      "#3B82F6", "#60A5FA", "#34D399", "#6EE7B7", "#A78BFA", "#8B5CF6", # transposable (6)
      "#F59E0B", "#FBBF24", "#FDE68A", # tandem (3)
      "#EC4899", "#F9A8D4", "#FCA5A5", "#FED7AA", "#E9D5FF", "#DDD6FE", # RNA (6)
      "#D1D5DB", # unknown (1)
      "#808080" # other(1)
    ),
    c(repeat_order, "other")
  )

  # Load reference library
  if (ref_genome == "hg38") {
    repeat_library_file <- download_rds("RepeatMasker_hg38.rds")
  } else if (ref_genome == "mm10") {
    repeat_library_file <- download_rds("RepeatMasker_mm10.rds")
  } else {
    stop("ref_genome must be 'hg38' or 'mm10'")
  }
  repeat_library <- readRDS(repeat_library_file)
  message(glue("Using reference genome {ref_genome} with {length(repeat_library)} repeats"))
  style <- seqlevelsStyle(query)[1]
  if (seqlevelsStyle(repeat_library)[1] != style) {
    seqlevelsStyle(repeat_library) <- style
  }

  # Annotate query regions with repeat elements
  repeat_df <- annotate_by_overlap(queries = query, library = repeat_library, feature_col = "rep_class", mode = mode, feature_order = repeat_order, BPPARAM = BPPARAM)
  write.table(repeat_df, file = file.path(out_dir, "repeat_distribution.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

  # Create stacked barplot
  if (plot) {
    plot_stacked_barplot(df = repeat_df, out_dir = out_dir, prefix = "repeat", title = "Repeat Elements", palette = repeat_colors)
  }
}

# ChromHMM Annotation Function
# Description:
#   Annotates genomic regions with chromatin states using ChromHMM data,
#   computes overlap-based chromatin state percentages, and generates a summary
#   table along with an optional visualization.
chromhmm_distribution <- function(query, out_dir = "./", ref_genome = "hg38", mode = "nearest", plot = TRUE, BPPARAM = BiocParallel::SerialParam()) {
  hmm_order <- c(
    "Acet", "EnhWk", "EnhA", "PromF", "TSS", # active enhancers & promoters
    "TxWk", "TxEx", "Tx", # transcription continuum
    "DNase", "OpenC", "TxEnh", # accessible chromatin & enhancer-coupled transcription
    "ReprPCopenC",
    "BivProm", "ZNF", # developmental regulation & ZNF clusters
    "ReprPC", "HET", # Polycomb repression & constitutive heterochromatin
    "GapArtf", # gap artifacts
    "Quies" # quiescent background
  )

  hmm_colors <- setNames(
    c(
      "#EF4444", "#FB9A99", "#FF7F00", "#FDBF6F", "#CAB2D6", # active
      "#B2DF8A", "#22C55E", "#10B981", # transcription
      "#CFE8F3", "#A6CEE3", "#6A3D9A", # accessible
      "#FFFF99", # intermediate
      "#B15928", "#FCCDE5", # developmental
      "#3B82F6", "#1E40AF", # repressed
      "#E2D6C8", # gap artifacts
      "#A0A0A0", # quiescent
      "#808080"
    ),
    c(hmm_order, "other")
  )

  # Load reference library
  if (ref_genome == "hg38") {
    hmm_library_file <- download_rds("ChromHMM_hg38.rds")
  } else if (ref_genome == "mm10") {
    hmm_library_file <- download_rds("ChromHMM_mm10.rds")
  } else {
    stop("ref_genome must be 'hg38' or 'mm10'")
  }
  hmm_library <- readRDS(hmm_library_file)
  message(glue("Using reference genome {ref_genome} with {length(hmm_library)} ChromHMMs"))
  style <- seqlevelsStyle(query)[1]
  if (seqlevelsStyle(hmm_library)[1] != style) {
    seqlevelsStyle(hmm_library) <- style
  }

  # Annotate query regions with ChromHMM elements
  hmm_df <- annotate_by_overlap(queries = query, library = hmm_library, feature_col = "clean_class", mode = mode, feature_order = hmm_order, BPPARAM = BPPARAM)
  write.table(hmm_df, file = file.path(out_dir, "chromhmm_distribution.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

  # Create stacked barplot
  if (plot) {
    plot_stacked_barplot(df = hmm_df, out_dir = out_dir, prefix = "chromhmm", title = "ChromHMM States", palette = hmm_colors)
  }
}


# cCRE Element Annotation Function
# Description:
#   Annotates genomic regions with cCRE elements using ENCODE data,
#   computes overlap-based repeat class percentages, and generates a summary
#   table along with an optional visualization.
ccre_distribution <- function(query, out_dir = "./", ref_genome = "hg38", mode = "nearest", plot = TRUE, BPPARAM = BiocParallel::SerialParam()) {
  ccre_order <- c(
    "dELS", "pELS", "PLS",
    "CA-H3K4me3", "CA-CTCF", "CA-TF", "TF", "CA"
  )
  ccre_colors <- setNames(
    c(
      "#8DD3C7", "#80B1D3", "#BEBADA",
      "#D55E00", "#C77C2E", "#E69F00", "#F4A261", "#F6E8C3",
      "#D9D9D9"
    ),
    c(ccre_order, "other")
  )

  if (ref_genome == "hg38") {
    ccre_library_file <- download_rds("ENCODE_cCRE_hg38.rds")
  } else if (ref_genome == "mm10") {
    ccre_library_file <- download_rds("ENCODE_cCRE_mm10.rds")
  } else {
    stop("ref_genome must be 'hg38' or 'mm10'")
  }

  ccre_library <- readRDS(ccre_library_file)
  style <- seqlevelsStyle(query)[1]
  if (seqlevelsStyle(ccre_library)[1] != style) {
    seqlevelsStyle(ccre_library) <- style
  }

  # Annotate query regions with cCRE elements
  ccre_df <- annotate_by_overlap(queries = query, library = ccre_library, feature_col = "ccre_class", mode = mode, feature_order = ccre_order, BPPARAM = BPPARAM)
  write.table(ccre_df, file = file.path(out_dir, "ccre_distribution.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

  # Create stacked barplot
  if (plot) {
    plot_stacked_barplot(df = ccre_df, out_dir = out_dir, prefix = "ccre", title = "cCRE States", palette = ccre_colors)
  }
}

# Genic Structure Annotation Function
# Description:
#   Annotates genomic regions with gene structural components based on
#   transcript annotations, including 5' UTR, exon, intron, and 3' UTR.
#   This function summarizes the distribution of query regions across
#   genic structures and optionally generates a visualization.
genic_distribution <- function(query, out_dir = "./", ref_genome = "hg38", ref_source = "knownGene", mode = "nearest", plot = TRUE, BPPARAM = BiocParallel::SerialParam()) {
  gene_order <- c(
    "Promoter", "5' UTR", "Exon", "Intron", "3' UTR"
  )
  gene_colors <- setNames(
    c(
      "#FC8D62", "#FFD92F", "#8DA0CB", "#A6D854", "#E78AC3", "#E0E0E0"
    ),
    c(gene_order, "other")
  )

  if (ref_genome == "hg38") {
    if (ref_source == "knownGene") {
      gene_library_file <- download_rds("knownGene_genic_hg38.rds")
    } else if (ref_source == "GENCODE") {
      gene_library_file <- download_rds("GENCODE_genic_hg38.rds")
    } else {
      stop("ref_source must be 'knownGene' or 'GENCODE'")
    }
  } else if (ref_genome == "mm10") {
    if (ref_source == "knownGene") {
      gene_library_file <- download_rds("knownGene_genic_mm10.rds")
    } else if (ref_source == "GENCODE") {
      gene_library_file <- download_rds("GENCODE_genic_mm10.rds")
    } else {
      stop("ref_source must be 'knownGene' or 'GENCODE'")
    }
  } else {
    stop("ref_genome must be 'hg38' or 'mm10'")
  }

  gene_library <- readRDS(gene_library_file)
  style <- seqlevelsStyle(query)[1]
  if (seqlevelsStyle(gene_library)[1] != style) {
    seqlevelsStyle(gene_library) <- style
  }

  gene_df <- annotate_by_overlap(queries = query, library = gene_library, feature_col = "class", mode = mode, feature_order = gene_order, BPPARAM = BPPARAM)
  write.table(gene_df, file = file.path(out_dir, "genic_distribution.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

  # Create stacked barplot
  if (plot) {
    plot_stacked_barplot(df = gene_df, out_dir = out_dir, prefix = "genic", title = "Genic States", palette = gene_colors)
  }
}

# Genomic Distribution Annotation Wrapper
# Description:
#   Comprehensive genomic annotation function that performs multiple types of
#   genomic feature distributions on query regions. Supports genic structures,
#   cCRE elements, ChromHMM states, and repeat elements.
#
# Parameters:
#   query: GRangesList or GRanges object containing query genomic regions
#          (e.g., peak sets, biclustering-derived clusters).
#   out_dir: Directory to save output files (tables and plots).
#   distributions: Character vector specifying which annotation types to perform.
#                Options: "genic", "ccre", "chromhmm", "repeat"
#                Default: c("genic", "ccre")
#   ref_genome: Reference genome version. Must be either "hg38" or "mm10".
#               Default: "hg38"
#   ref_source: Gene annotation source for genic annotation only.
#               Options: "knownGene" (UCSC) or "GENCODE"
#               Default: "knownGene"
#   mode: Annotation mode for assigning features to query regions.
#         Default: "nearest"
#         - "nearest": Assigns each query region to the closest overlapping feature
#                      based on distance from query center to feature boundaries.
#         - "weighted": Assigns features proportionally based on overlap length,
#                       weighted by (overlap_length / query_length).
#   plot: Boolean, whether to generate stacked barplots for each annotation.
#         Default: TRUE
#
# Output files (saved to out_dir for each selected annotation):
#   - genic_distribution.tsv & genic_distribution.pdf
#   - ccre_distribution.tsv & ccre_distribution.pdf
#   - chromhmm_distribution.tsv & chromhmm_distribution.pdf
#   - repeat_distribution.tsv & repeat_distribution.pdf
genomic_distribution <- function(query, out_dir, distributions = c("genic", "ccre"), ref_genome = "hg38", ref_source = "knownGene", mode = "nearest", plot = TRUE) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Set up parallel processing
  BPPARAM <- get_BPPARAM(backend = "multicore")
  
  # Perform selected distributions
  if ("genic" %in% distributions) {
    message("\n========== Running Genic annotation ==========")
    genic_distribution(
      query = query,
      out_dir = out_dir,
      ref_genome = ref_genome,
      ref_source = ref_source,
      mode = mode,
      plot = plot,
      BPPARAM = BPPARAM
    )
    message("Genic annotation complete. Results saved to: ", out_dir)
  }

  if ("ccre" %in% distributions) {
    message("\n========== Running cCRE annotation ==========")
    ccre_distribution(
      query = query,
      out_dir = out_dir,
      ref_genome = ref_genome,
      mode = mode,
      plot = plot,
      BPPARAM = BPPARAM
    )
    message("cCRE annotation complete. Results saved to: ", out_dir)
  }

  if ("chromhmm" %in% distributions) {
    message("\n========== Running ChromHMM annotation ==========")
    chromhmm_distribution(
      query = query,
      out_dir = out_dir,
      ref_genome = ref_genome,
      mode = mode,
      plot = plot,
      BPPARAM = BPPARAM
    )
    message("ChromHMM annotation complete. Results saved to: ", out_dir)
  }

  if ("repeat" %in% distributions) {
    message("\n========== Running Repeat annotation ==========")
    repeat_distribution(
      query = query,
      out_dir = out_dir,
      ref_genome = ref_genome,
      mode = mode,
      plot = plot,
      BPPARAM = BPPARAM
    )
    message("Repeat annotation complete. Results saved to: ", out_dir)
  }
}
