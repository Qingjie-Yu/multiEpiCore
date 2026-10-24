# Identify the nearest gene for each genomic interval.
#
# This function takes a set of query genomic regions and a reference set of gene
# annotations (both as GRanges objects), harmonizes their chromosome naming
# styles, and for each query region determines the nearest gene on the same
# chromosome. If the query midpoint overlaps one or more genes, the shortest
# overlapping gene is selected. Otherwise, the function evaluates the nearest
# upstream and downstream genes and chooses the closest gene body; ties are
# resolved by selecting the shorter gene. The output returns, for each query
# region, the gene name, gene length, and directed distance to the gene's TSS.
#
# Parameters:
#   query_gr: GRanges object of query regions.
#   genes_gr: GRanges object of gene annotations, with gene names stored in names().
#
# Returns:
#   A data.frame with one row per query region containing:
#     - chr: chromosome of the query region
#     - start, end: genomic coordinates of the query region
#     - nearest_gene: name of the nearest gene
#     - gene_length: width of the selected gene
#     - tss_distance: signed distance to the gene TSS
#                     (negative = upstream, positive = downstream)

get_nearest_gene <- function(query_gr, genes_gr) {
  suppressPackageStartupMessages({
    library(dplyr)
  })

  query_mid <- resize(query_gr, width = 1, fix = "center")
  overlaps <- findOverlaps(query_mid, genes_gr)
  n_query <- length(query_gr)
  nearest_gene <- rep(NA_character_, n_query)
  gene_length <- rep(NA_integer_, n_query)
  tss_distance <- rep(NA_integer_, n_query)

  # Handle overlapping genes - select shortest gene
  if (length(overlaps) > 0) {
    query_hits <- queryHits(overlaps)
    subject_hits <- subjectHits(overlaps)

    overlap_df <- data.frame(
      query_idx = query_hits,
      gene_idx = subject_hits,
      gene_width = width(genes_gr[subject_hits]),
      stringsAsFactors = FALSE
    )

    best_overlaps <- overlap_df |>
      dplyr::group_by(query_idx) |>
      dplyr::slice_min(gene_width, n = 1, with_ties = FALSE) |>
      dplyr::ungroup()
    best_overlaps <- as.data.frame(best_overlaps)

    for (i in seq_len(nrow(best_overlaps))) {
      q_idx <- best_overlaps$query_idx[i]
      g_idx <- best_overlaps$gene_idx[i]
      best_gene <- genes_gr[g_idx]
      q_mid_pos <- start(query_mid[q_idx])

      nearest_gene[q_idx] <- names(best_gene)
      gene_length[q_idx] <- width(best_gene)

      # TSS distance: negative = upstream, positive = downstream
      if (as.character(strand(best_gene)) == "-") {
        tss_distance[q_idx] <- q_mid_pos - end(best_gene)
      } else {
        tss_distance[q_idx] <- q_mid_pos - start(best_gene)
      }
    }
  }

  no_overlap_idx <- which(is.na(nearest_gene))

  # Handle non-overlapping genes
  if (length(no_overlap_idx) > 0) {
    query_chr <- as.character(seqnames(query_mid))

    for (chr in unique(query_chr[no_overlap_idx])) {
      chr_query_idx <- no_overlap_idx[query_chr[no_overlap_idx] == chr]
      chr_genes <- genes_gr[seqnames(genes_gr) == chr]

      if (length(chr_genes) == 0) next

      chr_query_mid <- query_mid[chr_query_idx]
      nearest_info <- distanceToNearest(chr_query_mid, chr_genes)

      if (length(nearest_info) > 0) {
        # Handle ties by selecting shorter gene
        dist_df <- data.frame(
          q_local_idx = queryHits(nearest_info),
          g_local_idx = subjectHits(nearest_info),
          distance = mcols(nearest_info)$distance,
          gene_width = width(chr_genes[subjectHits(nearest_info)])
        )

        best_matches <- dist_df |>
          dplyr::group_by(q_local_idx) |>
          dplyr::filter(distance == min(distance)) |>
          dplyr::slice_min(gene_width, n = 1, with_ties = FALSE)
        best_matches <- as.data.frame(best_matches)

        for (i in seq_len(nrow(best_matches))) {
          q_local_idx <- best_matches$q_local_idx[i]
          g_local_idx <- best_matches$g_local_idx[i]
          q_global_idx <- chr_query_idx[q_local_idx]

          best_gene <- chr_genes[g_local_idx]
          q_mid_pos <- start(chr_query_mid[q_local_idx])

          nearest_gene[q_global_idx] <- names(best_gene)
          gene_length[q_global_idx] <- width(best_gene)

          # TSS distance: negative = upstream, positive = downstream
          if (as.character(strand(best_gene)) == "-") {
            tss_distance[q_global_idx] <- q_mid_pos - end(best_gene)
          } else {
            tss_distance[q_global_idx] <- q_mid_pos - start(best_gene)
          }
        }
      }
    }
  }

  result_final <- data.frame(
    chr = as.character(seqnames(query_gr)),
    start = start(query_gr),
    end = end(query_gr),
    nearest_gene = nearest_gene,
    gene_length = gene_length,
    tss_distance = tss_distance,
    stringsAsFactors = FALSE
  )
  result_final <- result_final[!is.na(result_final$nearest_gene), ]
  return(result_final)
}

# Generate matched control regions for a single GRanges object
#
# Internal function that creates control regions by matching gene lengths
# and preserving TSS distances.
control_regions_single <- function(query_gr, genes_gr, chr_sizes, n_rep = 1, regions = 800, seed = 42, length_tolerance = 0.2, BPPARAM = BiocParallel::SerialParam()) {
  gene_lengths <- width(genes_gr)
  gene_names <- names(genes_gr)
  gene_chrs <- as.character(seqnames(genes_gr))
  gene_strands <- as.character(strand(genes_gr))
  gene_starts <- start(genes_gr)
  gene_ends <- end(genes_gr)

  target_with_genes <- get_nearest_gene(query_gr, genes_gr)
  n_targets <- nrow(target_with_genes)
  half_width <- floor(regions / 2)

  control_list <- BiocParallel::bplapply(seq_len(n_targets), function(i) {
    set.seed(seed + i)

    orig_gene <- target_with_genes$nearest_gene[i]
    orig_gene_length <- target_with_genes$gene_length[i]
    tss_distance <- target_with_genes$tss_distance[i]

    # Find candidate genes with similar length
    tol <- length_tolerance
    lower_bound <- orig_gene_length * (1 - tol)
    upper_bound <- orig_gene_length * (1 + tol)
    candidate_mask <- (gene_lengths >= lower_bound) &
      (gene_lengths <= upper_bound) &
      (gene_names != orig_gene)

    # Relax tolerance if insufficient candidates
    while (sum(candidate_mask) < n_rep * 2 && tol <= 0.9) {
      tol <- tol + 0.1
      lower_bound <- orig_gene_length * (1 - tol)
      upper_bound <- orig_gene_length * (1 + tol)
      candidate_mask <- (gene_lengths >= lower_bound) &
        (gene_lengths <= upper_bound) &
        (gene_names != orig_gene)
    }

    if (sum(candidate_mask) < n_rep) {
      candidate_mask <- (gene_names != orig_gene)
    }

    control_regions <- list()
    n_valid <- 0

    candidate_indices <- which(candidate_mask)
    if (length(candidate_indices) == 0) {
      attr(control_regions, "warning") <- sprintf("    No candidate genes found for target %d", i)
      return(control_regions)
    }
    max_attempts <- min(length(candidate_indices), n_rep * 5)
    sampled_indices <- sample(candidate_indices, size = max_attempts, replace = FALSE)

    for (idx in sampled_indices) {
      if (n_valid >= n_rep) break

      # Calculate control region position preserving TSS distance
      if (gene_strands[idx] == "-") {
        control_center <- gene_ends[idx] + tss_distance
      } else {
        control_center <- gene_starts[idx] + tss_distance
      }

      control_start <- control_center - half_width
      control_end <- control_center + half_width - 1
      control_chr <- gene_chrs[idx]

      # Validate chromosome boundaries
      if (control_start >= 1 && control_end <= chr_sizes[control_chr]) {
        control_regions[[n_valid + 1]] <- GRanges(
          seqnames = control_chr,
          ranges = IRanges(start = control_start, end = control_end)
        )
        n_valid <- n_valid + 1
      }
    }

    if (n_valid < n_rep) {
      attr(control_regions, "warning") <- sprintf("    Only found %d/%d valid control regions for target %d", n_valid, n_rep, i)
    }
    return(control_regions)
  }, BPPARAM = BPPARAM)

  warnings_msgs <- sapply(control_list, function(x) attr(x, "warning"))
  warnings_msgs <- warnings_msgs[!sapply(warnings_msgs, is.null)]
  if (length(warnings_msgs) > 0) {
    for (msg in warnings_msgs) warning(msg, call. = FALSE)
  }

  control_list <- control_list[sapply(control_list, length) > 0]
  if (length(control_list) == 0) {
    warning("No valid control regions found")
    return(GRanges())
  }
  control_list <- unlist(control_list, recursive = FALSE)
  control_gr <- unlist(GRangesList(control_list))
  return(control_gr)
}


# Generate matched control regions
#
# This function generates background control regions by randomly selecting genes
# with similar length to the nearest gene of each target peak, then placing
# control regions at the same TSS distance on the selected genes. This approach
# preserves the relationship between peaks and gene structure while avoiding
# biases from using the same genomic neighborhoods.
#
# Parameters:
#   query: GRanges or GRangesList object of target peak regions. If GRangesList,
#          control regions are generated separately for each element and then
#          combined into a single GRanges object.
#   ref_genome: Reference genome version ("hg38" or "mm10")
#   ref_source: Gene annotation source ("knownGene" or "GENCODE")
#               - "knownGene": UCSC knownGene annotation from TxDb packages;
#                 gene span = union of all transcripts
#               - "GENCODE": GENCODE basic annotation (v49 for hg38, vM25 for mm10);
#                 gene span = TSS–TES of one representative transcript per gene
#   style: Chromosome naming style (e.g., "UCSC" or "NCBI"). If NULL, inferred
#          from query object. Default: NULL
#   n_rep: Number of control regions to generate per target peak. Default: 1
#   regions: Width of control regions in base pairs. Default: 800
#   seed: Random seed for reproducibility. Default: 42
#   length_tolerance: Initial tolerance for gene length matching (±20% default).
#                     Automatically relaxed up to ±100% if insufficient candidates.
#
# Returns:
#   A GRanges object containing all control regions. If input is GRangesList,
#   control regions from all elements are combined into a single GRanges.

get_matched_control <- function(query, ref_genome = "hg38", ref_source = "knownGene", style = NULL, n_rep = 1, regions = 800, seed = 42, length_tolerance = 0.2) {
  # Set up parallel processing
  BPPARAM <- get_BPPARAM(verbose = FALSE, backend = "multicore")

  # Validate input
  if (!inherits(query, "GRanges") && !inherits(query, "GRangesList")) {
    stop("query must be a GRanges or GRangesList object")
  }
  txdb_pkgs <- c(
    hg38 = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    mm10 = "TxDb.Mmusculus.UCSC.mm10.knownGene"
  )
  if (!ref_genome %in% names(txdb_pkgs)) {
    stop("Unsupported genome. Please use 'hg38' or 'mm10'.")
  }
  if (!ref_source %in% c("knownGene", "GENCODE")) {
    stop("ref_source must be 'knownGene' or 'GENCODE'")
  }

  # Infer chromosome style from query if not provided
  if (is.null(style)) {
    style <- seqlevelsStyle(query)[1]
  }

  # Chromosome sizes from TxDb (standard chromosomes, no chrM)
  txdb <- getExportedValue(txdb_pkgs[[ref_genome]], txdb_pkgs[[ref_genome]])
  si <- seqinfo(txdb)
  si <- si[grep("^chr([0-9]+|X|Y)$", seqnames(si), value = TRUE)]
  seqlevelsStyle(si) <- style
  chr_list <- seqnames(si)
  chr_sizes <- seqlengths(si)

  # Load gene annotations
  if (ref_source == "knownGene") {
    genes_gr <- suppressMessages(GenomicFeatures::genes(txdb))
  } else {
    gencode_file <- download_rds(sprintf("GENCODE_rep_%s.rds", ref_genome))
    genes_gr <- readRDS(gencode_file)
    # One representative transcript per gene: start/end = TSS/TES
    genes_gr <- genes_gr[genes_gr$type == "transcript"]
    # Keep only genes whose representative transcript has supporting evidence
    genes_gr <- genes_gr[genes_gr$rep_rank <= 7]
    names(genes_gr) <- genes_gr$gene_name
    mcols(genes_gr) <- NULL
  }

  # Filter and sort gene annotations
  seqlevelsStyle(genes_gr) <- style
  genes_gr <- genes_gr[seqnames(genes_gr) %in% chr_list]
  seqlevels(genes_gr) <- chr_list
  genes_gr <- sort(genes_gr)

  # Generate control regions
  if (inherits(query, "GRanges")) {
    control_gr <- control_regions_single(
      query_gr = query,
      genes_gr = genes_gr,
      chr_sizes = chr_sizes,
      n_rep = n_rep,
      regions = regions,
      seed = seed,
      length_tolerance = length_tolerance,
      BPPARAM = BPPARAM
    )
  } else {
    # Process each element of GRangesList separately
    control_grl <- lapply(names(query), function(label) {
      cat("\nProcessing cluster:", label, "with", length(query[[label]]), "regions\n")
      control_regions_single(
        query_gr = query[[label]],
        genes_gr = genes_gr,
        chr_sizes = chr_sizes,
        n_rep = n_rep,
        regions = regions,
        seed = seed,
        length_tolerance = length_tolerance,
        BPPARAM = BPPARAM
      )
    })
    names(control_grl) <- names(query)
    control_gr <- unlist(GRangesList(control_grl), use.names = FALSE)
  }

  control_gr
}