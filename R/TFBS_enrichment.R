# Write region-level TFBS overlap records
#
# Output: gzipped BED-like TSV (#chrom, start, end, TF), 0-based half-open;
# one row per region-TF pair where the region overlaps >= 1 motif site.

write_TFBS_region <- function(region_gr, TFBS_library, out_path) {
  # GRangesList as query: hits are unique (TF, region) pairs,
  # same counting rule as countOverlaps(TFBS_library, region_gr)
  hits <- findOverlaps(TFBS_library, region_gr)
  idx <- subjectHits(hits)

  dt <- unique(data.table(
    chrom = as.character(seqnames(region_gr))[idx],
    start = start(region_gr)[idx] - 1L,
    end   = end(region_gr)[idx],
    TF    = names(TFBS_library)[queryHits(hits)]
  ))
  setorder(dt, chrom, start, end, TF)
  setnames(dt, "chrom", "#chrom")
  fwrite(dt, out_path, sep = "\t")

  message(glue("Region-level TFBS records ({nrow(dt)} region-TF pairs) saved to {out_path}"))
  invisible(out_path)
}

# TFBS Enrichment Analysis for Single Region Set
#
# Internal function that performs TFBS enrichment analysis for a single set of
# target regions against control regions using Fisher's exact test.

TFBS_enrichment_single <- function(query_gr, control_gr, TFBS_library, regions = NULL, out_path = "./TFBS_enrichment.tsv", region_out_path = NULL, BPPARAM = BiocParallel::SerialParam()) {
  # Resize regions if specified
  if (!is.null(regions)) {
    query_gr <- resize(query_gr, width = regions, fix = "center")
    control_gr <- resize(control_gr, width = regions, fix = "center")
  }

  # Count overlaps between motif sites and regions
  target_overlap <- countOverlaps(TFBS_library, query_gr)
  control_overlap <- countOverlaps(TFBS_library, control_gr)
  n_target <- length(query_gr)
  n_control <- length(control_gr)
  tf_names <- names(TFBS_library)

  # Build contingency table and perform Fisher's exact test
  test_result <- data.table(
    TF = tf_names,
    target_hit = target_overlap,
    control_hit = control_overlap,
    target_off = n_target - target_overlap,
    control_off = n_control - control_overlap
  )

  fisher_results <- BiocParallel::bplapply(seq_len(nrow(test_result)), function(i) {
    row <- test_result[i, ]
    contingency_table <- matrix(
      c(row$target_hit + 1, row$control_hit + 1,
        row$target_off + 1, row$control_off + 1),
      nrow = 2,
      dimnames = list(c("query", "control"), c("hit", "off"))
    )
    test_res <- fisher.test(contingency_table, alternative = "greater")
    list(
      odds_ratio = as.numeric(test_res$estimate),
      pvalue = test_res$p.value,
      odds_ratio_se = sqrt(sum(1 / contingency_table))
    )
  }, BPPARAM = BPPARAM)

  test_result[, `:=`(
    odds_ratio = vapply(fisher_results, `[[`, numeric(1), "odds_ratio"),
    pvalue = vapply(fisher_results, `[[`, numeric(1), "pvalue"),
    odds_ratio_se = vapply(fisher_results, `[[`, numeric(1), "odds_ratio_se")
  )]

  # Multiple testing correction
  test_result[, FDR := p.adjust(pvalue, method = "BH")]
  setorder(test_result, pvalue)

  # Save results
  res <- as.data.frame(test_result)
  rownames(res) <- res$TF
  res <- res[, -1]
  write.table(res, file = out_path, sep = "\t", quote = FALSE, col.names = NA)

  n_significant <- sum(res$FDR < 0.05)
  message(glue("Results saved to {out_path}"))
  message(glue("Found {n_significant} significant motifs at FDR < 0.05"))

  if (!is.null(region_out_path)) {
    write_TFBS_region(query_gr, TFBS_library, region_out_path)
  }

  invisible(out_path)
}


# TFBS Enrichment Calculation
#
# This function analyzes the enrichment of transcription factor binding motifs
# in target genomic regions compared to control regions using Fisher's exact test.
# The analysis identifies motifs that are significantly over-represented in target
# regions, accounting for background binding patterns observed in matched controls.
#
# Parameters:
#   query: GRanges or GRangesList object of target genomic regions. If GRangesList,
#          enrichment analysis is performed separately for each element and results
#          are saved to individual TSV files.
#   control: GRanges object containing control/background genomic regions.
#            Same control set is used for all elements if query is GRangesList.
#   regions: Width in base pairs to resize all regions around their center.
#            If NULL, uses original region sizes. Default: NULL
#   ref_genome: Reference genome version ("hg38" or "mm10"). Default: "hg38"
#   functional_region: Optional GRanges object of functional regions (e.g., ATAC-seq
#                      peaks, DNase hypersensitive sites). If provided, only motif
#                      sites overlapping these regions are considered in the analysis.
#                      Default: NULL (use all motif sites)
#   out_dir: Output directory path for enrichment results. Directory is created if
#            it doesn't exist. Default: "./"
#   style: Chromosome naming style (e.g., "UCSC" or "NCBI"). If NULL, inferred
#          from query object. Default: NULL
#   save_region_TFBS: If TRUE, also write region-level TF overlap records to
#                     "<out_dir>/TFBS_region/" as "TFBS_region_<label>.tsv.gz"
#                     (or "TFBS_region.tsv.gz" for GRanges input): BED-like
#                     columns #chrom, start, end, TF (0-based half-open), one row
#                     per region-TF pair. Coordinates are those of the tested
#                     regions (after resizing by `regions`, if set). Default: FALSE
#
# Returns:
#   Character vector of output file path(s), returned invisibly:
#   - For GRanges input: Single path to "TFBS_enrichment.tsv"
#   - For GRangesList input: Vector of paths to "TFBS_enrichment_<label>.tsv" files
#   Each TSV file contains motif enrichment statistics (odds ratio, p-value, FDR)

TFBS_enrichment <- function(query, control, regions = NULL, ref_genome = "hg38", functional_region = NULL, out_dir = "./", style = NULL, save_region_TFBS = FALSE) {
  # Load required packages
  suppressPackageStartupMessages({
    library(IRanges)
    library(GenomicRanges)
    library(data.table)
    library(glue)
  })

  # Set up parallel processing
  BPPARAM <- get_BPPARAM(backend = "multicore")

  # Parameter validation
  if (!inherits(query, "GRanges") && !inherits(query, "GRangesList")) {
    stop("Query must be a GRanges or GRangesList object")
  }
  if (!inherits(control, "GRanges")) {
    stop("control must be a GRanges object")
  }
  if (!ref_genome %in% c("hg38", "mm10")) {
    stop("Unsupported genome. Please use 'hg38' or 'mm10'.")
  }

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  region_dir <- file.path(out_dir, "TFBS_region")
  if (save_region_TFBS && !dir.exists(region_dir)) {
    dir.create(region_dir, recursive = TRUE)
  }

  # Load TFBS library
  TFBS_library_file <- switch(ref_genome,
    hg38 = "TFBS_lib_hg38.rds",
    mm10 = "TFBS_lib_mm10.rds"
  )
  TFBS_library_file <- download_rds(TFBS_library_file)
  TFBS_library <- readRDS(TFBS_library_file)
  message(glue("Using reference genome {ref_genome} with {length(TFBS_library)} TFs"))

  # Harmonize chromosome styles
  if (is.null(style)) {
    style <- seqlevelsStyle(query)[1]
  } else {
    seqlevelsStyle(query) <- style
  }
  seqlevelsStyle(control) <- style
  seqlevelsStyle(TFBS_library) <- style

  # Filter by functional regions if provided
  if (!is.null(functional_region)) {
    if (!inherits(functional_region, "GRanges")) {
      stop("functional_region must be a GRanges object")
    }
    seqlevelsStyle(functional_region) <- style
    message(glue("Filtering motif sites using {length(functional_region)} functional regions"))
    TFBS_library <- endoapply(TFBS_library, subsetByOverlaps, functional_region)
    TFBS_library <- TFBS_library[lengths(TFBS_library) > 0]

    if (sum(lengths(TFBS_library)) == 0) {
      warning("No motif sites overlap with the provided functional regions.")
      return(invisible(NULL))
    }
    message(glue("Retained {length(TFBS_library)} TFs after filtering"))
  }

  # Perform enrichment analysis
  if (inherits(query, "GRanges")) {
    # Single GRanges input
    message("Counting overlaps and performing enrichment analysis...")
    out_path <- file.path(out_dir, "TFBS_enrichment.tsv")

    result_path <- TFBS_enrichment_single(
      query_gr = query,
      control_gr = control,
      TFBS_library = TFBS_library,
      regions = regions,
      out_path = out_path,
      region_out_path = if (save_region_TFBS) file.path(region_dir, "TFBS_region.tsv.gz") else NULL,
      BPPARAM = BPPARAM
    )

    return(invisible(result_path))
  } else {
    # GRangesList input - batch processing
    message(glue("Processing {length(query)} target region sets..."))

    result_paths <- sapply(names(query), function(label) {
      message(glue("\nProcessing: {label} ({length(query[[label]])} regions)"))

      out_path <- file.path(out_dir, paste0("TFBS_enrichment_", label, ".tsv"))

      TFBS_enrichment_single(
        query_gr = query[[label]],
        control_gr = control,
        TFBS_library = TFBS_library,
        regions = regions,
        out_path = out_path,
        region_out_path = if (save_region_TFBS) file.path(region_dir, paste0("TFBS_region_", label, ".tsv.gz")) else NULL,
        BPPARAM = BPPARAM
      )
    })

    message(glue("\n{length(result_paths)} enrichment analyses completed"))
    return(invisible(result_paths))
  }
}
