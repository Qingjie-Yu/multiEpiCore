is_bam_paired <- function(bam_path) {
  param <- ScanBamParam(what = "flag")
  flags <- scanBam(bam_path[1], param = param)[[1]]$flag
  if (length(flags) == 0) return(FALSE)
  any(bitwAnd(head(flags, 1000), 0x1) != 0)
}

# TRUE if any of the first 1000 records of ONE BAM carries the "read paired" flag (0x1); NA if the
# BAM has no records. Reads only those records (yieldSize, no index). Used to check every BAM
# against the PE/SE decision; the decision itself still comes from is_bam_paired() under "auto".
.bam_paired_head <- function(bam) {
  bf <- BamFile(bam, index = character(0), yieldSize = 1000)
  open(bf)
  on.exit(close(bf))
  flags <- scanBam(bf, param = ScanBamParam(what = "flag"))[[1]]$flag
  if (length(flags) == 0) return(NA)
  any(bitwAnd(flags, 0x1) != 0)
}

# ---- SE extension resolver (per-BAM 'ext' helper for build_count_matrix) ----
# Resolve 'ext' into an integer vector aligned to bam_path in the order given.
#   single number            -> recycled to every BAM
#   unnamed vector / list    -> positional; length must equal length(bam_path)
#   named vector / list      -> keyed lookup; a key is a BAM's full path (as given), its basename
#                               ("x.bam"), or its output column name ("x"). Every BAM must match
#                               exactly one key (no default fallback); unused keys are allowed.
#   single string            -> path to a table with columns 'bam' (key) and 'ext'; read and used
#                               as the named lookup above (other columns ignored)
.resolve_ext <- function(ext, bam_path) {
  src <- ""
  if (is.character(ext)) {
    if (length(ext) != 1)
      stop("A character 'ext' must be a single path to a table with columns 'bam' and 'ext'.")
    if (!file.exists(ext))
      stop("'ext' is a string but no such file: ", ext)
    tab <- data.table::fread(ext)
    if (!all(c("bam", "ext") %in% names(tab)))
      stop("The 'ext' table must have columns 'bam' and 'ext': ", ext)
    src <- paste0(" from ", ext)
    ext <- setNames(tab$ext, as.character(tab$bam))
  }
  if (is.list(ext)) {
    if (!all(vapply(ext, function(e) is.numeric(e) && length(e) == 1, logical(1))))
      stop("Every element of an 'ext' list must be a single number.")
    nms <- names(ext)
    ext <- unlist(ext, use.names = FALSE)
    names(ext) <- nms
  }
  if (!is.numeric(ext) || length(ext) == 0)
    stop("'ext' must be a number, a vector/list of numbers, or a path to an ext table.")
  # value check applies only to the values actually used (a lookup table may hold unused rows)
  check_values <- function(v)
    if (any(!is.finite(v)) || any(v != round(v)))
      stop("'ext' values must be finite whole numbers (bp).")

  n   <- length(bam_path)
  nms <- names(ext)
  if (is.null(nms)) {
    check_values(ext)
    if (length(ext) == 1) return(rep(as.integer(ext), n))
    if (length(ext) != n)
      stop(sprintf(paste0("Unnamed 'ext' has length %d but there are %d BAMs. Give one number, ",
                          "one per BAM (in bam_path order), or a named vector/list."),
                   length(ext), n))
    return(as.integer(ext))
  }
  if (anyNA(nms) || any(nms == ""))
    stop("'ext' is partially named; name every element or none.")

  base <- basename(bam_path)
  col  <- tools::file_path_sans_ext(base)
  # hit[i, j]: key i refers to BAM j
  hit <- matrix(vapply(seq_len(n), function(j) nms %in% c(bam_path[j], base[j], col[j]),
                       logical(length(nms))),
                nrow = length(nms))
  problems <- c(
    if (any(rowSums(hit) > 1))
      paste0("keys matching several BAMs: ", paste(nms[rowSums(hit) > 1], collapse = ", ")),
    if (any(colSums(hit) == 0))
      paste0("BAMs with no 'ext' key: ", paste(bam_path[colSums(hit) == 0], collapse = ", ")),
    if (any(colSums(hit) > 1))
      paste0("BAMs given more than one 'ext' key: ", paste(bam_path[colSums(hit) > 1], collapse = ", "))
  )
  if (length(problems))
    stop("Named 'ext' must give every BAM exactly one key:\n  ", paste(problems, collapse = "\n  "))
  used <- ext[vapply(seq_len(n), function(j) which(hit[, j]), integer(1))]
  check_values(used)
  cat(sprintf("ext: used %d of %d keys%s\n", n, length(nms), src))
  as.integer(used)
}

# ---- hg38 blacklist helpers ('blacklist' = TRUE) ----
# The ENCODE hg38 blacklist v2 (Boyle Lab; Amemiya et al. 2019) ships with the package as
# inst/extdata/hg38-blacklist.v2.bed.gz.
.hg38_blacklist_file <- function() {
  system.file("extdata", "hg38-blacklist.v2.bed.gz", package = "multiEpiCore", mustWork = TRUE)
}

.read_hg38_blacklist <- function() {
  dt <- data.table::fread(.hg38_blacklist_file(), header = FALSE, select = 1:3)
  GRanges(seqnames = dt$V1, ranges = IRanges::IRanges(start = dt$V2 + 1L, end = dt$V3))
}

# 'bl' renamed to a seqlevels style ("UCSC" chr1, "NCBI" 1, ...) unless it already has it
.bl_in_style <- function(bl, style) {
  if (!is.na(style) && !style %in% seqlevelsStyle(bl)) seqlevelsStyle(bl) <- style
  bl
}

# TRUE for each range of 'gr' that does not overlap 'bl' by >= 1 bp. Whole ranges are kept or dropped,
# never trimmed (as bedtools intersect -v). 'bl' is first renamed to the naming style of 'gr' (chr1 <-> 1).
.not_blacklisted <- function(gr, bl) {
  style <- tryCatch(seqlevelsStyle(gr)[1], error = function(e) NA_character_)
  !overlapsAny(gr, .bl_in_style(bl, style), ignore.strand = TRUE)
}

# The hg38 blacklist is for hg38 data only: every BAM's chromosome lengths must be hg38 (names in the
# chr1 or the 1 style). Names alone cannot tell: mm10 chromosomes are also called chr1, chr2, ...
.check_bams_hg38 <- function(bam_path) {
  hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  len_hg <- seqlengths(hg38)[GenomeInfoDb::standardChromosomes(hg38)]
  len_hg <- c(len_hg, setNames(len_hg, sub("^chr", "", names(len_hg))))
  for (b in bam_path) {
    len <- scanBamHeader(b)[[1]]$targets
    common <- intersect(names(len), names(len_hg))
    if (length(common) == 0 || any(len[common] != len_hg[common]))
      stop("blacklist = TRUE uses the hg38 blacklist, but the chromosome lengths of this BAM are not hg38: ", b)
  }
}

# ---- Peak-directory helper ('regions' = a directory) ----
# Union every peak file (.bed / .narrowPeak / .broadPeak, optionally .gz) under 'dir' (recursive,
# sorted, excluding anything under out_dir) with build_peak_set(min_support = 1) into
# <out_dir>/<basename(dir)>.bed (natural chromosome order, then start) and return that path. Each BAM must have exactly one same-named
# peak file in its own folder under 'dir' (H3K27ac.bam <-> H3K27ac.bed.gz). Identical copies of a
# peak file (e.g. one per replicate) do not change a min_support = 1 union.
# With a blacklist (GRanges), blacklisted peaks are dropped file by file BEFORE the union, as the
# bl_filtered peak sets are made; a peak overlapping it would otherwise take its merged neighbours along.
.union_peak_dir <- function(dir, bam_path, out_dir, ref_genome, blacklist = NULL) {
  peak_re <- "\\.(bed|narrowPeak|broadPeak)(\\.gz)?$"
  peaks <- sort(list.files(dir, pattern = peak_re, recursive = TRUE, full.names = TRUE))
  out_abs <- normalizePath(out_dir, mustWork = FALSE)
  peaks <- peaks[!startsWith(normalizePath(peaks), paste0(out_abs, "/"))]
  if (length(peaks) == 0)
    stop("No peak files (.bed/.narrowPeak/.broadPeak, optionally .gz) found under: ", dir)

  # BAM <-> peak pairing: same file stem, same folder
  peak_stem <- sub(peak_re, "", basename(peaks))
  peak_dir  <- normalizePath(dirname(peaks))
  problems <- character(0)
  for (b in bam_path) {
    hit <- which(peak_stem == tools::file_path_sans_ext(basename(b)) &
                 peak_dir == normalizePath(dirname(b)))
    if (length(hit) != 1)
      problems <- c(problems, sprintf("%s: %d same-named peak file(s) next to it under %s",
                                      b, length(hit), dir))
  }
  if (length(problems))
    stop("With 'regions' = a directory, every BAM needs exactly one same-named peak file ",
         "in its own folder:\n  ", paste(problems, collapse = "\n  "))

  cat(sprintf("regions: union of %d peak files under %s (build_peak_set, min_support = 1)\n",
              length(peaks), dir))

  # Blacklist: build_peak_set reads filtered 3-column copies (it uses only the first 3 columns).
  # A file with no peak left is skipped (read_peak_file stops on an empty file).
  union_input <- peaks
  if (!is.null(blacklist)) {
    td <- tempfile("blacklist_peaks_")
    dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)
    union_input <- character(0)
    n_in <- n_kept <- 0
    for (i in seq_along(peaks)) {
      gr <- read_peak_file(peaks[i])
      keep <- .not_blacklisted(gr, blacklist)
      n_in <- n_in + length(gr)
      n_kept <- n_kept + sum(keep)
      if (any(keep)) {
        f <- file.path(td, sprintf("%04d.bed", i))
        data.table::fwrite(data.table::data.table(as.character(seqnames(gr))[keep], start(gr)[keep] - 1L,
                                                  end(gr)[keep]), f, sep = "\t", col.names = FALSE)
        union_input <- c(union_input, f)
      }
    }
    cat(sprintf("blacklist: removed %d of %d peaks (hg38-blacklist.v2, overlap >= 1 bp)\n",
                n_in - n_kept, n_in))
    if (length(union_input) == 0)
      stop("Every peak under ", dir, " overlaps the blacklist.")
  }

  bed <- build_peak_set(peak_path = union_input, pair = basename(normalizePath(dir)),
                        conditions = rep("all", length(union_input)), ref_genome = ref_genome,
                        out_dir = out_dir, window_size = NULL, min_support = 1L)

  # build_peak_set orders chromosomes by first appearance in the peak files (ENCODE: chr1, chr10, ...);
  # rewrite the union in natural order (chr1, chr2, ..., chrX, chrY), then by start
  dt <- data.table::fread(bed, header = FALSE)
  gr <- sort(GenomeInfoDb::sortSeqlevels(
    GRanges(seqnames = dt$V1, ranges = IRanges::IRanges(start = dt$V2 + 1L, end = dt$V3))))
  data.table::fwrite(data.table::data.table(as.character(seqnames(gr)), start(gr) - 1L, end(gr)),
                     bed, sep = "\t", col.names = FALSE)
  bed
}

# Build Count Matrix Function
# Post: Build a fragment-overlap count matrix from paired-end BAM files over user-specified genomic regions and save it as a Feather file.
# Parameter:
#   bam_path  : Character vector of BAM file paths. Each BAM file becomes one column in the output matrix.
#   regions   : Either
#               - a single integer (bin size in bp), e.g. regions = 5000, which tiles the genome into fixed-size bins; or
#               - a single file path to a BED / TSV / TXT / CSV file containing custom genomic regions; or
#               - a single directory: every peak file under it (.bed / .narrowPeak / .broadPeak, optionally
#                 .gz; recursive) is unioned with build_peak_set(min_support = 1) into
#                 <out_dir>/<dir name>.bed (rows in natural chromosome order: chr1, chr2, ..., chrX, chrY,
#                 then by start), which is then used as the region file. Each BAM must have
#                 exactly one same-named peak file in its own folder (H3K27ac.bam <-> H3K27ac.bed.gz).
#                 Calls that pass the same directory get the same rows.
#   out_dir  : Directory where the output Feather file will be written. Default "./".
#   ref_genome: Reference genome used when regions is numeric. One of "hg38" or "mm10". Ignored when custom regions are provided.
#   sample_name: Optional character. If provided, it is prepended to the output filename.
#   force_chr_coord: When TRUE, region IDs ("pos" column) are always  "CHR_start_end", even if gene_id is available. When FALSE and a non-empty gene_id column exists, gene_id is used as the region identifier.
#   seq_type  : "auto" (default), "PE" or "SE". "auto" decides from the first BAM (sorted by name), as
#               before, and then stops if any other BAM's flags say otherwise (a mixed PE/SE batch
#               would otherwise be miscounted). "PE"/"SE" declare the type; every BAM is checked
#               against it and a mismatch stops with an error.
#   ext       : SE read extension in bp (default 0 = no extension, the original behaviour). Ignored for
#               PE. > 0 extends each read from its 5' end toward 3' to 'ext' bp. One of: a single
#               number (every BAM); an unnamed vector / list (one per BAM, in bam_path order); a named
#               vector / list keyed by the BAM's full path, basename or column name (every BAM must
#               match exactly one key; unused keys are allowed); or a path to a table with columns
#               'bam' and 'ext' used as that keyed lookup.
#   blacklist : FALSE (default) or TRUE. TRUE drops every region (or bin) that overlaps the ENCODE hg38
#               blacklist v2 by >= 1 bp; whole regions are dropped, never trimmed. With 'regions' = a
#               directory the peaks are filtered file by file before the union. hg38 only: needs
#               ref_genome = "hg38" and BAMs whose chromosome lengths are hg38. File names do not change.
# Output: Writes a Feather file whose first column is 'pos' and remaining columns are fragment-overlap counts per BAM file. Returns the full output file path (character).

build_count_matrix <- function(bam_path, regions, out_dir = "./", ref_genome = "hg38", sample_name = NULL, force_chr_coord = FALSE,
                               seq_type = c("auto", "PE", "SE"), ext = 0L, blacklist = FALSE) {
  # Load Libraries
  suppressPackageStartupMessages({
    library(GenomicAlignments)
    library(Rsamtools)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(BSgenome.Mmusculus.UCSC.mm10)
    library(arrow)
    library(matrixStats)
    library(rtracklayer)
  })

  seq_type    <- match.arg(seq_type)
  ext_arg     <- ext        # captured here: the region block below reuses the name 'ext'
  bam_path_in <- bam_path   # order as given, for a positional 'ext'
  if (!is.logical(blacklist) || length(blacklist) != 1 || is.na(blacklist))
    stop("'blacklist' must be TRUE or FALSE.")

  # Set up parallel processing
  BPPARAM <- get_BPPARAM()
  
  # Create folder
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Optional hg38 blacklist, checked before any work is done
  bl <- NULL
  if (blacklist) {
    if (ref_genome != "hg38")
      stop("blacklist = TRUE is only available for hg38 (ref_genome = \"hg38\").")
    .check_bams_hg38(bam_path)
    bl <- .read_hg38_blacklist()
  }
  bl_regions <- bl  # applied to regions / bins below; NULL in directory mode, where the peaks are filtered

  # A directory of paired peak files -> their union as the region file (see .union_peak_dir)
  if (is.character(regions) && length(regions) == 1 && dir.exists(regions)) {
    regions <- .union_peak_dir(regions, bam_path, out_dir, ref_genome, blacklist = bl)
    bl_regions <- NULL
  }

  # Define Regions
  if (is.numeric(regions)) {
    if (length(regions) != 1) {
      stop("Error: numeric 'regions' must be a single bin size.")
    }
    BINSIZE <- regions
    use_custom_region <- FALSE
  } else if (is.character(regions) && all(file.exists(regions))) {
    if (length(regions) != 1) {
      stop("Error: Only one region file is allowed. Please provide a merged file.")
    }
    region_path <- regions
    ext <- tools::file_ext(region_path)
    if (all(ext == "tsv") || all(ext == "txt")) {
      region_df <- data.table::fread(region_path)
    } else if (all(ext == "csv")) {
      region_df <- read.csv(region_path, row.names = NULL)
    } else if (all(ext == "bed")) {
      region_df <- NULL
    } else {
      stop("Error: Invalid region file format. Only support .csv, .tsv, .txt, .bed")
    }
    use_custom_region <- TRUE
  } else {
    stop("Error: Custom regions must be provided. 'regions' argument is missing or invalid.")
  }

  # Detect chromosome naming style
  bam_path <- sort(bam_path) # sort by name
  bam_header <- scanBamHeader(bam_path[1])
  bam_seqinfo <- Seqinfo(seqnames = names(bam_header[[1]]$targets), seqlengths = bam_header[[1]]$targets)
  bam_style <- seqlevelsStyle(bam_seqinfo)[1]

  # Detect single or paired
  if (seq_type == "auto") {
    is_paired <- is_bam_paired(bam_path)
  } else {
    is_paired <- (seq_type == "PE")
  }
  # Every BAM must agree with that decision (first 1000 records each; BAMs with no records skipped)
  detected <- vapply(bam_path, .bam_paired_head, logical(1), USE.NAMES = FALSE)
  bad <- which(!is.na(detected) & detected != is_paired)
  if (length(bad))
    stop(sprintf("seq_type = '%s' treats this batch as %s, but the flags of these BAMs say %s:\n  %s\n%s",
                 seq_type, if (is_paired) "paired-end" else "single-end",
                 if (is_paired) "single-end" else "paired-end",
                 paste(bam_path[bad], collapse = "\n  "),
                 "A batch must be all PE or all SE; split it, or set seq_type to match the data."))

  # Per-BAM SE extension, aligned to the sorted bam_path (sort(x) == x[order(x)])
  ext_vec <- if (is_paired) rep(NA_integer_, length(bam_path)) else
    .resolve_ext(ext_arg, bam_path_in)[order(bam_path_in)]
  if (!is_paired && any(ext_vec > 0))
    cat("SE read extension (per BAM):\n",
        paste0("  ", tools::file_path_sans_ext(basename(bam_path)), " = ", ext_vec, "\n"), sep = "")

  # Generate Count Matrix for each chr
  # Process custom regions (no chromosome loop needed)
  if (use_custom_region) {
    # Handle CSV/GTF file
    if (!is.null(region_df)) {
      region_df <- fix_region_colnames(region_df)
      bin <- GRanges(seqnames = region_df$seqnames, ranges = IRanges(start = region_df$start, end = region_df$end))
      seqlevelsStyle(bin) <- bam_style
      binChriDataframe <- as.data.frame(bin)[, c("seqnames", "start", "end")]
      colnames(binChriDataframe)[1] <- "CHR"
      binChriDataframe$CHR <- as.character(binChriDataframe$CHR)

      if ("gene_id" %in% colnames(region_df)) {
        mcols(bin)$gene_id <- region_df$gene_id
        binChriDataframe$gene_id <- region_df$gene_id
      }
    }
    # Handle BED file
    else {
      bed_data <- read.table(region_path, sep = "\t", stringsAsFactors = FALSE)
      bin <- GRanges(seqnames = bed_data$V1, ranges = IRanges(start = bed_data$V2 + 1, end = bed_data$V3))
      if (ncol(bed_data) >= 4) {
        mcols(bin)$gene_id <- bed_data$V4
      }
      seqlevelsStyle(bin) <- bam_style
      binChriDataframe <- as.data.frame(bin)[, c("seqnames", "start", "end")]
      colnames(binChriDataframe)[1] <- "CHR"
      binChriDataframe$CHR <- as.character(binChriDataframe$CHR)

      if (!is.null(mcols(bin)) && "gene_id" %in% colnames(mcols(bin))) {
        binChriDataframe$gene_id <- mcols(bin)$gene_id
      }
    }

    # Blacklist: drop whole regions that overlap it (the other rows keep their order)
    if (!is.null(bl_regions)) {
      keep <- .not_blacklisted(bin, bl_regions)
      cat(sprintf("blacklist: removed %d of %d regions (hg38-blacklist.v2, overlap >= 1 bp)\n",
                  sum(!keep), length(keep)))
      if (!any(keep))
        stop("Every region overlaps the blacklist.")
      bin <- bin[keep]
      binChriDataframe <- binChriDataframe[keep, , drop = FALSE]
      rownames(binChriDataframe) <- NULL
    }

    # Process BAM files for custom regions (parallelized by chromosome)
    bin_chr <- as.character(seqnames(bin))
    chr_list <- unique(bin_chr)
    
    result_list <- bplapply(chr_list, function(chr_i) {
      bin_idx <- which(bin_chr == chr_i)
      bin_i <- bin[bin_idx]
      chrSizei <- seqlengths(bam_seqinfo)[chr_i]
      param <- ScanBamParam(which = GRanges(chr_i, IRanges(1, chrSizei)))

      sub_df <- binChriDataframe[bin_idx, , drop = FALSE]

      for (k in seq_along(bam_path)) {
        bam <- bam_path[k]
        overlapCount <- numeric(length(bin_i))

        if (is_paired) {
          temp <- readGAlignmentPairs(bam, param = param)
          if (length(temp) == 0) {
            bamContent <- GRanges()
          } else {
            locus <- data.frame(
              first_start = start(temp@first),
              first_end   = end(temp@first),
              last_start  = start(temp@last),
              last_end    = end(temp@last)
            )
            frag_start <- rowMin(as.matrix(locus))
            frag_end   <- rowMax(as.matrix(locus))
            bamContent <- makeGRangesFromDataFrame(data.frame(
              seqnames = as.vector(seqnames(temp)), strand = "*",
              start = frag_start, end = frag_end
            ))
          }
        } else {
          temp <- readGAlignments(bam, param = param)
          if (length(temp) == 0) {
            bamContent <- GRanges()
          } else if (ext_vec[k] > 0) {
            # SE read extension: keep the 5' end, extend toward 3' to ext bp, clip at chromosome ends
            bamContent <- trim(resize(as(temp, "GRanges"), width = ext_vec[k], fix = "start"))
            strand(bamContent) <- "*"
          } else {
            bamContent <- GRanges(
              seqnames = seqnames(temp),
              ranges   = IRanges(start = start(temp), end = end(temp)),
              strand   = "*"
            )
          }
        }

        if (length(bamContent) > 0) {
          overlaps <- findOverlaps(bamContent, bin_i, ignore.strand = TRUE)
          qh <- queryHits(overlaps)
          sh <- subjectHits(overlaps)

          fragment_starts  <- start(bamContent)[qh]
          fragment_ends    <- end(bamContent)[qh]
          fragment_lengths <- fragment_ends - fragment_starts + 1

          bin_starts <- start(bin_i)[sh]
          bin_ends   <- end(bin_i)[sh]

          overlap_starts  <- pmax(fragment_starts, bin_starts)
          overlap_ends    <- pmin(fragment_ends, bin_ends)
          overlap_lengths <- pmax(0, overlap_ends - overlap_starts + 1)
          proportions     <- overlap_lengths / fragment_lengths

          if (length(sh) > 0) {
            sums <- rowsum(proportions, group = sh)[, 1]
            overlapCount[as.integer(names(sums))] <- sums
          }
        }

        bamName <- tools::file_path_sans_ext(basename(bam))
        sub_df[[bamName]] <- overlapCount
      }
      return(sub_df)
    }, BPPARAM = BPPARAM)
    binChriDataframe <- as.data.frame(data.table::rbindlist(result_list))

    # Format output
    has_gene_id <- "gene_id" %in% colnames(binChriDataframe) &&
      all(!is.na(binChriDataframe$gene_id)) &&
      all(binChriDataframe$gene_id != "")
    if (has_gene_id && !force_chr_coord) {
      pos_vec <- binChriDataframe$gene_id
      drop_cols <- c("CHR", "start", "end", "gene_id")
    } else {
      pos_vec <- paste0(binChriDataframe$CHR, "_", binChriDataframe$start, "_", binChriDataframe$end)
      drop_cols <- c("CHR", "start", "end")
    }
    pos_df <- data.frame(pos = pos_vec)
    tmp_wgc <- binChriDataframe[, !(names(binChriDataframe) %in% drop_cols), drop = FALSE]
    binChriDataframe_full <- cbind(pos_df, tmp_wgc)
    
  } else {
    # Get reference genome size
    if (ref_genome == "hg38") {
      chrSizes0 <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
    } else if (ref_genome == "mm10") {
      chrSizes0 <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
    } else {
      stop("Error: 'ref_genome' must be either 'hg38' or 'mm10'.")
    }

    ref_seqinfo <- Seqinfo(
      seqnames = names(chrSizes0),
      seqlengths = chrSizes0
    )
    seqlevelsStyle(ref_seqinfo) <- bam_style

    chr_list <- GenomeInfoDb::standardChromosomes(ref_seqinfo)
    chr_list <- chr_list[!tolower(chr_list) %in% c("mt", "chrm", "m", "mito")]
    chrSizes <- seqlengths(ref_seqinfo)[chr_list]

    # Blacklist in the BAM's naming style, ready for the workers (NULL = no blacklist)
    bl_bins <- if (!is.null(bl_regions)) .bl_in_style(bl_regions, bam_style)

    # Process fixed bins with parallel chromosome processing
    binChriDataframe_list <- bplapply(chr_list, function(chr_i) {
      chrSizei <- chrSizes[chr_i]
      bin <- tileGenome(chrSizei, tilewidth = BINSIZE, cut.last.tile.in.chrom = TRUE)
      if (!is.null(bl_bins)) {
        bin <- bin[!overlapsAny(bin, bl_bins, ignore.strand = TRUE)]
        if (length(bin) == 0) return(NULL)  # every bin of this chromosome is blacklisted
      }
      binChriDataframe <- as.data.frame(bin)[, c("start", "end")]
      chr_df <- data.frame(CHR = names(chrSizei), stringsAsFactors = FALSE)
      binChriDataframe <- cbind(chr_df, binChriDataframe)

      for (k in seq_along(bam_path)) {
        bam <- bam_path[k]
        param <- ScanBamParam(which = GRanges(chr_i, IRanges(1, chrSizei)))
        overlapCount <- numeric(length(bin))

        if (is_paired) {
          temp <- readGAlignmentPairs(bam, param = param)
          if (length(temp) == 0) {
            bamContent <- GRanges()
          } else {
            locus <- data.frame(
              first_start = start(temp@first),
              first_end   = end(temp@first),
              last_start  = start(temp@last),
              last_end    = end(temp@last)
            )
            frag_start <- rowMin(as.matrix(locus))
            frag_end   <- rowMax(as.matrix(locus))
            bamContent <- makeGRangesFromDataFrame(data.frame(
              seqnames = as.vector(seqnames(temp)), strand = "*",
              start = frag_start, end = frag_end
            ))
          }
        } else {
          temp <- readGAlignments(bam, param = param)
          if (length(temp) == 0) {
            bamContent <- GRanges()
          } else if (ext_vec[k] > 0) {
            # SE read extension: keep the 5' end, extend toward 3' to ext bp, clip at chromosome ends
            bamContent <- trim(resize(as(temp, "GRanges"), width = ext_vec[k], fix = "start"))
            strand(bamContent) <- "*"
          } else {
            bamContent <- GRanges(
              seqnames = seqnames(temp),
              ranges   = IRanges(start = start(temp), end = end(temp)),
              strand   = "*"
            )
          }
        }

        if (length(bamContent) > 0) {
          overlaps <- findOverlaps(bamContent, bin, ignore.strand = TRUE)
          qh <- queryHits(overlaps)
          sh <- subjectHits(overlaps)

          fragment_starts  <- start(bamContent)[qh]
          fragment_ends    <- end(bamContent)[qh]
          fragment_lengths <- fragment_ends - fragment_starts + 1

          bin_starts <- start(bin)[sh]
          bin_ends   <- end(bin)[sh]

          overlap_starts  <- pmax(fragment_starts, bin_starts)
          overlap_ends    <- pmin(fragment_ends, bin_ends)
          overlap_lengths <- pmax(0, overlap_ends - overlap_starts + 1)
          proportions     <- overlap_lengths / fragment_lengths

          if (length(sh) > 0) {
            summed <- aggregate(proportions, by = list(bin_id = sh), FUN = sum)
            overlapCount[summed$bin_id] <- summed$x
          }
        }

        bamName <- tools::file_path_sans_ext(basename(bam))
        binChriDataframe[[bamName]] <- overlapCount
      }

      tmp_pos <- binChriDataframe[, c("CHR", "start", "end")]
      tmp_pos$pos <- paste0(tmp_pos$CHR, "_", tmp_pos$start, "_", tmp_pos$end)
      pos_df <- data.frame(pos = tmp_pos$pos)
      tmp_wgc <- binChriDataframe[, !(names(binChriDataframe) %in% c("CHR", "start", "end"))]
      binChriDataframe_final <- cbind(pos_df, tmp_wgc)

      return(binChriDataframe_final)
    }, BPPARAM = BPPARAM)

    binChriDataframe_full <- as.data.frame(data.table::rbindlist(binChriDataframe_list))
    if (!is.null(bl_regions)) {
      n_bins <- sum(ceiling(chrSizes / BINSIZE))
      cat(sprintf("blacklist: removed %d of %d bins (hg38-blacklist.v2, overlap >= 1 bp)\n",
                  n_bins - nrow(binChriDataframe_full), n_bins))
    }
  }

  if (is.numeric(regions)) {
    filename <- paste0("Count_Matrix_", BINSIZE)
  } else if (is.character(regions)) {
    prefix <- basename(tools::file_path_sans_ext(regions))
    filename <- paste0("Count_Matrix_", prefix[1])
  }

  if (is.null(sample_name)) {
    output_filename <- paste0(filename, ".feather")
  } else {
    output_filename <- paste0(sample_name, "_", filename, ".feather")
  }

  # Report
  output_path <- file.path(out_dir, output_filename)
  write_feather(binChriDataframe_full, output_path)
  cat("Successfully saved to: ", output_path, "\n")

  # Return
  # mat <- as.matrix(binChriDataframe_full[, -1])
  # rownames(mat) <- binChriDataframe_full$pos
  # invisible(list(mat = mat, path = output_path))
  return(output_path)
}