# Link each peak to every gene whose window it overlaps, and write one
# peak x gene BED file per sample: <out_dir>/<pair>_peak_gene.bed
#
# Gene library: GENCODE basic, one representative transcript per gene
#   (download_rds("GENCODE_gene_<ref_genome>.rds")), filtered by gene_types
#   (NULL or "all" keeps every biotype). Chromosome naming is converted to
#   the style of the peaks (chr1 vs 1).
#
# Gene window (1-based closed, strand-aware; left end clamped to 1):
#   + strand (TSS = start): [start - upstream,   end + downstream]
#   - strand (TSS = end):   [start - downstream, end + upstream]
#   A peak is linked to a gene if it shares >= 1 bp with the window
#   (ignore.strand = TRUE). A peak overlapping k windows gives k rows;
#   peaks overlapping no window are not written.
#
# Output columns (tab-separated, header line starting with "#"):
#   #chrom       peak chromosome
#   start        peak start, 0-based (= original narrowPeak chromStart)
#   end          peak end (= original narrowPeak chromEnd)
#   name         original narrowPeak column 4; "chr:start-end" if absent
#   score        original narrowPeak column 5; "." if absent
#   strand       peak strand; GRanges "*" written as BED "."
#   gene_id      GENCODE gene ID (versioned)
#   gene_name    gene symbol
#   gene_type    GENCODE biotype
#   dist_to_tss  signed distance from peak summit s to TSS along transcription:
#                  + strand: s - TSS;  - strand: TSS - s
#                  < 0 = upstream of TSS, 0 = at TSS, > 0 = into gene body
#                s = start(gr) + narrowPeak column 10 offset (1-based);
#                peak midpoint if the offset is absent, -1 or outside the peak.
#                Linking uses the whole peak but the distance uses the summit,
#                so |dist_to_tss| can slightly exceed the window size.
#
# Row order: peaks in genome order; within a peak, genes by |dist_to_tss|.
#
# query: GRanges or named GRangesList (names become output file prefixes),
#        1-based coordinates as returned by read_peak_to_grl(); optional
#        mcols name / score / peak are used as described above.
# Returns (invisibly) a named list of output file paths.
gene_annotation <- function(query, out_dir = "./", ref_genome = "hg38",
                            gene_types = "protein_coding",
                            upstream = 20000, downstream = 0) {
  # ---- Parameter validation ----
  if (!ref_genome %in% c("hg38", "mm10")) {
    stop("ref_genome must be 'hg38' or 'mm10'")
  }
  check_window <- function(x, nm) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0) {
      stop(nm, " must be a single non-negative number (bp).")
    }
    as.integer(round(x))
  }
  upstream   <- check_window(upstream, "upstream")
  downstream <- check_window(downstream, "downstream")
 
  if (methods::is(query, "GRanges")) {
    query <- GenomicRanges::GRangesList(query = query)
  }
  if (!methods::is(query, "GRangesList")) {
    stop("query must be a GRanges or a named GRangesList.")
  }
  if (is.null(names(query)) || anyNA(names(query)) || any(names(query) == "") ||
      anyDuplicated(names(query))) {
    stop("query must be a GRangesList with unique, non-empty names (used as output file prefixes).")
  }
 
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
 
  # ---- Gene library ----
  gene_library_file <- switch(
    ref_genome,
    "hg38" = download_rds("GENCODE_gene_hg38.rds"),
    "mm10" = download_rds("GENCODE_gene_mm10.rds")
  )
  genes <- readRDS(gene_library_file)
 
  # gene_types: NULL or "all" keeps every biotype
  if (!is.null(gene_types) && !identical(gene_types, "all")) {
    unknown <- setdiff(gene_types, unique(genes$gene_type))
    if (length(unknown) > 0L) {
      warning("gene_types not found in the gene library: ", paste(unknown, collapse = ", "))
    }
    genes <- genes[genes$gene_type %in% gene_types]
    if (length(genes) == 0L) {
      stop("No genes left after filtering by gene_types: ", paste(gene_types, collapse = ", "))
    }
  }
 
  # Match chromosome naming style to the peaks (chr1 vs 1)
  q_style <- tryCatch(GenomeInfoDb::seqlevelsStyle(query)[1L], error = function(e) NA_character_)
  if (!is.na(q_style)) {
    tryCatch(
      GenomeInfoDb::seqlevelsStyle(genes) <- q_style,
      error = function(e) warning("Could not convert gene library to '", q_style,
                                  "' chromosome style: ", conditionMessage(e))
    )
  }
 
  # ---- Gene windows (1-based closed, strand-aware) ----
  g_start  <- BiocGenerics::start(genes)
  g_end    <- BiocGenerics::end(genes)
  is_minus <- as.character(BiocGenerics::strand(genes)) == "-"   # "*" treated as "+"
 
  tss <- ifelse(is_minus, g_end, g_start)
  win_start <- as.integer(pmax(1L, ifelse(is_minus, g_start - downstream, g_start - upstream)))
  win_end   <- as.integer(ifelse(is_minus, g_end + upstream, g_end + downstream))
 
  windows <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(genes),
    ranges   = IRanges::IRanges(start = win_start, end = win_end),
    strand   = "*"
  )
 
  out_cols <- c("#chrom", "start", "end", "name", "score", "strand",
                "gene_id", "gene_name", "gene_type", "dist_to_tss")
 
  write_bed <- function(dt, pair) {
    output_file <- file.path(out_dir, paste0(pair, "_peak_gene.bed"))
    data.table::fwrite(dt, output_file, sep = "\t", quote = FALSE,
                       col.names = TRUE, na = ".")
    output_file
  }
 
  empty_dt <- function() {
    dt <- data.table::data.table(
      chrom = character(0), start = integer(0), end = integer(0),
      name = character(0), score = numeric(0), strand = character(0),
      gene_id = character(0), gene_name = character(0), gene_type = character(0),
      dist_to_tss = integer(0)
    )
    data.table::setnames(dt, out_cols)
    dt
  }
 
  # ---- Annotate one sample ----
  annotate_one <- function(peaks, pair) {
    if (length(peaks) == 0L) {
      message(sprintf("[%s] no peaks", pair))
      return(write_bed(empty_dt(), pair))
    }
 
    p_start <- BiocGenerics::start(peaks)
    p_end   <- BiocGenerics::end(peaks)
    p_chr   <- as.character(GenomicRanges::seqnames(peaks))
    mc      <- S4Vectors::mcols(peaks)
 
    # Summit (1-based): narrowPeak chromStart(0-based) + offset + 1 = start(gr) + offset
    p_mid <- as.integer((p_start + p_end) %/% 2L)
    summit <- p_mid
    if ("peak" %in% names(mc) && is.numeric(mc$peak)) {
      off <- as.integer(mc$peak)
      ok  <- !is.na(off) & off >= 0L & (p_start + off) <= p_end
      summit[ok] <- p_start[ok] + off[ok]
    }
 
    # name / score: keep the original narrowPeak columns 4 and 5 as-is.
    # Fallbacks only for BED3-style input: name -> "chr:start0-end", score -> "." (NA)
    p_name <- if ("name" %in% names(mc)) as.character(mc$name) else rep(NA_character_, length(peaks))
    miss <- is.na(p_name)
    p_name[miss] <- paste0(p_chr[miss], ":", p_start[miss] - 1L, "-", p_end[miss])
    p_score <- if ("score" %in% names(mc)) mc$score else rep(NA_real_, length(peaks))
    p_strand <- as.character(BiocGenerics::strand(peaks))
    p_strand[p_strand == "*"] <- "."
 
    hits <- GenomicRanges::findOverlaps(peaks, windows, ignore.strand = TRUE)
    if (length(hits) == 0L) {
      message(sprintf("[%s] 0 of %d peaks linked to a gene", pair, length(peaks)))
      return(write_bed(empty_dt(), pair))
    }
    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)
 
    ps <- p_start[qh]; pe <- p_end[qh]
    gm <- is_minus[sh]; t0 <- tss[sh]
 
    dist_to_tss <- as.integer(ifelse(gm, t0 - summit[qh], summit[qh] - t0))
 
    dt <- data.table::data.table(
      chrom       = p_chr[qh],
      start       = ps - 1L,          # back to 0-based BED
      end         = pe,
      name        = p_name[qh],
      score       = p_score[qh],
      strand      = p_strand[qh],
      gene_id     = genes$gene_id[sh],
      gene_name   = genes$gene_name[sh],
      gene_type   = genes$gene_type[sh],
      dist_to_tss = dist_to_tss
    )
 
    # Order: genome order of peaks, then nearest gene first within a peak
    chrom_rank <- match(dt$chrom, GenomeInfoDb::seqlevels(peaks))
    dt <- dt[order(chrom_rank, dt$start, dt$end, abs(dt$dist_to_tss))]
    data.table::setnames(dt, out_cols)
 
    n_linked <- length(unique(qh))
    message(sprintf("[%s] %d of %d peaks linked to %d gene%s (%d rows)",
                    pair, n_linked, length(peaks), length(unique(sh)),
                    if (length(unique(sh)) == 1L) "" else "s", nrow(dt)))
 
    write_bed(dt, pair)
  }
 
  result_list <- lapply(names(query), function(pair) annotate_one(query[[pair]], pair))
  names(result_list) <- names(query)
  invisible(result_list)
}