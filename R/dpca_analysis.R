# ---- dPCA analysis ----
# One call: dPCA from the per-sample count matrices, then every figure and table drawn from it, written
# into one folder. The steps live in their own environment, .dpca_tools, so that only dpca_analysis()
# is exported; the dPCA engine is .dpca_engine (dpca_engine.R).
.dpca_tools <- local({

# Step 1 (dPCA) then step 2 (downstream); returns tables, plots (pdf + png) and plots_png (png only)
dpca_comprehensive <- function(engine, counts, de, gene_windows, config) {
  dpca <- dpca_from_counts(engine, counts, config$marks, config$groups)
  res <- downstream(dpca, de, gene_windows, config)
  res$tables <- c(dpca, res$tables)
  res
}

# ---- Step 1: dPCA ----
# Every library (one mark of one sample) is scaled to the smallest total count over the windows,
# then the engine runs at the published defaults (log2(x + 1), columns of D centred, pooled
# variance, |z| capped at 7.5) in double precision. group_id 1 = groups[1], so D = groups[1] - groups[2].
dpca_from_counts <- function(engine, counts, marks, groups) {
  N <- colSums(counts$x)
  fit <- engine$dpca_run(x = counts$x * rep(min(N) / N, each = nrow(counts$x)),
                         group_id = match(counts$group, groups), dataset_id = match(counts$mark, marks),
                         transform = TRUE, single_precision = FALSE)
  V <- fit$V
  dimnames(V) <- list(marks, sprintf("dPC%d", seq_len(ncol(V))))
  D <- fit$D
  colnames(D) <- marks
  pw <- data.frame(counts$loci, stats::setNames(as.data.frame(D), paste0("d_", marks)), check.names = FALSE)
  for (j in seq_len(ncol(V))) {
    pw[[sprintf("dPC%d_beta", j)]] <- fit$beta[, j]
    pw[[sprintf("dPC%d_z", j)]]    <- fit$t[, j]
    pw[[sprintf("dPC%d_q", j)]]    <- fit$qvalue[, j]
  }
  ut <- which(upper.tri(diag(length(marks))), arr.ind = TRUE)
  cr <- list(spearman = stats::cor(D, method = "spearman"), pearson = stats::cor(D))
  list(per_window       = pw,
       loadings         = data.frame(mark = marks, V, check.names = FALSE, row.names = NULL),
       components       = data.frame(dPC = colnames(V), eigenvalue = fit$eigenvalue,
                                     var_explained = fit$percent_var, SNR = unname(fit$SNR),
                                     sig_prc = 1 - unname(fit$pi0)),
       mark_correlation = data.frame(mark1 = marks[ut[, 1]], mark2 = marks[ut[, 2]],
                                     lapply(cr, function(m) m[ut])))
}

# ---- Reading the inputs (inputs are checked here, where they are read) ----
# Count matrices: one feather per sample with 'pos' (chr_start_end, 1-based) and one column per mark;
# every file must list the same windows in the same order. marks = NULL takes every column but
# 'pos' of the first file, in file order. Columns of x: sample by sample, marks within a sample.
read_counts <- function(cm_path, conditions, marks, groups) {
  if (length(cm_path) != length(conditions))
    stop("cm_path and conditions must have the same length")
  if (length(groups) != 2L || !all(conditions %in% groups) ||
      any(table(factor(conditions, groups)) < 2L))
    stop("conditions must fall into the two groups of `contrast`, with at least two samples per group")
  tabs <- lapply(cm_path, function(f) as.data.frame(arrow::read_feather(f)))
  if (is.null(marks)) marks <- setdiff(names(tabs[[1]]), "pos")
  for (i in seq_along(tabs)) {
    miss <- setdiff(c("pos", marks), names(tabs[[i]]))
    if (length(miss)) stop("columns missing in ", cm_path[i], ": ", paste(miss, collapse = ", "))
  }
  pos <- tabs[[1]]$pos
  for (tab in tabs[-1]) if (!identical(tab$pos, pos)) stop("the count files do not list the same windows in the same order")
  x <- unname(do.call(cbind, lapply(tabs, function(tab) as.matrix(tab[marks]))))
  if (anyNA(x) || any(x < 0)) stop("NA or negative counts")
  m <- regmatches(pos, regexec("^(.+)_([0-9]+)_([0-9]+)$", pos))
  bad <- lengths(m) != 4L
  if (any(bad)) stop(sum(bad), " window ids are not of the form chr_start_end; first: ", pos[bad][1])
  list(x = x,
       marks = marks,
       mark  = rep(marks, length(cm_path)),
       group = rep(conditions, each = length(marks)),
       loci  = data.frame(chr   = vapply(m, `[`, "", 2L),
                          start = as.integer(vapply(m, `[`, "", 3L)),
                          end   = as.integer(vapply(m, `[`, "", 4L))))
}

# DESeq2 table: rows whose gene_id has the species' Ensembl prefix (ENSG / ENSMUSG), ids without version;
# log2FC is turned into groups[1] - groups[2] by de_sign (+1 or -1)
read_de <- function(de_file, de_sign, id_prefix) {
  de <- utils::read.delim(de_file)
  miss <- setdiff(c("gene_id", "log2FoldChange", "lfcSE", "pvalue", "padj", "baseMean"), names(de))
  if (length(miss)) stop("columns missing in ", de_file, ": ", paste(miss, collapse = ", "))
  de <- de[grepl(id_prefix, de$gene_id), ]
  id <- normalize_gene_id(de$gene_id)
  if (anyDuplicated(id))
    stop("gene ids repeated in ", de_file, ": ", paste(utils::head(unique(id[duplicated(id)]), 5), collapse = ", "))
  data.frame(gene_id = id, log2FC = de_sign * de$log2FoldChange, lfcSE = de$lfcSE,
             pvalue = de$pvalue, padj = de$padj, baseMean = de$baseMean)
}

# Gene windows from a GENCODE GTF: protein-coding genes on the standard chromosomes (chrM dropped),
# window = [TSS - 2000, TES] on the gene's strand, trimmed to the chromosome (chrom_sizes: named
# lengths from the BSgenome of ref_genome). Same windows as winA of ENCODE/tools/gene_windows.R.
read_gene_windows <- function(gtf_file, chrom_sizes) {
  suppressPackageStartupMessages({
    library(rtracklayer); library(GenomicRanges); library(GenomeInfoDb); library(S4Vectors)
  })

  gtf   <- rtracklayer::import(gtf_file)
  genes <- gtf[gtf$type == "gene"]
  genes <- genes[genes$gene_type == "protein_coding"]
  if (length(genes) == 0L) stop("no genes left after biotype filter: protein_coding")

  genes <- keepStandardChromosomes(genes, pruning.mode = "coarse")
  if ("chrM" %in% seqlevels(genes)) genes <- dropSeqlevels(genes, "chrM", pruning.mode = "coarse")

  sl <- chrom_sizes
  common <- intersect(seqlevels(genes), names(sl))
  if (!length(common)) stop("genes and chrom.sizes share no seqlevels - check chr naming.")
  genes <- keepSeqlevels(genes, common, pruning.mode = "coarse")
  seqlengths(genes) <- sl[seqlevels(genes)]

  tss  <- resize(genes, width = 1L, fix = "start")
  prom <- promoters(genes, upstream = 2000L, downstream = 0L)
  w    <- trim(punion(granges(genes), granges(prom), fill.gap = TRUE))
  data.frame(chr = as.character(seqnames(w)), start = start(w), end = end(w),
             strand = as.character(strand(w)), gencode_id = genes$gene_id,
             gene_name = genes$gene_name, tss = start(tss))
}

# Write every table (tsv; per_window as feather) and figure (pdf + png, 7 x 7 in; RNA plots png only)
write_outputs <- function(res, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (n in names(res$tables)) {
    if (n == "per_window") arrow::write_feather(res$tables[[n]], file.path(out_dir, "per_window.feather"))
    else utils::write.table(res$tables[[n]], file.path(out_dir, paste0(n, ".tsv")),
                            sep = "\t", quote = FALSE, row.names = FALSE)
  }
  for (n in names(res$plots)) for (f in c("pdf", "png"))
    ggplot2::ggsave(file.path(out_dir, paste0(n, ".", f)), res$plots[[n]], width = 7, height = 7)
  for (n in names(res$plots_png))
    ggplot2::ggsave(file.path(out_dir, paste0(n, ".png")), res$plots_png[[n]], width = 7, height = 7)
}

# ---- Step 2: figures and tables from the dPCA result ----
# Always: 1A, 1B, 1C, mark correlation heatmap, dPC scatter plot and their tables.
# With gene_windows (GTF): window_gene. With de as well: de, the gene lists and the RNA plots.
downstream <- function(dpca, de, gene_windows, config) {
  pw <- dpca$per_window
  retained <- retained_components(dpca$components, config$retain_min_var)

  scaled <- scaled_loadings(dpca$loadings, dpca$components)
  variance <- variance_table(dpca$components, config$retain_min_var)
  cl <- cluster_marks(scaled, retained, config$linkage)

  axes <- c(config$scatter_x, config$scatter_y)
  beta <- paste0(axes, "_beta")
  grid <- density_grid(pw[[beta[1]]], pw[[beta[2]]], config$density_bins)
  stem <- paste0(axes[1], "_vs_", axes[2])

  tables <- list(scaled_loadings           = scaled,
                 variance_explained        = variance,
                 silhouette_by_k           = cl$silhouette_by_k,
                 mark_clusters             = cl$mark_clusters,
                 denoised_correlation      = denoised_correlation(dpca$loadings, dpca$components,
                                                                  dpca$mark_correlation[c("mark1", "mark2")]))

  # genes of the windows (needs the GTF); gene lists and RNA plots (need the DE table as well)
  if (!is.null(gene_windows)) {
    wg <- assign_genes(pw[c("chr", "start", "end")], gene_windows)
    if (!is.null(de)) tables$de <- de
    tables$window_gene <- wg
  }
  if (!is.null(de)) {
    walk <- walk_windows(pw, wg, de, config$rank_dpc, config$copy_number_chr,
                         config$padj_max, config$log2fc_min, config$basemean_min)
    tw <- c("window_rank", "chr", "start", "end", grep("_beta$", names(walk), value = TRUE),
            "copy_like", "gene_id", "gene_name", "log2FC", "padj", "baseMean", "reason")

    rna <- rna_regions(pw, wg, de, config$marks, retained, config$q_max, config$copy_number_chr)
    thr <- rna_beta_thresholds(pw, rna$regions, retained, config$q_max)
    sp <- rna_spearman(rna$regions, config$marks, retained)

    tables <- c(tables,
                list(track_genes_copy_gate_on  = track_genes(walk, TRUE),
                     track_genes_copy_gate_off = track_genes(walk, FALSE),
                     track_windows             = walk[seq_len(min(config$report_windows, nrow(walk))), tw],
                     rna_regions               = rna$regions,
                     rna_spearman              = sp,
                     rna_beta_thresholds       = thr,
                     rna_counts                = rna$counts))
  }
  tables[[paste0(stem, "_density")]] <- stats::setNames(grid, c(beta, "count"))

  plots <- list(loadings_clustered_heatmap = plot_loadings_heatmap(scaled, cl$mark_clusters,
                                                                   dendrogram_segments(cl$hc), retained),
                optimal_number_of_dpcs     = plot_variance(variance, config$retain_min_var),
                optimal_number_of_clusters = plot_silhouette(cl$silhouette_by_k),
                marks_correlation_heatmap  = plot_mark_correlation(dpca$mark_correlation,
                                                                   mark_order(dpca$loadings, dpca$components$dPC[1]),
                                                                   config$correlation))
  plots[[paste0(stem, "_scatter")]] <- plot_density(grid, config$density_bins, axes)
  plots_png <- list()
  if (!is.null(de)) {
    plots_png <- lapply(retained, function(k)
      plot_rna_projection(rna$regions, k, thr$beta_threshold[thr$dPC == k], sp, config$q_max))
    names(plots_png) <- sprintf("rna_l2fc_vs_%s_projection", retained)
  }

  list(tables = tables, plots = plots, plots_png = plots_png)
}

# Retained dPCs: var_explained >= min_var
retained_components <- function(components, min_var) {
  components$dPC[components$var_explained >= min_var]
}

scaled_loadings <- function(loadings, components) {
  for (k in names(loadings)[-1])
    loadings[[k]] <- loadings[[k]] * sqrt(components$eigenvalue[components$dPC == k])
  loadings
}

variance_table <- function(components, min_var) {
  data.frame(dPC = components$dPC, var_explained = components$var_explained,
             percent = 100 * components$var_explained,
             retained = components$dPC %in% retained_components(components, min_var))
}

# Marks clustered on their scaled loadings (loading x sqrt(eigenvalue)) of the retained dPCs;
# k_best = largest mean silhouette over k = 2 .. marks - 1 (the smaller k on a tie)
cluster_marks <- function(scaled, dpcs, linkage) {
  x <- as.matrix(scaled[dpcs])
  rownames(x) <- scaled$mark
  d <- stats::dist(x)
  hc <- stats::hclust(d, method = linkage)
  k <- 2:(nrow(x) - 1)
  cut <- lapply(k, stats::cutree, tree = hc)
  sil <- lapply(cut, function(cl) unname(cluster::silhouette(cl, d)[, "sil_width"]))
  m <- vapply(sil, mean, 0)
  best <- which.max(m)
  cl <- cut[[best]][hc$order]
  list(hc = hc,
       silhouette_by_k = data.frame(
         k = k, mean_silhouette = m,
         sizes = vapply(cut, function(c) paste(tabulate(c), collapse = "/"), ""),
         best = seq_along(k) == best),
       mark_clusters = data.frame(
         mark = names(cl), cluster = unname(cl),
         silhouette = sil[[best]][hc$order], dendrogram_order = seq_along(cl)))
}

# Line segments of the hclust dendrogram, for drawing it next to the heatmap
dendrogram_segments <- function(hc) {
  pos <- order(hc$order)
  mid <- numeric(nrow(hc$merge))
  seg <- NULL
  for (i in seq_len(nrow(hc$merge))) {
    a <- hc$merge[i, ]
    x <- ifelse(a < 0, pos[abs(a)], mid[abs(a)])
    y <- ifelse(a < 0, 0, hc$height[abs(a)])
    mid[i] <- mean(x)
    seg <- rbind(seg, data.frame(x = c(x, x[1]), y = c(y, hc$height[i]),
                                 xend = c(x, x[2]), yend = hc$height[i]))
  }
  seg
}

mark_order <- function(loadings, dpc) loadings$mark[order(-loadings[[dpc]])]

# Correlation of the marks implied by the dPCs: cov2cor(V diag(lambda) V')
denoised_correlation <- function(loadings, components, pairs) {
  v <- as.matrix(loadings[-1])
  rownames(v) <- loadings$mark
  e <- stats::cov2cor(v %*% (components$eigenvalue[match(colnames(v), components$dPC)] * t(v)))
  data.frame(mark1 = pairs$mark1, mark2 = pairs$mark2, cor = e[cbind(pairs$mark1, pairs$mark2)])
}

# 2-D histogram of two betas: bins x bins cells over their ranges; non-empty cells only
density_grid <- function(x, y, bins) {
  bx <- seq(min(x), max(x), length.out = bins + 1)
  by <- seq(min(y), max(y), length.out = bins + 1)
  n <- tabulate((findInterval(x, bx, rightmost.closed = TRUE, all.inside = TRUE) - 1) * bins +
                findInterval(y, by, rightmost.closed = TRUE, all.inside = TRUE), bins^2)
  j <- which(n > 0) - 1
  data.frame(x = (bx[-1] + bx[-(bins + 1)])[j %/% bins + 1] / 2,
             y = (by[-1] + by[-(bins + 1)])[j %% bins + 1] / 2,
             count = n[j + 1])
}

normalize_gene_id <- function(id) sub("_PAR_Y$", "", sub("\\..*$", "", id))

# All (window, gene window) pairs that overlap by >= 1 bp, per chromosome
window_overlaps <- function(qchr, qstart, qend, schr, sstart, send) {
  qi <- split(seq_along(qchr), qchr)
  si <- split(seq_along(schr), schr)
  hits <- lapply(intersect(names(qi), names(si)), function(ch) {
    q <- qi[[ch]][order(qstart[qi[[ch]]])]
    s <- si[[ch]]
    lo <- findInterval(sstart[s] - max(qend[q] - qstart[q]) - 1, qstart[q]) + 1L
    n <- pmax(findInterval(send[s], qstart[q]) - lo + 1L, 0L)
    hq <- q[sequence(n, lo)]
    hs <- rep(s, n)
    k <- qend[hq] >= sstart[hs]
    cbind(hq[k], hs[k])
  })
  h <- do.call(rbind, c(list(matrix(integer(0), 0L, 2L)), hits))
  list(query = h[, 1], subject = h[, 2])
}

# Gene of each window: among the overlapping gene windows, the TSS nearest the window midpoint
# (ties: smaller TSS, then smaller GENCODE id); NA when none overlaps
assign_genes <- function(loci, gene_windows) {
  h <- window_overlaps(loci$chr, loci$start, loci$end,
                       gene_windows$chr, gene_windows$start, gene_windows$end)
  mid <- floor((as.numeric(loci$start) + loci$end) / 2)
  tss <- gene_windows$tss[h$subject]
  dist <- abs(mid[h$query] - tss)
  o <- order(h$query, dist, tss, gene_windows$gencode_id[h$subject], method = "radix")
  o <- o[!duplicated(h$query[o])]
  g <- rep(NA_integer_, nrow(loci))
  g[h$query[o]] <- h$subject[o]
  d <- rep(NA_real_, nrow(loci))
  d[h$query[o]] <- dist[o]
  tes <- ifelse(gene_windows$strand == "+", gene_windows$end, gene_windows$start)
  data.frame(chr = loci$chr, start = loci$start, end = loci$end,
             gene_id = normalize_gene_id(gene_windows$gencode_id[g]),
             gencode_id = gene_windows$gencode_id[g], gene_name = gene_windows$gene_name[g],
             strand = gene_windows$strand[g], tss = gene_windows$tss[g], tes = tes[g], dist = d)
}

# Gene walk: windows ordered by the rank_dpc beta, largest first, each marked with its gates:
# copy number (window on copy_number_chr), gene (has one) and RNA (padj, log2FC, baseMean)
walk_windows <- function(per_window, window_gene, de, rank_dpc, copy_number_chr,
                         padj_max, log2fc_min, basemean_min) {
  bcol <- paste0(rank_dpc, "_beta")
  k <- match(window_gene$gene_id, de$gene_id)
  fc <- de$log2FC[k]
  padj <- de$padj[k]
  bm <- de$baseMean[k]
  rna <- !is.na(fc) & !is.na(padj) & !is.na(bm) &
    padj <= padj_max & fc >= log2fc_min & bm >= basemean_min
  chr <- per_window$chr
  beta <- per_window[[bcol]]
  copy_like <- chr %in% copy_number_chr
  reason <- ifelse(is.na(window_gene$gene_id), "no gene", ifelse(rna, "passes", "RNA gate"))
  o <- order(-beta)
  walk <- data.frame(window_rank = seq_along(o), chr = chr[o],
                     start = per_window$start[o], end = per_window$end[o],
                     beta = beta[o], copy_like = copy_like[o],
                     window_gene[o, c("gene_id", "gene_name", "strand", "tss", "tes")],
                     log2FC = fc[o], padj = padj[o], baseMean = bm[o], reason = reason[o])
  names(walk)[5] <- bcol
  rownames(walk) <- NULL
  walk
}

# Genes for signal tracks: each gene's first passing window in walk order (copy_gate = TRUE also
# skips copy-number-like windows)
track_genes <- function(walk, copy_gate) {
  w <- walk[walk$reason == "passes" & !(copy_gate & walk$copy_like), , drop = FALSE]
  w <- w[!duplicated(w$gene_id), , drop = FALSE]
  cols <- c("gene_id", "gene_name", "strand", "tss", "tes", "window_rank", "chr", "start", "end",
            grep("_beta$", names(walk), value = TRUE), "copy_like", "log2FC", "padj", "baseMean")
  out <- data.frame(gene_rank = seq_len(nrow(w)), w[cols], check.names = FALSE)
  rownames(out) <- NULL
  out
}

# RNA plot points: windows left after dropping, in turn, those on copy_number_chr, without a gene,
# with a gene not in the DE table, or with padj NA; counts of each step
rna_regions <- function(per_window, window_gene, de, marks, dpcs, q_max, copy_number_chr) {
  wg <- window_gene[c("gene_id", "gene_name")]
  g <- match(wg$gene_id, de$gene_id)
  drop <- list(copy_number_chr = per_window$chr %in% copy_number_chr,
               no_gene         = is.na(wg$gene_id),
               gene_not_in_de  = is.na(g),
               padj_na         = is.na(de$padj[g]))
  keep <- rep(TRUE, nrow(per_window))
  n <- integer(length(drop))
  for (j in seq_along(drop)) {
    n[j] <- sum(keep & drop[[j]])
    keep <- keep & !drop[[j]]
  }
  i <- which(keep)
  regions <- data.frame(per_window[i, c("chr", "start", "end")],
                        gene_id = wg$gene_id[i], gene_name = wg$gene_name[i],
                        de[g[i], c("log2FC", "padj", "baseMean")],
                        per_window[i, paste0("d_", marks), drop = FALSE],
                        check.names = FALSE)
  for (k in dpcs) {
    regions[[paste0(k, "_beta")]] <- per_window[[paste0(k, "_beta")]][i]
    regions[[paste0(k, "_q")]] <- per_window[[paste0(k, "_q")]][i]
    regions[[paste0(k, "_significant")]] <- regions[[paste0(k, "_q")]] <= q_max
  }
  rownames(regions) <- NULL
  counts <- data.frame(step = c("windows", names(drop), "regions", "genes"),
                       n = c(nrow(per_window), n, length(i), length(unique(regions$gene_id))))
  list(regions = regions, counts = counts)
}

# Per retained dPC: the smallest |beta| with q <= q_max, and the numbers of windows past it
rna_beta_thresholds <- function(per_window, regions, dpcs, q_max) {
  sig <- lapply(dpcs, function(k) which(per_window[[paste0(k, "_q")]] <= q_max))
  thr <- vapply(seq_along(dpcs), function(j)
    if (length(sig[[j]])) min(abs(per_window[[paste0(dpcs[j], "_beta")]][sig[[j]]])) else NA_real_, 0)
  data.frame(dPC = as.character(dpcs), q_max = rep(q_max, length(dpcs)), beta_threshold = thr,
             n_windows_significant = lengths(sig),
             n_regions_significant = vapply(dpcs, function(k) sum(regions[[paste0(k, "_q")]] <= q_max),
                                            0L, USE.NAMES = FALSE))
}

# Spearman rho of RNA log2FC against each retained dPC's beta (all / significant / not) and each mark's d
rna_spearman <- function(regions, marks, dpcs) {
  y <- regions$log2FC
  one <- function(axis, subset, x, y)
    data.frame(axis = axis, subset = subset, n = length(x),
               rho = if (length(x) >= 3L) stats::cor(x, y, method = "spearman") else NA_real_)
  rows <- c(lapply(dpcs, function(k) {
              x <- regions[[paste0(k, "_beta")]]
              s <- regions[[paste0(k, "_significant")]]
              rbind(one(k, "all", x, y), one(k, "significant", x[s], y[s]),
                    one(k, "not_significant", x[!s], y[!s]))
            }),
            lapply(marks, function(m) one(paste0("d_", m), "all", regions[[paste0("d_", m)]], y)))
  out <- do.call(rbind, c(list(data.frame(axis = character(0), subset = character(0), n = integer(0),
                                          rho = numeric(0))), rows))
  rownames(out) <- NULL
  out
}

# ---- Figures (ggplot2 default theme and colours; diverging scales centred at 0) ----
# 1A: scaled loadings of the retained dPCs, marks in dendrogram order, k_best clusters boxed
plot_loadings_heatmap <- function(scaled, clusters, segments, dpcs) {
  y <- seq_along(clusters$mark) + 0.1 * cumsum(c(0, diff(clusters$cluster) != 0))
  end <- cumsum(rle(clusters$cluster)$lengths)
  blocks <- data.frame(xmin = 0.5, xmax = length(dpcs) + 0.5,
                       ymin = y[c(1, utils::head(end, -1) + 1)] - 0.5, ymax = y[end] + 0.5)
  v <- as.matrix(scaled[match(clusters$mark, scaled$mark), dpcs])
  tiles <- data.frame(x = c(col(v)), y = y[c(row(v))], value = c(v))
  tiles$label <- sub("^-(0\\.00)$", "\\1", sprintf("%.2f", tiles$value))
  at <- function(u) stats::approx(seq_along(y), y, u)$y
  tree <- data.frame(x = 0.4 - segments$y / max(segments$y), xend = 0.4 - segments$yend / max(segments$y),
                     y = at(segments$x), yend = at(segments$xend))
  lim <- max(abs(v))
  ggplot2::ggplot() +
    ggplot2::geom_tile(ggplot2::aes(x, y, fill = value), tiles) +
    ggplot2::geom_text(ggplot2::aes(x, y, label = label), tiles) +
    ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), blocks,
                       fill = NA, colour = "black") +
    ggplot2::geom_segment(ggplot2::aes(x, y, xend = xend, yend = yend), tree) +
    ggplot2::scale_fill_gradient2("scaled loading", limits = c(-lim, lim),
                                  breaks = lim * c(-1, -0.5, 0, 0.5, 1),
                                  labels = function(b) sprintf("%.2f", b)) +
    ggplot2::scale_x_continuous(NULL, breaks = seq_along(dpcs), labels = dpcs, position = "top") +
    ggplot2::scale_y_reverse(NULL, breaks = y, labels = clusters$mark, position = "right")
}

# 1B: % variance per dPC, dashed line at retain_min_var
plot_variance <- function(variance, min_var) {
  variance$dPC <- factor(variance$dPC, variance$dPC)
  ggplot2::ggplot(variance, ggplot2::aes(dPC, percent)) +
    ggplot2::geom_col(ggplot2::aes(fill = retained)) +
    ggplot2::geom_hline(yintercept = 100 * min_var, linetype = "dashed") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", percent)), variance[variance$retained, ],
                       vjust = -0.5) +
    ggplot2::annotate("text", x = nrow(variance) + 0.62, y = 100 * min_var,
                      label = paste0(format(100 * min_var), "%"), hjust = 0, vjust = -0.45) +
    ggplot2::scale_fill_discrete(guide = "none") +
    ggplot2::scale_x_discrete(NULL, expand = ggplot2::expansion(add = c(0.6, 1.5))) +
    ggplot2::scale_y_continuous("% variance", expand = ggplot2::expansion(mult = c(0, 0.09)))
}

# 1C: mean silhouette by k, k_best marked
plot_silhouette <- function(silhouette) {
  ggplot2::ggplot(silhouette, ggplot2::aes(k, mean_silhouette)) +
    ggplot2::geom_vline(xintercept = silhouette$k[silhouette$best], linetype = "dashed") +
    ggplot2::geom_line() +
    ggplot2::geom_point(ggplot2::aes(fill = best), shape = 21, size = 3) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", mean_silhouette)), vjust = -1.2) +
    ggplot2::scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "black"), guide = "none") +
    ggplot2::scale_x_continuous("Number of Clusters", breaks = silhouette$k) +
    ggplot2::scale_y_continuous("Mean Silhouette", expand = ggplot2::expansion(mult = c(0.085, 0.14)))
}

# Lower triangle of the marks' d correlation, marks ordered by dPC1 loading
plot_mark_correlation <- function(mark_correlation, marks, column) {
  i <- match(mark_correlation$mark1, marks)
  j <- match(mark_correlation$mark2, marks)
  d <- data.frame(col = factor(marks[pmin(i, j)], marks[-length(marks)]),
                  row = factor(marks[pmax(i, j)], rev(marks[-1])),
                  value = mark_correlation[[column]])
  ggplot2::ggplot(d, ggplot2::aes(col, row, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = sub("^(-?)0\\.", "\\1.", sub("^-(0\\.00)$", "\\1", sprintf("%.2f", value))))) +
    ggplot2::scale_fill_gradient2(bquote(rho ~ (.(column))), limits = c(-1, 1)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
                   legend.position = "inside", legend.position.inside = c(1, 1),
                   legend.justification = c(1, 1), legend.direction = "horizontal",
                   legend.title.position = "top")
}

# dPC scatter plot drawn from density_grid() (log10 count per cell)
plot_density <- function(grid, bins, dpcs) {
  ggplot2::ggplot(grid, ggplot2::aes(x, y, fill = log10(count))) +
    ggplot2::geom_tile(width = diff(range(grid$x)) / (bins - 1),
                       height = diff(range(grid$y)) / (bins - 1)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::scale_fill_continuous(name = NULL, breaks = range(log10(grid$count)),
                                   labels = c("Low", "High")) +
    ggplot2::labs(x = paste(dpcs[1], "projection"), y = paste(dpcs[2], "projection")) +
    ggplot2::theme(legend.position = "bottom")
}

# RNA log2FC of each window's gene against its beta on one dPC; significant = q <= q_max
plot_rna_projection <- function(regions, dpc, threshold, spearman, q_max) {
  lev <- c("not_significant", "significant")
  d <- data.frame(x = regions[[paste0(dpc, "_beta")]], y = regions$log2FC,
                  tier = factor(lev[1L + regions[[paste0(dpc, "_significant")]]], levels = lev))
  d <- d[order(d$tier), ]
  row <- function(s) spearman[spearman$axis == dpc & spearman$subset == s, ]
  num <- function(n) formatC(n, format = "d", big.mark = ",")
  fq <- format(q_max)
  ft <- function(op) if (is.finite(threshold)) sprintf(", |x| %s %.3f", op, threshold) else ""
  si <- row("significant")
  ns <- row("not_significant")
  lab <- c(sprintf("significant: rho = %+.3f, n = %s", si$rho, num(si$n)),
           sprintf("not significant: rho = %+.3f, n = %s", ns$rho, num(ns$n)))
  ggplot2::ggplot(d, ggplot2::aes(x = x, y = y, colour = tier)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_point(size = 0.3, alpha = 0.3) +
    (if (is.finite(threshold))
       ggplot2::geom_vline(xintercept = c(-threshold, threshold), linetype = "dashed")) +
    ggplot2::annotate("text", x = min(d$x), y = max(d$y), hjust = 0, vjust = 1,
                      label = paste(lab, collapse = "\n")) +
    ggplot2::scale_colour_discrete(
      name = NULL, breaks = lev, drop = FALSE,
      labels = c(sprintf("not significant (q > %s%s)", fq, ft("<")),
                 sprintf("significant (q <= %s%s)", fq, ft(">=")))) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1))) +
    ggplot2::labs(title = sprintf("RNA l2fc vs %s projection", dpc),
                  subtitle = sprintf("All %s assigned regions plotted.", num(nrow(regions))),
                  x = sprintf("%s projection", dpc),
                  y = "RNA log2FC") +
    ggplot2::theme(legend.position = "bottom", legend.direction = "vertical")
}

environment()
})

# dPCA Analysis Function
# Post: Differential principal component analysis (dPCA; Ji et al., PNAS 2013) of per-sample count matrices
#       comparing two groups, then every figure and table drawn from it, written into out_dir.
# Parameter:
#   cm_path    : Character vector of feather count matrices, one per sample, as written by build_count_matrix:
#                'pos' (chr_start_end, 1-based) and one column per mark; same windows in the same order.
#   conditions : Group of each file, same order as cm_path; exactly the two groups of contrast, >= 2 each.
#   contrast   : c(A, B): the direction of dPCA, D = A - B (positive = higher in A).
#   ref_genome : "hg38" or "mm10". Chromosome lengths (from BSgenome) and the Ensembl id prefix of the DE table.
#   out_dir    : Output folder; it must not exist or be empty. Default "./".
#   de_file    : Optional DESeq2 results table (tab-separated; gene_id, log2FoldChange, lfcSE, pvalue, padj,
#                baseMean). NULL (default) skips de, the gene lists and the RNA plots.
#   de_contrast: Required with de_file: c(X, Y), the table's log2FoldChange is X relative to Y. Must name the two
#                groups of contrast, in the same order (sign +1) or reversed (sign -1).
#   gtf_file   : GENCODE GTF (.gtf or .gtf.gz) for the genes of the windows. Required with de_file; given alone it
#                adds window_gene only.
#   marks      : Mark columns to use, in output order. NULL (default) = every column but 'pos' of the first file.
#   copy_number_chr: Chromosomes whose windows differ between the groups by copy number (e.g. "chrY"); such
#                windows are copy-number-like in the gene lists and dropped from the RNA plots. Default none.
#   retain_min_var: Retained dPCs: var_explained >= this. Default 0.10.
#   linkage    : hclust() method for clustering the marks. Default "average".
#   correlation: "spearman" (default) or "pearson", for the mark correlation heatmap.
#   scatter_x, scatter_y: dPCs of the scatter plot. Default "dPC1", "dPC2".
#   density_bins: Bins per axis of the scatter plot. Default 200.
#   q_max      : RNA plots: a window is significant along a dPC if its q <= this. Default 1e-3.
#   rank_dpc   : Gene lists: windows are walked by this dPC's beta, largest first. Default "dPC1".
#   padj_max, log2fc_min, basemean_min: Gene lists, RNA gate: padj <=, log2FC (A - B) >=, baseMean >=.
#                Default 0.05, 1, 50.
#   report_windows: Number of walk windows listed in track_windows. Default 3000.
# Output: Writes per_window.feather and the tables (.tsv) and figures (.pdf and .png; RNA plots .png only) into
#         out_dir. Returns out_dir (character), invisibly. No random numbers: the same inputs give the same outputs.

dpca_analysis <- function(cm_path, conditions, contrast, ref_genome, out_dir = "./",
                          de_file = NULL, de_contrast = NULL, gtf_file = NULL, marks = NULL,
                          copy_number_chr = character(0), retain_min_var = 0.10, linkage = "average",
                          correlation = c("spearman", "pearson"), scatter_x = "dPC1", scatter_y = "dPC2",
                          density_bins = 200, q_max = 1e-3, rank_dpc = "dPC1",
                          padj_max = 0.05, log2fc_min = 1, basemean_min = 50, report_windows = 3000) {
  if (missing(cm_path) || missing(conditions) || missing(contrast) || missing(ref_genome))
    stop("cm_path, conditions, contrast and ref_genome are required.")
  ref_genome  <- match.arg(ref_genome, c("hg38", "mm10"))
  correlation <- match.arg(correlation)
  if (length(contrast) != 2L || contrast[1] == contrast[2])
    stop("contrast must name two different groups, c(A, B) for D = A - B.")

  # DE table: its direction must be stated, and genes need the GTF
  de_sign <- NULL
  if (!is.null(de_file)) {
    if (is.null(de_contrast) || is.null(gtf_file))
      stop("With de_file, de_contrast and gtf_file are required.")
    if (identical(as.character(de_contrast), as.character(contrast))) {
      de_sign <- 1
    } else if (identical(as.character(de_contrast), rev(as.character(contrast)))) {
      de_sign <- -1
    } else {
      stop("de_contrast (", paste(de_contrast, collapse = " vs "), ") must name the two groups of contrast (",
           paste(contrast, collapse = " vs "), "), in either order.")
    }
  }

  if (length(list.files(out_dir, all.files = TRUE, no.. = TRUE)))
    stop("out_dir is not empty: ", out_dir)

  # Read the inputs
  counts <- .dpca_tools$read_counts(cm_path, conditions, marks, contrast)
  de <- if (!is.null(de_file)) .dpca_tools$read_de(de_file, de_sign, c(hg38 = "^ENSG", mm10 = "^ENSMUSG")[[ref_genome]])
  gene_windows <- NULL
  if (!is.null(gtf_file)) {
    genome <- if (ref_genome == "hg38") BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38 else
      BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    chrom_sizes <- GenomeInfoDb::seqlengths(genome)
    gene_windows <- .dpca_tools$read_gene_windows(gtf_file, setNames(as.numeric(chrom_sizes), names(chrom_sizes)))
  }

  config <- list(marks = counts$marks, groups = contrast, retain_min_var = retain_min_var,
                 copy_number_chr = copy_number_chr, linkage = linkage, correlation = correlation,
                 scatter_x = scatter_x, scatter_y = scatter_y, density_bins = density_bins, q_max = q_max,
                 rank_dpc = rank_dpc, padj_max = padj_max, log2fc_min = log2fc_min,
                 basemean_min = basemean_min, report_windows = report_windows)

  # Compute everything, then write (an error during the analysis writes nothing)
  res <- .dpca_tools$dpca_comprehensive(.dpca_engine, counts, de, gene_windows, config)
  .dpca_tools$write_outputs(res, out_dir)
  cat(sprintf("dPCA: %d windows x %d marks, %s vs %s; %d tables and %d figures written to %s\n",
              nrow(counts$loci), length(counts$marks), contrast[1], contrast[2], length(res$tables),
              length(res$plots) + length(res$plots_png), out_dir))
  invisible(out_dir)
}
