# ---- dPCA engine ----
# Base-R port of the CisGenome dPCA program (Ji, Li, Wang, Ning. PNAS 2013;110(17):6789-6794).
# The functions live in their own environment, .dpca_engine, so that none of them is exported;
# dpca_analysis() calls .dpca_engine$dpca_run(). Only the non-paired, pooled-variance, no-peak-call
# path is used by dpca_analysis(); the other branches are kept from the port as they are.
.dpca_engine <- local({

# Precision helpers: the C program stores signals as float; single_precision = TRUE emulates that,
# FALSE (used by dpca_analysis) keeps double precision
as_float32 <- function(x) {
  r <- readBin(writeBin(as.double(x), raw(), size = 4L), "double",
               n = length(x), size = 4L)
  dim(r) <- dim(x)
  dimnames(r) <- dimnames(x)
  r
}

as_float_identity <- function(x) x

precision_fun <- function(single_precision) {
  if (isTRUE(single_precision)) as_float32 else as_float_identity
}

# Standard normal CDF with |z| capped at 7.5, as CisGenome normcdf() (smallest two-sided p ~ 6.4e-14)
normcdf_dpca <- function(z, clamp = TRUE) {
  if (isTRUE(clamp)) z <- pmax(pmin(z, 7.5), -7.5)
  stats::pnorm(z)
}

# Two-sided p = 2 * (1 - Phi(|t|)), the same arithmetic as the C code
pval_two_sided <- function(t, clamp = TRUE) {
  p <- 2 * (1 - normcdf_dpca(abs(t), clamp = clamp))
  dim(p) <- dim(t)
  dimnames(p) <- dimnames(t)
  p
}

# Replicates per group (rows) and dataset, i.e. mark (columns)
build_rep_num <- function(group_id, dataset_id, n_group = 2L, n_dataset = max(dataset_id)) {
  rep_num <- matrix(0L, n_group, n_dataset)
  for (s in seq_along(group_id)) {
    rep_num[group_id[s], dataset_id[s]] <- rep_num[group_id[s], dataset_id[s]] + 1L
  }
  rep_num
}

# Omega: variance multiplier of each dataset's difference, 1/K1 + 1/K2 (1/K1 when paired)
build_omega <- function(rep_num, paired = FALSE) {
  if (isTRUE(paired)) 1 / rep_num[1, ] else 1 / rep_num[1, ] + 1 / rep_num[2, ]
}


# Variance models: return Var[d] for the covariance correction, the SNR and the test statistic.
# pooled = sigma^2 * Omega, the same for every window (the published method)
var_model_pooled <- function() {
  function(pp) {
    list(var_d = pp$sigma2 * pp$omega, per_locus = FALSE, sigma2 = pp$sigma2,
         label = "pooled sigma^2 * Omega (original)")
  }
}

# custom = caller-supplied variances, per dataset or per window x dataset
var_model_custom <- function(var_d, label = "custom") {
  force(var_d)
  function(pp) {
    per_locus <- is.matrix(var_d) && nrow(var_d) > 1L
    if (per_locus && (nrow(var_d) != pp$n_loci || ncol(var_d) != length(pp$omega))) {
      stop("var_model_custom: var_d must be n_loci x M")
    }
    if (!per_locus && length(var_d) != length(pp$omega)) {
      stop("var_model_custom: var_d must have length M")
    }
    list(var_d = var_d, per_locus = per_locus, sigma2 = pp$sigma2, label = label)
  }
}

# Evaluate a variance model and check its output
apply_var_model <- function(var_model, pp) {
  vm <- var_model(pp)
  if (!is.list(vm) || is.null(vm$var_d)) {
    stop("apply_var_model: a variance model must return a list with 'var_d'")
  }
  if (is.null(vm$per_locus)) vm$per_locus <- is.matrix(vm$var_d)
  if (any(vm$var_d < 0)) stop("apply_var_model: negative variances")
  vm
}

# Preprocessing (port of dPCA_Preprocess_NonPaired / _Paired): log2 transform, optional per-sample
# centring/scaling, replicate means, D = X1 - X2, column centring/scaling, pooled sigma^2.
# x: windows x samples; group_id in {1, 2}; dataset_id in 1..M (one dataset per mark)
dpca_preprocess <- function(x, group_id, dataset_id, rep_num = NULL,
                            transform = FALSE, sample_center = FALSE,
                            sample_scale = FALSE, col_center = TRUE,
                            col_scale = FALSE, paired = FALSE,
                            peak_prob = NULL, use_z = FALSE,
                            peak_fdr_cut = 0.5, single_precision = TRUE) {
  flt <- precision_fun(single_precision)
  x <- flt(as.matrix(x))
  G <- nrow(x); S <- ncol(x)
  M <- max(dataset_id)
  if (is.null(rep_num)) rep_num <- build_rep_num(group_id, dataset_id, 2L, M)
  paired <- isTRUE(paired) || identical(paired, 1L)
  use_z  <- isTRUE(use_z)  || identical(use_z, 1L)
  if (use_z && is.null(peak_prob)) {
    stop("dpca_preprocess: use_z = TRUE requires peak_prob")
  }
  peak_cut <- if (is.null(peak_prob)) NULL else flt(1 - peak_fdr_cut)
  if (!is.null(peak_prob)) peak_prob <- flt(as.matrix(peak_prob))

  # log2(x + 1) per sample, before averaging
  if (isTRUE(transform) || identical(transform, 1L)) x <- flt(log(x + 1) / log(2))

  if (paired) {
    if (S %% 2L != 0L) stop("dpca_preprocess: paired mode needs an even sample count")
    H <- S %/% 2L
    if (any(dataset_id[seq_len(H)] != dataset_id[H + seq_len(H)])) {
      stop("dpca_preprocess: paired samples must share a dataset id")
    }
    dD <- flt(x[, seq_len(H), drop = FALSE] - x[, H + seq_len(H), drop = FALSE])

    if (sample_center || sample_scale) {
      for (h in seq_len(H)) {
        dm <- sum(dD[, h]) / G
        ds <- if (sample_center) {
          if (G > 1) sqrt(sum((dD[, h] - dm)^2) / (G - 1)) else 0
        } else {
          if (G > 0) sqrt(sum(dD[, h]^2) / G) else 0
        }
        if (sample_center) {
          dD[, h] <- flt(dD[, h] - dm)
          x[, h]  <- flt(x[, h]  - dm)  # the C code does not centre x[, h + H]
        }
        if (sample_scale) {
          dD[, h]     <- flt(dD[, h]     / (ds + 1e-6))
          x[, h]      <- flt(x[, h]      / (ds + 1e-6))
          x[, h + H]  <- flt(x[, h + H]  / (ds + 1e-6))
        }
      }
    }

    D <- X1 <- X2 <- matrix(0, G, M)
    for (m in seq_len(M)) {
      idx <- which(dataset_id[seq_len(H)] == m)
      if (!length(idx)) next
      D[,  m] <- rowSums(dD[, idx, drop = FALSE])
      X1[, m] <- rowSums(x[,  idx,     drop = FALSE])
      X2[, m] <- rowSums(x[,  idx + H, drop = FALSE])
      if (rep_num[1L, m] > 0) {
        D[,  m] <- D[,  m] / rep_num[1L, m]
        X1[, m] <- X1[, m] / rep_num[1L, m]
        X2[, m] <- X2[, m] / rep_num[1L, m]
      }
    }
    A <- X1 + X2

    SS <- matrix(0, G, M)
    for (m in seq_len(M)) {
      idx <- which(dataset_id[seq_len(H)] == m)
      if (!length(idx)) next
      SS[, m] <- rowSums((dD[, idx, drop = FALSE] - D[, m])^2)
    }
    df <- sum(pmax(rep_num[1L, ] - 1L, 0L) * (rep_num[1L, ] > 0L))
  } else {
    if (sample_center || sample_scale) {
      for (s in seq_len(S)) {
        dm <- sum(x[, s]) / G
        ds <- if (sample_center) {
          if (G > 1) sqrt(sum((x[, s] - dm)^2) / (G - 1)) else 0
        } else {
          if (G > 0) sqrt(sum(x[, s]^2) / G) else 0
        }
        if (sample_center) x[, s] <- flt(x[, s] - dm)
        if (sample_scale)  x[, s] <- flt(x[, s] / (ds + 1e-6))
      }
    }

    # replicate means per group and dataset; D = X1 - X2
    X1 <- X2 <- matrix(0, G, M)
    for (m in seq_len(M)) {
      i1 <- which(group_id == 1L & dataset_id == m)
      i2 <- which(group_id == 2L & dataset_id == m)
      if (length(i1)) X1[, m] <- rowSums(x[, i1, drop = FALSE])
      if (length(i2)) X2[, m] <- rowSums(x[, i2, drop = FALSE])
      if (rep_num[1L, m] > 0) X1[, m] <- X1[, m] / rep_num[1L, m]
      if (rep_num[2L, m] > 0) X2[, m] <- X2[, m] / rep_num[2L, m]
    }
    D <- X1 - X2
    A <- X1 + X2

    # residual sum of squares around the replicate means, one term per group
    SS <- matrix(0, G, M)
    for (m in seq_len(M)) {
      i1 <- which(group_id == 1L & dataset_id == m)
      i2 <- which(group_id == 2L & dataset_id == m)
      if (length(i1)) SS[, m] <- SS[, m] + rowSums((x[, i1, drop = FALSE] - X1[, m])^2)
      if (length(i2)) SS[, m] <- SS[, m] + rowSums((x[, i2, drop = FALSE] - X2[, m])^2)
    }
    df <- sum(pmax(rep_num[1L, ] - 1L, 0L) * (rep_num[1L, ] > 0L)) +
          sum(pmax(rep_num[2L, ] - 1L, 0L) * (rep_num[2L, ] > 0L))
  }

  # dPCA-Z (only with peak probabilities): zero d where the mark has no peak; applied here and after centring
  if (use_z) D[peak_prob <= peak_cut] <- 0

  col_mean <- colSums(D) / G
  col_sd <- if (G > 1) {
    sqrt(colSums((D - rep(col_mean, each = G))^2) / (G - 1))
  } else {
    colSums((D - rep(col_mean, each = G))^2)
  }

  # centre (and optionally scale) every column of D over the windows; X2 is not shifted, as in the C code
  if (col_center) {
    D  <- D  - rep(col_mean, each = G)
    A  <- A  - rep(col_mean, each = G)
    X1 <- X1 - rep(col_mean, each = G)
  }
  if (col_scale) {
    D  <- D  / rep(col_sd + 1e-6, each = G)
    A  <- A  / rep(col_sd + 1e-6, each = G)
    X1 <- X1 / rep(col_sd + 1e-6, each = G)
    X2 <- X2 / rep(col_sd + 1e-6, each = G)
  }

  if (use_z) D[peak_prob <= peak_cut] <- 0

  # pooled noise variance over all windows, marks and both groups
  # (with col_scale the residuals are divided by sd^2 + 1e-6 but D by sd + 1e-6, as in the C code)
  SP <- if (col_scale) {
    rowSums(SS / rep(col_sd^2 + 1e-6, each = G))
  } else {
    rowSums(SS)
  }

  sigma2 <- if (df > 0) {
    sum(SP) / (df * G)
  } else {
    dpca_norep_var_est(D, col_mean, col_sd, paired = paired,
                       col_center = col_center, col_scale = col_scale)
  }

  list(D = D, X1 = X1, X2 = X2, A = A, sigma2 = sigma2, df = df,
       col_mean = col_mean, col_sd = col_sd, rep_num = rep_num,
       omega = build_omega(rep_num, paired), paired = paired, n_loci = G)
}

# Accessors onto dpca_preprocess()
build_D <- function(x, group_id, dataset_id, ...) {
  dpca_preprocess(x, group_id, dataset_id, ...)$D
}

estimate_sigma2 <- function(x, group_id, dataset_id, ...) {
  pp <- dpca_preprocess(x, group_id, dataset_id, ...)
  list(sigma2 = pp$sigma2, df = pp$df, total_df = pp$df * pp$n_loci)
}

# Fit (port of dPCA_SVD): Ehat = D'D / G minus Var[d] on the diagonal, eigendecomposition,
# projections beta = D V and SNR per dPC. A dataset whose corrected variance is negative is dropped.
dpca_fit <- function(pp, var_model = var_model_pooled(), snr_cutoff = 5,
                     sign_convention = c("max_abs", "none")) {
  sign_convention <- match.arg(sign_convention)
  D <- pp$D
  G <- nrow(D); M <- ncol(D)
  vm <- apply_var_model(var_model, pp)

  var_diag <- if (vm$per_locus) colMeans(vm$var_d) else vm$var_d

  Ehat <- crossprod(D) / G
  corrected <- diag(Ehat) - var_diag
  dropped <- which(corrected < 0)
  diag(Ehat) <- corrected
  if (length(dropped)) {
    Ehat[dropped, ] <- 0
    Ehat[, dropped] <- 0
  }

  # eigenvalues come sorted descending; negative ones are set to 0
  e <- eigen(Ehat, symmetric = TRUE)
  lambda <- pmax(e$values, 0)
  esum <- sum(lambda)
  if (esum < 1e-6) {
    stop("dpca_fit: eigenvalues too small. System is unstable! ",
         "(matches CisGenome dPCA_SVD:26852)")
  }
  percent_var <- lambda / esum

  # sign each dPC so that its largest |loading| is positive
  V <- e$vectors
  if (sign_convention == "max_abs") {
    flip <- apply(V, 2L, function(v) if (v[which.max(abs(v))] < 0) -1 else 1)
    V <- V * rep(flip, each = M)
  }
  dimnames(V) <- list(sprintf("dataset%d", seq_len(M)), sprintf("dPC%d", seq_len(M)))

  beta <- D %*% V
  colnames(beta) <- colnames(V)

  if (vm$per_locus) {
    var_denom <- vm$var_d %*% V^2
    den_mean <- colMeans(var_denom)
  } else {
    den_mean <- colSums(vm$var_d * V^2)
    var_denom <- den_mean
  }
  snr <- colSums(beta^2) / G / (den_mean + 1e-6)
  names(snr) <- colnames(V)

  # the C code reports the last dPC passing the SNR cutoff (not the count), 0 when none passes
  passing <- which(snr >= snr_cutoff)
  rank <- if (length(passing)) max(passing) else 0L

  list(V = V, eigenvalue = lambda, percent_var = percent_var, SNR = snr,
       beta = beta, rank = rank, retained = passing,
       dropped_datasets = dropped, Ehat = Ehat, var_denom = var_denom,
       var_model = vm, snr_cutoff = snr_cutoff, n_loci = G, n_dataset = M)
}

# Storey q-values with lambda fixed at 0.5, as dPCA_FDR (not the qvalue package, which smooths pi0);
# windows ordered by rank_by (dPCA passes -|t|) with a stable sort
storey_qvalue <- function(p, rank_by = NULL, lambda = 0.5) {
  G <- length(p)
  if (is.null(rank_by)) rank_by <- p
  ord <- order(rank_by, method = "radix")
  ps <- p[ord]

  pi0 <- (G - sum(p <= lambda)) / G / (1 - lambda)
  if (pi0 > 1) pi0 <- 1

  qraw <- pi0 * ps * G / seq_len(G)
  q <- pmin(rev(cummin(rev(qraw))), 1)

  out <- numeric(G)
  out[ord] <- q
  list(qvalue = out, pi0 = pi0, order = ord)
}

# Per-window test of every dPC: t = beta / SE against the standard normal, then q-values per dPC
dpca_test <- function(fit, clamp_z = TRUE, lambda = 0.5) {
  beta <- fit$beta
  G <- nrow(beta); M <- ncol(beta)
  vd <- fit$var_denom
  denom <- if (is.matrix(vd)) sqrt(vd) + 1e-10 else rep(sqrt(vd) + 1e-10, each = G)

  tstat <- beta / denom
  pval <- pval_two_sided(tstat, clamp = clamp_z)

  qval <- matrix(0, G, M, dimnames = dimnames(beta))
  sort_id <- matrix(0L, G, M, dimnames = dimnames(beta))
  pi0 <- numeric(M)
  for (j in seq_len(M)) {
    sq <- storey_qvalue(pval[, j], rank_by = -abs(tstat[, j]), lambda = lambda)
    qval[, j] <- sq$qvalue
    pi0[j] <- sq$pi0
    sort_id[, j] <- sq$order
  }
  names(pi0) <- colnames(beta)

  list(t = tstat, p = pval, qvalue = qval, pi0 = pi0, sig_prc = 1 - pi0,
       sort_id = sort_id)
}

# R^B statistic (port of dPCA_RB): share of a window's projected signal from marks with a called peak.
# Needs peak probabilities; not used by dpca_analysis()
rb_statistic <- function(D, V, peak_prob, peak_fdr_cut = 0.5) {
  G <- nrow(D); M <- ncol(D)
  cut <- 1 - peak_fdr_cut
  keep <- peak_prob > cut
  rb <- matrix(0, G, M, dimnames = list(NULL, colnames(V)))
  for (j in seq_len(M)) {
    contrib <- (D * rep(V[, j], each = G))^2
    denom <- rowSums(contrib)
    numer <- rowSums(contrib * peak_prob * keep)
    rb[, j] <- (numer + 1e-6) / (denom + 1e-6)
  }
  rb
}

# dPCA-Z zeroing rule on its own (dpca_preprocess applies it internally)
dpca_z <- function(D, peak_prob, peak_fdr_cut = 0.5, single_precision = TRUE) {
  flt <- precision_fun(single_precision)
  D[flt(peak_prob) <= flt(1 - peak_fdr_cut)] <- 0
  D
}

# Run dPCA end to end (mirrors the dpca command-line tool and its defaults): preprocess, fit, test.
# dpca_analysis() passes x, group_id and dataset_id directly; input files, peak calls and
# output_folder belong to the original program (their readers and writers are not part of this port)
dpca_run <- function(input = NULL, x = NULL, group_id = NULL, dataset_id = NULL,
                     loci = NULL, output_folder = NULL, output_title = "dpca",
                     transform = FALSE, col_center = TRUE, col_scale = FALSE,
                     sample_center = FALSE, sample_scale = FALSE,
                     snr_cutoff = 5, use_z = FALSE, peak_call_file = NULL,
                     peak_prob = NULL, peak_fdr_cut = 0.5, paired = FALSE,
                     var_model = var_model_pooled(), single_precision = TRUE,
                     sign_convention = c("max_abs", "none"), clamp_z = TRUE) {
  sign_convention <- match.arg(sign_convention)

  if (!is.null(input)) {
    dat <- io_read_dpca_input(input, single_precision = single_precision)
    x <- dat$x; group_id <- dat$group_id; dataset_id <- dat$dataset_id
    loci <- dat$loci; paired <- dat$paired == 1L
    rep_num <- dat$rep_num
  } else {
    if (is.null(x) || is.null(group_id) || is.null(dataset_id)) {
      stop("dpca_run: supply either `input` or all of `x`, `group_id`, `dataset_id`")
    }
    x <- as.matrix(x)
    rep_num <- build_rep_num(group_id, dataset_id, 2L, max(dataset_id))
    if (is.null(loci)) {
      loci <- data.frame(name = sprintf("L%d", seq_len(nrow(x))), chr = "chr1",
                         start = (seq_len(nrow(x)) - 1L) * 1000L + 1L,
                         end = seq_len(nrow(x)) * 1000L, strand = "+",
                         stringsAsFactors = FALSE)
    }
  }
  M <- max(dataset_id)

  use_rb <- !is.null(peak_call_file) || !is.null(peak_prob)
  if ((use_z || use_rb) && is.null(peak_prob)) {
    if (is.null(peak_call_file)) {
      stop("dpca_run: use_z = TRUE needs `peak_call_file` or `peak_prob`")
    }
    peak_prob <- io_read_peakprob(peak_call_file, M, loci,
                                  single_precision = single_precision)
  }

  pp <- dpca_preprocess(x, group_id, dataset_id, rep_num = rep_num,
                        transform = transform, sample_center = sample_center,
                        sample_scale = sample_scale, col_center = col_center,
                        col_scale = col_scale, paired = paired,
                        peak_prob = peak_prob, use_z = use_z,
                        peak_fdr_cut = peak_fdr_cut,
                        single_precision = single_precision)

  fit <- dpca_fit(pp, var_model = var_model, snr_cutoff = snr_cutoff,
                  sign_convention = sign_convention)
  tst <- dpca_test(fit, clamp_z = clamp_z)

  rb <- if (use_z || use_rb) {
    rb_statistic(pp$D, fit$V, peak_prob, peak_fdr_cut = peak_fdr_cut)
  } else NULL

  files <- character(0)
  if (!is.null(output_folder)) {
    dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
    p <- function(suffix) file.path(output_folder, paste0(output_title, suffix))
    files <- c(
      write_eigvec(p("_eigvec.txt"), fit$V),
      write_aobs(p("_Aobs.txt"), pp$X1, pp$X2, pp$A, loci),
      write_aproj(p("_Aproj.txt"), pp$X1, pp$X2, pp$A, fit$V, loci),
      write_dobs(p("_Dobs.txt"), pp$D, pp$sigma2, loci),
      write_proj(p("_proj.txt"), fit, tst, pp$sigma2, loci, rb),
      write_rank_files(output_folder, output_title, fit, tst, loci, rb)
    )
  }

  structure(list(preprocess = pp, fit = fit, test = tst, rb = rb, loci = loci,
                 sigma2 = pp$sigma2, df = pp$df, D = pp$D, V = fit$V,
                 beta = fit$beta, t = tst$t, p = tst$p, qvalue = tst$qvalue,
                 pi0 = tst$pi0, SNR = fit$SNR, eigenvalue = fit$eigenvalue,
                 percent_var = fit$percent_var, rank = fit$rank,
                 retained = fit$retained, sort_id = tst$sort_id,
                 dropped_datasets = fit$dropped_datasets, files = files),
            class = "dpca")
}

# Summary printer for a dpca_run() result (inside this environment it is not registered as an S3 method)
print.dpca <- function(x, ...) {
  cat(sprintf("dPCA: %d loci x %d datasets\n", nrow(x$D), ncol(x$D)))
  cat(sprintf("  sigma2 = %.6g  (residual df = %d per locus, %d total)\n",
              x$sigma2, x$df, x$df * nrow(x$D)))
  if (length(x$dropped_datasets)) {
    cat(sprintf("  datasets dropped (negative corrected variance): %s\n",
                paste(x$dropped_datasets, collapse = ", ")))
  }
  tab <- data.frame(eigenvalue = x$eigenvalue, percent_var = x$percent_var,
                    SNR = x$SNR, sig_prc = 1 - x$pi0,
                    row.names = colnames(x$V))
  print(round(tab, 6L))
  cat(sprintf("  rank = %d (dPCs with SNR >= %g: %s)\n", x$rank,
              x$fit$snr_cutoff,
              if (length(x$retained)) paste(x$retained, collapse = ", ") else "none"))
  invisible(x)
}

environment()
})
