## =============================================================================
## tail_functions.R -- function definitions for the tail-dependent synchrony
## analysis (Walter et al. 2022). Sourced by tail_dependent_synchrony_RASY.R and
## by test_tail_functions.R, so the tests exercise exactly the code the analysis
## runs. Only base R + MASS are needed here.
## =============================================================================

#' Partial Spearman correlation over a quantile-bound window (verbatim from
#' partialSpearman.R)
partialSpearman <- function(u, v, bounds) {
  if (!is.numeric(u) | !is.numeric(v)) stop("u and v must be numeric")
  if (length(u) != length(v)) stop("u and v have different length")
  if (min(bounds) < 0 | max(bounds) > 1) stop("bounds must be between 0 and 1")
  
  n <- length(u)
  x <- rank(u) / (n + 1)
  y <- rank(v) / (n + 1)
  
  x_b <- x[x + y > (2 * min(bounds)) & x + y < (2 * max(bounds))]
  y_b <- y[x + y > (2 * min(bounds)) & x + y < (2 * max(bounds))]
  
  sum((x_b - mean(x)) * (y_b - mean(y))) / ((n - 1) * sqrt(var(x) * var(y)))
}

#' Align two time series (common years, drop NA pairs) and reflect the
#' second series if the pair is significantly negatively correlated --
#' the vivj_matrix.R step that has to happen before partialSpearman()/
#' Corbds() can be applied, since both assume a positively-correlated pair.
align_pair <- function(x, y, year_x = NULL, year_y = NULL, level = 0.05) {
  if (!is.null(year_x) && !is.null(year_y)) {
    common <- intersect(year_x, year_y)
    x <- x[match(common, year_x)]
    y <- y[match(common, year_y)]
  }
  keep <- complete.cases(x, y)
  x <- x[keep]; y <- y[keep]
  
  ct <- suppressWarnings(cor.test(x, y, alternative = "two.sided",
                                  method = "spearman", exact = FALSE))
  if (ct$p.value < level && unname(ct$estimate) < 0) y <- -y  # reflect, per vivj_matrix.R
  
  list(x = x, y = y)
}

#' Tail dependence strength = upper-half minus lower-half partial Spearman.
#' Same return signature as before, so nothing downstream (surrogate test,
#' distance decay, GLS regression) needs to change.
tail_dependence <- function(x, y, min_n = 6, ...) {
  al <- align_pair(x, y, ...)
  if (length(al$x) < min_n) return(c(upper = NA_real_, lower = NA_real_, tail_dep = NA_real_))
  lower <- partialSpearman(al$x, al$y, c(0, 0.5))
  upper <- partialSpearman(al$x, al$y, c(0.5, 1))
  c(upper = upper, lower = lower, tail_dep = upper - lower)
}

#' Pairwise upper-tail, lower-tail, and tail-dependence matrices for a
#' SITE x YEAR matrix. tail_dep is the pairwise function (here
#' tail_dependence(), returning c(upper, lower, tail_dep)), passed in so the
#' same code builds the observed matrices and the surrogate ones.
#'
#' j runs from i (not i + 1), so the diagonal is filled by the same call: a
#' pond paired with itself has full Spearman rho = 1, split across the lower
#' and upper halves (Corl + Coru = 1), so those entries land near 0.5 each
#' rather than 1. We loop over i <= j and mirror, rather than over all i, j,
#' because align_pair() reflects the SECOND series when a pair is
#' significantly negatively correlated, and reflecting x instead of y swaps
#' the pair's lower- and upper-tail values. Computing (i, j) and (j, i)
#' separately could therefore give a non-symmetric matrix for negatively
#' correlated pairs.
tail_matrices <- function(x, tail_dep) {
  n   <- nrow(x)
  ids <- rownames(x)
  M_upper <- matrix(NA_real_, n, n, dimnames = list(ids, ids))
  M_lower <- M_upper
  for (i in seq_len(n)) {
    for (j in i:n) {
      td <- tail_dep(x[i, ], x[j, ])
      M_upper[i, j] <- M_upper[j, i] <- td["upper"]
      M_lower[i, j] <- M_lower[j, i] <- td["lower"]
    }
  }
  # positive = upper-tail dominant (synchronous booms)
  # negative = lower-tail dominant (synchronous crashes)
  list(upper = M_upper, lower = M_lower, taildep = M_upper - M_lower)
}

#' The two test statistics from Walter et al. (2022, Methods), computed on the
#' off-diagonal (upper-triangle) entries of the tail matrices:
#'   sum_taildep : sum over site pairs of tail-dependence strength (global
#'                 test of net upper- vs lower-tail bias)
#'   struct      : 1 - cor(M_upper, M_lower); large when the upper- and
#'                 lower-tail synchrony matrices have different geographies
#' Used for both the observed data and every surrogate, so they are computed
#' identically.
tail_stats <- function(M_upper, M_lower) {
  ut <- upper.tri(M_upper)
  c(sum_taildep = sum(M_upper[ut] - M_lower[ut], na.rm = TRUE),
    struct      = 1 - cor(M_upper[ut], M_lower[ut], use = "complete.obs"))
}

#' Nearest valid (positive semi-definite) correlation matrix. Uses
#' Matrix::nearPD(corr = TRUE, maxit = 1000), the same repair used in the
#' kelp paper's ncsurrog.R. If the Matrix package is not installed, falls back
#' to raising negative eigenvalues to 1e-8 and rescaling to a correlation matrix.
nearest_corr <- function(R) {
  if (requireNamespace("Matrix", quietly = TRUE)) {
    return(as.matrix(Matrix::nearPD(R, corr = TRUE, maxit = 1000)$mat))
  }
  e <- eigen(R, symmetric = TRUE)
  e$values <- pmax(e$values, 1e-8)
  cov2cor(e$vectors %*% diag(e$values) %*% t(e$vectors))
}

#' Gaussian-copula surrogate test for tail dependence.
#'
#' @param n_surr   number of surrogate datasets
#' @param corr_mat site x site Spearman correlation matrix of the observed data
#'   (converted internally to the Pearson matrix of the underlying normal)
#' @param tail_dep pairwise function returning c(upper, lower, tail_dep),
#'   e.g. tail_dependence
#' @param n_years  number of years (time points) per surrogate series
#' @param seed     optional RNG seed for reproducibility
#' @return n_surr x 2 matrix of surrogate sum_taildep and struct values
#'
#' RUNTIME: each surrogate computes choose(n, 2) + n pairwise tail_dep calls
#' (~1,400 for 52 ponds). Start with a small n_surr (e.g. 199) to confirm the
#' pipeline runs, then scale up to 999-1999; consider parallelizing over
#' surrogates (e.g. future.apply::future_lapply) if it gets slow.
surrogate_tail_test <- function(n_surr, corr_mat, tail_dep, n_years, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n   <- nrow(corr_mat)
  ids <- rownames(corr_mat)
  
  # Gaussian copula: convert target Spearman rho -> Pearson rho for the
  # underlying normal (rho_pearson = 2*sin(pi*rho_s/6))
  R <- 2 * sin(pi * corr_mat / 6)
  R[is.na(R)] <- 0
  diag(R) <- 1
  # ensure positive semi-definite (nearest valid correlation matrix if needed)
  if (min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) < 1e-8) R <- nearest_corr(R)
  
  out <- matrix(NA_real_, n_surr, 2, dimnames = list(NULL, c("sum_taildep", "struct")))
  for (s in seq_len(n_surr)) {
    # mvrnorm returns years x sites; transpose to site x year to match dat
    Z <- t(MASS::mvrnorm(n_years, mu = rep(0, n), Sigma = R))
    rownames(Z) <- ids
    sm <- tail_matrices(Z, tail_dep)
    out[s, ] <- tail_stats(sm$upper, sm$lower)
  }
  out
}

#' Faithful base-R port of ncsurrog() (Spearman case) from the Ghosh et al.
#' supplement (S9): draw from a Gaussian copula matched to the sample Spearman
#' matrix, then give each site's draw the SORTED OBSERVED values of that site,
#' so marginals and ties match the data exactly, and put NAs where the original
#' has them.
#'
#' @param m         years x sites matrix (same orientation as ncsurrog)
#' @param numsurrog number of surrogate datasets
#' @param seed      optional RNG seed
#' @return years x sites x numsurrog array. attr(, "psd") is FALSE if the
#'   transformed correlation matrix was not positive semidefinite. As in the
#'   kelp repo's ncsurrog.R, a warning is given and the nearest valid
#'   correlation matrix (Matrix::nearPD) is used so sampling can continue.
#'   Unlike the kelp version, NAs in m are allowed and kept in place.
ncsurrog_r <- function(m, numsurrog, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  N <- nrow(m); n <- ncol(m)
  scor <- cor(m, use = "pairwise.complete.obs", method = "spearman")
  ncov <- 2 * sin(pi * scor / 6)   # = copula::iRho(normalCopula(0, 2), scor)
  psd  <- min(eigen(ncov, symmetric = TRUE, only.values = TRUE)$values) >= -1e-10
  if (!psd) {
    warning("ncov is not positive semi-definite; replacing it with the nearest valid correlation matrix")
    ncov <- nearest_corr(ncov)
  }
  out <- array(NA_real_, c(N, n, numsurrog))
  for (b in seq_len(numsurrog)) {
    Z <- MASS::mvrnorm(N, rep(0, n), ncov)
    for (ca in seq_len(n)) {
      good <- is.finite(m[, ca]); mg <- m[good, ca]; sg <- Z[good, ca]
      sg[order(sg)] <- mg[order(mg)]
      Z[good, ca] <- sg; Z[!good, ca] <- NA
    }
    out[, , b] <- Z
  }
  attr(out, "psd") <- psd
  out
}