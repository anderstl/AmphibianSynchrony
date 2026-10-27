## =============================================================================
## examples_for_jon_dan.R -- three small examples of the issues described in the
## notes document. Run from the folder holding tail_functions.R:
##   source("examples_for_jon_dan.R")
## Examples 1 and 3 run in seconds. Example 2 takes a few minutes with the
## defaults below; the numbers in the notes came from a larger run (see there).
## =============================================================================
source("code/tail_functions.R")

## -----------------------------------------------------------------------------
## EXAMPLE 1: boundary years. Upper-half + lower-half does not equal Spearman.
## Nine years for two sites. Years 4, 5 and 6 have ranks that add up to exactly
## n + 1 = 10 (x + y = 1 on the pseudo-observation scale), so the strict
## inequalities in partialSpearman() leave them out of BOTH halves.
## -----------------------------------------------------------------------------
cat("\n==== Example 1: boundary years ====\n")
site_a <- c(12, 18, 25, 31, 40, 47, 52, 66, 80)
site_b <- c( 9,  5, 14, 30, 22, 17, 41, 55, 48)
n  <- length(site_a)
rk <- rbind(rank_a = rank(site_a), rank_b = rank(site_b))
rk <- rbind(rk, rank_sum = colSums(rk))
print(rk)
cat("Years on the dividing line (rank sum = n + 1 =", n + 1, "):", which(rk["rank_sum", ] == n + 1), "\n")
lower <- partialSpearman(site_a, site_b, c(0, 0.5))
upper <- partialSpearman(site_a, site_b, c(0.5, 1))
rho   <- cor(site_a, site_b, method = "spearman")
cat(sprintf("lower = %.4f  upper = %.4f  lower + upper = %.4f  Spearman = %.4f\n", lower, upper, lower + upper, rho))
cat("The gap is the contribution of the left-out years. With 17 years, about 43% of\n",
    "pairs of real-looking series had at least one such year in my simulations.\n", sep = "")

## -----------------------------------------------------------------------------
## EXAMPLE 2: false-positive rate of the two surrogate tests when there is NO
## tail dependence. Data are drawn from a Gaussian copula (symmetric tails by
## construction), so a calibrated test should reject at the 5% level about 5% of
## the time. Surrogates come from ncsurrog_r(), our port of ncsurrog().
## The structure statistic is 1 - cor(upper matrix, lower matrix) over pairs;
## the global statistic is the sum of upper - lower over pairs.
## -----------------------------------------------------------------------------
null_calibration <- function(n_sites = 32, n_years = 17, n_data = 100, n_surr = 100,
                             design = c("equal", "distance"), seed = 1, cores = 1) {
  design <- match.arg(design)
  idx <- seq_len(n_sites)
  R <- if (design == "equal") { M <- matrix(0.3, n_sites, n_sites); diag(M) <- 1; M
  } else 0.6 * exp(-abs(outer(idx, idx, "-")) / 8) + diag(0.4, n_sites)
  one <- function(i) {
    set.seed(seed * 100000 + i)
    m   <- MASS::mvrnorm(n_years, rep(0, n_sites), R)          # years x sites
    obs <- tail_matrices(t(m), tail_dependence)
    ut  <- upper.tri(obs$taildep)
    so  <- tail_stats(obs$upper, obs$lower)
    sa  <- suppressWarnings(ncsurrog_r(m, n_surr))
    ss  <- t(sapply(seq_len(n_surr), function(b) {
      S  <- t(sa[, , b]); tm <- tail_matrices(S, tail_dependence)
      c(tail_stats(tm$upper, tm$lower),
        sd_rho = sd(cor(t(S), method = "spearman")[ut]), sd_td = sd(tm$taildep[ut]))
    }))
    c(p_global = min(1, 2 * min(mean(ss[, "sum_taildep"] <= so["sum_taildep"]),
                                mean(ss[, "sum_taildep"] >= so["sum_taildep"]))),
      p_struct = (1 + sum(ss[, "struct"] >= so["struct"])) / (1 + n_surr),
      sd_rho_obs = sd(cor(m, method = "spearman")[ut]), sd_rho_surr = mean(ss[, "sd_rho"]),
      sd_td_obs  = sd(obs$taildep[ut]),                  sd_td_surr  = mean(ss[, "sd_td"]))
  }
  res <- if (cores > 1) parallel::mclapply(seq_len(n_data), one, mc.cores = cores) else lapply(seq_len(n_data), one)
  do.call(rbind, res)
}

cat("\n==== Example 2: false-positive rates with no tail dependence (small run) ====\n")
cat("32 sites x 17 years, 10 null datasets x 50 surrogates each. Expect a few minutes.\n")
cat("Increase n_data and n_surr (and cores) for steadier rates.\n")
res <- null_calibration(n_sites = 32, n_years = 17, n_data = 10, n_surr = 50, design = "equal")
cat(sprintf("Global test rejected at 5%%:    %.0f%% of datasets\n", 100 * mean(res[, "p_global"] <= 0.05)))
cat(sprintf("Structure test rejected at 5%%: %.0f%% of datasets\n", 100 * mean(res[, "p_struct"] <= 0.05)))
cat(sprintf("SD of pairwise Spearman rho:  data %.3f, surrogates %.3f\n", mean(res[, "sd_rho_obs"]), mean(res[, "sd_rho_surr"])))
cat(sprintf("SD of pairwise tail dep:      data %.3f, surrogates %.3f\n", mean(res[, "sd_td_obs"]), mean(res[, "sd_td_surr"])))
cat("The surrogates' pairwise correlations are more spread out than the data's, while\n",
    "the spread of tail dependence matches. That is the likely reason the structure\n",
    "statistic (which compares two matrices that both track overall correlation)\n",
    "rejects too often.\n", sep = "")

## -----------------------------------------------------------------------------
## EXAMPLE 3: more sites than years makes the surrogate correlation matrix invalid.
## Random data, 32 sites x 17 years. After converting the Spearman matrix to the
## underlying-normal matrix (2 sin(pi rho / 6)), it has negative eigenvalues, so
## it is not a valid correlation matrix. ncsurrog warns and uses Matrix::nearPD;
## ncsurrog_r() does the same.
## -----------------------------------------------------------------------------
cat("\n==== Example 3: invalid correlation matrix ====\n")
set.seed(5)
m  <- matrix(rnorm(17 * 32), 17, 32)                 # years x sites
ev <- eigen(2 * sin(pi * cor(m, method = "spearman") / 6), symmetric = TRUE, only.values = TRUE)$values
cat(sprintf("Smallest eigenvalue: %.3f; negative eigenvalues: %d of %d\n", min(ev), sum(ev < 0), length(ev)))
invisible(ncsurrog_r(m, 1))                          # prints the warning, repairs the matrix with nearPD
