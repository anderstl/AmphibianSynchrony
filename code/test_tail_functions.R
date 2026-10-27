## =============================================================================
## test_tail_functions.R -- checks the functions in tail_functions.R against
## cases with known answers. Run from the folder holding tail_functions.R:
##   source("test_tail_functions.R")
## Each test prints PASS or stops with an error describing what failed.
## Uses base R only (plus MASS for the copula check).
## =============================================================================
source("code/tail_functions.R")

check <- function(label, ok) {
  if (!isTRUE(ok)) stop("FAIL: ", label, call. = FALSE)
  cat("PASS:", label, "\n")
}
near <- function(a, b, tol = 1e-8) isTRUE(all.equal(unname(a), unname(b), tolerance = tol))

## -----------------------------------------------------------------------
## Known two-site dataset (20 years): the two sites crash together (lower
## half tracks closely, ranks differ only by two adjacent swaps) while their
## upper halves are scrambled relative to each other. Site B is NOT a copy of
## A: different scale, and 3 different upper-half orderings of ranks 11-20.
## Expected values below were computed independently in Python.
## -----------------------------------------------------------------------
A <- c(27, 41, 58, 60, 65, 71, 50, 66, 21, 31, 15, 37, 30, 78, 31, 101, 73, 22, 48, 51)
B <- c(28, 46, 84, 47, 101, 76, 48, 71, 22, 31, 18, 44, 30, 63, 42, 48, 148, 21, 47, 75)

exp_lower <- 0.49793003150227794
exp_upper <- 0.38125708383356577
exp_rho   <- 0.8791871153358436

## 1. partialSpearman reproduces independently computed values
lo <- partialSpearman(A, B, c(0, 0.5))
up <- partialSpearman(A, B, c(0.5, 1))
check("partialSpearman lower half matches independent value", near(lo, exp_lower))
check("partialSpearman upper half matches independent value", near(up, exp_upper))

## 2. Decomposition: lower + upper = full Spearman rho (holds here because no
##    pair of years has ranks summing to exactly n + 1, which would sit on the
##    x + y = 1 boundary and be excluded from both bands; see note at bottom)
check("lower + upper equals cor(A, B, method = 'spearman')",
      near(lo + up, cor(A, B, method = "spearman")))
check("... and equals the independent rho", near(lo + up, exp_rho))
check("partialSpearman(0, 1) equals full Spearman rho",
      near(partialSpearman(A, B, c(0, 1)), cor(A, B, method = "spearman")))

## 3. tail_dependence returns the right pieces and sign (lower-tail dominant)
td <- tail_dependence(A, B)
check("tail_dependence names", identical(names(td), c("upper", "lower", "tail_dep")))
check("tail_dependence upper/lower match", near(td["upper"], exp_upper) && near(td["lower"], exp_lower))
check("tail_dep = upper - lower", near(td["tail_dep"], exp_upper - exp_lower))
check("tail_dep is NEGATIVE for lower-tail-dependent data", td["tail_dep"] < 0)

## 4. Reflecting both series swaps the tails (upper-tail-dependent version)
tdr <- tail_dependence(-A, -B)
check("reflecting both series swaps lower and upper",
      near(tdr["lower"], td["upper"]) && near(tdr["upper"], td["lower"]))
check("tail_dep is POSITIVE after reflecting both", tdr["tail_dep"] > 0)

## 5. Order of arguments does not matter for a positively correlated pair
check("symmetric in its two arguments",
      near(tail_dependence(A, B), tail_dependence(B, A)))

## 6. Invariant to monotone transforms of either series (rank-based)
check("invariant to monotone transforms",
      near(tail_dependence(A, B), tail_dependence(log(A), B^3)))

## -----------------------------------------------------------------------
## Exact analytic cases
## -----------------------------------------------------------------------
x <- 1:20

## x paired with itself: by symmetry of the ranks about 0.5 each half carries
## exactly half of the (perfect) correlation
self <- tail_dependence(x, x)
check("self-pair: lower = 0.5 and upper = 0.5 exactly",
      near(self["lower"], 0.5) && near(self["upper"], 0.5))
check("self-pair: tail_dep = 0", near(self["tail_dep"], 0))
check("self-pair works for odd n too",
      near(tail_dependence(1:21, 1:21)[c("lower", "upper")], c(0.5, 0.5)))

## perfectly NEGATIVE pair: raw partialSpearman sees x + y = 1 for every year,
## so the strict inequalities exclude all of them -> both halves 0 (this is why
## the sign reflection in align_pair() has to happen first)
check("raw partialSpearman on a perfect negative pair is 0 in both halves",
      near(partialSpearman(x, rev(x), c(0, 0.5)), 0) &&
        near(partialSpearman(x, rev(x), c(0.5, 1)), 0))
al <- align_pair(x, rev(x))
check("align_pair reflects a significantly negative pair", cor(al$x, al$y) > 0.99)
neg <- tail_dependence(x, rev(x))
check("tail_dependence on a perfect negative pair behaves like a perfect positive pair",
      near(neg[c("lower", "upper")], c(0.5, 0.5)))

## -----------------------------------------------------------------------
## align_pair behaviour
## -----------------------------------------------------------------------
check("align_pair leaves a positive pair untouched",
      identical(align_pair(A, B)$y, B))
xna <- c(A, NA); yb <- c(B, 5)
check("NA years are dropped pairwise (result equals the no-NA result)",
      near(tail_dependence(xna, yb), tail_dependence(A, B)))
check("short series return NA", all(is.na(tail_dependence(1:4, c(2, 1, 4, 3)))))

## -----------------------------------------------------------------------
## tail_matrices
## -----------------------------------------------------------------------
set.seed(11)
X <- rbind(S1 = A, S2 = B, S3 = A * 1.5 + rnorm(20, 0, 8), S4 = rev(B))
tm <- tail_matrices(X, tail_dependence)
check("tail_matrices: upper matrix symmetric", isSymmetric(tm$upper))
check("tail_matrices: lower matrix symmetric", isSymmetric(tm$lower))
check("tail_matrices: diag(upper) + diag(lower) = 1",
      near(diag(tm$upper) + diag(tm$lower), rep(1, 4)))
# Ties nudge a series' two halves slightly off 0.5/0.5 (the sum stays exactly 1),
# so allow a small tolerance; the tie-free series S3 should be exactly 0.5.
check("tail_matrices: diagonals are ~0.5 each (within 0.001)",
      all(abs(diag(tm$upper) - 0.5) < 1e-3))
check("tail_matrices: tie-free series has diagonal exactly 0.5",
      near(tm$upper["S3", "S3"], 0.5) && near(tm$lower["S3", "S3"], 0.5))
check("tail_matrices: dimnames kept", identical(rownames(tm$upper), rownames(X)))
check("tail_matrices: entries match direct tail_dependence calls",
      near(tm$upper["S1", "S2"], tail_dependence(A, B)["upper"]) &&
        near(tm$lower["S1", "S3"], tail_dependence(X["S1", ], X["S3", ])["lower"]) &&
        near(tm$lower["S2", "S4"], tail_dependence(X["S2", ], X["S4", ])["lower"]))
check("tail_matrices: taildep = upper - lower", near(tm$taildep, tm$upper - tm$lower))

## -----------------------------------------------------------------------
## tail_stats on hand-built matrices
## -----------------------------------------------------------------------
Mu <- matrix(0, 3, 3); Ml <- Mu
Mu[upper.tri(Mu)] <- c(0.5, 0.4, 0.3)   # entries (1,2), (1,3), (2,3)
Ml[upper.tri(Ml)] <- c(0.3, 0.4, 0.5)
st <- tail_stats(Mu, Ml)
check("tail_stats: sum_taildep cancels when biases offset", near(st["sum_taildep"], 0))
check("tail_stats: struct = 1 - (-1) = 2 for perfectly opposed geographies", near(st["struct"], 2))
st2 <- tail_stats(Mu, Mu)
check("tail_stats: identical matrices -> sum 0, struct 0",
      near(st2["sum_taildep"], 0) && near(st2["struct"], 0))
Ml2 <- Mu; Ml2[upper.tri(Ml2)] <- Mu[upper.tri(Mu)] - 0.1
st3 <- tail_stats(Mu, Ml2)
check("tail_stats: uniform 0.1 upper bias -> sum = 0.3, struct = 0 (same geography)",
      near(st3["sum_taildep"], 0.3) && near(st3["struct"], 0))
check("tail_stats ignores the diagonal",
      near(tail_stats({m <- Mu; diag(m) <- 99; m}, Ml), st))

## -----------------------------------------------------------------------
## Gaussian-copula conversion used by surrogate_tail_test:
## rho_pearson = 2 sin(pi rho_s / 6) should give the target Spearman rho
## -----------------------------------------------------------------------
target <- 0.6
R2 <- matrix(c(1, 2 * sin(pi * target / 6), 2 * sin(pi * target / 6), 1), 2)
set.seed(3)
Zbig <- MASS::mvrnorm(20000, c(0, 0), R2)
check("copula conversion reproduces the target Spearman rho (within 0.01)",
      abs(cor(Zbig[, 1], Zbig[, 2], method = "spearman") - target) < 0.01)

## -----------------------------------------------------------------------
## surrogate_tail_test
## -----------------------------------------------------------------------
Rtrue <- matrix(0.5, 6, 6, dimnames = list(paste0("S", 1:6), NULL)); diag(Rtrue) <- 1
colnames(Rtrue) <- rownames(Rtrue)
s1 <- surrogate_tail_test(20, Rtrue, tail_dependence, 20, seed = 1)
s2 <- surrogate_tail_test(20, Rtrue, tail_dependence, 20, seed = 1)
check("surrogate_tail_test: output is n_surr x 2 with named columns",
      identical(dim(s1), c(20L, 2L)) && identical(colnames(s1), c("sum_taildep", "struct")))
check("surrogate_tail_test: same seed gives identical output", identical(s1, s2))
check("surrogate_tail_test: no NAs", !anyNA(s1))

## Calibration: with symmetric Gaussian-copula 'observed' data (no tail
## dependence by construction), the statistics should fall inside the
## surrogate distribution. Surrogate mean of sum_taildep should sit near 0.
sbig <- surrogate_tail_test(200, Rtrue, tail_dependence, 20, seed = 2)
check("null surrogates: mean sum_taildep is near zero (|mean| < 0.5 over 15 pairs)",
      abs(mean(sbig[, "sum_taildep"])) < 0.5)

set.seed(5)
Zobs <- t(MASS::mvrnorm(20, rep(0, 6), 2 * sin(pi * Rtrue / 6)))
rownames(Zobs) <- rownames(Rtrue)
ob <- tm0 <- tail_matrices(Zobs, tail_dependence)
obs_stats <- tail_stats(ob$upper, ob$lower)
p_global <- min(1, 2 * min(mean(sbig[, "sum_taildep"] <= obs_stats["sum_taildep"]),
                           mean(sbig[, "sum_taildep"] >= obs_stats["sum_taildep"])))
p_struct <- (1 + sum(sbig[, "struct"] >= obs_stats["struct"])) / (1 + 200)
check("null data: global test not significant (p > 0.05)", p_global > 0.05)
check("null data: structure test not significant (p > 0.05)", p_struct > 0.05)

## -----------------------------------------------------------------------
## Checks based on Jon's figure code (kelpsynch-tailassoc repo)
## -----------------------------------------------------------------------

## (a) kelp_theory.R (the theory figure, paper Fig. 2): the same simulation, copied
## without the plotting. Location 1 sits below the sigmoid threshold (stronger
## UPPER-tail association, blue panel d); location 2 sits above it (stronger
## LOWER-tail association, red panel e). 50 years per pair, as in the figure.
sim_theory <- function(seed) {
  burn <- 100; tmax <- burn + 50; sd.reg <- 1; mu.reg1 <- -2; mu.reg2 <- 2
  n.each <- 2; r <- 0.9; K <- 100
  set.seed(seed)
  sigmoid <- function(x) 1 / (1 + exp(-x))
  ricker  <- function(N, r, K) N * exp(r * (1 - (N / K)))
  rescale <- function(x, newsd) { y <- x - mean(x); y / (sd(x) / newsd) }
  Nt.1 <- matrix(NA, n.each, tmax); Nt.2 <- Nt.1
  Et.1 <- rnorm(tmax, mu.reg1, sd.reg); Et.2 <- rnorm(tmax, mu.reg2, sd.reg)
  Nt.1[, 1] <- rnorm(n.each, K, K / 10); Nt.2[, 1] <- rnorm(n.each, K, K / 10)
  for (tt in 2:tmax) {
    Nt.1[, tt] <- ricker(Nt.1[, tt - 1], r, K) + rescale(sigmoid(Et.1), K / 10)[tt] + rnorm(n.each, sd = K / 20)
    Nt.2[, tt] <- ricker(Nt.2[, tt - 1], r, K) + rescale(sigmoid(Et.2), K / 10)[tt] + rnorm(n.each, sd = K / 20)
  }
  keep <- (burn + 1):tmax
  list(blue = Nt.1[, keep], red = Nt.2[, keep])
}
fig <- sim_theory(19)   # the seed Jon used for the figure
td_blue <- tail_dependence(fig$blue[1, ], fig$blue[2, ])
td_red  <- tail_dependence(fig$red[1, ],  fig$red[2, ])
check("theory figure (seed 19): blue pair has upper-tail dependence (tail_dep > 0)", td_blue["tail_dep"] > 0)
check("theory figure (seed 19): red pair has lower-tail dependence (tail_dep < 0)",  td_red["tail_dep"]  < 0)
tds <- t(sapply(1:200, function(s) {
  z <- sim_theory(s)
  c(tail_dependence(z$blue[1, ], z$blue[2, ])["tail_dep"], tail_dependence(z$red[1, ], z$red[2, ])["tail_dep"])
}))
check("theory figure over 200 seeds: blue > 0 and red < 0 in at least 85% of runs",
      mean(tds[, 1] > 0) >= 0.85 && mean(tds[, 2] < 0) >= 0.85)

## (b) pedagogicalFig.R (paper Fig. 1): lower-tail (Clayton), no tail dependence
## (Gaussian) and upper-tail (flipped Clayton), all with Spearman rho = 0.9,
## 500 points. Clayton samples are drawn by the conditional method so no extra
## package is needed.
rclayton <- function(n, th) { u <- runif(n); w <- runif(n); cbind(u, ((w^(-th / (1 + th)) - 1) * u^(-th) + 1)^(-1 / th)) }
set.seed(1)
th <- uniroot(function(t) cor(rclayton(2e5, t), method = "spearman")[1, 2] - 0.9, c(0.5, 20))$root
set.seed(11)
Lw <- rclayton(500, th)
Up <- 1 - Lw
Gs <- MASS::mvrnorm(500, c(0, 0), matrix(c(1, 0.9, 0.9, 1), 2))
tdL <- tail_dependence(Lw[, 1], Lw[, 2]); tdU <- tail_dependence(Up[, 1], Up[, 2]); tdG <- tail_dependence(Gs[, 1], Gs[, 2])
check("Fig 1 style: Clayton (lower-tail) gives tail_dep < 0", tdL["tail_dep"] < -0.03)
check("Fig 1 style: flipped Clayton (upper-tail) gives tail_dep > 0", tdU["tail_dep"] > 0.03)
check("Fig 1 style: Gaussian (no tail dependence) gives tail_dep near 0", abs(tdG["tail_dep"]) < 0.03)
check("Fig 1 style: flipping the Clayton swaps its tails exactly",
      near(tdL[c("upper", "lower")], tdU[c("lower", "upper")]))
check("Fig 1 style: upper + lower is within 0.01 of the full Spearman rho",
      abs(tdL["upper"] + tdL["lower"] - cor(Lw[, 1], Lw[, 2], method = "spearman")) < 0.01)

cat("\nAll checks passed.\n")
cat("Note: lower + upper = full Spearman rho only when no pair of years has\n",
    "ranks summing to exactly n + 1 (x + y = 1 sits on the strict boundary and\n",
    "is dropped from both bands). Real data can have such points, so expect\n",
    "lower + upper to differ from rho by a small amount there.\n")

