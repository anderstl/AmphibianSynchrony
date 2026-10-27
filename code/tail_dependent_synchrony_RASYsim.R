## =============================================================================
## Tail-dependent spatial synchrony in a pond-breeding amphibian metapopulation
## Framework follows Walter, Castorani, Bell, Sheppard, Cavanaugh & Reuman (2022)
## "Tail-dependent spatial synchrony arises from nonlinear driver-response
## relationships" Ecology Letters 25:1189-1201, which builds on the "partial
## Spearman correlation" tail-association method of Ghosh, Sheppard, Reid &
## Reuman (2020, Ecol. Evol.) and Ghosh, Sheppard & Reuman (2020, Ecosphere).
##
## This script is written to mirror your 32-pond wood frog (Rana sylvatica)
## dataset's structure (pond x year abundance/egg-mass matrix + coordinates +
## a synchronous climate driver such as SPEI). Because I don't have your raw
## data in this session, Part 1 simulates a dataset with the same shape and
## with a built-in threshold-type driver-response relationship so you can
## confirm the pipeline recovers a known tail-dependence signal before you
## swap in real RASY egg-mass counts and your SPEI series.
##
## UPDATE: Part 2 now uses Jon's partialSpearman() (partialSpearman.R) in
## place of Corbds()/CorlCoru(). The two are numerically equivalent -- same
## pseudo-observations (rank/(n+1)), same anti-diagonal bounds filter on
## vi+vj, same subset-cross-products-over-full-sample-moments formula -- but
## partialSpearman() computes its own ranks internally, so it needs raw (not
## pre-ranked) vectors and drops the VineCopula dependency entirely. The
## vivj_pair()-style alignment/NA-handling/negative-correlation reflection
## from vivj_matrix.R is kept (that logic lives upstream of the correlation
## calculation itself, so it still applies regardless of which of the two
## equivalent functions computes the correlation).
## =============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(nlme)
library(ncf)

# Function definitions live in tail_functions.R (partialSpearman, align_pair,
# tail_dependence, tail_matrices, tail_stats, surrogate_tail_test) so that
# test_tail_functions.R can test exactly the code used here.
source("tail_functions.R")

set.seed(42)

## -----------------------------------------------------------------------
## PART 1. Simulate a 32-pond, 17-year dataset with a threshold-type
## driver-response relationship, analogous to the kelp/wave-intensity case.
## Swap this block out for your real RASY egg-mass matrix + coordinates +
## SPEI series -- everything downstream just needs these objects:
##   coords  : data.frame(pond, lon, lat)
##   dat     : SITE x YEAR numeric matrix of abundance (rows = ponds in the
##             same order as coords$pond, columns = years)
##   drv     : SITE x YEAR numeric matrix of the driver (e.g., SPEI), same
##             dimensions and row order as dat
## (The simulation below builds long-format abund/driver first, then converts
## them to dat/drv at the end of Part 1.)
## -----------------------------------------------------------------------

n_ponds <- 32
n_years <- 17
years   <- seq_len(n_years)

# Fake coordinates spread along a ~150 km transect with some scatter,
# analogous to a set of ponds spanning a gradient in mean hydro-regime.
coords <- data.frame(
  pond = paste0("P", sprintf("%02d", 1:n_ponds)),
  lon  = seq(-90.5, -89.0, length.out = n_ponds) + rnorm(n_ponds, 0, 0.02),
  lat  = seq(38.5, 39.6, length.out = n_ponds) + rnorm(n_ponds, 0, 0.02)
)

# Each pond has a mean "hydroperiod favorability" position relative to a
# threshold (mirrors mean wave calmness in the kelp study). Ponds toward one
# end of the transect run chronically drier (below threshold, on the lower
# shoulder of the response curve); ponds toward the other end are chronically
# wetter (above threshold, upper shoulder). This spatial gradient in mean
# driver level is what generates a *geography* of tail dependence.
pond_mean_driver <- seq(-1.5, 1.5, length.out = n_ponds)[rank(coords$lon)]

# A single spatially synchronous regional driver (e.g., a standardized
# spring SPEI/PDSI-type index) plus pond-level local noise, matching the
# general model structure in Walter et al. (2022): e_s(t) shared, e_l,i(t)
# local. Regional driver has AR(1)-ish persistence.
regional_driver <- as.numeric(arima.sim(list(ar = 0.4), n = n_years, sd = 1))
regional_driver <- scale(regional_driver)[, 1]

driver <- expand.grid(pond = coords$pond, year = years) %>%
  left_join(data.frame(pond = coords$pond, mean_driver = pond_mean_driver), by = "pond") %>%
  mutate(
    driver = mean_driver + regional_driver[year] + rnorm(n(), 0, 0.4)
  ) %>%
  select(pond, year, driver)

# Sigmoid (threshold-like) response of recruitment/growth to the driver --
# the direct analogue of kelp's damage-vs-wave-height relationship. f() is
# flat in both extremes and steep near the threshold (here, driver = 0).
sigmoid_response <- function(x, steepness = 2.5) {
  1 / (1 + exp(-steepness * x))
}

abund <- driver %>%
  group_by(pond) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    growth_effect = sigmoid_response(driver),
    local_noise   = rnorm(n(), 0, 0.25),
    log_abund     = 3 + 2.2 * growth_effect + local_noise
  ) %>%
  # remove linear temporal trend per pond, as in Walter et al. (2022)
  group_by(pond) %>%
  mutate(log_abund = residuals(lm(log_abund ~ year)) + mean(log_abund)) %>%
  ungroup() %>%
  transmute(pond, year, abund = exp(log_abund))

# Convert to SITE x YEAR matrices (rows = ponds, columns = years), the layout
# the rest of the script expects. If your real RASY data is already in this
# shape, skip the conversion and just define `dat` and `drv` directly:
#   dat : numeric matrix, rownames = pond IDs, colnames = years (abundance)
#   drv : numeric matrix, same dimensions/order as dat (e.g., SPEI)
# Row order must match coords$pond.
to_site_by_year <- function(long_df, value_col) {
  long_df %>%
    pivot_wider(names_from = year, values_from = all_of(value_col)) %>%
    arrange(match(pond, coords$pond)) %>%
    tibble::column_to_rownames("pond") %>%
    as.matrix()
}
dat <- to_site_by_year(abund,  "abund")
drv <- to_site_by_year(driver, "driver")
stopifnot(identical(rownames(dat), coords$pond),
          identical(dim(dat), dim(drv)))

## -----------------------------------------------------------------------
## PART 2. Partial Spearman correlation (tail-restricted association)
## (functions defined in tail_functions.R, sourced above)
## partialSpearman() copied verbatim from partialSpearman.R (Jon's version).
## align_pair() below reproduces the alignment/NA-handling/sign-reflection
## that vivj_matrix.R does upstream of the correlation calculation (common
## years, drop NA pairs, flip the second series if the pair is significantly
## negatively correlated, since the bounds-filter derivation assumes a
## positively-correlated pair).
## -----------------------------------------------------------------------


## -----------------------------------------------------------------------
## PART 3. Pairwise tail-dependent synchrony matrices across all ponds
## (tail_matrices() and tail_stats() are defined in tail_functions.R)
## -----------------------------------------------------------------------


pond_ids <- rownames(dat)
n        <- nrow(dat)
n_years  <- ncol(dat)

obs      <- tail_matrices(dat, tail_dependence)
M_upper  <- obs$upper
M_lower  <- obs$lower
M_taildep <- obs$taildep

## -----------------------------------------------------------------------
## PART 4. Surrogate test for overall / spatially structured tail asymmetry
## (surrogate_tail_test() is defined in tail_functions.R)
## Null model: Gaussian-copula surrogates that preserve each pair's overall
## Spearman correlation but have NO tail asymmetry by construction. This
## follows the logic of Walter et al.'s surrogate test (their Methods,
## following Ghosh et al. 2020a Appendix S9), implemented via a multivariate
## normal draw matched to the empirical rank-correlation structure.
## -----------------------------------------------------------------------


# cor() correlates columns, so transpose the site x year matrix to get a
# site x site Spearman matrix (pairwise-complete to tolerate missing years)
emp_spearman <- suppressWarnings(
  cor(t(dat), method = "spearman", use = "pairwise.complete.obs")
)

n_surr    <- 199
obs_stats <- tail_stats(M_upper, M_lower)
surr      <- surrogate_tail_test(n_surr, emp_spearman, tail_dependence,
                                 n_years, seed = 123)

# Global test (two-sided): is net upper- vs lower-tail bias unusual?
p_global <- min(1, 2 * min(mean(surr[, "sum_taildep"] <= obs_stats["sum_taildep"]),
                           mean(surr[, "sum_taildep"] >= obs_stats["sum_taildep"])))
# Structure test (one-sided, large = different geographies of upper vs lower
# tail synchrony); (1 + count) / (1 + n_surr) avoids reporting p = 0
p_struct <- (1 + sum(surr[, "struct"] >= obs_stats["struct"])) / (1 + n_surr)

cat(sprintf(
  "Global tail-dependence test: observed sum = %.3f, surrogate mean = %.3f, p = %.4f\n",
  obs_stats["sum_taildep"], mean(surr[, "sum_taildep"]), p_global
))
cat(sprintf(
  "Spatial-structure test: observed 1 - cor(Mu, Ml) = %.4f, surrogate mean = %.4f, p = %.4f\n",
  obs_stats["struct"], mean(surr[, "struct"]), p_struct
))

## -----------------------------------------------------------------------
## PART 5. Distance decay of upper- vs lower-tail synchrony (cf. Fig 3a-c)
## -----------------------------------------------------------------------

dist_km <- as.matrix(dist(coords[, c("lon", "lat")])) * 111  # rough deg->km

decay_df <- expand.grid(i = 1:n, j = 1:n) %>%
  filter(i < j) %>%
  mutate(
    pond_i   = pond_ids[i],
    pond_j   = pond_ids[j],
    distance = dist_km[cbind(i, j)],
    upper    = M_upper[cbind(i, j)],
    lower    = M_lower[cbind(i, j)]
  )

# Loess-smoothed distance-decay curves with confidence bands, separately for
# each tail -- directly comparable to Fig. 3a-c of Walter et al. (2022). If
# you want the exact nonparametric spatial covariance function + bootstrap
# CI approach they used, swap this for ncf::Sncf() called twice (once on a
# "synchrony" variable built from M_upper, once from M_lower), matching the
# workflow you're already using for RASY correlograms.
plot_decay <- ggplot(decay_df %>%
                       pivot_longer(c(upper, lower), names_to = "tail", values_to = "sync")) +
  geom_smooth(aes(distance, sync, color = tail), method = "loess", se = TRUE) +
  scale_color_manual(values = c(upper = "steelblue", lower = "firebrick")) +
  labs(x = "Distance (km)", y = "Partial Spearman correlation",
       color = "Tail", title = "Distance decay in upper- vs. lower-tail synchrony") +
  theme_minimal()

## -----------------------------------------------------------------------
## PART 6. Tail dependence in the pond-driver relationship, and its spatial
## regression against tail dependence in synchrony (cf. Fig 6a, Table S2)
## -----------------------------------------------------------------------

pond_driver_taildep <- map_dfr(pond_ids, function(p) {
  td <- tail_dependence(dat[p, ], drv[p, ])
  tibble(pond = p, driver_tail_dep = unname(td["tail_dep"]))
})

# Per-pond summary of synchrony tail dependence: the mean tail dependence between
# each pond and ALL other ponds. Walter et al. (2022) averaged only over pairs
# within 25 km because kelp synchrony decays with distance and they had hundreds
# of segments. Here there is no distance window: the ponds are close together and
# synchrony shows no distance decay, so a local window would add nothing. We still
# need one value per pond, since the regression below is pond-level (one row per
# pond, matched to that pond's driver tail dependence).
M_td_nodiag <- M_taildep
diag(M_td_nodiag) <- NA   # drop each pond's pairing with itself
pond_sync_taildep <- tibble(pond = pond_ids,
                            sync_tail_dep = rowMeans(M_td_nodiag, na.rm = TRUE))

reg_df <- pond_sync_taildep %>%
  left_join(pond_driver_taildep, by = "pond") %>%
  left_join(coords, by = "pond") %>%
  filter(is.finite(sync_tail_dep), is.finite(driver_tail_dep)) %>%
  # coordinates in km (approximate), so the spatial range parameter is on a
  # sensible scale; in raw degrees gls() failed to converge on the simulated data
  mutate(x_km = (lon - min(lon)) * 111 * cos(mean(lat) * pi / 180),
         y_km = (lat - min(lat)) * 111)

# GLS with spatially correlated (exponential) errors, as in Walter et al. (2022)
gls_fit <- gls(
  sync_tail_dep ~ driver_tail_dep,
  data = reg_df,
  correlation = corExp(form = ~ x_km + y_km, nugget = TRUE)
)
summary(gls_fit)

plot_reg <- ggplot(reg_df, aes(driver_tail_dep, sync_tail_dep)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "darkgreen") +
  labs(x = "Tail dependence: abundance vs. driver",
       y = "Tail dependence: spatial synchrony (mean over all other ponds)",
       title = "Does driver-response tail dependence predict synchrony tail dependence?") +
  theme_minimal()

## -----------------------------------------------------------------------
## PART 7. Map tail dependence strength across ponds (cf. Fig 4c)
## -----------------------------------------------------------------------

map_df <- coords %>%
  left_join(pond_sync_taildep, by = "pond")

plot_map <- ggplot(map_df, aes(lon, lat, color = sync_tail_dep)) +
  geom_point(size = 3) +
  scale_color_gradient2(low = "firebrick", mid = "white", high = "steelblue",
                        midpoint = 0, name = "Tail dependence\n(upper - lower)") +
  labs(title = "Geography of tail-dependent synchrony across ponds") +
  theme_minimal()

## -----------------------------------------------------------------------
## Inspect results
## -----------------------------------------------------------------------
print(plot_decay)
print(plot_reg)
print(plot_map)