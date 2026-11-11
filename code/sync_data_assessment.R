## =============================================================================
## site_year_coverage.R
## Tools for deciding how many sites vs. how many years to keep from wide-format
## monitoring data (Site, Year1, Year2, ... with NA = no survey), for spatial
## synchrony analyses.
##
## Trade-off being assessed:
##   * More years  -> longer series -> longer timescales can be resolved
##                    (wavelet / timescale-specific synchrony)
##   * More sites  -> more site pairs -> better resolution of geographic structure
##
## Wavelet methods (e.g., wsyn) need complete, CONTIGUOUS series, so the main
## tool here (window_tradeoff) searches every contiguous block of years and
## counts how many sites qualify within it. A site qualifies if, within the
## window, it is mostly complete (default >= 90% of years observed; the rest to
## be imputed) and not too zero-inflated (max_zero_prop; you set the cutoff).
## All window functions share one set of rules (.qualify), so the trade-off
## table, the extracted data, and the plots always agree.
##
## All functions are base R. Orientation is site x year throughout.
## =============================================================================


## -----------------------------------------------------------------------------
## 1. prep_wide(): wide data frame -> site x year numeric matrix
## -----------------------------------------------------------------------------
## df        data frame with a site column and one column per year
## site_col  name of the site column
## year_cols optional: which columns are years (default = all other columns)
## Column names are ordered by the number they contain ("2004", "X2004",
## "Year4" all work). Columns are labeled with that number.
prep_wide <- function(df, site_col = "Site", year_cols = NULL) {
  df <- as.data.frame(df)
  if (is.null(year_cols)) year_cols <- setdiff(names(df), site_col)
  yr <- suppressWarnings(as.numeric(gsub("\\D", "", year_cols)))
  if (anyNA(yr)) stop("Could not find a number in these column names: ",
                      paste(year_cols[is.na(yr)], collapse = ", "))
  ord <- order(yr)
  m <- as.matrix(df[, year_cols[ord], drop = FALSE])
  storage.mode(m) <- "numeric"
  dimnames(m) <- list(as.character(df[[site_col]]), as.character(yr[ord]))
  if (any(duplicated(rownames(m)))) warning("Duplicated site names found.")
  if (any(diff(yr[ord]) != 1))
    warning("Year columns are not consecutive (gaps in the column sequence). ",
            "Window searches treat adjacent columns as adjacent years.")
  m
}


## -----------------------------------------------------------------------------
## Helper: longest run of TRUE in a logical vector
## -----------------------------------------------------------------------------
.longest_run <- function(x) {
  r <- rle(x)
  if (!any(r$values)) return(0L)
  max(r$lengths[r$values])
}


## -----------------------------------------------------------------------------
## 2. site_summary(): one row per site
## -----------------------------------------------------------------------------
## n_surveyed       years with a value
## n_missing        years with NA (across the full record)
## first, last      first and last surveyed year
## span             last - first + 1
## internal_gaps    NAs between first and last (gaps inside the record)
## longest_run      longest stretch of consecutive surveyed years
## prop_complete    n_surveyed / span (completeness within the site's own record)
## n_nonzero        surveyed years with a count > 0
## n_zero           surveyed years with a count of 0
## prop_zero        n_zero / n_surveyed (zero inflation over the full record)
## longest_zero_run longest stretch of consecutive surveyed zeros (NAs break a run)
site_summary <- function(m) {
  yrs <- as.numeric(colnames(m))
  obs <- !is.na(m)
  out <- data.frame(
    site        = rownames(m),
    n_surveyed  = rowSums(obs),
    n_missing   = rowSums(!obs),
    first       = apply(obs, 1, function(o) if (any(o)) min(yrs[o]) else NA),
    last        = apply(obs, 1, function(o) if (any(o)) max(yrs[o]) else NA),
    longest_run = apply(obs, 1, .longest_run),
    n_nonzero   = rowSums(m > 0, na.rm = TRUE),
    n_zero      = rowSums(m == 0, na.rm = TRUE),
    longest_zero_run = apply(!is.na(m) & m == 0, 1, .longest_run),
    row.names   = NULL, stringsAsFactors = FALSE
  )
  out$span          <- out$last - out$first + 1
  out$internal_gaps <- out$span - out$n_surveyed
  out$prop_complete <- round(out$n_surveyed / out$span, 3)
  out$prop_zero     <- round(out$n_zero / out$n_surveyed, 3)
  out[order(-out$n_surveyed, -out$longest_run),
      c("site", "n_surveyed", "n_missing", "first", "last", "span",
        "internal_gaps", "prop_complete", "longest_run", "n_nonzero",
        "n_zero", "prop_zero", "longest_zero_run")]
}


## -----------------------------------------------------------------------------
## 3. year_summary(): one row per year
## -----------------------------------------------------------------------------
year_summary <- function(m) {
  data.frame(
    year       = as.numeric(colnames(m)),
    n_sites    = colSums(!is.na(m)),
    n_missing  = colSums(is.na(m)),
    n_nonzero  = colSums(m > 0, na.rm = TRUE),
    n_zero     = colSums(m == 0, na.rm = TRUE),
    row.names  = NULL
  )
}


## -----------------------------------------------------------------------------
## 4. sites_by_years(): how many sites have at least X years?
## -----------------------------------------------------------------------------
## Counts sites with >= X surveyed years (anywhere in the record) and sites
## with a consecutive run of >= X years. The gap between the two columns shows
## how much scattered missingness is costing you.
sites_by_years <- function(m) {
  ss <- site_summary(m)
  X  <- seq_len(ncol(m))
  data.frame(
    min_years             = X,
    sites_any_years       = sapply(X, function(x) sum(ss$n_surveyed  >= x)),
    sites_consecutive_run = sapply(X, function(x) sum(ss$longest_run >= x))
  )
}


## -----------------------------------------------------------------------------
## Helper: which sites qualify in a window of columns idx?
## -----------------------------------------------------------------------------
## Rules (all within the window only):
##   min_complete   proportion of years that must be observed. 0.9 means up to
##                  floor(0.1 * window length) NAs are allowed (to be imputed):
##                  0 NAs for windows < 10 yrs, 1 for 10-19, 2 for 20-29, ...
##   max_zero_prop  largest allowed proportion of zeros among OBSERVED years
##                  (1 = no zero filter)
##   min_nonzero    minimum number of non-zero counts
##   max_gap        longest allowed run of consecutive NAs (Inf = no limit);
##                  e.g. 1 allows only isolated missing years
##   ends_observed  first and last year of the window must be observed, so
##                  imputation is always between two real values
## Returns a data frame of per-site window stats plus a logical 'ok'.
.allowed_NA <- function(L, min_complete) floor((1 - min_complete) * L + 1e-9)

.qualify <- function(m, idx, min_complete = 0.9, max_zero_prop = 1,
                     min_nonzero = 1, max_gap = Inf, ends_observed = TRUE) {
  sub <- m[, idx, drop = FALSE]
  obs <- !is.na(sub)
  n_obs  <- rowSums(obs)
  n_NA   <- ncol(sub) - n_obs
  n_zero <- rowSums(obs & sub == 0)
  n_nz   <- n_obs - n_zero
  pz     <- ifelse(n_obs > 0, n_zero / n_obs, NA)
  gap    <- if (is.finite(max_gap)) apply(!obs, 1, .longest_run) else rep(0, nrow(sub))
  ok <- n_NA <= .allowed_NA(ncol(sub), min_complete) &
    n_nz >= min_nonzero &
    !is.na(pz) & pz <= max_zero_prop + 1e-9 &
    gap <= max_gap
  if (ends_observed) ok <- ok & obs[, 1] & obs[, ncol(sub)]
  data.frame(site = rownames(sub), n_NA = n_NA, n_zero = n_zero,
             prop_zero = round(pz, 3), ok = ok, stringsAsFactors = FALSE)
}


## -----------------------------------------------------------------------------
## 5. window_tradeoff(): the main sites-vs-years table
## -----------------------------------------------------------------------------
## For every contiguous block of years of length >= min_len, counts how many
## sites qualify within that block (rules in .qualify above), then keeps the
## best block for each length.
##
## Output, one row per window length (n_years):
##   n_years      window length
##   allowed_NA   NAs per site allowed at this length (from min_complete)
##   n_sites      most sites that qualify in any window of that length
##   n_pairs      site pairs = n_sites choose 2
##   start, end   the best window (ties: the most recent window is reported;
##                n_ties says how many windows tied)
##   n_NA_total   NAs that would need imputing in that best window
##   max_timescale  rough longest timescale worth interpreting (n_years / tscale_div;
##                  default divisor 3, a common rule of thumb given edge effects)
##   pareto       TRUE if no longer window keeps as many sites
window_tradeoff <- function(m, min_len = 5, min_complete = 0.9,
                            max_zero_prop = 1, min_nonzero = 1, max_gap = Inf,
                            ends_observed = TRUE, tscale_div = 3) {
  nY  <- ncol(m)
  yrs <- as.numeric(colnames(m))
  res <- list()
  for (L in min_len:nY) {
    best <- -1; best_s <- NA; ties <- 0; best_na <- NA
    for (s in 1:(nY - L + 1)) {
      q <- .qualify(m, s:(s + L - 1), min_complete, max_zero_prop,
                    min_nonzero, max_gap, ends_observed)
      k <- sum(q$ok)
      if (k > best) { best <- k; best_s <- s; ties <- 1; best_na <- sum(q$n_NA[q$ok]) }
      else if (k == best) { best_s <- s; ties <- ties + 1; best_na <- sum(q$n_NA[q$ok]) }
    }
    res[[length(res) + 1]] <- data.frame(
      n_years = L, allowed_NA = .allowed_NA(L, min_complete),
      n_sites = best, n_pairs = choose(best, 2),
      start = yrs[best_s], end = yrs[best_s + L - 1], n_ties = ties,
      n_NA_total = best_na, max_timescale = round(L / tscale_div, 1))
  }
  out <- do.call(rbind, res)
  # Pareto: keep a length only if it has more sites than every longer length
  out$pareto <- out$n_sites > c(rev(cummax(rev(out$n_sites)))[-1], -1)
  out
}


## -----------------------------------------------------------------------------
## 6. window_sites(): pull the actual data for a chosen window
## -----------------------------------------------------------------------------
## Returns the site x year sub-matrix for years start:end, keeping sites that
## meet the same rules as window_tradeoff() (NAs left in place for imputation).
## Attributes:
##   site_stats      per-site NAs, zeros, prop_zero in the window, and whether kept
##   dropped_sites   sites that failed
##   n_NA_remaining  NAs left to impute
window_sites <- function(m, start, end, min_complete = 0.9, max_zero_prop = 1,
                         min_nonzero = 1, max_gap = Inf, ends_observed = TRUE) {
  yrs <- as.numeric(colnames(m))
  idx <- which(yrs >= start & yrs <= end)
  q   <- .qualify(m, idx, min_complete, max_zero_prop, min_nonzero, max_gap,
                  ends_observed)
  out <- m[q$ok, idx, drop = FALSE]
  attr(out, "site_stats")     <- q
  attr(out, "dropped_sites")  <- q$site[!q$ok]
  attr(out, "n_NA_remaining") <- sum(is.na(out))
  out
}


## -----------------------------------------------------------------------------
## 7. Plots (base graphics)
## -----------------------------------------------------------------------------
## plot_coverage(): site x year heatmap. Grey = not surveyed, light = surveyed
## with 0, dark = surveyed with count > 0. Sites sorted by number of years.
## Optionally outlines a window (start, end) and marks its qualifying sites
## (red triangles). Pass the same rules as window_tradeoff() through '...'.
plot_coverage <- function(m, start = NULL, end = NULL, main = "Survey coverage",
                          ...) {
  ord <- order(rowSums(!is.na(m)), rowSums(m > 0, na.rm = TRUE))
  m2  <- m[ord, , drop = FALSE]
  code <- ifelse(is.na(m2), 1, ifelse(m2 > 0, 3, 2))
  yrs <- as.numeric(colnames(m2))
  op <- par(mar = c(4, max(4, max(nchar(rownames(m2))) * 0.5), 4.5, 1))
  on.exit(par(op))
  image(x = yrs, y = seq_len(nrow(m2)), z = t(code), zlim = c(1, 3),
        col = c("grey88", "#9ecae1", "#08519c"), axes = FALSE,
        xlab = "Year", ylab = "", main = "")
  title(main, line = 2.8)
  axis(1); axis(2, at = seq_len(nrow(m2)), labels = rownames(m2),
                las = 1, cex.axis = min(1, 30 / nrow(m2)))
  box()
  if (!is.null(start) && !is.null(end)) {
    keep <- rownames(window_sites(m, start, end, ...))
    y_in <- which(rownames(m2) %in% keep)
    rect(start - 0.5, 0.5, end + 0.5, nrow(m2) + 0.5, border = "red", lwd = 2)
    points(rep(start - 0.5, length(y_in)), y_in, pch = 17, col = "red", cex = 0.7)
  }
  legend("topleft", fill = c("grey88", "#9ecae1", "#08519c"), bty = "n",
         legend = c("not surveyed", "zero", "> 0"), cex = 0.8, horiz = TRUE,
         inset = c(0, -0.06), xpd = TRUE)
  invisible(m2)
}

## plot_tradeoff(): number of sites (and pairs, right axis) against window
## length. Pareto-optimal lengths are filled; the top axis shows the rough
## longest interpretable timescale.
plot_tradeoff <- function(tr, main = "Sites vs. years (best contiguous window)") {
  op <- par(mar = c(4.5, 4.5, 4.5, 4.5)); on.exit(par(op))
  plot(tr$n_years, tr$n_sites, type = "b", pch = ifelse(tr$pareto, 19, 1),
       xlab = "Years in window (contiguous)", ylab = "Sites with qualifying data",
       main = "", ylim = c(0, max(tr$n_sites)))
  title(main, line = 3)
  axis(3, at = tr$n_years, labels = tr$max_timescale, cex.axis = 0.7)
  mtext("Rough longest timescale (years)", 3, line = 1.8, cex = 0.8)
  pr <- pretty(c(0, max(tr$n_pairs)))
  sc <- max(tr$n_sites) / max(1, max(tr$n_pairs))
  axis(4, at = pr * sc, labels = pr, col.axis = "grey40")
  mtext("Site pairs", 4, line = 3, col = "grey40")
  legend("topright", pch = c(19, 1), bty = "n", cex = 0.8,
         legend = c("Pareto-optimal", "dominated by a longer window"))
  invisible(tr)
}


## -----------------------------------------------------------------------------
## 8. compare_missing_tolerance(): how much does allowing imputation buy you?
## -----------------------------------------------------------------------------
## Runs window_tradeoff() at several completeness requirements and returns
## n_sites by window length, one column per requirement. min_complete = 1 is
## "no imputation". Other rules pass through '...'.
compare_missing_tolerance <- function(m, completeness = c(1, 0.9, 0.8), ...) {
  tabs <- lapply(completeness, function(k) window_tradeoff(m, min_complete = k, ...))
  out <- data.frame(n_years = tabs[[1]]$n_years)
  for (i in seq_along(completeness))
    out[[paste0("sites_complete", completeness[i] * 100)]] <- tabs[[i]]$n_sites
  out
}


## -----------------------------------------------------------------------------
## 9. Zero inflation
## -----------------------------------------------------------------------------
## zero_summary(): how many sites would survive each zero cutoff, using each
## site's full record (proportion of surveyed years that are zero). A quick
## first look before choosing a cutoff; window-specific results come from
## compare_zero_tolerance().
zero_summary <- function(m, cutoffs = c(0, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1)) {
  ss <- site_summary(m)
  pz <- ss$prop_zero[!is.na(ss$prop_zero)]
  data.frame(max_zero_prop = cutoffs,
             n_sites_kept  = sapply(cutoffs, function(z) sum(pz <= z + 1e-9)),
             prop_sites_kept = round(sapply(cutoffs, function(z) mean(pz <= z + 1e-9)), 3))
}

## compare_zero_tolerance(): the window trade-off at several zero cutoffs.
## Returns n_sites by window length, one column per cutoff (zeros counted
## within each window). Other rules (min_complete, min_len, ...) pass through.
compare_zero_tolerance <- function(m, cutoffs = c(0.25, 0.5, 0.75, 1), ...) {
  tabs <- lapply(cutoffs, function(z) window_tradeoff(m, max_zero_prop = z, ...))
  out <- data.frame(n_years = tabs[[1]]$n_years)
  for (i in seq_along(cutoffs))
    out[[paste0("sites_zero", cutoffs[i] * 100)]] <- tabs[[i]]$n_sites
  out
}

## plot_zeros(): each site's proportion of zeros (full record), sorted, with
## optional cutoff lines; point size shows years surveyed, so short records
## with misleadingly low or high zero proportions stand out.
plot_zeros <- function(m, cutoffs = c(0.25, 0.5), main = "Zero inflation by site") {
  ss <- site_summary(m)
  ss <- ss[!is.na(ss$prop_zero), ]
  ss <- ss[order(ss$prop_zero), ]
  plot(seq_len(nrow(ss)), ss$prop_zero, pch = 19, col = "#08519c",
       cex = 0.5 + 1.5 * ss$n_surveyed / max(ss$n_surveyed),
       ylim = c(0, 1), xaxt = "n", xlab = "", ylab = "Proportion of surveyed years = 0",
       main = main)
  axis(1, at = seq_len(nrow(ss)), labels = ss$site, las = 2,
       cex.axis = min(0.8, 30 / nrow(ss)))
  for (z in cutoffs) {
    abline(h = z, lty = 2, col = "red")
    text(1, z, paste0(sum(ss$prop_zero <= z + 1e-9), " sites <= ", z),
         adj = c(0, -0.4), col = "red", cex = 0.75)
  }
  legend("topleft", pch = 19, col = "#08519c", pt.cex = c(0.8, 2), bty = "n",
         legend = c("few years surveyed", "many years surveyed"), cex = 0.75)
  invisible(ss)
}


## =============================================================================
## Example workflow
## =============================================================================

# dat  <- read.csv("egg_masses_wide.csv")
# m    <- prep_wide(dat, site_col = "Location")
# 
# site_summary(m)            # per-site years, gaps, completeness, zeros
# year_summary(m)            # per-year sites surveyed, zeros
# sites_by_years(m)          # sites with >= X years (any vs consecutive)
# 
# # zero inflation: first look, then choose a cutoff
# zero_summary(m)
# plot_zeros(m, cutoffs = c(0.25, 0.5))
# compare_zero_tolerance(m, cutoffs = c(0.25, 0.5, 0.75, 1), min_len = 10)
# 
# # imputation tolerance (default is 90% complete)
# compare_missing_tolerance(m, completeness = c(1, 0.9), min_len = 10)
# 
# # main trade-off with chosen rules
# tr <- window_tradeoff(m, min_len = 10, min_complete = 0.9, max_zero_prop = 0.5)
# tr[tr$pareto, ]
# plot_tradeoff(tr)
# 
# # chosen option, e.g. 17 years (use the SAME rules throughout):
# best <- tr[tr$n_years == 21, ]
# plot_coverage(m, best$start, best$end, min_complete = 0.9, max_zero_prop = 0.5)
# dat_use <- window_sites(m, best$start, best$end, min_complete = 0.9,
#                         max_zero_prop = 0.5)
# attr(dat_use, "site_stats")       # NAs and zeros per site in the window
# attr(dat_use, "n_NA_remaining")   # values to impute
