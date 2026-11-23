###############################################################################
# Helper functions for geography-of-synchrony analyses
#
# Original functions written by Jon Walter; modified by Tom Anderson.
#
# Source this script before running the geography-of-synchrony analyses:
#   source("synchrony_geography_helpers.R")
#
# Packages used (called with :: where possible, so they only need to be
# installed, not loaded):
#   wsyn     - clust(), cleandat()          [ClustSignif]
#   MASS     - mvrnorm()                    [ClustSignif]
#   ggplot2  - plotting                     [plotmap_gg; must be loaded]
#   ecodist  - full()                       [spline.fit]
#
# Contents:
#   1. ClustSignif  - significance test for cluster (module) structure in a
#                     synchrony network
#   2. spline.fit   - smoothing spline of pairwise synchrony vs. distance,
#                     with bootstrap confidence bands
#   3. mod.in.out   - flags each pair of sites as same-module or
#                     different-module
#   4. plotmap_gg   - ggplot2 map of sites colored by module
###############################################################################

library(ggplot2)


#==============================================================================
# 1. ClustSignif
#==============================================================================
# Tests whether the module structure found by wsyn::clust is stronger than
# what clustering would find in data with no module structure.
#
# How it works:
#   - Runs wsyn::clust on the real data and records the modularity (totQ) of
#     the final partition.
#   - Simulates 'numreps' surrogate datasets as AR(1) time series in which
#       * every pair of sites has the same correlation (the mean pairwise
#         correlation of the real data), so there are no true modules, and
#       * every site has the same lag-1 autocorrelation (the mean across sites
#         in the real data).
#   - Runs the same clust analysis on each surrogate and records its modularity.
#   - The p-value is the fraction of surrogates with modularity >= the real
#     data's modularity.
#
# Note: clust nearly always splits sites into 2+ modules, even for surrogates,
# so the number of clusters alone is not evidence of structure. Modularity
# compared against the surrogates is the test.
#
# Arguments:
#   d        matrix, time (rows) x sites (columns). Each column must have unit
#            variance (e.g., cleaned with wsyn::cleandat, clev >= 3); the
#            function stops if not.
#   numreps  number of surrogate datasets
#   burnin   number of initial time steps discarded from each simulated series
#            so the surrogates are not affected by the starting value of 0
#   method   synchrony measure passed to wsyn::clust, e.g. "ReXWT", "coh"
#            (timescale-specific) or "pearson", "spearman", "kendall"
#            (all timescales combined)
#   tsranges named list of timescale bands, e.g.
#            list(short = c(2, 4), long = c(4, Inf)). Used only for the
#            timescale-specific methods; for correlation methods it is ignored
#            and one analysis labeled "all" is run.
#   clev     cleaning level applied to each surrogate with wsyn::cleandat
#   seed     optional random seed, for reproducible surrogates
#
# Returns a list with:
#   d         the input data
#   surrs     array of surrogate data, time x sites x numreps
#   method    synchrony method used
#   tsranges  timescale bands analyzed
#   datres    for each band, the real-data results: Q (modularity), nclust
#             (number of split levels in the clust object), nmod (number of
#             modules in the final split), membership (final module of each
#             site), clustres (full clust object)
#   surrmod   matrix of surrogate modularities, numreps x bands
#   nclust    matrix of surrogate split-level counts, numreps x bands.
#             This is length(clust$clusters), which counts the unsplit level
#             too: 1 = never split, 2 = one split (2 modules), and so on. It is
#             not always the number of modules; use nmod for that.
#   nmod      matrix of surrogate module counts in the final split,
#             numreps x bands
#   summary   one row per band:
#               Q_obs        modularity of the real data
#               nclust_obs   number of split levels for the real data
#               nmod_obs     number of modules in the real data's final split
#               Q_null_mean  mean surrogate modularity
#               Q_null_q025, Q_null_q975   2.5% and 97.5% surrogate quantiles
#               Q_ses        standardized effect size:
#                            (Q_obs - mean surrogate Q) / sd of surrogate Q
#               pval         fraction of surrogates with Q >= Q_obs
#------------------------------------------------------------------------------
ClustSignif <- function(d, numreps, burnin, method = "ReXWT",
                        tsranges = list(short = c(2, 4), long = c(4, Inf)),
                        clev = 2, seed = NULL)
{
  if (!is.null(seed)) set.seed(seed)
  
  # Decide whether the method is timescale-specific (needs tsrange) or a
  # correlation method (uses all timescales at once)
  timescale_method <- grepl("^(ReXWT|coh|ampl)", method)
  if (timescale_method) {
    # Give unnamed bands names like "2-4"
    if (is.null(names(tsranges)) || any(names(tsranges) == "")) {
      names(tsranges) <- sapply(tsranges, paste, collapse = "-")
    }
  } else {
    tsranges <- list(all = c(0, Inf))   # placeholder; clust ignores tsrange here
  }
  
  # Covariance matrix of the data. Because each column has unit variance, this
  # is also the correlation matrix. Stop if the data were not standardized.
  covd <- cov(d)
  if (max(abs(diag(covd) - 1)) > 1e-5) {
    stop("Error at loc 1: d should have a covariance matrix with 1s on the diagonal")
  }
  
  # Lag-1 autocorrelation of each site, estimated as the slope of a regression
  # of each value on the previous year's value
  numplots <- dim(d)[2]
  rhos <- NA * numeric(numplots)
  for (counter in 1:numplots) {
    h <- d[, counter]
    x <- unname(h[1:(length(h) - 1)])
    y <- unname(h[2:length(h)])
    rhos[counter] <- unname(coef(lm(y ~ x))[2])
  }
  # Surrogates use the average autocorrelation at every site.
  # (To give each site its own autocorrelation instead, use
  #  matrix(rep(rhos, times = numreps), numplots, numreps).)
  rhos <- matrix(mean(rhos), numplots, numreps)
  
  # --- Simulate the surrogate datasets ---
  lendat <- dim(d)[1]
  surrs <- array(NA, c(burnin + lendat, numplots, numreps))
  surrs[1, , ] <- 0
  
  # Noise correlation matrix: every off-diagonal entry is the mean pairwise
  # correlation of the real data, so the surrogates have no module structure
  mncov <- (sum(covd) - sum(diag(covd))) / (numplots^2 - numplots)
  sig <- matrix(mncov, numplots, numplots)
  diag(sig) <- 1
  
  # Draw correlated noise for all time steps and replicates, then arrange it
  # as time x sites x replicates
  noise <- MASS::mvrnorm((burnin + lendat - 1) * numreps,
                         mu = rep(0, numplots), Sigma = sig)
  noise <- array(noise, c(burnin + lendat - 1, numreps, numplots))
  noise <- aperm(noise, c(1, 3, 2))
  
  # AR(1) recursion. The sqrt(1 - rho^2) scaling keeps the variance at 1.
  for (counter in 2:(burnin + lendat)) {
    surrs[counter, , ] <- rhos * surrs[counter - 1, , ] +
      sqrt(1 - rhos^2) * noise[counter - 1, , ]
  }
  # Drop the burn-in so each surrogate has the same length as the real data
  surrs <- surrs[(burnin + 1):(dim(surrs)[1]), , ]
  
  # clust requires coordinates, but they don't affect modularity, so
  # placeholder coordinates are used
  fakecoords <- data.frame(X = 1:numplots, Y = 1:numplots)
  
  # Helper: run clust on a sites x time matrix and extract the modularity of
  # the final split, the number of split levels, the number of modules in the
  # final split, and the final membership
  run_clust <- function(dat, tsr) {
    args <- list(dat = dat, times = 1:lendat, coords = fakecoords, method = method)
    if (timescale_method) args$tsrange <- tsr
    res <- do.call(wsyn::clust, args)
    list(Q = res$modres[[length(res$modres)]]$totQ,
         nclust = length(res$clusters),
         nmod = max(unlist(res$clusters[[length(res$clusters)]])),
         membership = res$clusters[[length(res$clusters)]],
         clustres = res)
  }
  
  # --- Clustering of the real data (clust wants sites x time, hence t()) ---
  datres <- lapply(tsranges, function(tsr) run_clust(t(d), tsr))
  
  # --- Clustering of each surrogate ---
  nband <- length(tsranges)
  surrmod <- matrix(NA, numreps, nband, dimnames = list(NULL, names(tsranges)))
  nclust  <- matrix(NA, numreps, nband, dimnames = list(NULL, names(tsranges)))
  nmod    <- matrix(NA, numreps, nband, dimnames = list(NULL, names(tsranges)))
  for (counter in 1:numreps) {
    surrdat <- t(surrs[, , counter])
    surrdat <- wsyn::cleandat(surrdat, times = 1:lendat, clev = clev)$cdat
    for (b in seq_len(nband)) {
      h <- run_clust(surrdat, tsranges[[b]])
      surrmod[counter, b] <- h$Q
      nclust[counter, b]  <- h$nclust
      nmod[counter, b]    <- h$nmod
    }
  }
  
  # --- Summary table: one row per timescale band ---
  datmod <- sapply(datres, `[[`, "Q")
  summ <- data.frame(
    method      = method,
    timescales  = if (timescale_method)
      sapply(tsranges, paste, collapse = "-") else "all",
    Q_obs       = datmod,
    nclust_obs  = sapply(datres, `[[`, "nclust"),
    nmod_obs    = sapply(datres, `[[`, "nmod"),
    Q_null_mean = colMeans(surrmod),
    Q_null_q025 = apply(surrmod, 2, quantile, 0.025),
    Q_null_q975 = apply(surrmod, 2, quantile, 0.975),
    Q_ses       = (datmod - colMeans(surrmod)) / apply(surrmod, 2, sd),
    pval        = sapply(seq_len(nband),
                         function(b) sum(surrmod[, b] >= datmod[b]) / numreps),
    row.names   = names(tsranges)
  )
  
  return(list(d = d, surrs = surrs, method = method, tsranges = tsranges,
              datres = datres, surrmod = surrmod, nclust = nclust, nmod = nmod,
              summary = summ))
}


#==============================================================================
# 2. spline.fit
#==============================================================================
# Fits a smoothing spline of pairwise synchrony against pairwise distance
# (the solid line in a synchrony-vs-distance plot), and gets confidence bands
# by bootstrapping sites (the dashed lines).
#
# How it works:
#   - Pairs each site-pair distance with its synchrony value from the lower
#     triangle of zmat, and fits smooth.spline with df = sqrt(number of sites).
#   - Bootstrap: resamples sites with replacement, rebuilds the distance and
#     synchrony pairs for the resampled set, refits the spline with the same
#     df, and predicts it at the x values of the empirical spline. Resampling
#     sites (not pairs) keeps the non-independence of pairs that share a site.
#   - Pairs with NA synchrony are dropped before each bootstrap fit. When a
#     site is drawn twice it gets paired with itself, and because the
#     diagonal of zmat is NA, those self-pairs are dropped here.
#
# Arguments:
#   dists      vector of pairwise distances between sites, taken from the
#              lower triangle of a distance matrix, so its order matches
#              zmat[lower.tri(zmat)]. ecodist::full() rebuilds the full
#              distance matrix from it for the bootstrap.
#   zmat       sites x sites synchrony matrix with NA on the diagonal,
#              e.g. the $adj element of a wsyn::clust object
#   nresamp    number of bootstrap resamples (default 1000)
#   quantiles  bootstrap quantiles for the confidence band (default 95%)
#
# Example:
#   distmat <- ncf::gcdist(dat$Longitude, dat$Latitude)   # great-circle km
#   geog.dist <- distmat[lower.tri(distmat)]
#   sp <- spline.fit(geog.dist, clustobj$adj)
#   plot(sp$emp.spline$x, sp$emp.spline$y, type = "l")        # spline
#   lines(sp$emp.spline$x, sp$spline.quantiles[1, ], lty = 2) # lower band
#   lines(sp$emp.spline$x, sp$spline.quantiles[2, ], lty = 2) # upper band
#
# Returns a list with:
#   emp.spline        the smooth.spline fit to the real data
#                     (plot emp.spline$x vs. emp.spline$y for the line)
#   resamp.splines    matrix, nresamp x (number of x values), of bootstrap
#                     spline predictions
#   spline.quantiles  2 x (number of x values) matrix of the lower and upper
#                     bootstrap quantiles at each x (the confidence band)
#------------------------------------------------------------------------------
spline.fit <- function(dists, zmat, nresamp = 1000, quantiles = c(0.025, 0.975)) {
  # Empirical data: distance and synchrony for each pair of sites
  xemp <- c(dists)
  triang <- lower.tri(zmat)
  yemp <- zmat[triang]
  dfs <- sqrt(nrow(zmat))          # spline degrees of freedom
  distmat <- ecodist::full(dists)  # full sites x sites distance matrix
  out <- list()
  
  # Spline fit to the real data
  emp.spline <- smooth.spline(xemp, yemp, df = dfs)
  out$emp.spline <- emp.spline
  
  # Bootstrap: resample sites, refit, predict at the empirical x values
  resamp.splines <- matrix(NA, nrow = nresamp, ncol = length(emp.spline$y))
  for (ii in 1:nresamp) {
    shuffle <- sample(1:nrow(zmat), size = nrow(zmat), replace = TRUE)
    xres <- distmat[shuffle, shuffle][triang]
    yres <- zmat[shuffle, shuffle][triang]
    drop.NaNs <- !is.na(yres)       # drop pairs with missing synchrony
    xres <- xres[drop.NaNs]; yres <- yres[drop.NaNs]
    res.spline <- smooth.spline(xres, yres, df = dfs)
    resamp.splines[ii, ] <- predict(res.spline, x = emp.spline$x)$y
  }
  out$resamp.splines <- resamp.splines
  
  # Confidence band: quantiles of the bootstrap predictions at each x value
  out$spline.quantiles <- apply(resamp.splines, 2, quantile, probs = quantiles)
  return(out)
}


#==============================================================================
# 3. mod.in.out
#==============================================================================
# Labels every pair of sites as belonging to the same module or to different
# modules. Used to color within- vs. between-module pairs in a
# synchrony-vs-distance plot.
#
# Arguments:
#   mod.ids  vector of module memberships, one per site (e.g., the final
#            split from a clust object: clustobj$clusters[[length(clustobj$clusters)]])
#
# Returns a sites x sites matrix:
#   2 = both sites are in the same module (including the diagonal)
#   1 = the sites are in different modules
#
# Example (color within-module pairs blue, between-module pairs red, in the
# same order as the distances passed to spline.fit):
#   inout <- mod.in.out(clustobj$clusters[[length(clustobj$clusters)]])
#   colors <- c("red", "blue")            # 1 = between (red), 2 = within (blue)
#   points(geog.dist, ecodist::lower(clustobj$adj), pch = 20,
#          col = colors[ecodist::lower(inout)])
#   # Mean synchrony of within- and between-module pairs:
#   mean(ecodist::lower(clustobj$adj)[ecodist::lower(inout) == 2], na.rm = TRUE)
#   mean(ecodist::lower(clustobj$adj)[ecodist::lower(inout) == 1], na.rm = TRUE)
#------------------------------------------------------------------------------
mod.in.out <- function(mod.ids) {
  xx <- matrix(1, nrow = length(mod.ids), ncol = length(mod.ids))
  for (i in 1:length(mod.ids)) {
    for (j in 1:length(mod.ids)) {
      if (mod.ids[i] == mod.ids[j]) { xx[i, j] <- 2 }
    }
  }
  return(xx)
}


#==============================================================================
# 4. plotmap_gg
#==============================================================================
# A ggplot2 version of wsyn::plotmap. Maps the sites, colored by module, with
# point size scaled by each site's contribution to its module's modularity
# (nodeQ; larger = more strongly a member of its module). Returns a ggplot
# object, so layers, themes, scales, basemaps, etc. can be added afterward.
#
# Arguments:
#   inclust        a clust object from wsyn::clust
#   spltlvl        which split level to plot (index into inclust$clusters);
#                  default is the final split
#   nodesize       c(min, max) point size in ggplot units; default c(1.5, 4.5),
#                  which roughly matches plotmap's cex = c(1, 3)
#   palette        colors for the modules; default RColorBrewer Set1 (as in
#                  plotmap). Supply a longer vector for more than 9 modules.
#   module_labels  optional legend labels, one per module; default
#                  "module 1", "module 2", ...
#   shape          point shape (default 16, solid circle, as in plotmap)
#   alpha          point transparency (default 1)
#   size_breaks    nodeQ values shown in the size legend; default min, mean, max
#   coord_map      if TRUE and coordinates are lat/lon, adds coord_quickmap()
#                  for an approximately correct aspect ratio
#   filename       file name without extension; if given, saves a .pdf with
#                  ggsave
#   ...            additional arguments passed to ggsave (e.g., width, height)
#
# Coordinates are read from inclust$coords, which must have columns X/Y,
# lat/lon, or latitude/longitude.
#
# Returns a ggplot object; the plotted data (X, Y, module, nodeQ) are in $data.
#------------------------------------------------------------------------------
plotmap_gg <- function(inclust, spltlvl = length(inclust$clusters),
                       nodesize = c(1.5, 4.5),
                       palette = NULL,
                       module_labels = NULL,
                       shape = 16, alpha = 1,
                       size_breaks = NULL,
                       coord_map = FALSE,
                       filename = NA, ...)
{
  # ---- Input checks and defaults ----
  if (!inherits(inclust, "clust")) {
    stop("Error in plotmap_gg: inclust must be a clust object")
  }
  memb <- unlist(inclust$clusters[[spltlvl]])   # module of each site
  nmod <- max(memb)                              # number of modules
  
  if (is.null(palette)) {
    # Set1 from RColorBrewer, same as plotmap
    palette <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
                 "#FFFF33","#A65628","#F781BF","#999999")
  }
  if (nmod > length(palette)) {
    stop(paste0("Error in plotmap_gg: ", nmod, " modules but only ", length(palette),
                " colors in palette; supply a longer palette"))
  }
  if (is.null(module_labels)) {
    module_labels <- paste0("module ", 1:nmod)
  }
  if (length(module_labels) != nmod) {
    stop("Error in plotmap_gg: module_labels must have one entry per module")
  }
  
  # ---- Put coordinates into a common X/Y format ----
  cc <- inclust$coords
  if (all(c("X","Y") %in% names(cc))) {
    X <- cc$X; Y <- cc$Y; x_label <- "X"; y_label <- "Y"; geo <- FALSE
  } else if (all(c("lat","lon") %in% names(cc))) {
    X <- cc$lon; Y <- cc$lat; x_label <- "longitude"; y_label <- "latitude"; geo <- TRUE
  } else if (all(c("latitude","longitude") %in% names(cc))) {
    X <- cc$longitude; Y <- cc$latitude; x_label <- "longitude"; y_label <- "latitude"; geo <- TRUE
  } else {
    stop("Error in plotmap_gg: coords must have X/Y, lat/lon, or latitude/longitude columns")
  }
  
  # ---- Plotting data: one row per site ----
  nodeQ <- inclust$modres[[spltlvl]]$nodeQ
  df <- data.frame(X = X, Y = Y,
                   module = factor(memb, levels = 1:nmod, labels = module_labels),
                   nodeQ = nodeQ)
  
  if (is.null(size_breaks)) {
    size_breaks <- c(min(nodeQ, na.rm = TRUE), mean(nodeQ, na.rm = TRUE),
                     max(nodeQ, na.rm = TRUE))
  }
  
  # ---- Plot ----
  # scale_radius maps nodeQ linearly to point radius, matching plotmap's
  # linear cex scaling
  p <- ggplot(df, aes(x = X, y = Y)) +
    geom_point(aes(color = module, size = nodeQ), shape = shape, alpha = alpha) +
    scale_color_manual(values = setNames(palette[1:nmod], module_labels),
                       drop = FALSE, name = "Membership") +
    scale_radius(range = nodesize, breaks = size_breaks,
                 labels = round(size_breaks, 3), name = "Node weight") +
    guides(color = guide_legend(order = 1, override.aes = list(size = 3)),
           size  = guide_legend(order = 2, override.aes = list(shape = 1))) +
    labs(x = x_label, y = y_label) +
    theme_bw()
  
  if (coord_map && geo) {
    p <- p + coord_quickmap()
  }
  
  # Optionally save to PDF
  if (!is.na(filename)) {
    ggsave(paste0(filename, ".pdf"), plot = p, ...)
  }
  
  p
}