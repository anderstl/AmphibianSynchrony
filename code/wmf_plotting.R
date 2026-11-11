## Custom plotting for wsyn wmf and wpmf objects
## Based on wsyn::plotmag.tts / plotmag.wpmf
## Requires: wsyn, fields (base version), ggplot2 (ggplot version)

library(wsyn)

# ---- shared helper: pull out and prepare everything both versions need ----
.prep_wmf <- function(object, zlims = NULL, neat = TRUE, sigthresh = 0.95) {
  is_wpmf <- inherits(object, "wpmf")
  
  wav        <- Mod(get_values(object))
  times      <- get_times(object)
  timescales <- get_timescales(object)
  signif     <- if (is_wpmf) get_signif(object) else NA
  
  if (is_wpmf && any(sigthresh >= 1 | sigthresh <= 0)) {
    stop("sigthresh must be strictly between 0 and 1")
  }
  
  rg <- range(wav, na.rm = TRUE)
  if (is.null(zlims)) {
    zlims <- rg
  } else if (rg[1] < zlims[1] || rg[2] > zlims[2]) {
    stop("zlims must encompass the range of the values being plotted")
  }
  
  if (neat) {
    inds       <- which(!is.na(colMeans(wav, na.rm = TRUE)))
    wav        <- wav[, inds, drop = FALSE]
    timescales <- timescales[inds]
    if (!identical(signif, NA) && signif[[1]] %in% c("fft", "aaft")) {
      signif[[3]] <- signif[[3]][, inds, drop = FALSE]
    }
  }
  
  list(wav = wav, times = times, timescales = timescales,
       signif = signif, zlims = zlims, is_wpmf = is_wpmf)
}

.jet <- function() {
  grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
}


# =====================================================================
# 1. BASE R VERSION
# =====================================================================
plotwmf <- function(object,
                    zlims = NULL, neat = TRUE, colorfill = NULL,
                    sigthresh = 0.95, colorbar = TRUE,
                    title = NULL, xlab = "Time", ylab = "Timescale",
                    n_xticks = 8, n_yticks = 8,
                    contour_col = "black", contour_lwd = 2, contour_lty = 1,
                    axis_cex = 1,
                    filename = NA, ...) {
  
  p <- .prep_wmf(object, zlims, neat, sigthresh)
  if (is.null(colorfill)) colorfill <- .jet()
  
  ylocs <- pretty(p$timescales, n = n_yticks)
  xlocs <- pretty(p$times, n = n_xticks)
  logts <- log2(p$timescales)
  
  if (!is.na(filename)) grDevices::pdf(paste0(filename, ".pdf"))
  
  if (colorbar) {
    fields::image.plot(x = p$times, y = logts, z = p$wav, xlab = xlab, ylab = ylab,
                       zlim = p$zlims, axes = FALSE, col = colorfill(100),
                       main = title, ...)
  } else {
    graphics::image(x = p$times, y = logts, z = p$wav, xlab = xlab, ylab = ylab,
                    zlim = p$zlims, axes = FALSE, col = colorfill(100),
                    main = title, ...)
  }
  graphics::axis(1, at = xlocs,        labels = xlocs, cex.axis = axis_cex)
  graphics::axis(2, at = log2(ylocs),  labels = ylocs, cex.axis = axis_cex)
  graphics::box()
  
  # significance contours (wpmf only)
  if (p$is_wpmf && !all(is.na(p$signif))) {
    sg <- p$signif
    if (sg[[1]] == "quick") {
      q <- stats::quantile(sg[[2]], sigthresh)
      graphics::contour(x = p$times, y = logts, z = p$wav, levels = q,
                        drawlabels = FALSE, lwd = contour_lwd, lty = contour_lty,
                        col = contour_col, add = TRUE)
    } else if (sg[[1]] %in% c("fft", "aaft")) {
      graphics::contour(x = p$times, y = logts, z = sg[[3]], levels = sigthresh,
                        drawlabels = FALSE, lwd = contour_lwd, lty = contour_lty,
                        col = contour_col, add = TRUE)
    }
  }
  
  if (!is.na(filename)) grDevices::dev.off()
  invisible(NULL)
}


# =====================================================================
# 2. GGPLOT2 VERSION
# =====================================================================
plotwmf_gg <- function(object,
                       zlims = NULL, neat = TRUE, colorfill = NULL,
                       sigthresh = 0.95, colorbar = TRUE,
                       title = NULL, xlab = "Time", ylab = "Timescale",
                       n_xticks = 8, n_yticks = 8,
                       contour_col = "black", contour_linewidth = 0.8,
                       contour_linetype = "solid",
                       filename = NA, width = 7, height = 5, ...) {
  
  p <- .prep_wmf(object, zlims, neat, sigthresh)
  if (is.null(colorfill)) colorfill <- .jet()
  
  ylocs <- pretty(p$timescales, n = n_yticks)
  xlocs <- pretty(p$times, n = n_xticks)
  logts <- log2(p$timescales)
  
  # matrix is [time, timescale], column-major -> times vary fastest
  grid <- expand.grid(time = p$times, logts = logts)
  grid$value <- as.vector(p$wav)
  
  g <- ggplot2::ggplot(grid, ggplot2::aes(x = time, y = logts, fill = value)) +
    ggplot2::geom_raster(interpolate = FALSE) +   # use geom_tile() if log2 spacing is irregular
    ggplot2::scale_fill_gradientn(colours = colorfill(100), limits = p$zlims,
                                  na.value = "white", name = NULL) +
    ggplot2::scale_x_continuous(breaks = xlocs, labels = xlocs,
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = log2(ylocs), labels = ylocs,
                                expand = c(0, 0)) +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::theme_classic() +
    ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA, colour = "black"))
  
  if (!colorbar) g <- g + ggplot2::guides(fill = "none")
  
  # significance contours (wpmf only)
  if (p$is_wpmf && !all(is.na(p$signif))) {
    sg <- p$signif
    if (sg[[1]] == "quick") {
      cdat   <- grid
      brks   <- stats::quantile(sg[[2]], sigthresh)
    } else if (sg[[1]] %in% c("fft", "aaft")) {
      cdat       <- grid
      cdat$value <- as.vector(sg[[3]])
      brks       <- sigthresh
    }
    g <- g + ggplot2::geom_contour(
      data = cdat, ggplot2::aes(x = time, y = logts, z = value),
      breaks = brks, colour = contour_col, linewidth = contour_linewidth,
      linetype = contour_linetype, inherit.aes = FALSE,na.rm=T)
  }
  
  if (!is.na(filename)) {
    ggplot2::ggsave(paste0(filename, ".pdf"), g, width = width, height = height)
  }
  g
}


## ---- usage ----
# plotwmf(res_wmf)                                  # wmf, base R
# plotwmf(res_wpmf, sigthresh = c(0.95, 0.99))      # wpmf with contours
# g <- plotwmf_gg(res_wpmf, title = "Wood frogs")   # returns a ggplot you can keep adding to
# g + ggplot2::theme(text = ggplot2::element_text(size = 14))