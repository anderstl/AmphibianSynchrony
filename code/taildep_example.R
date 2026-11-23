library(ggplot2)

col_lower <- "#D55E00"   # crashes (lower half)
col_upper <- "#0072B2"   # booms   (upper half)

A <- c(27, 41, 58, 60, 65, 71, 50, 66, 21, 31, 15, 37, 30, 78, 31, 101, 73, 22, 48, 51)
B <- c(28, 46, 84, 47, 101, 76, 48, 71, 22, 31, 18, 44, 30, 63, 42, 48, 148, 21, 47, 75)

rng  <- range(c(A, B))
flip <- function(x) sum(rng) - x

type_levels <- c("Synchronized crashes", "Synchronized booms")

toydat <- data.frame(count = c(A, B, flip(A), flip(B)),
                     type  = factor(rep(type_levels, each = 40), levels = type_levels),
                     site  = rep(rep(c("A", "B"), each = 20), 2),
                     years = 1:20)

# Top row: time series
toy1 <- ggplot(toydat, aes(years, count, color = site, shape = site)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("black", "coral4")) +
  scale_shape_manual(values = c(17, 19)) +
  labs(x = "Year", y = "Count", color = "Site", shape = "Site") +
  facet_grid(~type) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

# Bottom row: tail dependence
panel_df <- function(x, y, label) {
  n <- length(x); u <- rank(x) / (n + 1); v <- rank(y) / (n + 1)
  td <- tail_dependence(x, y)
  data.frame(u = u, v = v,
             band = ifelse(u + v < 1, "Low years (lower half)",
                           ifelse(u + v > 1, "High years (upper half)",
                                  "On the dividing line (left out)")),
             panel = sprintf("%s\nlower half %.2f | upper half %.2f",
                             label, td["lower"], td["upper"]))
}

toy <- rbind(panel_df(A, B, "Crash together"),
             panel_df(flip(A), flip(B), "Boom together"))
toy$panel <- factor(toy$panel, levels = unique(toy$panel))  # crash left, boom right

toy2 <- ggplot(toy, aes(u, v, colour = band, shape = band)) +
  geom_abline(intercept = 1, slope = -1, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3.2) +
  scale_colour_manual(values = c("Low years (lower half)"          = col_lower,
                                 "High years (upper half)"         = col_upper,
                                 "On the dividing line (left out)" = "grey40")) +
  scale_shape_manual(values  = c("Low years (lower half)"          = 16,
                                 "High years (upper half)"         = 17,
                                 "On the dividing line (left out)" = 4)) +
  facet_wrap(~panel) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Site A (rank, low to high)", y = "Site B (rank, low to high)",
       colour = NULL, shape = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

toytaildep <- cowplot::plot_grid(toy1, toy2, labels = c("A", "B"), ncol = 1,
                           align = "v", axis = "lr", rel_heights = c(1, 1.15))

ggsave("results/taildep_concept.png", toytaildep, width = 9, height = 9, dpi = 600, bg = "white")
