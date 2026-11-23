###Author: Tom Anderson
###Last Updated: November 19, 2026

library(ggplot2)
library(wsyn)
library(plyr)
library(dplyr)
library(tidyr)
library(Matrix)

#Load and combine data
rasy <- read.csv("data/rasy.csv",fileEncoding="latin1")

#Data Cleaning

#replace 2014 with time series median due to shit surveys that year
rasy<-rasy %>%
  rowwise() %>%
  mutate(X2014 = if_else(
    is.na(X2014),
    X2014,
    as.integer(median(c_across(X2006:X2026), na.rm = TRUE))
  )) %>%
  ungroup()

# Assess data completeness ------------------------------------------------

#assess missing values and 0 inflation
source("code/sync_data_assessment.R")
rasy.wide<-rasy%>%
  select(pondname,X2006:X2026)
m <- prep_wide(rasy.wide, site_col = "pondname")

site_summary(m)            # per-site years, gaps, completeness, zeros
year_summary(m)            # per-year sites surveyed, zeros
sites_by_years(m)          # sites with >= X years (any vs consecutive)

# zero inflation: first look, then choose a cutoff
zero_summary(m)
plot_zeros(m, cutoffs = c(0.25, 0.5))
compare_zero_tolerance(m, cutoffs = c(0.25, 0.5, 0.75, 1), min_len = 10)

# imputation tolerance (default is 90% complete)
compare_missing_tolerance(m, completeness = c(1, 0.9), min_len = 10)

# main trade-off with chosen rules
tr <- window_tradeoff(m, min_len = 10, min_complete = 0.9, max_zero_prop = 0.5)
tr[tr$pareto, ]
plot_tradeoff(tr)

# chosen option, e.g. 17 years (use the SAME rules throughout):
best <- tr[tr$n_years == 21, ]
plot_coverage(m, best$start, best$end, min_complete = 0.9, max_zero_prop = 0.5, label_cex=0.3)
dat_use <- window_sites(m, best$start, best$end, min_complete = 0.9,
                        max_zero_prop = 0.5)
attr(dat_use, "site_stats")       # NAs and zeros per site in the window
attr(dat_use, "n_NA_remaining")   # values to impute

rules <- list(min_complete = 0.9, max_zero_prop = 0.5, min_nonzero = 8)
start <- 2006; end <- 2026          # your chosen window

sel <- do.call(window_sites, c(list(m, start, end), rules))

keep_sites <- rownames(sel)
keep_years <- start:end

eggs_use <- rasy|>
  select(pondname,X2006:X2026)|>
  pivot_longer(X2006:X2026,names_to="Year",values_to="Count")|>
  mutate(Year=as.numeric(sub('X', '', Year))) |>
  filter(as.character(pondname) %in% keep_sites,
         Year >= start, Year <= end)

#best option is 21 years for 21 sites (2001-2026)

rasy_comp <- sel %>%
  as.data.frame()%>%
  mutate(row_mean = round(rowMeans(select(., `2006`:`2026`), na.rm = TRUE),0),
         across(`2006`:`2026`, ~ coalesce(., row_mean)))%>%
  mutate(pondname=row.names(sel))%>%
  left_join(rasy%>%select(pondname,Latitude,Longitude),by="pondname")

rasy_mat<-rasy_comp%>%
  select(`2006`:`2026`)%>%
  as.matrix()

#standardize data sets using Box-Cox transformation
rasy_dt<-cleandat(as.matrix(rasy_mat),clev=5,times=keep_years)

#seems like I must need spatiotemporal variable for driver? single weather time series wouldn't work? must be positively correlated?

## -----------------------------------------------------------------------
## PART 3. Pairwise tail-dependent synchrony matrices across all ponds
## -----------------------------------------------------------------------
#source('Fns_taildependence.R')
source('code/tail_functions.R')

pond_ids <- rasy_comp$pondname
n        <- nrow(rasy_mat)
n_years  <- ncol(rasy_mat)

obs      <- tail_matrices(rasy_dt$cdat, tail_dependence)
M_upper  <- obs$upper
M_lower  <- obs$lower
M_taildep <- obs$taildep

#M_upper+M_lower= should equal total spearman correlation but don't
uppervec<-as.vector(M_upper)
lowervec<-as.vector(M_lower)
allvec<-uppervec+lowervec
rasy.spear<-clust(dat = rasy_dt$cdat,
                  times = keep_years,
                  coords = data.frame(lon=rasy_comp$Longitude,
                                      lat=rasy_comp$Latitude),
                  method="spearman",
                  weighted = T)
spearcor<-as.vector(rasy.spear$adj)
spearcor[is.na(spearcor)]<-1
cbind(spearcor,allvec)
sum(spearcor==allvec)/length(spearcor)
#double check all functions that things in correct matrix set up (years vs rows in columns)

dat <- rasy_dt$cdat                       # site x year, detrended
# row names are needed by tail_matrices(); this assumes rows of cdat are in the
# same order as rasy_comp
if (is.null(rownames(dat))) rownames(dat) <- pond_ids
ut     <- upper.tri(M_taildep)            # reuse the observed matrices from above
obs_td <- M_taildep[ut]
obs_stats <- tail_stats(M_upper, M_lower)

## ---- surrogates: one set of 999, each analysed once ----
n_surr <- 999
sa <- ncsurrog_r(t(dat), n_surr, seed = 123)    # years x sites x surrogates
if (!attr(sa, "psd")) message("Note: correlation matrix was not positive semidefinite; nudged (see ncsurrog_r).")

surr      <- matrix(NA_real_, n_surr, 2, dimnames = list(NULL, c("sum_taildep", "struct")))
surr_mats <- vector("list", n_surr)             # kept for poster_checks.R
for (b in seq_len(n_surr)) {
  S <- t(sa[, , b]); rownames(S) <- rownames(dat)
  sm <- tail_matrices(S, tail_dependence)
  surr[b, ]      <- tail_stats(sm$upper, sm$lower)
  surr_mats[[b]] <- sm$taildep
}
rm(sa)                                          # large; no longer needed

## ---- tests ----
p_global <- min(1, 2 * min(mean(surr[, "sum_taildep"] <= obs_stats["sum_taildep"]),
                           mean(surr[, "sum_taildep"] >= obs_stats["sum_taildep"])))
p_struct <- (1 + sum(surr[, "struct"] >= obs_stats["struct"])) / (1 + n_surr)
cat(sprintf("Global tail-dependence test: observed sum = %.3f, surrogate mean = %.3f, p = %.4f\n",
            obs_stats["sum_taildep"], mean(surr[, "sum_taildep"]), p_global))
cat(sprintf("Spatial-structure test (exploratory): observed 1 - cor(Mu, Ml) = %.4f, surrogate mean = %.4f, p = %.4f\n",
            obs_stats["struct"], mean(surr[, "struct"]), p_struct))

## ---- numbers behind the figure ----
surr_td  <- unlist(lapply(surr_mats, function(M) M[ut]))
per_surr <- t(sapply(surr_mats, function(M) { td <- M[ut]; c(mean = mean(td, na.rm = TRUE), sd = sd(td, na.rm = TRUE)) }))
cat(sprintf("Spread: observed SD %.3f is larger than the SD of %d of %d surrogate datasets\n",
            sd(obs_td, na.rm = TRUE), sum(per_surr[, "sd"] < sd(obs_td, na.rm = TRUE)), n_surr))
cat(sprintf("Mean: observed %.3f; %d of %d surrogate datasets have a mean above it, %d below\n",
            mean(obs_td, na.rm = TRUE), sum(per_surr[, "mean"] > mean(obs_td, na.rm = TRUE)),
            n_surr, sum(per_surr[, "mean"] < mean(obs_td, na.rm = TRUE))))

idx <- which(ut, arr.ind = TRUE)
ord <- order(obs_td, decreasing = TRUE)[1:15]
top <- data.frame(pond_i = rownames(dat)[idx[ord, 1]], pond_j = rownames(dat)[idx[ord, 2]],
                  tail_dep = round(obs_td[ord], 2))
print(top)
cat("\nHow often each pond appears among those 15 pairs:\n")
print(sort(table(c(top$pond_i, top$pond_j)), decreasing = TRUE))

## ---- Figure 5 ----
## Density curve pools all 999 surrogates (smooth, free). If you ever want
## fewer, subsample: surr_td <- unlist(lapply(surr_mats[1:200], function(M) M[ut]))
txt <- sprintf("Mean: observed %.3f, surrogates %.3f\nSpread (SD): observed %.3f, surrogates %.3f",
               mean(obs_td, na.rm = TRUE), mean(surr_td, na.rm = TRUE),
               sd(obs_td, na.rm = TRUE), sd(surr_td, na.rm = TRUE))

library(ggplot2)
col_upper <- "#0072B2"

fig5a <- ggplot() +
  geom_histogram(data = data.frame(td = obs_td), aes(td, after_stat(density), fill = "Observed ponds"),
                 binwidth = 0.05, colour = "white", alpha = 0.85) +
  geom_density(data = data.frame(td = surr_td), aes(td, colour = "Surrogates (no tail asymmetry)"),
               linewidth = 1.2, bw = 0.04, key_glyph = "path") +
  geom_vline(xintercept = 0, linetype = "dashed", show.legend = FALSE) +
  scale_fill_manual(values = c("Observed ponds" = "grey60")) +
  scale_colour_manual(values = c("Surrogates (no tail asymmetry)" = col_upper)) +
  labs(x = "Tail dependence \n<- crash together more          boom together more ->",
       y = "Density", fill = NULL, colour = NULL) +
  theme_classic(base_size = 14) + theme(legend.position = "bottom")


# Example data (matching fields::image.plot(x, y, z))
x <- 1:nrow(M_taildep)
y <- 1:ncol(M_taildep)
z <- volcano

# Reshape matrix to a long data frame
library(viridis)
df <- expand.grid(x = x, y = y)
df$z <- as.vector(M_taildep)

png('results/taildep_map.png',res=600,height=5,width=5.5,units="in")
fig5b<-ggplot(df,aes(x=x,y=y,fill=z))+
  geom_raster()+
  scale_fill_viridis(name="Tail\nDependence")+
  theme_classic()+
  labs(y="Site 2", x="Site 1")+
  scale_x_continuous(breaks=0:20, labels = NULL, expand = c(0, 0))+
  scale_y_continuous(breaks=0:20, labels = NULL, expand = c(0, 0))
dev.off()

fig5<-cowplot::plot_grid(fig5a,fig5b,ncol=1,labels=c("A","B"))
ggsave("results/fig5_taildep.png", fig5, width = 5, height = 9, dpi = 600, bg = "white")


#other tests, not shown since significance was not apparent
library(ncf)
test1<-spline.correlog.simmat(coords=cbind(rasy_comp$Longitude,rasy_comp$Latitude),M_upper)
test2<-spline.correlog.simmat(coords=cbind(rasy_comp$Longitude,rasy_comp$Latitude),M_lower)
test3<-spline.correlog.simmat(coords=cbind(rasy_comp$Longitude,rasy_comp$Latitude),emp_spearman)
par(mfrow=c(1,3))
plot(test1)
plot(test2)
plot(test3)
#no patterns present

#png('tail_assoc_distance_with diag.png',res=500,width=5,height=5,units="in")
diag(M_lower)<-0.5;diag(M_upper)<-0.5 #diag was skewing the tail up

ggplot()+
  geom_point(aes(x=distmat17yr,y=rasy.spear$adj,color="Spearman"))+
  geom_point(aes(x=distmat17yr,y=M_lower,color="Lower",))+
  geom_point(aes(x=distmat17yr,y=M_upper,color="Upper"))+
  geom_smooth(aes(x=distmat17yr,y=rasy.spear$adj,color="Spearman"),se=F)+
  geom_smooth(aes(x=distmat17yr,y=M_lower,color="Lower",),se=F)+
  geom_smooth(aes(x=distmat17yr,y=M_upper,color="Upper"),se=F)+
  labs(x="Distance (km)",y="Partial Spearman Correlation")+
  theme_classic()+
  scale_x_continuous(breaks=seq(0,3.5,0.5))+
  scale_color_manual(name="Tail",values=c("Upper"="red","Lower"="blue","Spearman"="green"))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=17),legend.position = c(0.8,0.9))
dev.off()

#Test tail dependence in the pond-driver relationship, and its spatial
## regression against tail dependence in synchrony
rownames(rasy_dt$cdat) <- keep_sites
colnames(rasy_dt$cdat) <- keep_years       # your vector of years, 2010-2026 or the 17 you use

#load climate predictors
clim<-runif(21,5,10)
drv <- matrix(rep(clim, each = nrow(rasy_dt$cdat)), nrow = nrow(rasy_dt$cdat),
              dimnames = list(rownames(rasy_dt$cdat), colnames(rasy_dt$cdat)))

pond_driver_taildep <- map_dfr(pond_ids, function(p) {
  td <- tail_dependence(rasy_dt$cdat[p, ], drv[p, ])
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
coords<-data.frame(pond=rasy_comp$pondname,Longitude=rasy_comp$Longitude,Latitude=rasy_comp$Latitude)
reg_df <- pond_sync_taildep %>%
  left_join(pond_driver_taildep, by = "pond") %>%
  left_join(coords, by = "pond") %>%
  filter(is.finite(sync_tail_dep), is.finite(driver_tail_dep)) %>%
  # coordinates in km (approximate), so the spatial range parameter is on a
  # sensible scale; in raw degrees gls() failed to converge on the simulated data
  mutate(x_km = (Longitude - min(Longitude)) * 111 * cos(mean(Latitude) * pi / 180),
         y_km = (Latitude - min(Latitude)) * 111)

# GLS with spatially correlated (exponential) errors, as in Walter et al. (2022)
library(nlme)
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
plot_reg
