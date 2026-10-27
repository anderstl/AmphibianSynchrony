########################################################
#########Wood Frog Synchrony###########
########################################################

###Author: Tom Anderson
###Last Updated: March 9, 2026
library(ncf)
library(ggplot2)
library(wsyn)
library(plyr)
library(dplyr)
library(tidyr)

####Daniel Boone Conservation Area Analysis####

#Load data
rasy <- read.csv("data/rasy.csv",fileEncoding="latin1")
pond_info <- read.csv("data/pond_attributes.csv",na.strings=c(""))
rasy<-merge(rasy,pond_info,by="pondname")

#replace 2014 with time series median due to shit surveys that year
rasy<-rasy %>%
  rowwise() %>%
  mutate(X2014 = if_else(
    is.na(X2014),
    X2014,
    as.integer(median(c_across(X2006:X2026), na.rm = TRUE))
  )) %>%
  ungroup()

#Convert to long format
rasy_long<-rasy%>%
  select(pondname,X2006:X2026)%>%
  pivot_longer(X2006:X2026,names_to="Year",values_to="Count")
rasy_long$Year<-as.numeric(sub('X', '', rasy_long$Year))

#Plot time series
png("rasy.timeseries.png",height=5,width=5,res=500,units="in")
#pdf('woodfrog_synchrony.pdf')
ggplot(rasy_long,aes(Year,Count,color=pondname))+
  geom_point(show.legend = F)+
  geom_line(aes(color=pondname),show.legend = F)+
  labs(y="Number of Egg Masses")+
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=17))+
  scale_x_continuous(breaks=seq(2006,2026,1),labels=c("2006",rep("",3),"2010",rep("",3),"2014",rep("",3),"2018",rep("",3),"2022",rep("",3),"2026"))
dev.off()

#Scale data- includes all sites and keeps NAs
rasy_mat <- rasy %>%
  select(X2006:X2026)%>%
  rowwise() %>%
  mutate(
    # Store row mean and standard deviation
    row_mean = mean(c_across(everything()),na.rm=T),
    row_sd = sd(c_across(everything()),na.rm=T),
    # Scale all original columns
    across(X2006:X2026, ~ (.x - row_mean) / row_sd)
  ) %>%
  ungroup()%>%
  select(X2006:X2026)

# Standard Synchrony Analysis ---------------------------------------------

#perform cross correlation analysis
rasy.syn<-Sncf(x=rasy$Longitude,
               y=rasy$Latitude,na.rm=T,
               z=rasy_mat,resamp = 1000,latlon = T,quiet=T)

#find degree of synchrony across study area
summary(rasy.syn)$Regional.sync #0.12

#png("rasy_sncf.png",res=500,height=5,width=5,units="in")
par(mar=c(4,5,2,1))
plot(rasy.syn,xlab="Distance (km)",cex.lab=1.5,cex.axis=1.25,las=1,ylim=c(-1,1))
title(main="DBCA Wood Frog Egg Mass Counts, 2006-2026",cex=0.5)
#dev.off()

#manually plot geographic distance vs correlation
distmat<-gcdist(rasy$Longitude,rasy$Latitude)
cormat <- cor(t(as.matrix(rasy_mat)), method="pearson")
diag(cormat) <- NA

# plot(as.vector(distmat),as.vector(cormat),yaxt="n",
#      xlab="Distance (km)",ylab="Correlation",ylim=c(-0.5,1.0))
# abline(h=0)
# axis(2,at = c(-.5,-.25,0,.25,.5,.75,1),las=1,cex.lab=1,cex.axis=1.5)
# lines(rasy.syn$real$predicted$x,rasy.syn$real$predicted$y,lwd=3)
# title(main="DBCA Wood Frog Egg Mass Counts\n2006-2026")
# text(x = 0.1,y=0.9,labels="Regional Avg: 0.12",adj=0)
#dev.off()

png("rasy_dist.png",res=500,height=5,width=5,units="in")
ggplot()+
  geom_ribbon(aes(x=as.vector(rasy.syn$boot$boot.summary$predicted$x),
                  ymin = as.vector(rasy.syn$boot$boot.summary$predicted$y[2,]),
                  ymax = as.vector(rasy.syn$boot$boot.summary$predicted$y[10,])),fill="gray")+
  geom_point(aes(as.vector(distmat),as.vector(cormat)),shape=21,fill="cornflowerblue")+
  labs(y="Pearson Correlation",x="Distance (km)")+
  geom_line(aes(as.vector(rasy.syn$real$predicted$x),as.vector(rasy.syn$real$predicted$y)),size=1)+
  theme_classic()+
  scale_x_continuous(breaks=seq(0,3,0.5),limits=c(0,3.1))+
  scale_y_continuous(breaks=seq(-.5,1,.25))+
  geom_hline(yintercept=0)+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=17))
dev.off()

# Geography of Synchrony --------------------------------------------------

#Truncate data to 2010-present to include more ponds (n =32 ponds)
rasy_17yr <- rasy %>%
  filter(if_all(X2010:X2026, ~ !is.na(.)))%>%
  filter(if_any(X2010:X2026, ~ . != 0))
rasy_17yr_mat<-rasy_17yr%>%
  select(X2010:X2026)

#Scale data
rasy_17yr_dt<-cleandat(as.matrix(rasy_17yr_mat),times=2010:2026,clev=3)
rasy.syn17yr<-Sncf(x=rasy_17yr$Longitude,
                   y=rasy_17yr$Latitude,
                   z=rasy_17yr_dt$cdat,resamp = 100,latlon = T,quiet=T)
plot(rasy.syn17yr)
title(main="DBCA Wood Frog Egg Mass Counts\n2010-2026")
summary(rasy.syn17yr)$Regional.sync

# Matrix Regression to Examine Geography ----------------------------------

#create density dependence matrix
# Function to calculate DD slope for a single site time series
dd_slope <- function(x) {
  xt  <- x[-length(x)]
  xt1 <- x[-1]
  keep <- !is.na(xt) & !is.na(xt1)  # zeros are real, only drop NAs
  xt  <- xt[keep]
  xt1 <- xt1[keep]
  if (sum(keep) < 2) return(NA)
  coef(lm(log1p(xt1) ~ log1p(xt)))[2]
}

# Function to calculate direct (t-1) and delayed (t-2) DD slopes for one site
dd_slope_lag2 <- function(x) {
  n <- length(x)
  if (n < 4) return(c(beta1 = NA, beta2 = NA))  # need enough points to fit 2 predictors
  
  xt  <- x[3:n]        # response, time t
  xt1 <- x[2:(n - 1)]  # lag 1 (t-1) -- direct DD
  xt2 <- x[1:(n - 2)]  # lag 2 (t-2) -- delayed DD
  
  keep <- !is.na(xt) & !is.na(xt1) & !is.na(xt2)
  xt  <- xt[keep]; xt1 <- xt1[keep]; xt2 <- xt2[keep]
  
  if (sum(keep) < 4) return(c(beta1 = NA, beta2 = NA))  # need df for 2 slopes + intercept
  
  fit <- lm(log1p(xt) ~ log1p(xt1) + log1p(xt2))
  co  <- coef(fit)
  c(beta1 = unname(co[2]), beta2 = unname(co[3]))
}

# Master function: takes a matrix (sites x years), can be of multiple species stored in a list
# and returns a pairwise Euclidean distance matrix 
# calc_eucldd <- function(sp_list) {
#   n <- nrow(sp_list[[1]])
#   
#   # --- site means and DD slopes ---
#   site.means <- lapply(sp_list, function(mat) rowMeans(log1p(mat), na.rm = TRUE))
#   beta       <- lapply(sp_list, function(mat) apply(mat, 1, dd_slope))
#   
#   # --- scale by pooled SD across species ---
#   sdM <- mean(sapply(site.means, sd, na.rm = TRUE))
#   sdS <- mean(sapply(beta,       sd, na.rm = TRUE))
#   
#   site.means.sc <- lapply(site.means, function(x) x / sdM)
#   beta.sc        <- lapply(beta,       function(x) x / sdS)
#   
#   # --- pairwise Euclidean distance per species ---
#   eucldd <- lapply(names(sp_list), function(i) {
#     tmp <- matrix(NA, n, n)
#     for (k in 1:(n - 1)) {
#       for (l in (k + 1):n) {
#         tmp[l, k] <- sqrt(
#           (site.means.sc[[i]][k] - site.means.sc[[i]][l])^2 +
#             (beta.sc[[i]][k] - beta.sc[[i]][l])^2
#         )
#       }
#     }
#     as.matrix(as.dist(tmp, upper = TRUE, diag = TRUE))
#   })
#   names(eucldd) <- names(sp_list)
#   
#   list(eucldd = eucldd, beta = beta, site.means = site.means)  # return intermediates too
# }

calc_dd_dissim <- function(sp_list, lag = 1) {
  n <- nrow(sp_list[[1]])
  
  site.means <- lapply(sp_list, function(mat) rowMeans(log1p(mat), na.rm = TRUE))
  
  if (lag == 1) {
    beta <- lapply(sp_list, function(mat) apply(mat, 1, dd_slope))
    beta_list <- list(beta1 = beta)
  } else if (lag == 2) {
    betas <- lapply(sp_list, function(mat) apply(mat, 1, dd_slope_lag2))
    beta_list <- list(
      beta1 = lapply(betas, function(b) b["beta1", ]),
      beta2 = lapply(betas, function(b) b["beta2", ])
    )
  }
  
  # simple 1-D pairwise distance helper (|xi - xj|), scaled by pooled SD
  dist1d <- function(x_list) {
    sdX <- mean(sapply(x_list, sd, na.rm = TRUE))
    lapply(x_list, function(x) {
      xs <- x / sdX
      as.matrix(dist(xs))
    })
  }
  
  list(
    mean_dist = dist1d(site.means), #matrix of mean site abundance (log scale)
    beta_dist = lapply(beta_list, dist1d) #matrices of mean site abundance (log scale)
  )
}

#rasy_dd<-calc_eucldd(list(RASY=rasy_17yr_mat))
rasy_dd<-calc_dd_dissim(list(RASY=rasy_17yr_mat),lag=2)

#turn mean abundance and lag1 and lag 2 Pollard estimates into dissimilarity matrices
meanabun_dist<-dist(rasy_dd$mean_dist$RASY)
dd1_dist<-dist(rasy_dd$beta_dist$beta1$RASY)
dd2_dist<-dist(rasy_dd$beta_dist$beta2$RASY)

#make synchrony matrix
cor17yr_mat<-synmat(rasy_17yr_dt$cdat,times = 2010:2026,method = "pearson")
cor17yr_xwtmat<-synmat(rasy_17yr_dt$cdat,times = 2010:2026,method = "ReXWT",tsrange=c(4,6))

#make predictors as distance objects
dbca_topo<-readRDS('data/dbca_topo.rds')
dbca_topo17yr<-dbca_topo%>%
  sf::st_drop_geometry() %>%
  dplyr::filter(pondname%in%c(rasy_17yr$pondname))
east<-dist(as.numeric(dbca_topo17yr$eastness))
north<-dist(as.numeric(dbca_topo17yr$northness))
slope<-dist(as.numeric(dbca_topo17yr$slope))
hydrodist<-dist(as.numeric(as.factor(rasy_17yr$Hydroperiod))-1)
areadist<-dist(as.numeric(rasy_17yr$Area))

distmat17yr<-gcdist(rasy_17yr$Longitude,rasy_17yr$Latitude)
geog.dist<-distmat17yr[lower.tri(distmat17yr)]

#run models with distance
#predictors- aspect, elevation?, area, hydroperiod, distance, density dependence in egg stage, density dependence in terrestrial environment
library(ecodist)
#with distance as covariate
areamrm<- MRM(as.dist(cor17yr_mat)~areadist+as.dist(distmat17yr)) #has a lot of NAs
areamrm
hydromrm<- MRM(as.dist(cor17yr_mat)~hydrodist+as.dist(distmat17yr))
hydromrm
meanabunmrm<- MRM(as.dist(cor17yr_mat)~meanabun_dist+as.dist(distmat17yr)) #mean site abundance
meanabunmrm
dd1mrm<- MRM(as.dist(cor17yr_mat)~dd1_dist+as.dist(distmat17yr)) #lag 1 density dependence
dd1mrm
dd2mrm<- MRM(as.dist(cor17yr_mat)~dd2_dist+as.dist(distmat17yr)) #lag 2 density dependence
dd2mrm



#without distance as covariate
meanabunmrm<- MRM(as.dist(cor17yr_mat)~meanabun_dist) #mean site abundance
meanabunmrm
dd1mrm<- MRM(as.dist(cor17yr_mat)+dd1_dist) #lag 1 density dependence
dd1mrm
dd2mrm<- MRM(as.dist(cor17yr_mat)~dd2_dist) #lag 2 density dependence
dd2mrm
areamrm2<- MRM(as.dist(cor17yr_mat)~areadist)
areamrm2
hydromrm2<- MRM(as.dist(cor17yr_mat)~hydrodist)
hydromrm2
northmrm2<- MRM(as.dist(cor17yr_mat)~north)
northmrm2
eastmrm2<- MRM(as.dist(cor17yr_mat)~east)
eastmrm2
slopemrm2<- MRM(as.dist(cor17yr_mat)~slope)
slopemrm2
topomrm2<- MRM(as.dist(cor17yr_mat)~slope+east+north)
topomrm2

#cross wavelet geography test
#with distance as covariate
areamrm<- MRM(as.dist(cor17yr_xwtmat)~areadist+as.dist(distmat17yr))
areamrm
hydromrm<- MRM(as.dist(cor17yr_xwtmat)~hydrodist+as.dist(distmat17yr))
hydromrm
meanabunmrm<- MRM(as.dist(cor17yr_xwtmat)~meanabun_dist+as.dist(distmat17yr)) #mean site abundance
meanabunmrm
dd1mrm<- MRM(as.dist(cor17yr_xwtmat)~dd1_dist+as.dist(distmat17yr)) #lag 1 density dependence
dd1mrm
dd2mrm<- MRM(as.dist(cor17yr_xwtmat)~dd2_dist+as.dist(distmat17yr)) #lag 2 density dependence
dd2mrm

#without distance as covariate
meanabunmrm<- MRM(as.dist(cor17yr_xwtmat)~meanabun_dist) #mean site abundance
meanabunmrm
dd1mrm<- MRM(as.dist(cor17yr_xwtmat)~dd1_dist) #lag 1 density dependence
dd1mrm
dd2mrm<- MRM(as.dist(cor17yr_xwtmat)~dd2_dist) #lag 2 density dependence
dd2mrm
areamrm2<- MRM(as.dist(cor17yr_xwtmat)~areadist)
areamrm2
hydromrm2<- MRM(as.dist(cor17yr_xwtmat)~hydrodist)
hydromrm2
northmrm2<- MRM(as.dist(cor17yr_xwtmat)~north)
northmrm2
eastmrm2<- MRM(as.dist(cor17yr_xwtmat)~east)
eastmrm2
slopemrm2<- MRM(as.dist(cor17yr_xwtmat)~slope)
slopemrm2
topomrm2<- MRM(as.dist(cor17yr_xwtmat)~slope+east+north)
topomrm2

#

# Tail Dependency ---------------------------------------------------------

#seems like I must need spatiotemporal variable for driver? single weather time series wouldn't work? must be positively correlated?

## -----------------------------------------------------------------------
## PART 3. Pairwise tail-dependent synchrony matrices across all ponds
## -----------------------------------------------------------------------
source('Fns_taildependence.R')
source('code/tail_functions.R')

pond_ids <- rasy_17yr$pondname
n        <- nrow(rasy_17yr_mat)
n_years  <- ncol(rasy_17yr_mat)

obs      <- tail_matrices(rasy_17yr_dt$cdat, tail_dependence)
M_upper  <- obs$upper
M_lower  <- obs$lower
M_taildep <- obs$taildep

#M_upper+M_lower= should equal total spearman correlation but don't
uppervec<-as.vector(M_upper)
lowervec<-as.vector(M_lower)
allvec<-uppervec+lowervec
rasy.spear<-clust(dat = rasy_17yr_dt$cdat,
                  times = 2010:2026,
                  coords = data.frame(lon=rasy_17yr$Longitude,
                                      lat=rasy_17yr$Latitude),
                  method="spearman",
                  weighted = T)
spearcor<-as.vector(rasy.spear$adj)
spearcor[is.na(spearcor)]<-1
cbind(spearcor,allvec)
sum(spearcor==allvec)/length(spearcor)
#double check all functions that things in correct matrix set up (years vs rows in columns)

# cor() correlates columns, so transpose the site x year matrix to get a
# site x site Spearman matrix (pairwise-complete to tolerate missing years)
emp_spearman <- suppressWarnings(
  cor(t(rasy_17yr_dt$cdat), method = "spearman", use = "pairwise.complete.obs")
)

n_surr    <- 999
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

library(ncf)
test1<-spline.correlog.simmat(coords=cbind(rasy_17yr$Longitude,rasy_17yr$Latitude),M_upper)
test2<-spline.correlog.simmat(coords=cbind(rasy_17yr$Longitude,rasy_17yr$Latitude),M_lower)
test3<-spline.correlog.simmat(coords=cbind(rasy_17yr$Longitude,rasy_17yr$Latitude),emp_spearman)
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
rownames(rasy_17yr_dt$cdat) <- rasy_17yr$pondname
colnames(rasy_17yr_dt$cdat) <- 2010:2026        # your vector of years, 2010-2026 or the 17 you use
drv <- matrix(rep(climind$enso_lag1[5:21], each = nrow(rasy_17yr_dt$cdat)), nrow = nrow(rasy_17yr_dt$cdat),
              dimnames = list(rownames(rasy_17yr_dt$cdat), colnames(rasy_17yr_dt$cdat)))

pond_driver_taildep <- map_dfr(pond_ids, function(p) {
  td <- tail_dependence(rasy_17yr_dt$cdat[p, ], drv[p, ])
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
coords<-data.frame(pond=rasy_17yr$pondname,Longitude=rasy_17yr$Longitude,Latitude=rasy_17yr$Latitude)
reg_df <- pond_sync_taildep %>%
  left_join(pond_driver_taildep, by = "pond") %>%
  left_join(coords, by = "pond") %>%
  filter(is.finite(sync_tail_dep), is.finite(driver_tail_dep)) %>%
  # coordinates in km (approximate), so the spatial range parameter is on a
  # sensible scale; in raw degrees gls() failed to converge on the simulated data
  mutate(x_km = (Longitude - min(Longitude)) * 111 * cos(mean(Latitude) * pi / 180),
         y_km = (Latitude - min(Latitude)) * 111)

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

# Test Modularity and Clustering ------------------------------------------
diag(cor17yr_mat)<-0
m<-cluseigen(cor17yr_mat)
modularity(cor17yr_mat,m[[length(m)]],decomp = T) #what does total Q mean

rasy.corclust<-clust(dat = rasy_17yr_dt$cdat,
                     times = 2010:2026,
                     coords = data.frame(lon=rasy_17yr$Longitude,
                                         lat=rasy_17yr$Latitude),
                     method="pearson",
                     weighted = T)
#figure out mean synchrony (pearson correlation) by cluster
mean(rasy.corclust$adj[rasy.corclust$clusters[[2]]==1,],na.rm=T)
mean(rasy.corclust$adj[rasy.corclust$clusters[[2]]==2,],na.rm=T)
plotmap(rasy.corclust)

rasy.xwtclust<-clust(dat = rasy_17yr_dt$cdat,
                     times = 2010:2026,
                     coords = data.frame(lon=rasy_17yr$Longitude,
                                         lat=rasy_17yr$Latitude),
                     method="ReXWT",
                     tsrange = c(4,6),
                     weighted = T)
#png('xwt_clusters.png',res=500,height=5,width=5,units="in")
plotmap(rasy.xwtclust)
#dev.off()

rasy.xwtclust23<-clust(dat = rasy_17yr_dt$cdat,
                     times = 2010:2026,
                     coords = data.frame(lon=rasy_17yr$Longitude,
                                         lat=rasy_17yr$Latitude),
                     method="ReXWT",
                     tsrange = c(2,3),
                     weighted = T)
plotmap(rasy.xwtclust23)

rasy.xwtclust<-addwpmfs(rasy.xwtclust)
par(mfrow=c(1,2))
#png('xwt_mf_clust2.png',res=500,width=5,height=4,units="in")
plotmag(get_wpmfs(rasy.xwtclust)[[2]][[1]],las=1)
plotmag(get_wpmfs(rasy.xwtclust)[[2]][[2]])

mean(rasy.xwtclust$adj[rasy.xwtclust$clusters[[2]]==1,],na.rm=T) #0.22
mean(rasy.xwtclust$adj[rasy.xwtclust$clusters[[2]]==2,],na.rm=T) #0.43

#cluster analysis following Jasper Ridge paper

distmat17yr<-gcdist(rasy_17yr$Longitude,rasy_17yr$Latitude)
geog.dist<-distmat17yr[lower.tri(distmat17yr)]

#define function for doing spling fitting
spline.fit<-function(dists,zmat,nresamp=1000,quantiles=c(0.025,0.975)){
  xemp<-c(dists)
  triang<-lower.tri(zmat)
  yemp<-zmat[triang]
  dfs=sqrt(nrow(zmat))
  distmat<-full(dists)
  out<-list()
  
  emp.spline<-smooth.spline(xemp,yemp,df=dfs)
  out$emp.spline<-emp.spline
  
  resamp.splines<-matrix(NA, nrow=nresamp, ncol=length(emp.spline$y))
  for(ii in 1:nresamp){
    shuffle<-sample(1:nrow(zmat), size=nrow(zmat), replace=TRUE)
    xres<-distmat[shuffle,shuffle][triang]
    yres<-zmat[shuffle,shuffle][triang]
    drop.NaNs<-!is.na(yres)
    xres<-xres[drop.NaNs]; yres<-yres[drop.NaNs]
    res.spline<-smooth.spline(xres,yres,df=dfs)
    resamp.splines[ii,]<-predict(res.spline, x=emp.spline$x)$y
  }
  out$resamp.splines<-resamp.splines
  out$spline.quantiles<-apply(resamp.splines,2,quantile,probs=quantiles)
  return(out)
}

#do spline fits
splinexwt.rasy<-spline.fit(geog.dist,rasy.xwtclust$adj) #3-5 yr timescales

colors=c("red","blue") #blue for in-module, red for out-of-module (or is this cluster 1 and cluster 2?)

mod.in.out<-function(mod.ids){
  xx<-matrix(1, nrow=length(mod.ids), ncol=length(mod.ids))
  for(i in 1:length(mod.ids)){
    for(j in 1:length(mod.ids)){
      if(mod.ids[i]==mod.ids[j]){xx[i,j]<-2}
    }
  }
  return(xx)
}
rasy.xwtinout<-mod.in.out(rasy.xwtclust$clusters[[length(rasy.xwtclust$clusters)]])

png('crosswavelet_clusters.png',res=600,height=5,width=5,units="in")
plot(splinexwt.rasy$emp.spline$x, splinexwt.rasy$emp.spline$y, type="l",# main="Clusters of ReXWT (4-6 yr)",
     xlab="Distance (km)", xaxt="n",las=1,ylab="Synchrony", ylim=c(-0.5,1), cex.main=1, lwd=2)
abline(h=0, col="grey")
axis(1,at=seq(0,4,0.5))
#abline(h=mean(lower(clust.mido.st$adj)),col="grey")
lines(splinexwt.rasy$emp.spline$x, splinexwt.rasy$emp.spline$y,lwd=2)
lines(splinexwt.rasy$emp.spline$x, splinexwt.rasy$spline.quantiles[1,], lty=2)
lines(splinexwt.rasy$emp.spline$x, splinexwt.rasy$spline.quantiles[2,], lty=2)
points(geog.dist, lower(rasy.xwtclust$adj), pch=20, col=colors[lower(rasy.xwtinout)],cex=0.8)
points(4,mean(lower(rasy.xwtclust$adj*ifelse(rasy.xwtinout==2,1,NA)),na.rm=T),col="blue",pch=18,cex=2)
points(4,mean(lower(rasy.xwtclust$adj*ifelse(rasy.xwtinout==1,1,NA)),na.rm=T),col="red",pch=18,cex=2)
#dev.off()

#using pearman correlation
splinecor.rasy<-spline.fit(geog.dist,rasy.corclust$adj)

colors=c("red","blue") #blue for in-module, red for out-of-module

rasy.corinout<-mod.in.out(rasy.corclust$clusters[[length(rasy.corclust$clusters)]])

#clusters using corman correlations
plot(splinecor.rasy$emp.spline$x, splinecor.rasy$emp.spline$y, type="l", main="Clusters of corman Correlation",
     xlab="Distance (km)", ylab="Synchrony", ylim=c(-0.5,1), cex.main=1, lwd=2)
abline(h=0, col="grey")
#abline(h=mean(lower(clust.mido.st$adj)),col="grey")
lines(splinecor.rasy$emp.spline$x, splinecor.rasy$emp.spline$y,lwd=2)
lines(splinecor.rasy$emp.spline$x, splinecor.rasy$spline.quantiles[1,], lty=2)
lines(splinecor.rasy$emp.spline$x, splinecor.rasy$spline.quantiles[2,], lty=2)
points(geog.dist, lower(rasy.corclust$adj), pch=20, col=colors[lower(rasy.corinout)],cex=0.8)
points(4,mean(lower(rasy.corclust$adj*ifelse(rasy.corinout==2,1,NA)),na.rm=T),col="blue",pch=18,cex=2)
points(4,mean(lower(rasy.corclust$adj*ifelse(rasy.corinout==1,1,NA)),na.rm=T),col="red",pch=18,cex=2)


ClustSignif<-function(d,numreps,burnin,shortts,longts)
{
  #get covariance matrix, which should be the same as the correlation matrix
  covd<-cov(d)
  if (max(abs(diag(covd)-1))>1e-5)
  {
    stop("Error at loc 1: d should have a covariance matrix with 1s on the diagonal")
  }
  
  #get autocorrelations, one for each plot
  numplots<-dim(d)[2]
  rhos<-NA*numeric(numplots)
  for (counter in 1:numplots)
  {
    h<-d[,counter]
    x<-unname(h[1:(length(h)-1)])
    y<-unname(h[2:length(h)])
    modl<-lm(y~x)
    rhos[counter]<-unname(coef(modl)[2])
  }
  
  #put autocorrelations in a convenient form for later
  #rhos<-matrix(rep(rhos,times=numreps),numplots,numreps) #makes sims have the same autocorrelations, in each plot, as the data
  rhos<-matrix(mean(rhos),numplots,numreps) #makes sims have the average autocorrelation in all plots
  
  #do the simulations
  lendat<-dim(d)[1]
  surrs<-array(NA,c(burnin+lendat,numplots,numreps))
  surrs[1,,]<-0
  mncov<-(sum(covd)-sum(diag(covd)))/(numplots^2-numplots)
  sig<-matrix(mncov,numplots,numplots)
  diag(sig)<-1
  noise<-MASS::mvrnorm((burnin+lendat-1)*numreps,mu=rep(0,numplots),Sigma=sig)
  noise<-array(noise,c(burnin+lendat-1,numreps,numplots))
  noise<-aperm(noise,c(1,3,2))
  for (counter in 2:(burnin+lendat))
  {
    surrs[counter,,]<-rhos*surrs[counter-1,,]+sqrt(1-rhos^2)*noise[counter-1,,]
  }
  surrs<-surrs[(burnin+1):(dim(surrs)[1]),,]
  
  #get adjacency matrices, perform clustering, and get modularity for the real data
  fakecoords<-data.frame(X=1:numplots,Y=1:numplots)
  clustres_short<-wsyn::clust(dat=t(d),times=1:lendat,coords=fakecoords,method="ReXWT",tsrange=shortts)
  datmod_short<-clustres_short$modres[[length(clustres_short$modres)]]$totQ
  clustres_long<-wsyn::clust(dat=t(d),times=1:lendat,coords=fakecoords,method="ReXWT",tsrange=longts)
  datmod_long<-clustres_long$modres[[length(clustres_long$modres)]]$totQ
  
  #get adjacency matrices and perform clustering for each surrogate dataset
  surrmod_short<-NA*numeric(numreps)
  surrmod_long<-NA*numeric(numreps)
  nclust_short<-rep(NA, numreps)
  nclust_long<-rep(NA, numreps)
  for (counter in 1:numreps)
  {
    #print(paste0("Working on surrogate ",counter," of ",numreps))
    surrdat<-t(surrs[,,counter])
    surrdat<-wsyn::cleandat(surrdat,times=1:lendat,clev=2)$cdat
    h_short<-wsyn::clust(dat=surrdat,times=1:lendat,coords=fakecoords,method="ReXWT",tsrange=shortts)
    surrmod_short[counter]<-h_short$modres[[length(h_short$modres)]]$totQ
    nclust_short[counter]<-length(h_short$clusters)
    h_long<-wsyn::clust(dat=surrdat,times=1:lendat,coords=fakecoords,method="ReXWT",tsrange=longts)
    surrmod_long[counter]<-h_long$modres[[length(h_long$modres)]]$totQ
    nclust_long[counter]<-length(h_long$clusters)
    
  }
  
  #get info for p values
  fracmore_short<-sum(surrmod_short>datmod_short)/numreps
  fracmore_long<-sum(surrmod_long>datmod_long)/numreps
  #fracmore_tab<-matrix(c(shortts,fracmore_short,longts,fracmore_long),nrow=2,ncol=3,byrow=T)
  fracmore_tab<-data.frame(timescales=c(paste(shortts[1],shortts[2],sep="-"),paste(longts[1],longts[2],sep="-")),
                           pvals=c(fracmore_short,fracmore_long))
  
  return(list(d=d,surrs=surrs,
              datmod_short=datmod_short,surrmod_short=surrmod_short,fracmore_short=fracmore_short,
              datmod_long=datmod_long,surrmod_long=surrmod_long,fracmore_long=fracmore_long, 
              pvals=fracmore_tab,nclust_short=nclust_short, nclust_long=nclust_long))
}

d_rasy<-t(rasy_17yr_dt$cdat) #why do we transpose this?
clustsig.rasy<-ClustSignif(d_rasy,numreps=1000,burnin=500,shortts=c(2,3),longts=c(4,6))

clustsig.rasy$pvals

hist(clustsig.rasy$surrmod_short,xlab="Surrogate Q Values",main="Real Q vs Surrogate Q")
abline(v=clustsig.rasy$datmod_short,col="red")
abline(v=clustsig.rasy$datmod_long,col="blue")

text(x=0.1,y=250,"P = 0.95")

source("../../jasperridge/plotClusterMap_JRG.R")
par(mfrow=c(1,1))
plotClusterMap(rasy.xwtclust, 
               nodewgt=NULL, 
               edgewgt=0.9, 
               nodesize=c(2,2), 
               edgesize = c(0.5,0.5), 
               addlegend = T,
               nodechar=c(1,2))

