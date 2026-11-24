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
png("results/mo_rasy.timeseries.png",height=5,width=5,res=500,units="in")
#pdf('woodfrog_synchrony.pdf')
ggplot(rasy_long,aes(Year,Count,color=pondname))+
  geom_point(show.legend = F)+
  geom_line(aes(color=pondname),show.legend = F)+
  labs(y="Number of Egg Masses")+
  theme_classic()+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))+
  scale_y_continuous(breaks=seq(0,200,25))+
  scale_x_continuous(breaks=seq(2006,2026,1),labels=c("2006",rep("",3),"2010",rep("",3),"2014",rep("",3),"2018",rep("",3),"2022",rep("",3),"2026"))
dev.off()

#Remove data sets that are > 50% 0 or all NA
rasymat <- rasy %>%
  select(X2006:X2026)

keep <- rowMeans(rasymat == 0, na.rm = TRUE) <= 0.5 & rowSums(!is.na(rasymat)) > 0
rasymat_filtered <- rasymat[keep, , drop = FALSE]

#scale matrix
rasymat_dt<-rasymat_filtered%>%
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
rasy_filtered <- rasy[keep, , drop = FALSE]

#perform cross correlation analysis
rasy.syn<-Sncf(x=rasy_filtered$Longitude,
               y=rasy_filtered$Latitude,na.rm=T,
               z=rasymat_dt,resamp = 1000,latlon = T,quiet=T)

#find degree of synchrony across study area
summary(rasy.syn)$Regional.sync #0.12

#png("rasy_sncf.png",res=500,height=5,width=5,units="in")
par(mar=c(4,5,2,1))
plot(rasy.syn,xlab="Distance (km)",cex.lab=1.5,cex.axis=1.25,las=1,ylim=c(-1,1))
title(main="DBCA Wood Frog Egg Mass Counts, 2006-2026",cex=0.5)
#dev.off()

#manually plot geographic distance vs correlation
distmat<-gcdist(rasy_filtered$Longitude,rasy_filtered$Latitude)
cormat <- cor(t(as.matrix(rasymat_dt)), method="pearson",use = "pairwise.complete.obs")
diag(cormat) <- NA

png("results/Fig2_distdecay.png",res=500,height=5,width=5,units="in")
ggplot()+
  geom_ribbon(aes(x=as.vector(rasy.syn$boot$boot.summary$predicted$x),
                  ymin = as.vector(rasy.syn$boot$boot.summary$predicted$y[2,]),
                  ymax = as.vector(rasy.syn$boot$boot.summary$predicted$y[10,])),fill="gray")+
  geom_point(aes(as.vector(distmat),as.vector(cormat)),shape=21,fill="cornflowerblue")+
  labs(y="Pearson Correlation",x="Distance (km)")+
  geom_line(aes(as.vector(rasy.syn$real$predicted$x),as.vector(rasy.syn$real$predicted$y)),size=1)+
  theme_classic()+
  scale_x_continuous(breaks=seq(0,5,0.5),limits=c(0,5),labels=c(0,"",1,"",2,"",3,"",4,"",5))+
  scale_y_continuous(breaks=seq(-1,1,.25),limits=c(-1,1))+
  geom_hline(yintercept=0)+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=17))
dev.off()

#perform cross correlation analysis on complete time series (rasy_comp from timescale code)
gps<-rasy%>%
  filter(pondname%in%c(unique(rasy_comp$pondname)))%>%
  select(pondname,Latitude,Longitude)%>%
  arrange(pondname)

rasy.syn<-Sncf(x=gps$Longitude,
               y=gps$Latitude,na.rm=T,
               z=rasy_compmat,resamp = 1000,latlon = T,quiet=T)

#find degree of synchrony across study area
summary(rasy.syn)$Regional.sync #0.12

#png("rasy_sncf.png",res=500,height=5,width=5,units="in")
par(mar=c(4,5,2,1))
plot(rasy.syn,xlab="Distance (km)",cex.lab=1.5,cex.axis=1.25,las=1,ylim=c(-1,1))
title(main="DBCA Wood Frog Egg Mass Counts, 2006-2026",cex=0.5)
#dev.off()

#manually plot geographic distance vs correlation
distmat<-gcdist(gps$Longitude,gps$Latitude)
cormat <- cor(t(as.matrix(rasy_compmat)), method="pearson")
diag(cormat) <- NA

png("results/rasycomp_dist.png",res=500,height=5,width=5,units="in")
mo_sncf<-ggplot()+
  geom_ribbon(aes(x=as.vector(rasy.syn$boot$boot.summary$predicted$x),
                  ymin = as.vector(rasy.syn$boot$boot.summary$predicted$y[2,]),
                  ymax = as.vector(rasy.syn$boot$boot.summary$predicted$y[10,])),fill="gray")+
  geom_point(aes(as.vector(distmat),as.vector(cormat)),shape=21,fill="cornflowerblue")+
  labs(y="Pearson Correlation",x="Distance (km)")+
  geom_line(aes(as.vector(rasy.syn$real$predicted$x),as.vector(rasy.syn$real$predicted$y)),size=1)+
  theme_classic()+
  scale_x_continuous(breaks=seq(0,3.5,0.5),limits=c(0,3.5))+
  scale_y_continuous(breaks=seq(-0.75,1,.25),limits=c(-0.75,1))+
  geom_hline(yintercept=0)+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=17))
mo_sncf
dev.off()
saveRDS(mo_sncf,'results/mo_sncf.rds')

# Geography of Synchrony --------------------------------------------------
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
best <- tr[tr$n_years == 17, ]
plot_coverage(m, best$start, best$end, min_complete = 0.9, max_zero_prop = 0.5)
dat_use <- window_sites(m, best$start, best$end, min_complete = 0.9,
                        max_zero_prop = 0.5)
attr(dat_use, "site_stats")       # NAs and zeros per site in the window
attr(dat_use, "n_NA_remaining")   # values to impute

rules <- list(min_complete = 0.9, max_zero_prop = 0.5, min_nonzero = 8)
start <- 2010; end <- 2026          # your chosen window

sel <- do.call(window_sites, c(list(m, start, end), rules))

keep_sites <- rownames(sel)
keep_years <- start:end

rasy_17yr <- sel %>%
  as.data.frame()%>%
  mutate(row_mean = round(rowMeans(select(., `2010`:`2026`), na.rm = TRUE),0),
         across(`2010`:`2026`, ~ coalesce(., row_mean)))%>%
  mutate(pondname=row.names(sel))%>%
  left_join(rasy%>%select(pondname,Latitude,Longitude),by="pondname")

rasy_17yr_mat<-rasy_17yr%>%
  select(`2010`:`2026`)%>%
  as.matrix()
#24 sites for 17 years which drops zero inflated ponds

#Scale data
rasy_17yr_dt<-cleandat(as.matrix(rasy_17yr_mat),times=2010:2026,clev=5)$cdat

#confirm timescale structure still exists in shorter time series
mo_wmf17<-wmf(rasy_17yr_dt,keep_years)
plotmag(mo_wmf17)
mo_wpmf17<-wpmf(rasy_17yr_dt,keep_years,sigmethod = "quick")
plotmag(mo_wpmf17,sigthresh=0.99)

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

# Master function: takes a matrix (sites x years), can be of multiple species (sites) stored in a list
# and returns a pairwise Euclidean distance matrix 

calc_dd_dissim <- function(d_list, lag = 1) {
  n <- nrow(d_list[[1]])
  
  site.means <- lapply(d_list, function(mat) rowMeans(log1p(mat), na.rm = TRUE))
  
  if (lag == 1) {
    beta <- lapply(d_list, function(mat) apply(mat, 1, dd_slope))
    beta_list <- list(beta1 = beta)
  } else if (lag == 2) {
    betas <- lapply(d_list, function(mat) apply(mat, 1, dd_slope_lag2))
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

rasy_dd<-calc_dd_dissim(list(RASY=rasy_17yr_mat),lag=2)

#turn mean abundance and lag1 and lag 2 Pollard estimates into dissimilarity matrices
meanabun_dist<-dist(rasy_dd$mean_dist$RASY)
dd1_dist<-dist(rasy_dd$beta_dist$beta1$RASY)
dd2_dist<-dist(rasy_dd$beta_dist$beta2$RASY)

#make predictors as distance objects
dbca_topo<-readRDS('data/dbca_topo.rds')
dbca_topo17yr<-dbca_topo%>%
  sf::st_drop_geometry() %>%
  dplyr::filter(pondname%in%c(rasy_17yr$pondname))
head(dbca_topo17yr)

#set up topographical variables as distance matrices (dissimilarity)
east<-dist(as.numeric(dbca_topo17yr$eastness))
north<-dist(as.numeric(dbca_topo17yr$northness))
slope<-dist(as.numeric(dbca_topo17yr$slope))
slopeavg<-dist(as.numeric(dbca_topo17yr$slope_mean_r300))
swnorth<-dist(as.numeric(dbca_topo17yr$sw_northness_r300))
sweast<-dist(as.numeric(dbca_topo17yr$sw_eastness_r300))
tpi<-dist(as.numeric(dbca_topo17yr$tpi_r300))
elevsd<-dist(as.numeric(dbca_topo17yr$elev_sd_r300))

#set up hydro and area distance matrices
rasy_17yr<-merge(rasy_17yr,rasy[,c("pondname","Hydroperiod","Area")],by="pondname")
hydrodist<-dist(as.numeric(as.factor(rasy_17yr$Hydroperiod))-1)
areadist<-dist(as.numeric(rasy_17yr$Area))

distmat17yr<-gcdist(rasy_17yr$Longitude,rasy_17yr$Latitude)
geog.dist<-distmat17yr[lower.tri(distmat17yr)]

#run models with distance
#predictors- aspect, elevation?, area, hydroperiod, distance, density dependence in egg stage, density dependence in terrestrial environment
library(ecodist)

#cross wavelet geography test
#make synchrony matrices
cor17yr_xwt4_6<-synmat(rasy_17yr_dt,times = 2010:2026,method = "ReXWT",tsrange=c(4,10))
cor17yr_xwt2_4<-synmat(rasy_17yr_dt,times = 2010:2026,method = "ReXWT",tsrange=c(2,4))
cor17yr_pear<-synmat(rasy_17yr_dt,times = 2010:2026,method = "pearson")

vars <- list(tpi=tpi, elevsd=elevsd, swnorth=swnorth, sweast=sweast, slopeavg=slopeavg,
             dd1_dist=dd1_dist, dd2_dist=dd2_dist, meanabun_dist=meanabun_dist,
             area_dist=areadist, hydro_dist=hydrodist)

y4_6 <- as.dist(cor17yr_xwt4_6)
y2_4 <- as.dist(cor17yr_xwt2_4)
#ycor <- as.dist(cor17yr_pear)

#mrmrescor <- list()
mrmres4_6 <- list()
mrmres2_4 <- list()

for (nm in names(vars)) {
  dat <- list(y4_6 = y4_6, y2_4 = y2_4, x = as.dist(vars[[nm]]))
  names(dat)[3] <- nm
  
  f4_6 <- as.formula(paste("y4_6 ~", nm))
  f2_4 <- as.formula(paste("y2_4 ~", nm))
  fcor <- as.formula(paste("ycor ~", nm))
  
  mrmres4_6[[nm]] <- MRM(f4_6, data = dat,nperm = 10000)$coef
  mrmres2_4[[nm]] <- MRM(f2_4, data = dat,nperm = 10000)$coef
  #mrmrescor[[nm]] <- MRM(fcor, data = dat,nperm = 10000)$coef
}

mrmrestab4_6<-apply(do.call(rbind,mrmres4_6),2,function(x){round(x,3)})
mrmrestab2_4<-apply(do.call(rbind,mrmres2_4),2,function(x){round(x,3)})
#mrmrestabcor<-apply(do.call(rbind,mrmrescor),2,function(x){round(x,3)})

tabs1<-cbind(mrmrestab2_4,mrmrestab4_6)
write.csv(tabs1,'results/mrm_results.csv',row.names = F)

# Test Modularity and Clustering ------------------------------------------
source('code/geogsync_helperfunctions.R') #SlustSignif, mod.in.out, spline_fit

diag(cor17yr_pear)<-0
m<-cluseigen(cor17yr_pear)
modularity(cor17yr_pear,m[[length(m)]],decomp = T) 

diag(cor17yr_xwt2_4)<-0
m<-cluseigen(cor17yr_xwt2_4)
modularity(cor17yr_xwt2_4,m[[length(m)]],decomp = T) 

diag(cor17yr_xwt4_6)<-0
m<-cluseigen(cor17yr_xwt4_6)
modularity(cor17yr_xwt4_6,m[[length(m)]],decomp = T) 

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

rasy.xwtclust46<-clust(dat = rasy_17yr_dt,
                     times = 2010:2026,
                     coords = data.frame(lon=rasy_17yr$Longitude,
                                         lat=rasy_17yr$Latitude),
                     method="ReXWT",
                     tsrange = c(4,10),
                     weighted = T)
#png('xwt_clusters.png',res=500,height=5,width=5,units="in")

#dev.off()

rasy.xwtclust24<-clust(dat = rasy_17yr_dt,
                     times = 2010:2026,
                     coords = data.frame(lon=rasy_17yr$Longitude,
                                         lat=rasy_17yr$Latitude),
                     method="ReXWT",
                     tsrange = c(2,4),
                     weighted = T)

png('results/timescale_clust.png',res=600,height=7,width=5,units="in")
clustmap46<-plotmap_gg(rasy.xwtclust46)+
  labs(y="Latitude",x="Longitude")
clustmap24<-plotmap_gg(rasy.xwtclust24)+
  labs(y="Latitude",x="")
cowplot::plot_grid(clustmap24,clustmap46,ncol=1,labels=c("A","B"))
dev.off()

rasy.xwtclust46<-addwpmfs(rasy.xwtclust46)
par(mfrow=c(1,2))
#png('xwt_mf_clust2.png',res=500,width=5,height=4,units="in")
plotmag(get_wpmfs(rasy.xwtclust46)[[2]][[1]],las=1)
plotmag(get_wpmfs(rasy.xwtclust46)[[2]][[2]])

mean(rasy.xwtclust46$adj[rasy.xwtclust46$clusters[[2]]==1,],na.rm=T) #0.22
mean(rasy.xwtclust46$adj[rasy.xwtclust46$clusters[[2]]==2,],na.rm=T) #0.43

rasy.xwtclust24<-addwpmfs(rasy.xwtclust24)
par(mfrow=c(1,2))
#png('xwt_mf_clust2.png',res=500,width=5,height=4,units="in")
plotmag(get_wpmfs(rasy.xwtclust24)[[2]][[1]],las=1)
plotmag(get_wpmfs(rasy.xwtclust24)[[2]][[2]])

mean(rasy.xwtclust24$adj[rasy.xwtclust24$clusters[[2]]==1,],na.rm=T) #0.22
mean(rasy.xwtclust24$adj[rasy.xwtclust24$clusters[[2]]==2,],na.rm=T) #0.43

#cluster analysis following Jasper Ridge paper

distmat17yr<-gcdist(rasy_17yr$Longitude,rasy_17yr$Latitude)
geog.dist<-distmat17yr[lower.tri(distmat17yr)]

#do spline fits
splinexwt46.rasy<-spline.fit(geog.dist,rasy.xwtclust46$adj) #4-6 yr timescales

colors=c("red","blue") #blue for in-module, red for out-of-module

rasy.xwtinout46<-mod.in.out(rasy.xwtclust46$clusters[[length(rasy.xwtclust46$clusters)]])

png('crosswavelet_clusters.png',res=600,height=5,width=5,units="in")
plot(splinexwt46.rasy$emp.spline$x, splinexwt46.rasy$emp.spline$y, type="l",# main="Clusters of ReXWT (4-6 yr)",
     xlab="Distance (km)", xaxt="n",las=1,ylab="Synchrony", ylim=c(-1,1), cex.main=1, lwd=2)
abline(h=0, col="grey")
axis(1,at=seq(0,4,0.5))
#abline(h=mean(lower(clust.mido.st$adj)),col="grey")
lines(splinexwt46.rasy$emp.spline$x, splinexwt46.rasy$emp.spline$y,lwd=2)
lines(splinexwt46.rasy$emp.spline$x, splinexwt46.rasy$spline.quantiles[1,], lty=2)
lines(splinexwt46.rasy$emp.spline$x, splinexwt46.rasy$spline.quantiles[2,], lty=2)
points(geog.dist, lower(rasy.xwtclust46$adj), pch=20, col=colors[lower(rasy.xwtinout46)],cex=0.8)
points(4,mean(lower(rasy.xwtclust46$adj*ifelse(rasy.xwtinout46==2,1,NA)),na.rm=T),col="blue",pch=18,cex=2)
points(4,mean(lower(rasy.xwtclust46$adj*ifelse(rasy.xwtinout46==1,1,NA)),na.rm=T),col="red",pch=18,cex=2)
#dev.off()

#do spline fits
splinexwt24.rasy<-spline.fit(geog.dist,rasy.xwtclust24$adj) #4-6 yr timescales

colors=c("red","blue") #blue for in-module, red for out-of-module

rasy.xwtinout24<-mod.in.out(rasy.xwtclust24$clusters[[length(rasy.xwtclust24$clusters)]])

png('crosswavelet_clusters.png',res=600,height=5,width=5,units="in")
plot(splinexwt24.rasy$emp.spline$x, splinexwt24.rasy$emp.spline$y, type="l",# main="Clusters of ReXWT (4-6 yr)",
     xlab="Distance (km)", xaxt="n",las=1,ylab="Synchrony", ylim=c(-1,1), cex.main=1, lwd=2)
abline(h=0, col="grey")
axis(1,at=seq(0,4,0.5))
#abline(h=mean(lower(clust.mido.st$adj)),col="grey")
lines(splinexwt24.rasy$emp.spline$x, splinexwt24.rasy$emp.spline$y,lwd=2)
lines(splinexwt24.rasy$emp.spline$x, splinexwt24.rasy$spline.quantiles[1,], lty=2)
lines(splinexwt24.rasy$emp.spline$x, splinexwt24.rasy$spline.quantiles[2,], lty=2)
points(geog.dist, lower(rasy.xwtclust24$adj), pch=20, col=colors[lower(rasy.xwtinout24)],cex=0.8)
points(4,mean(lower(rasy.xwtclust24$adj*ifelse(rasy.xwtinout24==2,1,NA)),na.rm=T),col="blue",pch=18,cex=2)
points(4,mean(lower(rasy.xwtclust24$adj*ifelse(rasy.xwtinout24==1,1,NA)),na.rm=T),col="red",pch=18,cex=2)
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

#run significance testing for modularity against surrogate datasets
d_rasy<-t(rasy_17yr_dt)
clustsig.rasy<-ClustSignif(d_rasy,numreps=1000,burnin=500,tsranges=list(short=c(2,4),long=c(4,10)))
clustsig.rasy$summary

clustsigcor.rasy<-ClustSignif(d_rasy,numreps=1000,burnin=500,method="pearson")
clustsigcor.rasy$summary

hist(clustsig.rasy$surrmod,xlab="Surrogate Q Values",main="Real Q vs Surrogate Q")
abline(v=clustsig.rasy$datmod_short,col="red")
abline(v=clustsig.rasy$datmod_long,col="blue")

text(x=0.1,y=250,"P = 0.95")

png('results/qvalues.png',res=500,height=3.5,width=3.5,units="in")
ggplot()+
  geom_histogram(aes(clustsig.rasy$surrmod),fill="lightgray",color="black")+
  geom_vline(aes(xintercept = clustsig.rasy$summary$Q_obs[1],color="Short"))+
  geom_vline(aes(xintercept = clustsig.rasy$summary$Q_obs[2],color="Long"))+
  theme_classic()+
  labs(y="Number of surrogates",x="Q values")+
  scale_color_manual(
    name = "Timescale:",
    values = c("Short" = "steelblue", "Long" = "firebrick")
  )+
  theme(legend.position="top")
dev.off()


