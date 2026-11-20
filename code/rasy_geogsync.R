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

eggs_use <- rasy_long |>
  filter(as.character(pondname) %in% keep_sites,
         Year >= start, Year <= end)
#24 sites for 17 years which drops zero inflated ponds and too many NA values

#Truncate data to 2010-present to include more ponds (n =32 ponds)
rasy_17yr<-eggs_use%>%
  pivot_wider(names_from = "Year",values_from = "Count",names_prefix = "X")%>%
  mutate(
    .row_mean = rowMeans(across(-pondname), na.rm = TRUE),
    across(-c(pondname, .row_mean), ~ coalesce(.x, .row_mean))
  ) %>%
  select(-.row_mean)%>%
  left_join(rasy %>% select(pondname, Latitude, Longitude),by="pondname")
rasy_17yr_mat <- rasy_17yr %>%
  select(X2010:X2026)

#Scale data
rasy_17yr_dt<-cleandat(as.matrix(rasy_17yr_mat),times=2010:2026,clev=5)
# rasy.syn17yr<-Sncf(x=rasy_17yr$Longitude,
#                    y=rasy_17yr$Latitude,
#                    z=rasy_17yr_dt$cdat,resamp = 100,latlon = T,quiet=T)
# plot(rasy.syn17yr)
# title(main="DBCA Wood Frog Egg Mass Counts\n2010-2026")
# summary(rasy.syn17yr)$Regional.sync

#confirm temporal structure still exists in shorter time series
mo_wmf17<-
mo_wpmf<-


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

#rasy_dd<-calc_eucldd(list(RASY=rasy_17yr_mat))
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
cor17yr_xwt4_6<-synmat(rasy_17yr_dt$cdat,times = 2010:2026,method = "ReXWT",tsrange=c(4,6))
cor17yr_xwt2_3<-synmat(rasy_17yr_dt$cdat,times = 2010:2026,method = "ReXWT",tsrange=c(2,3))
#cor17yr_pear<-synmat(rasy_17yr_dt$cdat,times = 2010:2026,method = "pearson")

vars <- list(tpi=tpi, elevsd=elevsd, swnorth=swnorth, sweast=sweast, slopeavg=slopeavg,
             dd1_dist=dd1_dist, dd2_dist=dd2_dist, meanabun_dist=meanabun_dist,
             area_dist=areadist, hydro_dist=hydrodist)

y4_6 <- as.dist(cor17yr_xwt4_6)
y2_3 <- as.dist(cor17yr_xwt2_3)
#ycor <- as.dist(cor17yr_pear)

#mrmrescor <- list()
mrmres4_6 <- list()
mrmres2_3 <- list()

for (nm in names(vars)) {
  dat <- list(y4_6 = y4_6, y2_3 = y2_3, x = as.dist(vars[[nm]]))
  names(dat)[3] <- nm
  
  f4_6 <- as.formula(paste("y4_6 ~", nm))
  f2_3 <- as.formula(paste("y2_3 ~", nm))
  fcor <- as.formula(paste("ycor ~", nm))
  
  mrmres4_6[[nm]] <- MRM(f4_6, data = dat,nperm = 10000)$coef
  mrmres2_3[[nm]] <- MRM(f2_3, data = dat,nperm = 10000)$coef
  #mrmrescor[[nm]] <- MRM(fcor, data = dat,nperm = 10000)$coef
}

mrmrestab4_6<-apply(do.call(rbind,mrmres4_6),2,function(x){round(x,3)})
mrmrestab2_3<-apply(do.call(rbind,mrmres2_3),2,function(x){round(x,3)})
#mrmrestabcor<-apply(do.call(rbind,mrmrescor),2,function(x){round(x,3)})

tabs1<-cbind(mrmrestab2_3,mrmrestab4_6)
write.csv(tabs1,'results/mrm_results.csv',row.names = F)

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

# source("../../jasperridge/plotClusterMap_JRG.R")
# par(mfrow=c(1,1))
# plotClusterMap(rasy.xwtclust, 
#                nodewgt=NULL, 
#                edgewgt=0.9, 
#                nodesize=c(2,2), 
#                edgesize = c(0.5,0.5), 
#                addlegend = T,
#                nodechar=c(1,2))

