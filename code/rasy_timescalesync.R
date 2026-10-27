###Author: Tom Anderson
###Last Updated: March 9, 2026
library(ggplot2)
library(wsyn)
library(plyr)
library(dplyr)
library(tidyr)
library(Matrix)

#Load and combine data
rasy <- read.csv("data/rasy.csv",fileEncoding="latin1")
pond_info <- read.csv("data/pond_attributes.csv",na.strings=c(""))
rasy<-merge(rasy,pond_info,by="pondname")

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

# Plot counts by year for all sites

#Convert to long format
rasy_long<-rasy%>%
  select(pondname,X2006:X2026)%>%
  pivot_longer(X2006:X2026,names_to="Year",values_to="Count")
rasy_long$Year<-as.numeric(sub('X', '', rasy_long$Year))

### Plot all time series ###
#png("rasy.timeseries.png",height=5,width=5,res=500,units="in")
#pdf('woodfrog_synchrony.pdf')
ggplot(rasy_long,aes(Year,Count,color=pondname))+
  geom_point(show.legend = F)+
  geom_line(aes(color=pondname),show.legend = F)+
  labs(y="Number of Egg Masses",subtitle="All time Series")+
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=17))+
  scale_x_continuous(breaks=seq(2006,2026,1),labels=c("2006",rep("",3),"2010",rep("",3),"2014",rep("",3),"2018",rep("",3),"2022",rep("",3),"2026"))
#dev.off()

# Timescale of Synchrony --------------------------------------------------
#extract time series into a matrix spanning all yrs, dropping ponds where no RASY were ever observed (N = 22 ponds)
rasy_comp <- rasy %>%
  filter(if_all(X2006:X2026, ~ !is.na(.)))%>%
  filter(if_any(X2006:X2026, ~ . != 0))
rasy_compmat <- rasy_comp %>%
  select(X2006:X2026)

#plot complete time series
rasy_comp%>%
  select(pondname,X2006:X2026)%>%
  pivot_longer(X2006:X2026,names_to="Year",values_to="Count")%>%
  mutate(Year=as.numeric(sub('X', '', Year)))%>%
  ggplot(aes(Year,Count,color=pondname))+
  geom_point(show.legend = F)+
  geom_line(aes(color=pondname),show.legend = F)+
  labs(y="Number of Egg Masses",subtitle="Complete Time Series")+
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=17))+
  scale_x_continuous(breaks=seq(2006,2026,1),labels=c("2006",rep("",3),"2010",rep("",3),"2014",rep("",3),"2018",rep("",3),"2022",rep("",3),"2026"))


#same as before but separating by hydroperiod (permanent vs ephemeral ponds)
#to examine whether wmfs differ by hydroperiod
rasy_matperm <- rasy_comp %>%
  filter(Hydroperiod=="perm") %>%
  select(X2006:X2026)
rasy_matephem <- rasy_comp %>%
  filter(Hydroperiod=="ephem") %>%
  select(X2006:X2026)

#standarize data sets using Box-Cox transformation
years<-2006:2026
years1<-2006:2025 #for use with predictors that only go to 2020
rasy_dt<-cleandat(as.matrix(rasy_compmat),clev=5,times=years)
rasyperm_dt<-cleandat(as.matrix(rasy_matperm),clev=5,times=years)
rasyephem_dt<-cleandat(as.matrix(rasy_matephem),clev=5,times=years)

#wavelet phasor mean fields- phase synchrony 
#png('rasy_wpmf.png',res=500,height=4,width=5,units="in")
rasy.wpmf<-wpmf(rasy_dt$cdat,times=years,sigmethod = "quick")
rasyperm.wpmf<-wpmf(rasyperm_dt$cdat,times=years,sigmethod = "quick")
rasyephem.wpmf<-wpmf(rasyephem_dt$cdat,times=years,sigmethod = "quick")
par(mfrow=c(1,3))
plotmag(rasy.wpmf,sigthresh = 0.99,las=1,title="all")
plotmag(rasyperm.wpmf,sigthresh = 0.99,las=1,title = "perm")
plotmag(rasyephem.wpmf,sigthresh = 0.99,las=1,title = "ephem")

#dev.off()

#wavelet mean fields
rasy.wmf<-wmf(rasy_dt$cdat,times=years)
rasyperm.wmf<-wmf(rasyperm_dt$cdat,times=years)
rasyephem.wmf<-wmf(rasyephem_dt$cdat,times=years)

png('rasy_wmf.png',res=500,height=4,width=5,units="in")
plotmag(rasy.wmf,title="all")
plotmag(rasyperm.wmf,title="perm")
plotmag(rasyephem.wmf,title="ephem")

dev.off()

#Coherence test
#to get lag: delay = phase / 2π × period

#load climate values from climind (rasy_regression.R)
#turn into a list of repeated entry matrices

climlist <- lapply(climind, function(x) {
  matrix(x, nrow = 22, ncol = 21, byrow = TRUE)
})

#february coherences (2006-2025)
febtmin<-cleandat(climlist$feb_tmin[,1:20],clev=5,times=2006:2025)
febprcp<-cleandat(climlist$feb_prcp[,1:20],clev=5,times=2006:2025)
febtmax<-cleandat(climlist$feb_tmax[,1:20],clev=5,times=2006:2025)
rasy_dt1<-cleandat(as.matrix(rasy_compmat[,1:20]),clev=5,times=2006:2025) #cuts off 2026 data for SPEI

res<-coh(dat1=febtmin$cdat,dat2=rasy_dt1$cdat,times=years1,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3)) #significant, but barely- phase would indicate almost a 1 year lag based on phase
res<-bandtest(res,c(4,6))
get_bandp(res)
plotphase(res)
plotrank(res)
plotmag(wt(febtmin$cdat[1,],times = 2006:2025)) #would really say it looks like the wmf of count data

res<-coh(dat1=febprcp$cdat,dat2=rasy_dt1$cdat,times=years1,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3))  #significant, phase indicates 0.75-1 yr lag
res<-bandtest(res,c(4,6))
get_bandp(res)
plotmag(wt(febprcp$cdat[1,],times = 2006:2025)) #wouldn't really say it looks like the wmf of count data

res<-coh(dat1=febtmax$cdat,dat2=rasy_dt1$cdat,times=years1,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3))  #marginal significance but phase would be 0.75-1 year lag
res<-bandtest(res,c(4,6))  
get_bandp(res)
plotmag(wt(febtmax$cdat[1,],times = 2006:2025)) #wouldn;t really say it looks like the wmf of count data

#winter effects
wintmin<-cleandat(climlist$win_tmin[,1:20],clev=5,times=2006:2025)
winprcp<-cleandat(climlist$win_prcp[,1:20],clev=5,times=2006:2025)
wintmax<-cleandat(climlist$win_tmax[,1:20],clev=5,times=2006:2025)
winthaw<-cleandat(climlist$winter_thaw[,1:20],clev=5,times=2006:2025)
winenso<-cleandat(climlist$winenso_avg,clev=5,times=2006:2026)

res<-coh(dat1=wintmin$cdat,dat2=rasy_dt1$cdat,times=years1,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3)) 
res<-bandtest(res,c(4,6))
get_bandp(res) #not significant

res<-coh(dat1=winprcp$cdat,dat2=rasy_dt1$cdat,times=years1,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3)) #marginal significance
res<-bandtest(res,c(4,6))
get_bandp(res)

res<-coh(dat1=wintmax$cdat,dat2=rasy_dt1$cdat,times=years1,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3)) 
res<-bandtest(res,c(4,6))
get_bandp(res) #not significant

res<-coh(dat1=winthaw$cdat,dat2=rasy_dt1$cdat,times=years1,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3)) 
res<-bandtest(res,c(4,6))
get_bandp(res) #not significant

res<-coh(dat1=winenso$cdat,dat2=rasy_dt$cdat,times=years,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3)) 
res<-bandtest(res,c(4,6)) #p = 0.11
get_bandp(res)
plotmag(wt(winenso$cdat[1,],times = 2006:2026)) #would really say it looks like the wmf of count data


#two year lagged effects
ensosumlag2<-cleandat(climlist$enso_lag2,clev=5,times=2006:2026)
speisumlag2<-cleandat(climlist$spei_lag2,clev=5,times=2006:2026)
ensowinlag1<-cleandat(climlist$winenso_lag,clev=5,times=2006:2026)
ensowin<-cleandat(climlist$winenso_avg,clev=5,times=2006:2026)
pdsisumlag2<-cleandat(climlist$pdsi_lag2,clev=5,times=2006:2026)
pdsisumlag1<-cleandat(climlist$pdsi_lag1,clev=5,times=2006:2026)

res<-coh(dat1=ensosumlag2$cdat,dat2=rasy_dt$cdat,times=years,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3)) #marginal significance
res<-bandtest(res,c(4,6))
get_bandp(res)

res<-coh(dat1=ensowinlag1$cdat,dat2=rasy_dt$cdat,times=years,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3)) 
res<-bandtest(res,c(4,6)) #significant, phase = 2.22; delay = 1.4-2.1 yrs but already a lag so 3.4-4.1
get_bandp(res)
plotphase(res)
phasemean(res = res,tsrange=c(4,6))
plotmag(wt(ensowinlag1$cdat[1,],times = 2006:2026)) #would really say it looks like the wmf of count data

res<-coh(dat1=ensowin$cdat,dat2=rasy_dt$cdat,times=years,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)
res<-bandtest(res,c(2,3)) 
res<-bandtest(res,c(4,6)) #not significant? should it be but with just a longer lag period in the phase? Or lack of power maybe?

res<-coh(dat1=pdsisumlag2$cdat,dat2=rasy_dt$cdat,times=years,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3)) #significant,in phase
res<-bandtest(res,c(4,6))
get_bandp(res)
phasemean(res,tsrange=c(2,3))
plotphase(res)
plotmag(wt(pdsisumlag2$cdat[1,],times = 2004:2024)) #wavelet power of summer pdsi, lagged two years

res<-coh(dat1=speisumlag2$cdat,dat2=rasy_dt$cdat,times=years,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3)) 
res<-bandtest(res,c(4,6))
get_bandp(res) #not significant


#build neighborhood competition model and calculate growth rate (as in Rowland paper)
library(tidyverse)
library(sf)

max_dist <- 500   # meters; theta = 1 / max_dist

# 1. Pond coordinates and pairwise distances (meters) -----------------------
ponds <- rasy_comp |> distinct(pondname, Latitude, Longitude)

pts <- st_as_sf(ponds, coords = c("Longitude", "Latitude"), crs = 4326)
d <- st_distance(pts)                        # geodesic distances in meters
d <- matrix(as.numeric(d), nrow = nrow(ponds),
            dimnames = list(ponds$pondname, ponds$pondname))

# 2. Weight matrix: exponential decay, zero beyond 500 m, zero on diagonal ---
W <- exp(-d / max_dist)
W[d > max_dist] <- 0
diag(W) <- 0

# 3. Complete pond x year grid so lags are correct --------------------------
years <- seq(min(rasy_long$Year), max(rasy_long$Year))

full <- rasy_long |>
  filter(pondname%in%unique(rasy_comp$pondname))|>
  select(pondname, Year, Count) |>
  complete(pondname = ponds$pondname, Year = years) |>
  arrange(pondname, Year)

# pond x year count matrix (rows in same order as W)
C <- full |>
  pivot_wider(names_from = Year, values_from = Count) |>
  slice(match(ponds$pondname, pondname))
C <- as.matrix(select(C, -pondname))
rownames(C) <- ponds$pondname

C0 <- C
C0[is.na(C0)] <- 0                            # unsurveyed ponds contribute nothing

# 4. Neighborhood sum for every pond-year: sum_j exp(-d_ij/500) * count_j ---
N <- W %*% C0

nbr <- as_tibble(N, rownames = "pondname") |>
  pivot_longer(-pondname, names_to = "Year", values_to = "neighbor_sum") |>
  mutate(Year = as.integer(Year))

nbr_wide<-nbr|>
  pivot_wider(names_from=Year,values_from = neighbor_sum)
nbr_mat<-as.matrix(nbr_wide|>
                     select(!pondname))
nbr_dt<-cleandat(nbr_mat,times=2006:2026,clev=1)#clev 3 and 5 don't work because one time series is all 0

plotmag(wpmf(nbr_dt$cdat,times=2006:2026,sigmethod = "quick"),sigthresh = 0.99)
nbr_res<-coh(dat1=nbr_dt$cdat,dat2=rasy_dt$cdat,times=years,norm="powall",
             sigmethod="fast",nrand=10000,
             f0=1)

nbr_res<-bandtest(nbr_res,c(4,6))#marginally significant, close to in phase relationship (-0.37)
nbr_res<-bandtest(nbr_res,c(2,3))
get_bandp(nbr_res)
phasemean(nbr_res,tsrange=c(4,6)) #phase means counts lead competition by 1/3 of year?



#load list of ENSO values- windows taking the average ENSO over the time period (from dbca_wx.R)
# current list is avg ENSO in the previous 1,2,6,12 and 24 months, using march as the anchor month

ensomat_list<-readRDS('data/enso.rds')
enso_dt_list<-lapply(ensomat_list,function(x){cleandat(x,years,clev=5)})

#plot enso wavelet transform
ensopost_wmf<-wt(enso_dt_list$postbreed$cdat[1,],years) 
ensoprebreed_wmf<-wt(enso_dt_list$prebreed6$cdat[1,],years) 

par(mfrow=c(1,2))
plotmag(ensopost_wmf,title="ENSO post breeding") #looks like it should have same periodicity (4-6 years)
plotmag(ensoprebreed_wmf,title="ENSO pre breeding") #looks like it should have same periodicity (4-6 years)

#test counts with enso during egg/tadpole development through juveniles entering hiberation (sept)
res<-coh(dat1=enso_dt_list$postbreed$cdat,dat2=rasy_dt$cdat,times=years,norm="powall",
         sigmethod="fast",nrand=10000,
         f0=1,scale.max.input=28)

res<-bandtest(res,c(2,3)) 
res<-bandtest(res,c(4,6))
get_bandp(res)

#test counts with enso during the previous 6 months
res<-coh(dat1=enso_dt_list$prebreed6$cdat,dat2=rasy_dt$cdat,times=years,norm="powall",
         sigmethod="fast",nrand=10000,
         f0=1,scale.max.input=28)

res<-bandtest(res,c(2,3))
res<-bandtest(res,c(4,6))
get_bandp(res)

#test counts with enso during the previous 12 months
res<-coh(dat1=enso_dt_list$prebreed12$cdat,dat2=rasy_dt$cdat,times=years,norm="powall",
         sigmethod="fast",nrand=10000,
         f0=1,scale.max.input=28)

res<-bandtest(res,c(2,3))
res<-bandtest(res,c(4,6))
get_bandp(res)

#test counts with winter temperatures from the previous two months (Jan & Feb of current year)
#load NOAA weather data (from dbca_wx.R)
noaa_wx<-readRDS('data/noaa_wx.rds')
tmin_dt<-cleandat(noaa_wx$tmin,years,clev=5)
res<-coh(dat1=tmin_dt$cdat,dat2=rasy_dt$cdat,times=years,norm="powall",
         sigmethod="fast",nrand=10000,
         f0=1,scale.max.input=28)

res<-bandtest(res,c(3,5))
res<-bandtest(res,c(4,6))
get_bandp(res)

tmax_dt<-cleandat(noaa_wx$tmax,years,clev=5)
res<-coh(dat1=tmax_dt$cdat,dat2=rasy_dt$cdat,times=years,norm="powall",
         sigmethod="fast",nrand=10000,
         f0=1,scale.max.input=28)

res<-bandtest(res,c(3,5))
res<-bandtest(res,c(4,6))
get_bandp(res)

prcp_dt<-cleandat(noaa_wx$prcp,years,clev=5)
res<-coh(dat1=prcp_dt$cdat,dat2=rasy_dt$cdat,times=years,norm="powall",
         sigmethod="fast",nrand=10000,
         f0=1,scale.max.input=28)

res<-bandtest(res,c(3,5))
res<-bandtest(res,c(4,6))
get_bandp(res)

#tests with SPEI, cutting off 2026 count data as Daymet/SPEI isn't available yet
spei_list<-readRDS('data/spei_list.rds')

rasy_dt1<-cleandat(as.matrix(rasy_compmat[,1:20]),clev=5,times=2006:2025) #cuts off 2026 data for SPEI
rasy_permdt1<-cleandat(as.matrix(rasy_matperm[,1:20]),clev=5,times=2006:2025) #cuts off 2026 data for SPEI
rasy_ephemdt1<-cleandat(as.matrix(rasy_matephem[,1:20]),clev=5,times=2006:2025) #cuts off 2026 data for SPEI

#complete time series, all hydroperiods
#test whether SPEI during months prior to breeding impact synchrony
spei3_dt<-cleandat(spei_list$prebreed3,2006:2025,clev=5)
res<-coh(dat1=spei3_dt$cdat,dat2=rasy_dt1$cdat,times=2006:2025,norm="powall",
         sigmethod="fast",nrand=10000,
         f0=1,scale.max.input=28)

res<-bandtest(res,c(2,3))
res<-bandtest(res,c(4,6))
get_bandp(res)

spei6_dt<-cleandat(spei_list$prebreed6,2006:2025,clev=5)
res<-coh(dat1=spei6_dt$cdat,dat2=rasy_dt1$cdat,times=2006:2025,norm="powall",
         sigmethod="fast",nrand=10000,
         f0=1,scale.max.input=28)

res<-bandtest(res,c(2,3))
res<-bandtest(res,c(4,6))
get_bandp(res)

spei12_dt<-cleandat(spei_list$breed12,2006:2025,clev=5)
res<-coh(dat1=spei12_dt$cdat,dat2=rasy_dt1$cdat,times=2006:2025,norm="powall",
         sigmethod="fast",nrand=10000,
         f0=1,scale.max.input=28)

res<-bandtest(res,c(2,3))
res<-bandtest(res,c(4,6))
get_bandp(res)

#SPEI during 7 months post breeding
spei7_dt<-cleandat(spei_list$postbreed7,2006:2025,clev=5)
res<-coh(dat1=spei7_dt$cdat,dat2=rasy_dt1$cdat,times=2006:2025,norm="powall",
         sigmethod="fast",nrand=10000,
         f0=1,scale.max.input=28)

res<-bandtest(res,c(2,3))
res<-bandtest(res,c(4,6))
get_bandp(res)

#complete time series, ephemeral hydroperiods
#test whether SPEI during months prior to breeding impact synchrony
res<-coh(dat1=spei3_dt$cdat[1:8,],dat2=rasy_ephemdt1$cdat,times=2006:2025,norm="powall",
         sigmethod="fast",nrand=10000,
         f0=1,scale.max.input=28)

res<-bandtest(res,c(2,3))
res<-bandtest(res,c(4,6))
get_bandp(res)

res<-coh(dat1=spei6_dt$cdat[1:8,],dat2=rasy_ephemdt1$cdat,times=2006:2025,norm="powall",
         sigmethod="fast",nrand=10000,
         f0=1,scale.max.input=28)

res<-bandtest(res,c(2,3))
res<-bandtest(res,c(4,6))
get_bandp(res)

res<-coh(dat1=spei12_dt$cdat[1:8,],dat2=rasy_ephemdt1$cdat,times=2006:2025,norm="powall",
         sigmethod="fast",nrand=10000,
         f0=1,scale.max.input=28)

res<-bandtest(res,c(2,3))
res<-bandtest(res,c(4,6))
get_bandp(res)

#SPEI during 7 months post breeding
res<-coh(dat1=spei7_dt$cdat[1:8,],dat2=rasy_ephemdt1$cdat,times=2006:2025,norm="powall",
         sigmethod="fast",nrand=10000,
         f0=1,scale.max.input=28)

res<-bandtest(res,c(2,3))
res<-bandtest(res,c(4,6))
get_bandp(res)

#complete time series, permanent hydroperiods
#test whether SPEI during months prior to breeding impact synchrony
res<-coh(dat1=spei3_dt$cdat[1:14,],dat2=rasy_permdt1$cdat,times=2006:2025,norm="powall",
         sigmethod="fast",nrand=10000,
         f0=1,scale.max.input=28)

res<-bandtest(res,c(2,3))
res<-bandtest(res,c(4,6))
get_bandp(res)

res<-coh(dat1=spei6_dt$cdat[1:14,],dat2=rasy_permdt1$cdat,times=2006:2025,norm="powall",
         sigmethod="fast",nrand=10000,
         f0=1,scale.max.input=28)

res<-bandtest(res,c(2,3))
res<-bandtest(res,c(4,6))
get_bandp(res)

res<-coh(dat1=spei12_dt$cdat[1:14,],dat2=rasy_permdt1$cdat,times=2006:2025,norm="powall",
         sigmethod="fast",nrand=10000,
         f0=1,scale.max.input=28)

res<-bandtest(res,c(2,3))
res<-bandtest(res,c(4,6))
get_bandp(res)

#SPEI during 7 months post breeding
res<-coh(dat1=spei7_dt$cdat[1:14,],dat2=rasy_permdt1$cdat,times=2006:2025,norm="powall",
         sigmethod="fast",nrand=10000,
         f0=1,scale.max.input=28)

res<-bandtest(res,c(2,3))
res<-bandtest(res,c(4,6))
get_bandp(res)

plotmag(wt(spei12_dt$cdat[1,],times=2006:2025))
plotmag(wt(spei6_dt$cdat[1,],times=2006:2025))
plotmag(wt(spei3_dt$cdat[1,],times=2006:2025))
plotmag(wt(spei7_dt$cdat[1,],times=2006:2025))

#Test significant effects of weather using Daymet data, cutting off 2026 as it's not available
daymet_wx<-readRDS('data/daymet_wx.rds')

#winter thaw from rowland- need to filter to only ponds with complete time series
daymet_wx<-lapply(daymet_wx,function(x){subset(x,site%in%rasy_comp$pondname)})
daymet_dt<-lapply(daymet_wx,function(x){cleandat(as.matrix(x[,2:21]),2006:2025,clev=5)})

res<-coh(dat1=daymet_dt$winthaw$cdat,dat2=rasy_dt1$cdat,times=2006:2025,norm="powall",
         sigmethod="fast",nrand=10000,
         f0=1,scale.max.input=10)

res<-bandtest(res,c(2,3))
res<-bandtest(res,c(4,6))
get_bandp(res)

par(mfrow=c(1,2))
plotmag(wmf(rasy_dt1$cdat,2006:2025),title="RASY")
plotmag(wt(daymet_dt$winthaw$cdat[1,],2006:2025),title="Winter Thaw") #seems decently similar to it?

# Wavelet coherence with regional abundance -------------------------------
#Wavelet coherence with regional abundance
totaleggs<-colSums(rasy_compmat)
totaleggs_perm<-colSums(rasy_compmat[rasy_comp$Hydroperiod=="perm",])
totaleggs_ephem<-colSums(rasy_compmat[rasy_comp$Hydroperiod=="ephem",])

totaleggs_dt<-cleandat(as.vector(totaleggs),clev=5,times=2006:2026)
totaleggs1_dt<-cleandat(as.vector(totaleggs[1:20]),clev=5,times=2006:2025)

totaleggs_permdt<-cleandat(as.vector(totaleggs_perm),clev=5,times=2006:2026)
totaleggs_ephemdt<-cleandat(as.vector(totaleggs_ephem),clev=5,times=2006:2026)

plot(2006:2026,totaleggs,type="b",xlab="Year",ylab="Total Eggs")
lines(2006:2026,totaleggs_ephem,type="b",col=c("green"))
lines(2006:2026,totaleggs_perm,type="b",col=c("blue"))

plotmag(wt(totaleggs_dt$cdat,times=years))
plotmag(wt(totaleggs_permdt$cdat,times=years))
plotmag(wt(totaleggs_ephemdt$cdat,times=years))

thaw_lag2<-cleandat(climind$winenso_lag,times=years,clev=5)

res1<-coh(dat1=thaw_lag2$cdat,dat2=totaleggs_dt$cdat,times=years,norm="powall",
          sigmethod="fast",nrand=10000,
          f0=1)
res1<-bandtest(res1,c(2,3))
res1<-bandtest(res1,c(4,5))
get_bandp(res1)

res2<-coh(dat1=tmax_dt$cdat[1,],dat2=totaleggs_dt$cdat,times=years,norm="powall",
          sigmethod="fast",nrand=10000,
          f0=1)
res2<-bandtest(res2,c(4,6))
res2<-bandtest(res2,c(5,6))
get_bandp(res2)
plotmag(wt(tmax_dt$cdat[1,],times=2006:2026))

totaleggs_dt<-cleandat(as.vector(totaleggs),clev=5,times=2006:2026)
res3<-coh(dat1=enso_dt_list$postbreed$cdat[1,],dat2=totaleggs_dt$cdat,times=years,norm="powall",
          sigmethod="fast",nrand=10000,
          f0=1,scale.max.input=10)
res3<-bandtest(res3,c(4,6))
res3<-bandtest(res3,c(2,3))
res3<-bandtest(res3,c(4,5))
res3<-bandtest(res3,c(2,6))
get_bandp(res3)

par(mfrow=c(1,2))
plotmag(wt(totaleggs_dt$cdat,2006:2026),title="RASY")
plotmag(wt(enso_dt_list$postbreed$cdat[1,],times=2006:2026))

#fit arma models of varying time lags using environmental predictors, compare sets up models
#examine frequncy information from the top models

#build neighborhood competition model
