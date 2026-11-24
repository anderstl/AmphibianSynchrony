###Author: Tom Anderson

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
plot_coverage(m, best$start, best$end, min_complete = 0.9, max_zero_prop = 0.5)
dat_use <- window_sites(m, best$start, best$end, min_complete = 0.9,
                        max_zero_prop = 0.5)
attr(dat_use, "site_stats")       # NAs and zeros per site in the window
attr(dat_use, "n_NA_remaining")   # values to impute

rules <- list(min_complete = 0.9, max_zero_prop = 0.5, min_nonzero = 8)
start <- 2006; end <- 2026          # your chosen window

sel <- do.call(window_sites, c(list(m, start, end), rules))

keep_sites <- rownames(sel)
keep_years <- start:end

#plot time series
eggs_use <- rasy|>
  select(pondname,X2006:X2026)|>
  pivot_longer(X2006:X2026,names_to="Year",values_to="Count")|>
  mutate(Year=as.numeric(sub('X', '', Year))) |>
  filter(as.character(pondname) %in% keep_sites,
         Year >= start, Year <= end)
ts_plot<-eggs_use%>%
  ggplot(aes(Year,Count,color=pondname))+
  geom_point(show.legend = F)+
  geom_line(aes(color=pondname),show.legend = F)+
  labs(y="Number of Egg Masses",subtitle="Complete Time Series")+
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=17))+
  scale_x_continuous(breaks=seq(2006,2026,1),labels=c("2006",rep("",3),"2010",rep("",3),"2014",rep("",3),"2018",rep("",3),"2022",rep("",3),"2026"))
ts_plot

#best option is 21 years for 21 sites (2001-2026)

rasy_comp <- sel %>%
  as.data.frame()%>%
  mutate(row_mean = round(rowMeans(select(., `2006`:`2026`), na.rm = TRUE),0),
         across(`2006`:`2026`, ~ coalesce(., row_mean)))%>%
  mutate(pondname=row.names(sel))%>%
  left_join(rasy%>%select(pondname,Latitude,Longitude,Hydroperiod),by="pondname")

rasy_compmat<-rasy_comp%>%
  select(`2006`:`2026`)%>%
  as.matrix()
#21 sites for 21 years which drops zero inflated ponds

#standardize data sets using Box-Cox transformation
rasy_dt<-cleandat(as.matrix(rasy_compmat),clev=5,times=keep_years)

#wavelet phasor mean fields- phase synchrony 
png('results/mo_wpmf.png',res=500,height=4,width=5,units="in")

#load custom plotting code
source('code/wmf_plotting.R')

rasy.wpmf<-wpmf(rasy_dt$cdat,times=keep_years,sigmethod = "quick")

#plotmag(rasy.wpmf,sigthresh = 0.99,las=1,title="all")
mo_wpmf<-plotwmf_gg(rasy.wpmf,sigthresh = 0.99)+
  labs(x="")+
  scale_x_continuous(breaks=seq(2006,2026,2),labels=c(2006,"",2010,"",2014,"",2018,"",2022,"",2026))+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))
mo_wpmf
dev.off()

#wavelet mean fields
rasy.wmf<-wmf(rasy_dt$cdat,times=keep_years)

png('rasy_wmf.png',res=500,height=4,width=5,units="in")
#plotmag(rasy.wmf,title="all")
mo_wmf<-plotwmf_gg(rasy.wmf)+
  labs(x="Year")+
  scale_x_continuous(breaks=seq(2006,2026,2),labels=c(2006,"",2010,"",2014,"",2018,"",2022,"",2026))+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))
mo_wmf
dev.off()

png('results/Fig3_wavs.png',res=600,width=4,height=7,units="in")
cowplot::plot_grid(mo_wpmf,mo_wmf,labels=c("A","B"),ncol=1)
dev.off()

#mean fields by hydroperiod
#convert to matrix
rasy_ephmat <- rasy_comp %>%
  filter(Hydroperiod=="ephem")%>%
  select(-c(pondname,Hydroperiod))
rasy_permmat <- rasy_comp %>%
  filter(Hydroperiod=="perm")%>%
  select(-c(pondname,Hydroperiod))
#standardize data sets using Box-Cox transformation
rasy_ephdt<-cleandat(as.matrix(rasy_ephmat),clev=5,times=keep_years)
rasy_permdt<-cleandat(as.matrix(rasy_permmat),clev=5,times=keep_years)

rasy.ephwpmf<-wpmf(rasy_ephdt$cdat,times=keep_years,sigmethod = "quick")
mo_ephwpmf<-plotwmf_gg(rasy.ephwpmf,sigthresh = 0.99)+
  labs(x="",subtitle = "Ephemeral")+
  scale_x_continuous(breaks=seq(2006,2026,2),labels=c(2006,"",2010,"",2014,"",2018,"",2022,"",2026))+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))
mo_ephwpmf
rasy.permwpmf<-wpmf(rasy_permdt$cdat,times=keep_years,sigmethod = "quick")
mo_permwpmf<-plotwmf_gg(rasy.permwpmf,sigthresh = 0.99)+
  labs(x="",subtitle = "Ephemeral")+
  scale_x_continuous(breaks=seq(2006,2026,2),labels=c(2006,"",2010,"",2014,"",2018,"",2022,"",2026))+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))
mo_permwpmf

rasy.ephwmf<-wmf(rasy_ephdt$cdat,times=keep_years)
mo_ephwmf<-plotwmf_gg(rasy.ephwmf)+
  labs(x="Year")+
  scale_x_continuous(breaks=seq(2006,2026,2),labels=c(2006,"",2010,"",2014,"",2018,"",2022,"",2026))+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))
mo_ephwmf
rasy.permwmf<-wmf(rasy_permdt$cdat,times=keep_years)
mo_permwmf<-plotwmf_gg(rasy.permwmf)+
  labs(x="Year")+
  scale_x_continuous(breaks=seq(2006,2026,2),labels=c(2006,"",2010,"",2014,"",2018,"",2022,"",2026))+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))
mo_permwmf

cowplot::plot_grid(mo_ephwpmf,mo_ephwmf,mo_permwpmf,mo_permwmf)


# Coherence Testing -------------------------------------------------------

#to get lag: delay = phase / 2π × period

#load 2-yr lagged climate and weather values
#turn into a list of repeated entry matrices
#allwx_lag2<-readRDS('data/allwx_lag2.rds') #daymet data- currently missing most recent year

#coherence with breeding conditions of the current year and 2 years prior
breed_wx<-readRDS('data/breedingwx_gridmet.rds')
mo_wx<-breed_wx%>%
  filter(year%in%keep_years & site=="mo")

breedlist <- lapply(mo_wx, function(x) {
  matrix(x, nrow = length(keep_sites), ncol = length(keep_years), byrow = TRUE)
})

#tmin
breedtmin<-cleandat(breedlist$tmin_mean,clev=5,times=keep_years)
plotmag(wt(breedtmin$cdat[1,],times = keep_years))

res<-coh(dat1=breedtmin$cdat,dat2=rasy_dt$cdat,times=keep_years,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,4))
res<-bandtest(res,c(4,6))
get_bandp(res)
plotphase(res)
plotrank(res)

#precip
breedprcp<-cleandat(breedlist$rain_total,clev=5,times=keep_years)
plotmag(wt(breedprcp$cdat[1,],times = keep_years))

res<-coh(dat1=breedprcp$cdat,dat2=rasy_dt$cdat,times=keep_years,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3))  
res<-bandtest(res,c(4,6))
get_bandp(res)

#mean temp
breedtmean<-cleandat(breedlist$tmean,clev=5,times=keep_years)
plotmag(wt(breedtmean$cdat[1,],times = keep_years))

res<-coh(dat1=breedtmean$cdat,dat2=rasy_dt$cdat,times=keep_years,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3))  
res<-bandtest(res,c(4,6))  
get_bandp(res)

#coherence with breeding conditions 2 years prior
breedtmin2<-cleandat(breedlist$tmin_mean2,clev=5,times=keep_years)
plotmag(wt(breedtmin2$cdat[1,],times = keep_years))

#tmin
res<-coh(dat1=breedtmin2$cdat,dat2=rasy_dt$cdat,times=keep_years,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3)) #significant, close to anti phase
res<-bandtest(res,c(4,6))
get_bandp(res)
plotphase(res)
plotrank(res) #not a strong signal

#precip
breedprcp2<-cleandat(breedlist$rain_total2,clev=5,times=keep_years)
plotmag(wt(breedprcp2$cdat[1,],times = keep_years))

res<-coh(dat1=breedprcp2$cdat,dat2=rasy_dt$cdat,times=keep_years,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3))  
res<-bandtest(res,c(4,6))
get_bandp(res)

#mean temp
breedtmean2<-cleandat(breedlist$tmean2,clev=5,times=keep_years)
plotmag(wt(breedtmean2$cdat[1,],times = keep_years))

res<-coh(dat1=breedtmean2$cdat,dat2=rasy_dt$cdat,times=keep_years,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3)) #significant  
res<-bandtest(res,c(4,6))  
get_bandp(res)
plotphase(res)
plotrank(res) #not a strong signal

#coherence with winter conditions 2 years prior
winter_wx<-readRDS('data/winterwx_gridmet.rds')
mowinter_wx<-winter_wx%>%
  filter(breed_year%in%keep_years & site=="mo")

winterlist <- lapply(mowinter_wx, function(x) {
  matrix(x, nrow = length(keep_sites), ncol = length(keep_years), byrow = TRUE)
})

#freeze thaw cycles
ftcycles<-cleandat(winterlist$ft_cycles,clev=5,times=keep_years)
plotmag(wt(ftcycles$cdat[1,],times = keep_years))

res<-coh(dat1=ftcycles$cdat,dat2=rasy_dt$cdat,times=keep_years,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3)) 
res<-bandtest(res,c(4,6))
get_bandp(res) 

#freezing degree days
fdd<-cleandat(winterlist$fdd,clev=5,times=keep_years)
plotmag(wt(fdd$cdat[1,],times = keep_years))

res<-coh(dat1=fdd$cdat,dat2=rasy_dt$cdat,times=keep_years,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3)) 
res<-bandtest(res,c(4,6))
get_bandp(res) 

#mean winter temp
tmean_djf<-cleandat(winterlist$tmean_djf,clev=5,times=keep_years)
plotmag(wt(tmean_djf$cdat[1,],times = keep_years))

res<-coh(dat1=tmean_djf$cdat,dat2=rasy_dt$cdat,times=keep_years,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3)) 
res<-bandtest(res,c(4,6))
get_bandp(res) 

#freeze thaw cycles 1 year lag
ftcycles1<-cleandat(winterlist$ft_cycles1,clev=5,times=keep_years)
plotmag(wt(ftcycles1$cdat[1,],times = keep_years))

res<-coh(dat1=ftcycles1$cdat,dat2=rasy_dt$cdat,times=keep_years,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3)) 
res<-bandtest(res,c(4,6))
get_bandp(res) 

#freezing degree days 1 year lag
fdd1<-cleandat(winterlist$fdd1,clev=5,times=keep_years)
plotmag(wt(fdd1$cdat[1,],times = keep_years))

res<-coh(dat1=fdd1$cdat,dat2=rasy_dt$cdat,times=keep_years,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3)) 
res<-bandtest(res,c(4,6))
get_bandp(res) 

#mean winter temp 1 year lag
tmean_djf1<-cleandat(winterlist$tmean_djf1,clev=5,times=keep_years)
plotmag(wt(tmean_djf1$cdat[1,],times = keep_years))

res<-coh(dat1=tmean_djf1$cdat,dat2=rasy_dt$cdat,times=keep_years,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3)) 
res<-bandtest(res,c(4,6))
get_bandp(res) 

#freeze thaw cycles 1 year lag
ftcycles2<-cleandat(winterlist$ft_cycles2,clev=5,times=keep_years)
plotmag(wt(ftcycles2$cdat[1,],times = keep_years))

res<-coh(dat1=ftcycles2$cdat,dat2=rasy_dt$cdat,times=keep_years,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3)) 
res<-bandtest(res,c(4,6))
get_bandp(res) 

#freezing degree days 1 year lag
fdd2<-cleandat(winterlist$fdd2,clev=5,times=keep_years)
plotmag(wt(fdd2$cdat[1,],times = keep_years))

res<-coh(dat1=fdd2$cdat,dat2=rasy_dt$cdat,times=keep_years,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3)) 
res<-bandtest(res,c(4,6))
get_bandp(res) 

#mean winter temp 2 year lag
tmean_djf2<-cleandat(winterlist$tmean2,clev=5,times=keep_years)
plotmag(wt(tmean_djf2$cdat[1,],times = keep_years))

res<-coh(dat1=tmean_djf2$cdat,dat2=rasy_dt$cdat,times=keep_years,
         norm="powall",
         sigmethod="fast",
         nrand=10000,
         f0=1,
         scale.max.input=28)

res<-bandtest(res,c(2,3)) 
res<-bandtest(res,c(4,6))
get_bandp(res) 

#coherence with drought and ENSO
ensosumlag2<-cleandat(climlist$enso_lag2,clev=5,times=2006:2026)
ensowin<-cleandat(climlist$winenso_avg,clev=5,times=2006:2026)
ensowinlag1<-cleandat(climlist$winenso_lag,clev=5,times=2006:2026)
ensowinlag2<-cleandat(climlist$winenso_lag,clev=5,times=2006:2026)

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

