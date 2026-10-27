###Author: Tom Anderson
###Last Updated: September 20, 2026
library(ggplot2)
library(wsyn)
library(plyr)
library(dplyr)
library(tidyr)
library(lubridate)
library(readxl)

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
  select(pondname,Latitude,Longitude,X2006:X2026)%>%
  pivot_longer(X2006:X2026,names_to="Year",values_to="Count")
rasy_long$Year<-as.numeric(sub('X', '', rasy_long$Year))

#png("rasy.timeseries.png",height=5,width=5,res=500,units="in")
#pdf('woodfrog_synchrony.pdf')
ggplot(rasy_long,aes(Year,Count,color=pondname))+
  geom_point(show.legend = F)+
  geom_line(aes(color=pondname),show.legend = F)+
  labs(y="Number of Egg Masses")+
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=17))+
  scale_x_continuous(breaks=seq(2006,2026,1),labels=c("2006",rep("",3),"2010",rep("",3),"2014",rep("",3),"2018",rep("",3),"2022",rep("",3),"2026"))
#dev.off()

#Regression analysis similar to Rowland et al.
#drop ponds with all zeros or nearly all zeros 
zero_ponds <- rasy %>%
  select(pondname, X2006:X2026) %>%
  rowwise() %>%
  mutate(
    n_years    = sum(!is.na(c_across(X2006:X2026))),
    n_zero     = sum(c_across(X2006:X2026) == 0, na.rm = TRUE),
    prop_zero  = n_zero / n_years
  ) %>%
  ungroup() %>%
  filter(prop_zero >= 0.9)   # adjust threshold: 1.0 = all zero, 0.9 = "nearly all"

zero_ponds %>% select(pondname, prop_zero) #7 ponds to be dropped

rasy_clean <- rasy %>%
  filter(!pondname %in% zero_ponds$pondname)%>%
  select(pondname,Latitude,Longitude,X2006:X2026)%>%
  pivot_longer(X2006:X2026,names_to="Year",values_to="Count")%>%
  mutate(Year=as.numeric(sub('X', '', Year)))

#load NOAA daily climate data
warren_wx<-read.csv('data/warren_county_wx_all.csv')

#parse dates and average across weather station in the county
# dailyavg<-warren_wx%>%
#   mutate(Date=mdy(DATE),
#          doy=yday(Date),
#          Year=year(Date),
#          Month=month(Date))%>%
#   group_by(Date,doy,Year,Month)%>%
#   summarise(tmax=mean(TMAX,na.rm=T),
#             tmin=mean(TMIN,na.rm=T),
#             prcp=mean(PRCP,na.rm=T))
daily_warren<-warren_wx%>%
  mutate(
    date = parse_date_time(DATE, orders = c("ymd", "mdy")),
    #date = ymd(DATE),
    yday = yday(date),
    year= year(date),
    #date = as.Date(yday - 1, origin = paste0(year, "-01-01")),
    month = lubridate::month(date),
    day   = lubridate::day(date),
    # Oct-Dec counts toward the *following* year's breeding season
    breeding_year = if_else(month >= 10, year + 1L, year),
    snow = SNOW/10,
    prcp  = PRCP/10,
    tmax  = TMAX/10,
    tmin  = TMIN /10
  )

feb_warren<-daily_warren%>%
  filter(yday%in%c(45:59) ) %>%
  group_by(breeding_year)%>%
  summarise(feb_tmax=mean(tmax,na.rm=T),
            feb_tmin=mean(tmin,na.rm=T),
            feb_prcp=sum(prcp,na.rm=T),
            feb_snow=sum(snow,na.rm=T))%>%
  rename(year = breeding_year)%>%
  filter(year>=2006)

#calculate winter (Jan-Feb) averages
win_avg<-dailyavg%>%
  filter(doy<60)%>%
  group_by(Year)%>%
  summarise(tmax=mean(tmax,na.rm=T),
            tmin=mean(tmin,na.rm=T),
            prcp=mean(prcp,na.rm=T))%>%
  

#download daymet data
library(daymetr)
dm <- download_daymet(
  site     = "DBCA",
  lat      = mean(rasy$Latitude),
  lon      = mean(rasy$Longitude) ,
  start    = 1980, # get all the data for examinig long term trends
  end      = 2025,
  #internal = TRUE    # keep in memory rather than writing to disk
  #simplify = TRUE     # return a flat data frame instead of a nested list
)$data

daymet_dbca <- dm %>%
  mutate(
    date = as.Date(yday - 1, origin = paste0(year, "-01-01")),
    month = lubridate::month(date),
    day   = lubridate::day(date),
    # Oct-Dec counts toward the *following* year's breeding season
    breeding_year = if_else(month >= 10, year + 1L, year)
  )

#weather during preceding feb t 
daym_feb_wx<-daymet_dbca%>%
  filter(yday%in%c(45:59) ) %>%
  group_by(breeding_year)%>%
  summarise(feb_tmax=mean(tmax..deg.c.),
            feb_tmin=mean(tmin..deg.c.),
            feb_prcp=sum(prcp..mm.day.),
            feb_snow=sum(swe..kg.m.2.))%>%
  rename(year = breeding_year)%>%
  filter(year>=2006)

#weather during preceding winter months (dec t-1 to feb t), and the previous winter 
daym_win_wx<-daymet_dbca%>%
  filter((month >= 12) | (month == 2 & day <= 28)) %>%
  group_by(breeding_year)%>%
  summarise(win_tmax=mean(tmax..deg.c.),
            win_tmin=mean(tmin..deg.c.),
            win_prcp=sum(prcp..mm.day.),
            win_snow=sum(swe..kg.m.2.))%>%
  rename(year = breeding_year)%>%
  mutate(win_tmax_lag2=lag(win_tmax, n = 2),
         win_tmin_lag2=lag(win_tmin, n = 2),
         win_prcp_lag2=lag(win_prcp, n =2),
         win_snow_lag2 = lag(win_snow, n =2))
daym_win_wx_lag<-daym_win_wx%>%
  select(year,win_tmax_lag2,win_tmin_lag2,win_prcp_lag2,win_snow_lag2)%>%
  filter(year%in%c(2006:2026))
daym_win_wx_cur<-daym_win_wx%>%
  select(year,win_tmax,win_tmin,win_prcp,win_snow)%>%
  filter(year%in%c(2006:2025)) # not climate data yet for 2026
  
#create winter thaw variable like in rowland paper but for our date windows of prebreeding
daym_win_thaw <- daymet_dbca %>%
  filter((month >= 12) | (month == 2 & day <= 28)) %>%
  group_by(breeding_year) %>%
  summarise(winter_thaw = sum(tmax..deg.c. > 0, na.rm = TRUE), .groups = "drop") %>%
  rename(year = breeding_year)%>%
  mutate(winter_thaw_lag2=lag(winter_thaw, n = 2),
         winter_thaw_lag1=lag(winter_thaw, n = 1))%>%
  filter(year%in%c(2006:2025))

#calculate SPEI
library(SPEI)
daily_daymet <- dm %>%
  mutate(
    date  = as.Date(paste(year, yday), format = "%Y %j"),
    month = month(date),
    tmean = (tmax..deg.c. + tmin..deg.c.) / 2,
    # Convert srad from W/m² to MJ/m²/day (needed by hargreaves())
    # dayl is in seconds; 1 W = 1 J/s
    Ra_MJ = srad..W.m.2. * dayl..s. / 1e6
  )

#aggregate to monthly
monthly_daymet <- daily_daymet %>%
  group_by(year, month) %>%
  summarise(
    prcp_mm  = sum(prcp..mm.day., na.rm = TRUE),  # total monthly precip (mm)
    tmin_mean = mean(tmin..deg.c., na.rm = TRUE),
    tmax_mean = mean(tmax..deg.c., na.rm = TRUE),
    tmean_mean = mean(tmean, na.rm = TRUE),
    Ra_mean   = mean(Ra_MJ, na.rm = TRUE),        # mean daily Ra for the month
    .groups = "drop"
  ) %>%
  arrange(year, month)

#compute PET using hargreaves method
monthly_daymet$PET_mm <- hargreaves(
  Tmin = monthly_daymet$tmin_mean,
  Tmax = monthly_daymet$tmax_mean,
  Ra   = monthly_daymet$Ra_mean,
  lat  = mean(rasy$Latitude)              # used internally only when Ra = NA; kept for reference
)

# calculate water balance  (Precip - PET)
monthly_daymet$BAL <- monthly_daymet$prcp_mm - monthly_daymet$PET_mm

# build monthly time series;ts() requires knowing the start year/month
start_ym <- c(monthly_daymet$year[1], monthly_daymet$month[1])

bal_daymet_ts <- ts(monthly_daymet$BAL, start = start_ym, frequency = 12)

# Compute SPEI at chosen scale
spei_out <- spei(bal_daymet_ts, scale = 7)

daymet_spei <- monthly_daymet %>%
  mutate(spei = as.numeric(spei_out$fitted)) %>%
  filter(month == 10)%>%
  mutate(spei_lag2=lag(spei,n=2))%>%
  filter(year%in%2006:2025)

#load ENSO data  
enso<-read.csv('data/enso.csv')
enso_long<-enso%>%
  rename(year = YEAR)%>%
  select(year:DEC)%>%
  pivot_longer(JAN:DEC,
               names_to= "Month",
               values_to="SOI")%>%
  mutate(
    MonthNum = match(toupper(Month), toupper(month.abb)),
    date = ymd(paste(year, MonthNum, 1, sep = "-"))
  ) %>%
  arrange(date)

enso_sum<-enso_long%>%
  filter(MonthNum%in%c(3:9))%>%
  group_by(year)%>%
  summarise(avg_enso=mean(SOI,na.rm=T))

enso_sum <- enso_sum %>%
  arrange(year) %>%
  mutate(enso_lag2 = lag(avg_enso, n = 2),
         enso_lag1 = lag(avg_enso, n = 1))%>%
  filter(year%in%c(2006:2026))

enso_win<-enso_long%>%
  mutate(breed_year = if_else(MonthNum >= 10, year + 1L, year))%>%
  filter((MonthNum >= 12) | (MonthNum < 3)) %>%
  group_by(breed_year)%>%
  summarise(winenso_avg=mean(SOI,na.rm=T))%>%
  rename(year=breed_year)%>%
  arrange(year) %>%
  mutate(winenso_lag=lag(winenso_avg,n = 1))%>%
  filter(year%in%c(2006:2026))

#load palmer drought data
# data from: https://www.ncei.noaa.gov/monitoring-content/data/us/climdiv/monthly/current/
pdsi<-read.csv('data/pdsi_climdiv2302.csv')
pdsi <- pdsi %>%
  mutate(
      state_code    = as.integer(substr(CODE, 1, 2)),
    division      = as.integer(substr(CODE, 3, 4)),
    element_code  = as.integer(substr(CODE, 5, 6)),
    year          = as.integer(substr(CODE, 7, 10))
  )
head(pdsi)

pdsi_long<-pdsi%>%
  select(year,JAN,FEB,MAR,APR,MAY,JUN,JUL,AUG,SEP,OCT,NOV,DEC)%>%
  pivot_longer(JAN:DEC,values_to="PDSI",names_to = "Month")%>%
  mutate(Month=as.factor(Month))%>%
  mutate(MonthNum = match(toupper(Month), toupper(month.abb)))

pdsi_avg<-pdsi_long%>%
  filter(MonthNum%in%c(3:9) & PDSI>=-4)%>%
  group_by(year)%>%
  summarise(avg_pdsi=mean(PDSI,na.rm=T))

pdsi_avg <- pdsi_avg %>%
  arrange(year) %>%
  mutate(pdsi_lag2 = lag(avg_pdsi, n = 2),
         pdsi_lag1 = lag(avg_pdsi, n = 1))%>%
  filter(year%in%c(2006:2026))

climind <- Reduce(function(x, y) merge(x, y, by = "year", all = TRUE), 
                    list(enso_win,enso_sum,
                         pdsi_avg,
                         daym_feb_wx,
                         daym_win_thaw,
                         daym_win_wx_cur,
                         daym_win_wx_lag,
                         daymet_spei[,c("year","spei","spei_lag2")]))

climind$spei_lag2[21]<-climind$spei[19]
head(climind)

#build neighborhood competition model and calculate growth rate (as in Rowland paper)
library(tidyverse)
library(sf)

max_dist <- 500   # meters; theta = 1 / max_dist

# 1. Pond coordinates and pairwise distances (meters) -----------------------
ponds <- rasy_clean |> distinct(pondname, Latitude, Longitude)

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

full <- rasy_clean |>
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

# 5. Lag by 2 years and compute growth rate ----------------------------------
rasy_clean2 <- full |>
  left_join(nbr, by = c("pondname", "Year")) |>
  group_by(pondname) |>
  arrange(Year, .by_group = TRUE) |>
  mutate(
    count_lag2    = lag(Count, 2),
    count_lag1    = lag(Count, 1),
    neighbor_lag2 = lag(neighbor_sum, 2),                 # neighborhood competition (t-2)
    neighbor_lag1 = lag(neighbor_sum, 1),                 # neighborhood competition (t-2)
    growth_rate   = log((Count + 1) / (count_lag2 + 1)) / 2
  ) |>
  ungroup()

alldat<-merge(rasy_clean2,climind,by.x="Year",by.y="year")
alldat<-merge(alldat,pond_info[,c("pondname","Area")])
alldat$Area[alldat$Area==0]<-mean(alldat$Area,na.rm=T)#impute missing areas
alldat$density_lag2<-alldat$count_lag2/alldat$Area 

# Regression analyses -----------------------------------------------------

#general trend analysis
library(glmmTMB)

nbmod<-glmmTMB(Count~scale(Year)+(scale(Year)|pondname),data=alldat,family=nbinom1,REML = T) #correlated random slopes and intercepts
nbmod1<-glmmTMB(Count~scale(Year)+(1|pondname),data=alldat,family=nbinom1,REML = T) #only site random effects
nbmod2<-glmmTMB(Count~scale(Year)+(1|pondname)+(0+scale(Year)|pondname),data=alldat,family=nbinom1,REML = T) #uncorrelated random slopes and intercepts

pmod<-glmmTMB(Count~scale(Year)+(scale(Year)|pondname),data=alldat,family=poisson,REML = T)
pmod1<-glmmTMB(Count~scale(Year)+(1|pondname),data=alldat,family=poisson,REML = T)
pmod2<-glmmTMB(Count~scale(Year)+(1|pondname)+(0+scale(Year)|pondname),data=alldat,family=poisson,REML = T)

AIC(nbmod,nbmod1,nbmod2,pmod,pmod1,pmod2)

nbmod2a<-update(nbmod2,REML=F) #best model
summary(nbmod2a) #trend not significant

png('yeartrends.png',res=600,height=3.5,width=3.5,units="in")
ggplot(alldat,aes(Year,Count))+
  geom_smooth(aes(color=pondname),size=0.5,se=F,method="lm",show.legend = F)+
  geom_smooth(color="black",size=2,se=F,method="lm",show.legend = F)+
  theme_classic()+
  labs(y="Number of Egg Masses")+
  scale_color_manual(values=gray.colors(n=length(unique(alldat$pondname))))
dev.off()

fixef(nbmod2)$cond[2]                       # β
VarCorr(nbmod2)                               # σ_b (SD of year_c by pond)
ps <- coef(nbmod2)$cond$pondname[,2]            # pond-specific slopes (β + b_p)
hist(ps); abline(v = fixef(nbmod2)$cond[2]  , col = "red")
  
#model evaluation
library(DHARMa)
nbfit<-simulateResiduals(nbmod2)
plot(nbfit)
testUniformity(nbfit)
testOverdispersion(nbfit)
testDispersion(nbfit)
testOutliers(nbfit)



#Run global models in glmmTMB with random intercept and uncorrelated random slopes
#SPEI and PDSI are correlated (r > 0.75) so only use one

#2-yr lag effects
cmod1<-glmmTMB(Count~scale(density_lag2)+scale(neighbor_lag2)+scale(pdsi_lag2)+scale(enso_lag1)+scale(enso_lag2)+(1|pondname),data=alldat,family=nbinom1)
summary(cmod1) #density, drought and enso significant
library(performance)
check_model(cmod1)

#winter effects
cmod2<-glmmTMB(Count~scale(win_tmin)+scale(win_prcp)+scale(win_snow)+(1|pondname),data=alldat,family=nbinom1)
summary(cmod2)
check_model(cmod2)

cmod3<-glmmTMB(Count~scale(feb_tmin)+scale(feb_prcp)+scale(feb_snow)+(1|pondname)+(0+Year|pondname),data=alldat,family=nbinom1)
summary(cmod3) #feb precip significant but negative- maybe spreads out the eggs more? 

#univariate models
cdensmod<-glmmTMB(Count~scale(density_lag2)+ (1 | pondname)+(0 + scale(density_lag2)|pondname),data=alldat,family="nbinom1")
censomod<-glmmTMB(Count~scale(enso_lag2)+ (1 | pondname)+(0 + scale(enso_lag2)|pondname),data=alldat,family="nbinom1")
cthawmod<-glmmTMB(Count~scale(winter_thaw)+ (1 | pondname) +(0+scale(winter_thaw)|pondname),data=alldat,family="nbinom1")
cpdsimod<-glmmTMB(Count~scale(pdsi_lag2)+ (1 | pondname) +(0 + scale(pdsi_lag2)|pondname),data=alldat,family="nbinom1")
cterrmod<-glmmTMB(Count~scale(neighbor_lag2)+ (1 | pondname)+(0+scale(neighbor_lag2)|pondname),data=alldat,family="nbinom1")
ctminmod<-glmmTMB(Count~scale(feb_tmin)+(1 | pondname) + (0 + scale(feb_tmin) | pondname),data=alldat,family="nbinom1")
cprcpmod<-glmmTMB(Count~scale(feb_prcp)+ (1 | pondname) + (0 + scale(feb_prcp) | pondname),data=alldat,family="nbinom1")
csnowmod<-glmmTMB(Count~scale(feb_snow)+ (1 | pondname) + (0 + scale(feb_snow) | pondname),data=alldat,family="nbinom1")
summary(cdensmod)
summary(cthawmod)
summary(censomod)
summary(cpdsimod)
summary(cterrmod)
summary(ctminmod)
summary(cprcpmod)
summary(csnowmod)

#population growth rate
gmod1<-glmmTMB(growth_rate~scale(density_lag2)+scale(neighbor_lag2)+scale(enso_lag2)+scale(enso_lag1)+scale(pdsi_lag2)+(1|pondname)+(0+Year|pondname),data=alldat)
summary(gmod1)
check_model(gmod1)

gmod2<-glmmTMB(growth_rate~scale(win_tmin)+scale(win_prcp)+scale(win_snow)+(1|pondname)+(0+Year|pondname),data=alldat)
summary(gmod2) #more snow, warmer temps = slower growth rates
check_model(gmod2)

gmod3<-glmmTMB(growth_rate~scale(feb_tmin)+scale(feb_prcp)+scale(feb_snow)+(1|pondname)+(0+Year|pondname),data=alldat)
summary(gmod3) #more snow, warmer temps = slower growth rates

#univariate growth models, controlling for pond level variability in responses
gdensmod<-glmmTMB(growth_rate~scale(density_lag2)+ (1 | pondname)+(0 + scale(density_lag2)|pondname),data=alldat)
gensomod<-glmmTMB(growth_rate~scale(enso_lag2)+ (1 | pondname)+(0 + scale(enso_lag2)|pondname),data=alldat)
genso1mod<-glmmTMB(growth_rate~scale(enso_lag1)+ (1 | pondname)+(0 + scale(enso_lag1)|pondname),data=alldat)
gthawmod<-glmmTMB(growth_rate~scale(winter_thaw)+ (1 | pondname) +(0+scale(winter_thaw)|pondname),data=alldat)
gpdsimod<-glmmTMB(growth_rate~scale(pdsi_lag2)+ (1 | pondname) +(0 + scale(pdsi_lag2)|pondname),data=alldat)
gterrmod<-glmmTMB(growth_rate~scale(neighbor_lag2)+ (1 | pondname)+(0+scale(neighbor_lag2)|pondname),data=alldat)
gtminmod<-glmmTMB(growth_rate~scale(feb_tmin)+(1 | pondname) + (0 + scale(feb_tmin) | pondname),data=alldat)
gprcpmod<-glmmTMB(growth_rate~scale(feb_prcp)+ (1 | pondname) + (0 + scale(feb_prcp) | pondname),data=alldat)
gsnowmod<-glmmTMB(growth_rate~scale(feb_snow)+ (1 | pondname) + (0 + scale(feb_snow) | pondname),data=alldat)
summary(gdensmod)
summary(gthawmod)
summary(gensomod)
summary(gpdsimod)
summary(gterrmod)
summary(gtminmod)
summary(gprcpmod)
summary(gsnowmod)

#run same models with log of counts and AR2 structure
library(nlme)
cmod1<-lme(log1p(Count)~scale(density_lag2)+scale(neighbor_lag2)+
             scale(pdsi_lag2)+scale(enso_lag2),
           random = ~1|pondname,
           data=alldat,
           na.action="na.omit",
           correlation = corARMA(form = ~ Year | pondname,p=2,q=0))
summary(cmod1)

gmod1<-lme(growth_rate~scale(density_lag2)+scale(neighbor_lag2)+
             scale(pdsi_lag2)+scale(enso_lag2),
           random = ~1|pondname,
           data=alldat,
           na.action="na.omit",
           correlation = corARMA(form = ~ Year | pondname,p=2,q=0))
summary(gmod1)

library(ggeffects)
pred_data <- ggpredict(cpdsimod, terms = "pdsi_lag2")

# 4. Plot using ggeffects' built-in ggplot wrapper
png('pdsi2_predicted.png',res=500,height=3.5,width=3.5,units="in")
plot(pred_data)+
  theme_classic() +
  labs(y="Predicted Count", x="Lagged Drought Severity",title = "")+
  lims(y=c(0,15)) 
dev.off()

pred_data <- ggpredict(cmod3, terms = "feb_prcp")
  plot(pred_data)+
  theme_classic() +
  labs(y="Predicted Count", x="Feb. Rain",title = "")+
  lims(y=c(0,15)) 

#amburgey covariates: precip during breeding season, pdsi or spei in summer that would affect recruitment
#enso effects are more likely in the winter in MO: https://ipm.missouri.edu/mpg/2026/8/el_nino-zl/

