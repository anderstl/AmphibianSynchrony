#load daymet weather for each study site for winter weather and SPEI
##download daymet data
Sys.setenv(CURL_SSL_BACKEND = "openssl")
library(daymetr)
library(tidyr)
library(dplyr)
#coords<-data.frame(site=c("ct","oh.bl","oh.sh","mo","ri","ny"),
#                   Latitude=c(41.95,39.94,40.11,38.78,41.59,41.77),
#                   Longitude=c(-72.14,-82.81,-82.96,-91.39,-71.72,-74.15)
#                   )
#write.csv(coords,'data/coords.csv',row.names=F)

dm <- download_daymet_batch("C:/Users/thander/Documents/GitRepos/AmphibianSynchrony/data/coords.csv",
  start    = 1980, # get all the data for examinig long term trends
  end      = 2025
  #internal = TRUE    # keep in memory rather than writing to disk
  #simplify = TRUE     # return a flat data frame instead of a nested list
)

wx_long <- bind_rows(
  lapply(dm, function(x) {
    df <- x$data
    df$site      <- x$site
    df$latitude  <- x$latitude
    df$longitude <- x$longitude
    df
  })
)

wx_long <- wx_long %>%
  mutate(
    date = as.Date(yday - 1, origin = paste0(year, "-01-01")),
    month = lubridate::month(date),
    day   = lubridate::day(date),
    # Oct-Dec counts toward the *following* year's breeding season
    breeding_year = if_else(month >= 10, year + 1L, year)
  )

win_wx<-wx_long%>%
  filter((month >= 12) | (month < 3) | (month == 3 & day <= 31)) %>%
  group_by(site,breeding_year)%>%
  summarise(tmax=mean(tmax..deg.c.),
            tmin=mean(tmin..deg.c.),
            prcp=mean(prcp..mm.day.),
            snow=mean(swe..kg.m.2.))%>%
  rename(year = breeding_year)

#weather during preceding feb 15-march 15 
breeding_wx<-wx_long%>%
  filter(yday%in%c(46:74)) %>%
  group_by(site,breeding_year)%>%
  summarise(breed_tmax=mean(tmax..deg.c.,na.rm=T),
            breed_tmin=mean(tmin..deg.c.,na.rm=T),
            breed_prcp=sum(prcp..mm.day.,na.rm=T),
            breed_snow=sum(swe..kg.m.2.,na.rm=T))%>%
  rename(year = breeding_year)%>%
  mutate(breed_tmax_lag2=lag(breed_tmax, n = 2),
         breed_tmin_lag2=lag(breed_tmin, n = 2),
         breed_prcp_lag2=lag(breed_prcp, n =2),
         breed_snow_lag2 = lag(breed_snow, n =2))%>%
  filter(year>=1988)

#weather during preceding winter months (dec t-1 to feb t), and the previous winter 
winter_wx<-wx_long%>%
  filter((month >= 12) | (month == 3 & day <= 31)) %>%
  group_by(site,breeding_year)%>%
  summarise(win_tmax=mean(tmax..deg.c.),
            win_tmin=mean(tmin..deg.c.),
            win_prcp=sum(prcp..mm.day.),
            win_snow=sum(swe..kg.m.2.))%>%
  rename(year = breeding_year)%>%
  mutate(win_tmax_lag2=lag(win_tmax, n = 2),
         win_tmin_lag2=lag(win_tmin, n = 2),
         win_prcp_lag2=lag(win_prcp, n =2),
         win_snow_lag2 = lag(win_snow, n =2))%>%
  filter(year>=1988)

#create winter thaw variable like in rowland paper but only for Dec-Feb
win_thaw <- wx_long %>%
  filter((month >= 12) | (month == 2 & day <= 28)) %>%
  group_by(site,breeding_year) %>%
  summarise(winter_thaw = sum(tmax..deg.c. > 0, na.rm = TRUE), .groups = "drop") %>%
  rename(year = breeding_year)%>%
  mutate(winter_thaw_lag2=lag(winter_thaw, n = 2))%>%
  filter(year>=1988)

winter_lag2<-winter_wx%>%
  select(site,year,win_tmax_lag2,win_tmin_lag2,win_prcp_lag2,win_snow_lag2)%>%
  left_join(win_thaw%>%select(site,year,winter_thaw_lag2),by=c('site',"year"))

allwx_lag2<-winter_lag2%>%
  left_join(breeding_wx%>%select(site,year,breed_tmax_lag2,breed_tmin_lag2,breed_prcp_lag2,breed_snow_lag2),by=c("site","year"))

saveRDS(allwx_lag2,"data/allwx_lag2.rds")

#calculate SPEI
#calculate SPEI
library(SPEI)
library(dplyr)
library(lubridate)

daily_daymet <- wx_long %>%
  mutate(
    date  = as.Date(paste(year, yday), format = "%Y %j"),
    month = month(date),
    tmean = (tmax..deg.c. + tmin..deg.c.) / 2,
    # Convert srad from W/m² to MJ/m²/day (needed by hargreaves())
    Ra_MJ = srad..W.m.2. * dayl..s. / 1e6
  )

# aggregate to monthly
monthly_daymet <- daily_daymet %>%
  group_by(site, year, month) %>%
  summarise(
    prcp_mm    = sum(prcp..mm.day., na.rm = TRUE),
    tmin_mean  = mean(tmin..deg.c., na.rm = TRUE),
    tmax_mean  = mean(tmax..deg.c., na.rm = TRUE),
    tmean_mean = mean(tmean, na.rm = TRUE),
    Ra_mean    = mean(Ra_MJ, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(site, year, month)          # <- site first, so each site's months are contiguous

# PET (Hargreaves) and water balance -- row-wise, so fine to do on the whole table
monthly_daymet <- monthly_daymet %>%
  mutate(
    PET_mm = as.numeric(hargreaves(Tmin = tmin_mean, Tmax = tmax_mean, Ra = Ra_mean)),
    BAL    = prcp_mm - PET_mm
  )

# SPEI per site: build a ts for each site with its own start year/month
spei_scale <- 7

monthly_spei <- monthly_daymet %>%
  group_by(site) %>%
  arrange(year, month, .by_group = TRUE) %>%
  group_modify(~ {
    start_ym <- c(.x$year[1], .x$month[1])
    bal_ts   <- ts(.x$BAL, start = start_ym, frequency = 12)
    spei_out <- spei(bal_ts, scale = spei_scale, verbose = FALSE)
    mutate(.x, spei = as.numeric(spei_out$fitted))
  }) %>%
  ungroup()

# October values, with a 2-year lag within each site
spei_lag2 <- monthly_spei %>%
  filter(month == 9) %>%
  group_by(site) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(spei_lag2 = lag(spei, n = 2)) %>%
  ungroup()%>%
  filter(year>=1988)%>%
  select(site,year,spei_lag2)
saveRDS(spei_lag2,"data/spei_lag2.rds")


