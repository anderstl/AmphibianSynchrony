library(climateR) # to get climate date
library(sf) #to make spatial points
library(purrr) #for data wrangling
library(memoise) #for caching data
library(cachem) #for caching dat
library(dplyr) #for data wrangling
library(lubridate) #for date wrangling
library(ggplot2) #for plotting

#examine what variable are availabe with climateR
catalog |>
  filter(id == "gridmet") |>
  select(variable, description, units) |>
  distinct()

#load data points
gps<-read.csv('data/coords.csv')
pts <- st_as_sf(gps, coords = c("Longitude", "Latitude"), crs = 4326)

#get data- only reloads everything if new year added. if new site added, only adds the new site
if (file.exists("wx_gridmet.rds")) {
  wx <- readRDS("wx_gridmet.rds")
} else {
  wx <- map(seq_len(nrow(pts)), function(i) {
    out <- get_gm(pts[i, ], c("tmmn", "tmmx", "pr"), "1987-12-01", "2026-09-30")
    out$site <- pts$site[i]
    out
  }) |> list_rbind()
  saveRDS(wx, "gridmet_sites.rds")
}

t_snow <- 0     # all snow at or below this daily mean (°C)
t_rain <- 2     # all rain at or above this

#load data and create winter averages
winter_wx <- wx |>
  mutate(tmin = tmmn - 273.15, tmax = tmmx - 273.15, #subract because temp data units are Kelvin
         tmean = (tmin + tmax) / 2,
         snow_frac = pmin(1, pmax(0, (t_rain - tmean) / (t_rain - t_snow))),
         snow = pr * snow_frac,
         rain = pr - snow,
         breed_year = if_else(month(date) == 12, year(date) + 1, year(date))) |>
  filter(month(date) %in% c(12, 1, 2)) |>
  group_by(site, breed_year) |>
  summarise(ft_cycles = sum(tmin < 0 & tmax > 0), #freeze-thaw cycles
            fdd       = sum(pmax(0, -tmean)), #freezing degree days
            tmean_djf = mean(tmean), #mean temp
            snow_total = sum(snow),
            rain_total = sum(rain),
            .groups = "drop")%>%
  left_join(gps,by="site")%>%
  mutate(ft_cycles1=lag(ft_cycles, n = 1),
         fdd1=lag(fdd, n = 1),
         tmean1 = lag(tmean_djf,n=1),
         snow_total1=lag(snow_total, n =1),
         rain_total1 = lag(rain_total, n =1),
         ft_cycles2=lag(ft_cycles, n = 2),
         fdd2=lag(fdd, n = 2),
         tmean2 = lag(tmean_djf,n=2),
         snow_total2=lag(snow_total, n =2),
         rain_total2 = lag(rain_total, n =2))
saveRDS(winter_wx,'data/winterwx_gridmet.rds')

ggplot(winter,aes(site,ft_cycles))+geom_boxplot()
ggplot(winter,aes(site,fdd))+geom_boxplot()
ggplot(winter,aes(site,tmean_djf))+geom_boxplot()

ggplot(winter,aes(Latitude,tmean_djf))+stat_summary()
ggplot(winter,aes(Latitude,ft_cycles))+stat_summary()
ggplot(winter,aes(Latitude,fdd))+stat_summary()

#breeding windows 
#ct = median day of april 6, first = march 15, last = apr 23 (arietta et al. 2020)
#mo = median day of march 9, first = feb 27, last = march 28 (unpublished)
#ri = "mid-March to mid April" so assigning value of april 1 as median (raithel et al. 2011)
#ny = median day of april 8, first = march, last = april 28 (unpublished)
#oh = 

window_days <- 28   # need to try 14 and 21 as a sensitivity check

# median breeding date by site (add both OH sites when available)
med <- tibble::tribble(
  ~site,   ~med_md,
  "ct",    "04-06",
  "mo",    "03-09",
  "ri",    "04-01",
  "ny",    "04-08"
  # "oh.bl", "MM-DD",
  # "oh.sh", "MM-DD"
)

breeding_wx <- wx |>
  mutate(date = as.Date(date),
         tmin = tmmn - 273.15, tmax = tmmx - 273.15,
         tmean = (tmin + tmax) / 2,
         snow_frac = pmin(1, pmax(0, (t_rain - tmean) / (t_rain - t_snow))),
         snow = pr * snow_frac,
         rain = pr - snow,
         year = year(date)) |>
  inner_join(med, by = "site") |>
  mutate(med_date = as.Date(paste(year, med_md, sep = "-")),
         rel_day  = as.numeric(date - med_date)) |>
  filter(rel_day > -window_days, rel_day <= 0) |>
  group_by(site, year) |>
  summarise(tmin_mean  = mean(tmin),
            tmax_mean = mean(tmax),
            tmean      = mean(tmean),
            snow_total = sum(snow),
            rain_total = sum(rain),
            .groups = "drop") |>
  mutate(tmean2 = lag(tmean,n=2),
         tmax_mean2 = lag(tmax_mean,n=2),
         tmin_mean2 = lag(tmin_mean,n=2),
         snow_total2 = lag(snow_total,n=2),
         rain_total2 = lag(rain_total,n=2))|>
  left_join(gps, by = "site")
saveRDS(breeding_wx,'data/breedingwx_gridmet.rds')




