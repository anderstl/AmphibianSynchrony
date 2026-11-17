library(lubridate)
library(tidyr)

#load NAO data
nao<-read.csv('data/nao.csv')
nao_long<-nao%>%
  rename(year = Year)%>%
  select(year:Dec)%>%
  pivot_longer(Jan:Dec,
               names_to= "Month",
               values_to="NAO")%>%
  mutate(
    MonthNum = match(toupper(Month), toupper(month.abb)),
    date = ymd(paste(year, MonthNum, 1, sep = "-"))
  ) %>%
  arrange(date)

#create average summer
nao_sum<-nao_long%>%
  filter(MonthNum%in%c(3:9))%>%
  group_by(year)%>%
  summarise(avg_nao=mean(NAO,na.rm=T))

#create lagged average summer
nao_sum <- nao_sum %>%
  arrange(year) %>%
  mutate(nao_lag2 = lag(avg_nao, n = 2),
         nao_lag1 = lag(avg_nao, n = 1))%>%
  filter(year>=1988)

#create average winter and one year lag in winter
nao_win<-nao_long%>%
  mutate(breed_year = if_else(MonthNum >= 10, year + 1L, year))%>%
  filter((MonthNum >= 12) | (MonthNum < 3)) %>%
  group_by(breed_year)%>%
  summarise(winnao_avg=mean(NAO,na.rm=T))%>%
  rename(year=breed_year)%>%
  arrange(year) %>%
  mutate(winnao_lag=lag(winnao_avg,n = 1))%>%
  filter(year>=1988)

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

#create average summer
enso_sum<-enso_long%>%
  filter(MonthNum%in%c(3:9))%>%
  group_by(year)%>%
  summarise(avg_enso=mean(SOI,na.rm=T))

#create lagged average summer
enso_sum <- enso_sum %>%
  arrange(year) %>%
  mutate(enso_lag2 = lag(avg_enso, n = 2),
         enso_lag1 = lag(avg_enso, n = 1))%>%
  filter(year>=1988)

#create average winter and one year lag in winter
enso_win<-enso_long%>%
  mutate(breed_year = if_else(MonthNum >= 10, year + 1L, year))%>%
  filter((MonthNum >= 12) | (MonthNum < 3)) %>%
  group_by(breed_year)%>%
  summarise(winenso_avg=mean(SOI,na.rm=T))%>%
  rename(year=breed_year)%>%
  arrange(year) %>%
  mutate(winenso_lag=lag(winenso_avg,n = 1))%>%
  filter(year>=1988)

#load palmer drought data
#pdsi codes: MO (2302); NY(3005?); RI (3701); OH(3310); CT (0602)
# url <- "https://www.ncei.noaa.gov/monitoring-content/data/us/climdiv/monthly/current/climdiv-pdsidv-v1.0.0-20261105"
# pdsi <- readr::read_table(url, col_names = c("code", month.abb),
#                           col_types = readr::cols(code = "c", .default = "d"))
# pdsi$state    <- substr(pdsi$code, 1, 2)
# pdsi$division <- substr(pdsi$code, 3, 4)
# pdsi$year     <- as.integer(substr(pdsi$code, 7, 10))

pdsi<-read.csv('data/pdsi.csv',colClasses = c(CODE = "character"))
pdsi <- pdsi %>%
  mutate(
    state_code    = as.integer(substr(CODE, 1, 2)),
    division      = as.integer(substr(CODE, 3, 4)),
    state_div    = as.character(substr(CODE, 1,4)),
    element_code  = as.integer(substr(CODE, 5, 6)),
    year          = as.integer(substr(CODE, 7, 10))
  )
head(pdsi)

pdsi_long<-pdsi%>%
  filter(state_div%in%c("2302","3005","3701","3310","0602"))%>%
  select(CODE,state_div,year,JAN,FEB,MAR,APR,MAY,JUN,JUL,AUG,SEP,OCT,NOV,DEC)%>%
  pivot_longer(JAN:DEC,values_to="PDSI",names_to = "Month")%>%
  mutate(Month=as.factor(Month))%>%
  mutate(MonthNum = match(toupper(Month), toupper(month.abb)))

pdsi_avg<-pdsi_long%>%
  filter(MonthNum%in%c(3:9) & PDSI>=-4)%>%
  group_by(year,state_div)%>%
  summarise(avg_pdsi=mean(PDSI,na.rm=T))

pdsi_avg <- pdsi_avg %>%
  group_by(state_div) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(pdsi_lag2 = lag(avg_pdsi, n = 2),
         pdsi_lag1 = lag(avg_pdsi, n = 1))%>%
  filter(year>=1988)

climind <- Reduce(function(x, y) merge(x, y, by = "year", all = TRUE), 
                  list(enso_win,enso_sum,
                       nao_sum,nao_win))
saveRDS(climind,"data/climind.rds")
saveRDS(pdsi_avg,"data/pdsi_avg.rds")
