library(dplyr)

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
  select(pondname,Latitude,Longitude,X2006:X2026)%>%
  pivot_longer(X2006:X2026,names_to="Year",values_to="Count")
rasy_long$Year<-as.numeric(sub('X', '', rasy_long$Year))


# --- 1. Yearly summaries ---------------------------------------------------

yearly_summary <- rasy_long %>%
  filter(!is.na(Count)) %>%
  group_by(Year) %>%
  summarise(
    n_sites_total    = n_distinct(pondname),
    n_sites_occupied = n_distinct(pondname[Count > 0]),
    prop_occupied    = n_sites_occupied / n_sites_total,
    mean_count_occ   = mean(Count[Count > 0], na.rm = TRUE),   # egg masses, occupied sites only
    sd_count_occ     = sd(Count[Count > 0], na.rm = TRUE),
    mean_count_all   = mean(Count, na.rm = TRUE),              # includes zeros, if useful
    .groups = "drop"
  )

# --- 2. Site-level (population) summaries -----------------------------------

site_summary <- rasy_long %>%
  filter(!is.na(Count)) %>%
  group_by(pondname) %>%
  summarise(
    n_years_sampled  = n_distinct(Year),
    n_years_occupied = sum(Count > 0, na.rm = TRUE),
    prop_years_occ   = n_years_occupied / n_years_sampled,
    mean_count_occ   = mean(Count[Count > 0], na.rm = TRUE),   # avg egg masses when occupied
    sd_count_occ     = sd(Count[Count > 0], na.rm = TRUE),
    mean_count_all   = mean(Count, na.rm = TRUE),              # avg egg masses across all years (incl. 0s)
    .groups = "drop"
  )

mean(site_summary$prop_years_occ)

# --- 3. Overall (dataset-wide) summary ---------------------------------------

overall_summary <- rasy_long %>%
  filter(!is.na(Count)) %>%summarise(
    mean_prop_sites_occ_per_year = mean(yearly_summary$prop_occupied, na.rm = TRUE),
    mean_egg_masses_when_occ     = mean(Count[Count > 0], na.rm = TRUE),
    sd_egg_masses_when_occ       = sd(Count[Count > 0], na.rm = TRUE),
    mean_prop_years_occ_per_site = mean(site_summary$prop_years_occ, na.rm = TRUE)
  )

