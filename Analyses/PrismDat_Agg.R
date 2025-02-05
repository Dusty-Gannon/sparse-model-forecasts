# This script contains different aggregations of prism data to split
# each variable into summer (June-Sept), winter(Oct-May), wateryea (Oct-Sept)

library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggcorrplot)

BMG_tmin <- read_csv("SparseTS_prismdata/BMG_tmin.csv")
BMG_tmax <- read_csv("SparseTS_prismdata/BMG_tmax.csv")
BMG_td_mean <- read_csv("SparseTS_prismdata/BMG_td_mean.csv")
BMG_ppt <- read_csv("SparseTS_prismdata/BMG_ppt.csv")
BMG_tmean <- read_csv("SparseTS_prismdata/BMG_tmean.csv")
BMG_vpdmin <- read_csv("SparseTS_prismdata/BMG_vpdmin.csv")
BMG_vpdmax <- read_csv("SparseTS_prismdata/BMG_vpdmax.csv")

## create summer season variables ##
summer_months <- c(6,7,8,9)

# tmin
summer_tmin <- BMG_tmin %>%
  filter(month %in% summer_months) %>%
  group_by(year) %>%
  summarise(summer_tmin = mean(tmin))

# tmax
summer_tmax <- BMG_tmax %>%
  filter(month %in% summer_months) %>%
  group_by(year) %>%
  summarise(summer_tmax = mean(tmax))

# tmean
summer_tmean <- BMG_tmean %>%
  filter(month %in% summer_months) %>%
  group_by(year) %>%
  summarise(summer_tmean = mean(tmean))

# tdew_mean
summer_td_mean <- BMG_td_mean %>%
  filter(month %in% summer_months) %>%
  group_by(year) %>%
  summarise(summer_td_mean = mean(td_mean))

# ppt
summer_ppt <- BMG_ppt %>%
  filter(month %in% summer_months) %>%
  group_by(year) %>%
  summarise(summer_ppt = mean(ppt))

# vpdmin
summer_vpdmin <- BMG_vpdmin %>%
  filter(month %in% summer_months) %>%
  group_by(year) %>%
  summarise(summer_vpdmin = mean(vpdmin))

# vpdmax
summer_vpdmax <- BMG_vpdmax %>%
  filter(month %in% summer_months) %>%
  group_by(year) %>%
  summarise(summer_vpdmax = mean(vpdmax))

## create winter season variables ##
winter_months <- c(1,2,3,4,5,10,11,12)

# ppt
BMG_ppt <- BMG_ppt %>%
  mutate(wateryear = case_when(month == "10" ~ year+1,
                               month == "11" ~ year+1,
                               month == "12" ~ year+1,
                               TRUE ~ as.numeric(year)))

winter_ppt <- BMG_ppt %>%
  filter(month %in% winter_months) %>%
  group_by(wateryear) %>%
  summarise(winter_ppt = mean(ppt))

# tmin
BMG_tmin <- BMG_tmin %>%
  mutate(wateryear = case_when(month == "10" ~ year+1,
                               month == "11" ~ year+1,
                               month == "12" ~ year+1,
                               TRUE ~ as.numeric(year)))

winter_tmin <- BMG_tmin %>%
  filter(month %in% winter_months) %>%
  group_by(wateryear) %>%
  summarise(winter_tmin = mean(tmin))

# tmax
BMG_tmax <- BMG_tmax %>%
  mutate(wateryear = case_when(month == "10" ~ year+1,
                               month == "11" ~ year+1,
                               month == "12" ~ year+1,
                               TRUE ~ as.numeric(year)))

winter_tmax <- BMG_tmax %>%
  filter(month %in% winter_months) %>%
  group_by(wateryear) %>%
  summarise(winter_tmax = mean(tmax))

# tmean
BMG_tmean <- BMG_tmean %>%
  mutate(wateryear = case_when(month == "10" ~ year+1,
                               month == "11" ~ year+1,
                               month == "12" ~ year+1,
                               TRUE ~ as.numeric(year)))

winter_tmean <- BMG_tmean %>%
  filter(month %in% winter_months) %>%
  group_by(wateryear) %>%
  summarise(winter_tmean = mean(tmean))


# tdew_mean
BMG_td_mean <- BMG_td_mean %>%
  mutate(wateryear = case_when(month == "10" ~ year+1,
                               month == "11" ~ year+1,
                               month == "12" ~ year+1,
                               TRUE ~ as.numeric(year)))

winter_td_mean <- BMG_td_mean %>%
  filter(month %in% winter_months) %>%
  group_by(wateryear) %>%
  summarise(winter_td_mean = mean(td_mean))


# vpdmin
BMG_vpdmin <- BMG_vpdmin %>%
  mutate(wateryear = case_when(month == "10" ~ year+1,
                               month == "11" ~ year+1,
                               month == "12" ~ year+1,
                               TRUE ~ as.numeric(year)))

winter_vpdmin <- BMG_vpdmin %>%
  filter(month %in% winter_months) %>%
  group_by(wateryear) %>%
  summarise(winter_vpdmin = mean(vpdmin))


# vpdmax
BMG_vpdmax <- BMG_vpdmax %>%
  mutate(wateryear = case_when(month == "10" ~ year+1,
                               month == "11" ~ year+1,
                               month == "12" ~ year+1,
                               TRUE ~ as.numeric(year)))

winter_vpdmax <- BMG_vpdmax %>%
  filter(month %in% winter_months) %>%
  group_by(wateryear) %>%
  summarise(winter_vpdmax = mean(vpdmax))


## create wateryear variables ##

# ppt
wateryear_ppt <- BMG_ppt %>%
  group_by(wateryear) %>%
  summarise(wy_ppt = mean(ppt))

# tmin
wateryear_tmin <- BMG_tmin %>%
  group_by(wateryear) %>%
  summarise(wy_tmin = mean(tmin))

# tmax
wateryear_tmax <- BMG_tmax %>%
  group_by(wateryear) %>%
  summarise(wy_tmax = mean(tmax))

# tmean
wateryear_tmean <- BMG_tmean %>%
  group_by(wateryear) %>%
  summarise(wy_tmean = mean(tmean))

# tdew_mean
wateryear_td_mean <- BMG_td_mean %>%
  group_by(wateryear) %>%
  summarise(wy_td_mean = mean(td_mean))

# vpdmin
wateryear_vpdmin <- BMG_vpdmin %>%
  group_by(wateryear) %>%
  summarise(wy_vpdmin = mean(vpdmin))

# vpdmax
wateryear_vpdmax <- BMG_vpdmax %>%
  group_by(wateryear) %>%
  summarise(wy_vpdmax = mean(vpdmax))


## save all variables ##
prism_summer <- left_join(summer_ppt, summer_tmin) %>%
  left_join(summer_tmax) %>%
  left_join(summer_tmean) %>%
  left_join(summer_td_mean) %>%
  left_join(summer_vpdmin) %>%
  left_join(summer_vpdmax)

prism_winter <- left_join(winter_ppt, winter_tmin) %>%
  left_join(winter_tmax) %>%
  left_join(winter_tmean) %>%
  left_join(winter_td_mean) %>%
  left_join(winter_vpdmin) %>%
  left_join(winter_vpdmax)%>%
  rename(year = "wateryear")

prism_wateryear <- left_join(wateryear_ppt, wateryear_tmin) %>%
  left_join(wateryear_tmax) %>%
  left_join(wateryear_tmean) %>%
  left_join(wateryear_td_mean) %>%
  left_join(wateryear_vpdmin) %>%
  left_join(wateryear_vpdmax)
