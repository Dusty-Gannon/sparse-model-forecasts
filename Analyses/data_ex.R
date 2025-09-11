# This script contains exploration, selection, and standardization of empirical
# data from the ITRDB

# load packages
library(tidyverse)
library(dplR)

# load Giant sequoia data downloaded from collection CA726 in ITRDB
ca726_se_raw <- read.table(url("https://www.ncei.noaa.gov/pub/data/paleo/treering/measurements/northamerica/usa/ca726-rwl-noaa.txt"), header = T, sep = "\t", stringsAsFactors = FALSE)

# pivot raw data into long format
ca726_se_long <- ca726_se_raw %>%
  rename(year = age_CE) %>%
  tidyr::pivot_longer(2:259, names_to = "tree_id", values_to = "trw")

# remove trees with missing data for a particular year (ie. tree was not alive that year)
ca726_se_long <- na.omit(ca726_se_long)

# determine the min and max year each tree has data
ca726_se_wide_ts <- ca726_se_long %>%
  group_by(tree_id) %>%
  summarise(min_year = min(year), max_year = max(year)) %>%
  filter(max_year > 1970) %>%
  filter(min_year < 1800)

# pivot longer
ca726_se_long_ts <- ca726_se_wide_ts %>%
  pivot_longer(2:3, names_to = "year_type", values_to = "year")

# create a loop to determine the number of trees with data for a particular timeseries
start_years <- seq(1600, 1800, by = 10)
end_years <- unique(ca726_se_wide_ts$max_year)
timeseries_container <- tibble()

for(i in start_years){
  for(j in end_years){
    temp <- ca726_se_wide_ts %>%
      select(min_year, max_year) %>%
      filter(min_year <= i) %>%
      filter(max_year == j)

    temp$start_year <- i
    temp$end_year <- j
    temp$num_trees <- dim(temp)[1]

    timeseries_container <- rbind(timeseries_container, temp)
  }
}

# select timeseries with at least 5 trees with data
timeseries_container <- timeseries_container %>%
  select(start_year, end_year, num_trees) %>%
  filter(num_trees >= 5) %>%
  distinct()

# extract the tree_ids for the trees with data for a particular timeseries
# 1600 - 2012
ca726_se_1600_2012 <- ca726_se_wide_ts %>%
  filter(max_year == 2012) %>%
  filter(min_year <= 1600)

dim(ca726_se_1600_2012)[1] # 5 trees with data

# collect rwl data for those particular trees with data for the full timeseries
trees <- ca726_se_1600_2012$tree_id
ca726_se_long_1600_2012_filtered <- ca726_se_long %>%
  filter(tree_id %in% trees)

# filter data for our chosen timeseries min
ca726_se_long_1600_2012_filtered <- ca726_se_long_1600_2012_filtered %>%
  filter(year >= 1600)

# sumarise avg growth across our timeseries
ca726_se_1600_2012_avg_growth <- ca726_se_long_1600_2012_filtered %>%
  group_by(year) %>%
  summarise(avg.trw = mean(trw), sd.trw = sd(trw),n.trw = n()) %>%
  mutate(se.trw = sd.trw / sqrt(n.trw),
         lower.ci.trw = avg.trw - qt(1 - (0.05 / 2), n.trw - 1) * se.trw,
         upper.ci.trw = avg.trw + qt(1 - (0.05 / 2), n.trw - 1) * se.trw)

ca726_se_1600_2012_avg_growth$tree <- 5

# extract the tree_ids for the trees with data for a particular timeseries
# 1700 - 2012
ca726_se_1700_2012 <- ca726_se_wide_ts %>%
  filter(max_year == 2012) %>%
  filter(min_year <= 1700)

dim(ca726_se_1700_2012)[1] # 20 trees with data

# collect rwl data for those particular trees with data for the full timeseries
trees <- ca726_se_1700_2012$tree_id
ca726_se_long_1700_2012_filtered <- ca726_se_long %>%
  filter(tree_id %in% trees)

# filter data for our chosen timeseries min
ca726_se_long_1700_2012_filtered <- ca726_se_long_1700_2012_filtered %>%
  filter(year >= 1700)

# sumarise avg growth across our timeseries
ca726_se_1700_2012_avg_growth <- ca726_se_long_1700_2012_filtered %>%
  group_by(year) %>%
  summarise(avg.trw = mean(trw), sd.trw = sd(trw),n.trw = n()) %>%
  mutate(se.trw = sd.trw / sqrt(n.trw),
         lower.ci.trw = avg.trw - qt(1 - (0.05 / 2), n.trw - 1) * se.trw,
         upper.ci.trw = avg.trw + qt(1 - (0.05 / 2), n.trw - 1) * se.trw)

ca726_se_1700_2012_avg_growth$tree <- 20


# extract the tree_ids for the trees with data for a particular timeseries
# 1800 - 2012
ca726_se_1800_2012 <- ca726_se_wide_ts %>%
  filter(max_year == 2012) %>%
  filter(min_year <= 1800)

dim(ca726_se_1800_2012)[1] # 41 trees with data

# collect rwl data for those particular trees with data for the full timeseries
trees <- ca726_se_1800_2012$tree_id
ca726_se_long_1800_2012_filtered <- ca726_se_long %>%
  filter(tree_id %in% trees)

# filter data for our chosen timeseries min
ca726_se_long_1800_2012_filtered <- ca726_se_long_1800_2012_filtered %>%
  filter(year >= 1800)

# sumarise avg growth across our timeseries
ca726_se_1800_2012_avg_growth <- ca726_se_long_1800_2012_filtered %>%
  group_by(year) %>%
  summarise(avg.trw = mean(trw), sd.trw = sd(trw),n.trw = n()) %>%
  mutate(se.trw = sd.trw / sqrt(n.trw),
         lower.ci.trw = avg.trw - qt(1 - (0.05 / 2), n.trw - 1) * se.trw,
         upper.ci.trw = avg.trw + qt(1 - (0.05 / 2), n.trw - 1) * se.trw)

ca726_se_1800_2012_avg_growth$tree <- 41


# extract the tree_ids for the trees with data for a particular timeseries
# 1760 - 2012
ca726_se_1760_2012 <- ca726_se_wide_ts %>%
  filter(max_year == 2012) %>%
  filter(min_year <= 1760)

dim(ca726_se_1760_2012)[1] # 20 trees with data

# collect rwl data for those particular trees with data for the full timeseries
trees <- ca726_se_1760_2012$tree_id
ca726_se_long_1760_2012_filtered <- ca726_se_long %>%
  filter(tree_id %in% trees)

# filter data for our chosen timeseries min
ca726_se_long_1760_2012_filtered <- ca726_se_long_1760_2012_filtered %>%
  filter(year >= 1760)

# sumarise avg growth across our timeseries
ca726_se_1760_2012_avg_growth <- ca726_se_long_1760_2012_filtered %>%
  group_by(year) %>%
  summarise(avg.trw = mean(trw), sd.trw = sd(trw),n.trw = n()) %>%
  mutate(se.trw = sd.trw / sqrt(n.trw),
         lower.ci.trw = avg.trw - qt(1 - (0.05 / 2), n.trw - 1) * se.trw,
         upper.ci.trw = avg.trw + qt(1 - (0.05 / 2), n.trw - 1) * se.trw)

ca726_se_1760_2012_avg_growth$tree <- 30


# combine datasets to facet figure
ca726_se_ts_compare_2012 <- rbind(ca726_se_1600_2012_avg_growth, ca726_se_1700_2012_avg_growth,
                             ca726_se_1760_2012_avg_growth, ca726_se_1800_2012_avg_growth)

ggplot(ca726_se_ts_compare_2012) +
  geom_line(aes(x = year, y=avg.trw), color = "black") +
  geom_ribbon(aes(x =year, y = avg.trw, ymin = lower.ci.trw, ymax = upper.ci.trw), alpha = 0.2) +
  xlab("Year") +
  ylab("Average TRW (mm)") +
  theme_minimal()+
  facet_grid(vars(tree), scales = "free")


# ca726_se_ts_compare$tree<- factor(ca726_se_ts_compare_2012$tree, levels = c(5, 20, 30, 41))
# ggplot(ca726_se_ts_compare, aes(x = tree, y = avg.trw)) +
#   geom_boxplot() +
#   ylab("Average TRW (mm)") +
#   theme_minimal()

# select trees from BMC 1800-2012 timeseries for final analyses
BMC_trees <- ca726_se_1800_2012$tree_id[1:22]

ca726_seq_BMC <- ca726_se_long %>%
  filter(tree_id %in% BMC_trees) %>%
  filter(year >= 1800)

# pull in data from BM site specifically CA719
ca719_seq_rwl <- read.table(url("https://www.ncei.noaa.gov/pub/data/paleo/treering/measurements/northamerica/usa/ca719-rwl-noaa.txt"), header = TRUE, sep = "", stringsAsFactors = FALSE)
# pivot longer and subset trees identified as alive for entire timeseries
ca719_seq_rwl_long <- ca719_seq_rwl %>%
  pivot_longer(2:121, names_to = "tree_id", values_to = "rwl") %>%
  filter(tree_id %in% BMC_trees)

# pivot wider and add row names so that tree_ids are the columns (representing a series)
# and year as rows
ca719_seq_rwl <- ca719_seq_rwl_long %>%
  pivot_wider(names_from = tree_id, values_from = rwl)
rownames(ca719_seq_rwl) <- ca719_seq_rwl$age_CE

# detrend the data using the spline method
ca719_seq_detrended <- detrend(ca719_seq_rwl, y.name=names(ca719_seq_rwl), make.plot = TRUE,
                               method = "AgeDepSpline", nyrs = NULL)

# make year a column again
ca719_seq_detrended$year <- rownames(ca719_seq_detrended)

# delete age_CE since it was also detrended (hense the negative warning)
ca719_seq_detrended <- ca719_seq_detrended %>%
  select(-age_CE)

# make year numeric and pivot data raw, subsetting after 1800
ca719_seq_detrended$year <- as.numeric(ca719_seq_detrended$year)
ca719_seq_detrended_long <- ca719_seq_detrended %>%
  pivot_longer(1:22, names_to = "tree_id", values_to = "rwi") %>%
  filter(year >= 1800)

# save data
# saveRDS("SparseTS_prismdata/ca719_BM_Seq_1800_2012.rds" )


