# This script contains investigation of trends, breakpoints, and correlations in empirical data.

# load necessary packages
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggcorrplot)
library(segmented)

# load data
treering_dat <- read_csv("SparseTS_prismdata/BMG_rwi.csv")
BMG_tmin <- read_csv("SparseTS_prismdata/BMG_tmin.csv")
BMG_tmax <- read_csv("SparseTS_prismdata/BMG_tmax.csv")
BMG_td_mean <- read_csv("SparseTS_prismdata/BMG_td_mean.csv")
BMG_ppt <- read_csv("SparseTS_prismdata/BMG_ppt.csv")
BMG_tmean <- read_csv("SparseTS_prismdata/BMG_tmean.csv")
BMG_vpdmin <- read_csv("SparseTS_prismdata/BMG_vpdmin.csv")
BMG_vpdmax <- read_csv("SparseTS_prismdata/BMG_vpdmax.csv")


#### data trends  ####
# tree-ring trends
treering_trend <- treering_dat%>%
  group_by(year)%>%
  summarize(rwi = mean(rwi))
treering_trend$tree_id <- "Trend"
treering_dat$year <- as.character(treering_dat$year)
treering_trend$year <- as.character(treering_trend$year)
raw_rwi <- ggplot()+
  geom_line(data = treering_dat, aes(x=year, y=rwi, group=tree_id, col=tree_id))+
  geom_line(data = treering_trend, mapping = aes(x=year, y=rwi, group = tree_id, col = tree_id), color = "black", linewidth = 1.25)+
  theme_bw()+
  scale_colour_grey(start = 1, end = 0.5)+
  scale_x_discrete(breaks = seq(1800,2012,10))+
  theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 45, hjust = 1, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, face = "plain"),
        legend.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(color = "grey20", size = 10,angle = 0, hjust = 0, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank()) +
  xlab("Year")+
  ylab("Ring Width Index")

# precipitaiton trends
BMG_ppt_trend <- BMG_ppt%>%
  group_by(year)%>%
  summarize(avg_ppt = mean(ppt))
raw_ppt <- ggplot(data = BMG_ppt_trend, aes(x=year, y=avg_ppt))+
  geom_line()+
  theme_bw()+
  xlab("Year")+
  ylab("Average Annual Precip")

# minimum temperature trends
BMG_tmin_trend <- BMG_tmin%>%
  group_by(year)%>%
  summarize(avg_tmin = mean(tmin))
raw_tmin <- ggplot(data = BMG_tmin_trend, aes(x=year, y=avg_tmin))+
  geom_line()+
  theme_bw()+
  xlab("Year")+
  ylab("Average Annual Tmin")

# maximum temperature trends
BMG_tmax_trend <- BMG_tmax%>%
  group_by(year)%>%
  summarize(avg_tmax = mean(tmax))
raw_tmax <- ggplot(data = BMG_tmax_trend, aes(x=year, y=avg_tmax))+
  geom_line()+
  theme_bw()+
  xlab("Year")+
  ylab("Average Annual Tmax")

# mean temperature trends
BMG_tmean_trend <- BMG_tmean%>%
  group_by(year)%>%
  summarize(avg_tmean = mean(tmean))
raw_tmean <- ggplot(data = BMG_tmean_trend, aes(x=year, y=avg_tmean))+
  geom_line()+
  theme_bw()+
  xlab("Year")+
  ylab("Average Annual Tmean")

# dewpoint temperature trends
BMG_td_mean_trend <- BMG_td_mean%>%
  group_by(year)%>%
  summarize(avg_td_mean = mean(td_mean))
raw_td_mean <- ggplot(data = BMG_td_mean_trend, aes(x=year, y=avg_td_mean))+
  geom_line()+
  theme_bw()+
  xlab("Year")+
  ylab("Average Annual Dew Point Temp")

# minimum vapor pressure deficit trends
BMG_vpdmin_trend <- BMG_vpdmin%>%
  group_by(year)%>%
  summarize(avg_vpdmin = mean(vpdmin))
raw_vpdmin <- ggplot(data = BMG_vpdmin_trend, aes(x=year, y=avg_vpdmin))+
  geom_line()+
  theme_bw()+
  xlab("Year")+
  ylab("Average Annual VPDmin")

# maximum vapor pressure deficit trends
BMG_vpdmax_trend <- BMG_vpdmax%>%
  group_by(year)%>%
  summarize(avg_vpdmax = mean(vpdmax))
raw_vpdmax <- ggplot(data = BMG_vpdmax_trend, aes(x=year, y=avg_vpdmax))+
  geom_line()+
  theme_bw()+
  xlab("Year")+
  ylab("Average Annual VPDmax")

#### breakpoint analyses ####
# tree ring break point
treering_trend$year <- as.numeric(treering_trend$year)
treering.p <- ggplot(data = treering_trend, aes(x=year, y = rwi)) + geom_line()
treering.lm <- lm(rwi~year, data = treering_trend)
summary(treering.lm)
treering.coeff <- coef(treering.lm)
treering.p.line <- treering.p + geom_abline(intercept = treering.coeff[1],
                                            slope = treering.coeff[2])
treering.seg <- segmented(treering.lm,
                          seg.Z = ~ year,
                          psi = list(year = c(1970)))
summary(treering.seg)
treering.seg$psi
slope(treering.seg)

treering.fitted <- fitted(treering.seg)
treering.model <- data.frame(year = treering_trend$year, trend = treering.fitted)
treering.p <- treering.p + geom_line(data = treering.model, aes(x=year, y = trend), color = "tomato")
treering.lines <- treering.seg$psi[,2]
treering.p <- treering.p + geom_vline(xintercept = treering.lines, linetype = "dashed", color = "darkgrey") +
  labs(y = "avg treering")+
  theme_bw()

# precipitation break point
ppt.p <- ggplot(data = BMG_ppt_trend, aes(x=year, y = avg_ppt)) + geom_line()
ppt.lm <- lm(avg_ppt~year, data = BMG_ppt_trend)
summary(ppt.lm)
ppt.coeff <- coef(ppt.lm)
ppt.p.line <- ppt.p + geom_abline(intercept = ppt.coeff[1],
                                  slope = ppt.coeff[2])
ppt.seg <- segmented(ppt.lm,
                     seg.Z = ~ year,
                     psi = list(year = c(1970)))
summary(ppt.seg)
ppt.seg$psi
slope(ppt.seg)

ppt.fitted <- fitted(ppt.seg)
ppt.model <- data.frame(year = BMG_ppt_trend$year, trend = ppt.fitted)
ppt.p <- ppt.p + geom_line(data = ppt.model, aes(x=year, y = trend), color = "tomato")
ppt.lines <- ppt.seg$psi[,2]
ppt.p <- ppt.p + geom_vline(xintercept = ppt.lines, linetype = "dashed", color = "darkgrey") +
  labs(y = "avg ppt")+
  theme_bw()

# minimum temperature break point
tmin.p <- ggplot(data = BMG_tmin_trend, aes(x=year, y = avg_tmin)) + geom_line()
tmin.lm <- lm(avg_tmin~year, data = BMG_tmin_trend)
summary(tmin.lm)
tmin.coeff <- coef(tmin.lm)
tmin.p.line <- tmin.p + geom_abline(intercept = tmin.coeff[1],
                                    slope = tmin.coeff[2])
tmin.seg <- segmented(tmin.lm,
                      seg.Z = ~ year,
                      psi = list(year = c(1970)))
summary(tmin.seg)
tmin.seg$psi
slope(tmin.seg)

tmin.fitted <- fitted(tmin.seg)
tmin.model <- data.frame(year = BMG_tmin_trend$year, trend = tmin.fitted)
tmin.p <- tmin.p + geom_line(data = tmin.model, aes(x=year, y = trend), color = "tomato")
tmin.lines <- tmin.seg$psi[,2]
tmin.p <- tmin.p + geom_vline(xintercept = tmin.lines, linetype = "dashed", color = "darkgrey") +
  labs(y = "avg min temp")+
  theme_bw()

# maximum temperature break point
tmax.p <- ggplot(data = BMG_tmax_trend, aes(x=year, y = avg_tmax)) + geom_line()
tmax.lm <- lm(avg_tmax~year, data = BMG_tmax_trend)
summary(tmax.lm)
tmax.coeff <- coef(tmax.lm)
tmax.p.line <- tmax.p + geom_abline(intercept = tmax.coeff[1],
                                  slope = tmax.coeff[2])
tmax.seg <- segmented(tmax.lm,
                     seg.Z = ~ year,
                     psi = list(year = c(1970)))
summary(tmax.seg)
tmax.seg$psi
slope(tmax.seg)

tmax.fitted <- fitted(tmax.seg)
tmax.model <- data.frame(year = BMG_tmax_trend$year, trend = tmax.fitted)
tmax.p <- tmax.p + geom_line(data = tmax.model, aes(x=year, y = trend), color = "tomato")
tmax.lines <- tmax.seg$psi[,2]
tmax.p <- tmax.p + geom_vline(xintercept = tmax.lines, linetype = "dashed", color = "darkgrey") +
  labs(y = "avg tmax")+
  theme_bw()


# mean temperature break point
tmean.p <- ggplot(data = BMG_tmean_trend, aes(x=year, y = avg_tmean)) + geom_line()
tmean.lm <- lm(avg_tmean~year, data = BMG_tmean_trend)
summary(tmean.lm)
tmean.coeff <- coef(tmean.lm)
tmean.p.line <- tmean.p + geom_abline(intercept = tmean.coeff[1],
                                    slope = tmean.coeff[2])
tmean.seg <- segmented(tmean.lm,
                      seg.Z = ~ year,
                      psi = list(year = c(1970)))
summary(tmean.seg)
tmean.seg$psi
slope(tmean.seg)

tmean.fitted <- fitted(tmean.seg)
tmean.model <- data.frame(year = BMG_tmean_trend$year, trend = tmean.fitted)
tmean.p <- tmean.p + geom_line(data = tmean.model, aes(x=year, y = trend), color = "tomato")
tmean.lines <- tmean.seg$psi[,2]
tmean.p <- tmean.p + geom_vline(xintercept = tmean.lines, linetype = "dashed", color = "darkgrey") +
  labs(y = "avg tmean")+
  theme_bw()

# dew point temperature break point
td_mean.p <- ggplot(data = BMG_td_mean_trend, aes(x=year, y = avg_td_mean)) + geom_line()
td_mean.lm <- lm(avg_td_mean~year, data = BMG_td_mean_trend)
summary(td_mean.lm)
td_mean.coeff <- coef(td_mean.lm)
td_mean.p.line <- td_mean.p + geom_abline(intercept = td_mean.coeff[1],
                                      slope = td_mean.coeff[2])
td_mean.seg <- segmented(td_mean.lm,
                       seg.Z = ~ year,
                       psi = list(year = c(1970)))
summary(td_mean.seg)
td_mean.seg$psi
slope(td_mean.seg)

td_mean.fitted <- fitted(td_mean.seg)
td_mean.model <- data.frame(year = BMG_td_mean_trend$year, trend = td_mean.fitted)
td_mean.p <- td_mean.p + geom_line(data = td_mean.model, aes(x=year, y = trend), color = "tomato")
td_mean.lines <- td_mean.seg$psi[,2]
td_mean.p <- td_mean.p + geom_vline(xintercept = td_mean.lines, linetype = "dashed", color = "darkgrey") +
  labs(y = "avg td_mean")+
  theme_bw()

# minimum vapor pressure deficit break point
vpdmin.p <- ggplot(data = BMG_vpdmin_trend, aes(x=year, y = avg_vpdmin)) + geom_line()
vpdmin.lm <- lm(avg_vpdmin~year, data = BMG_vpdmin_trend)
summary(vpdmin.lm)
vpdmin.coeff <- coef(vpdmin.lm)
vpdmin.p.line <- vpdmin.p + geom_abline(intercept = vpdmin.coeff[1],
                                          slope = vpdmin.coeff[2])
vpdmin.seg <- segmented(vpdmin.lm,
                         seg.Z = ~ year,
                         psi = list(year = c(1970)))
summary(vpdmin.seg)
vpdmin.seg$psi
slope(vpdmin.seg)

vpdmin.fitted <- fitted(vpdmin.seg)
vpdmin.model <- data.frame(year = BMG_vpdmin_trend$year, trend = vpdmin.fitted)
vpdmin.p <- vpdmin.p + geom_line(data = vpdmin.model, aes(x=year, y = trend), color = "tomato")
vpdmin.lines <- vpdmin.seg$psi[,2]
vpdmin.p <- vpdmin.p + geom_vline(xintercept = vpdmin.lines, linetype = "dashed", color = "darkgrey") +
  labs(y = "avg vpdmin")+
  theme_bw()

# maximum vapor pressure deficit breakpoint
vpdmax.p <- ggplot(data = BMG_vpdmax_trend, aes(x=year, y = avg_vpdmax)) + geom_line()
vpdmax.lm <- lm(avg_vpdmax~year, data = BMG_vpdmax_trend)
summary(vpdmax.lm)
vpdmax.coeff <- coef(vpdmax.lm)
vpdmax.p.line <- vpdmax.p + geom_abline(intercept = vpdmax.coeff[1],
                                        slope = vpdmax.coeff[2])
vpdmax.seg <- segmented(vpdmax.lm,
                        seg.Z = ~ year,
                        psi = list(year = c(1970)))
summary(vpdmax.seg)
vpdmax.seg$psi
slope(vpdmax.seg)

vpdmax.fitted <- fitted(vpdmax.seg)
vpdmax.model <- data.frame(year = BMG_vpdmax_trend$year, trend = vpdmax.fitted)
vpdmax.p <- vpdmax.p + geom_line(data = vpdmax.model, aes(x=year, y = trend), color = "tomato")
vpdmax.lines <- vpdmax.seg$psi[,2]
vpdmax.p <- vpdmax.p + geom_vline(xintercept = vpdmax.lines, linetype = "dashed", color = "darkgrey") +
  labs(y = "avg vpdmax")+
  theme_bw()

#### correlations ####
covariate.dat <- left_join(BMG_ppt, BMG_tmin) %>%
  left_join(BMG_tmax) %>%
  left_join(BMG_tmean) %>%
  left_join(BMG_td_mean) %>%
  left_join(BMG_vpdmin) %>%
  left_join(BMG_vpdmax)

covariate.mx <- covariate.dat %>%
  dplyr::select(ppt, tmin, tmax, tmean, td_mean, vpdmin, vpdmax)

covariate.cor<- cor(covariate.mx, use = "complete.obs")
ggcorrplot(covariate.cor, method = "square", type = "lower")

#### moving window correlations ####
# annual trends
covariate.dat.annual <- left_join(BMG_ppt_trend, BMG_tmin_trend) %>%
  left_join(BMG_tmax_trend) %>%
  left_join(BMG_tmean_trend) %>%
  left_join(BMG_td_mean_trend) %>%
  left_join(BMG_vpdmin_trend) %>%
  left_join(BMG_vpdmax_trend)

window_length <- 20
window_delta <- window_length - 1
start_year <- min(covariate.dat.annual$year)
end_year <- max(covariate.dat.annual$year) - window_delta

# create an empty tibble to store correlations from all windows
correlation_trends <- tibble()

# loop through the years to calculate correlations and save results
for (i in start_year:end_year) {
  window_start <- i
  window_end <- i + window_delta
  select_years <- window_start:window_end

  # filter the data for the selected years
  annual_covariate_matrix <- covariate.dat.annual %>%
    filter(year %in% select_years) %>%
    dplyr::select(avg_ppt, avg_tmin, avg_tmax, avg_tmean, avg_td_mean, avg_vpdmin, avg_vpdmax)

  # calculate correlations
  annual_correlation_matrix <- cor(annual_covariate_matrix, use = "complete.obs")
  annual_correlation <- as.data.frame(annual_correlation_matrix)
  annual_correlation$var1 <- row.names(annual_correlation)
  annual_correlation <- annual_correlation %>%
    pivot_longer(1:7, names_to = "var2", values_to = "cor")


  annual_correlation$start_year <- window_start
  annual_correlation$end_year <- window_end

  # combine the results to the main data frame
  correlation_trends <- bind_rows(correlation_trends, annual_correlation)
}

# filter out duplicates and self-correlations (var1 == var2)
correlation_trends <- correlation_trends %>%
  filter(var1 != var2)

# open a PDF device to save all plots
pdf(file = "correlation_trends_plots.pdf", width = 8, height = 6)

# get unique variable names from var1
unique_var1 <- unique(correlation_trends$var1)

# loop through each unique variable in var1
for (variable in unique_var1) {
  # filter the dataset for the current variable
  var_trends <- correlation_trends %>%
    filter(var1 == variable)

  plot <- ggplot(data = var_trends, aes(x = start_year, y = cor,
                                        group = var2, col = var2)) +
    geom_point(size = 1) +
    geom_line() +
    labs(title = paste("Correlation Trends for", variable),
         x = "Start Year", y = "Correlation",
         color = "Variable 2") +
    theme_minimal()

  # print the plot to the PDF
  print(plot)
}


dev.off()


# monthly trends
covariate.dat.mo <- left_join(BMG_ppt, BMG_tmin) %>%
  left_join(BMG_tmax) %>%
  left_join(BMG_tmean) %>%
  left_join(BMG_td_mean) %>%
  left_join(BMG_vpdmin) %>%
  left_join(BMG_vpdmax)

window_length <- 20
window_delta <- window_length - 1
start_year <- min(covariate.dat.mo$year)
end_year <- max(covariate.dat.mo$year) - window_delta

# create an empty tibble to store correlations from all windows
correlation_trends <- tibble()

# loop through the years to calculate correlations and save results
for (i in start_year:end_year) {
  window_start <- i
  window_end <- i + window_delta
  select_years <- window_start:window_end

  # filter the data for the selected years
  mo_covariate_matrix <- covariate.dat.mo %>%
    filter(year %in% select_years) %>%
    dplyr::select(ppt, tmin, tmax, tmean, td_mean, vpdmin, vpdmax)


  mo_correlation_matrix <- cor(mo_covariate_matrix, use = "complete.obs")
  mo_correlation <- as.data.frame(mo_correlation_matrix)
  mo_correlation$var1 <- row.names(mo_correlation)
  mo_correlation <- mo_correlation %>%
    pivot_longer(1:7, names_to = "var2", values_to = "cor")


  mo_correlation$start_year <- window_start
  mo_correlation$end_year <- window_end

  # combine the results to the main data frame
  correlation_trends <- bind_rows(correlation_trends, mo_correlation)
}

# filter out duplicates and self-correlations (var1 == var2)
correlation_trends <- correlation_trends %>%
  filter(var1 != var2)

# open a PDF device to save all plots
pdf(file = "correlation_trends_plots.pdf", width = 8, height = 6)

# get unique variable names from var1
unique_var1 <- unique(correlation_trends$var1)

# loop through each unique variable in var1
for (variable in unique_var1) {
  # filter the dataset for the current variable
  var_trends <- correlation_trends %>%
    filter(var1 == variable)


  plot <- ggplot(data = var_trends, aes(x = start_year, y = cor,
                                        group = var2, col = var2)) +
    geom_point(size = 1) +
    geom_line() +
    labs(title = paste("Correlation Trends for", variable),
         x = "Start Year", y = "Correlation",
         color = "Variable 2") +
    theme_minimal()

  # print the plot to the PDF
  print(plot)
}


dev.off()

