# constructing conceptual figures

library(tidyverse)
library(patchwork)
devtools::load_all()
set.seed(3400)

n <- 2
freq = 100

N <- n * freq

ts <-  sin(0.05 * 1:N) + 0.8 * sin(0.2 * 1:N) + 0.5 * cos(0.5 * 1:N) +
  0.2 * sin(0.7 * (1:N)) + arima.sim(model = list(ar = c(0.6, 0.2)), n = N, sd = 0.5)

plot(ts, type = "l")

# make y a time series object
y <- stats::ts(ts, frequency = 120)

# do the fourier decomposition
decomp <- forecast::fourier(
  y, K = 50
)

# combine into dataframe
df <- data.frame(
  y = ts,
  time = 1:N
)

# get weights using sparse regression model
hs_reg <- rstan::stan_model(here::here("Stan/sparse_reg_FHS.stan"))

datlist <- list(
  N = N,
  P0 = 1,
  P = ncol(decomp) + 1,
  y = ts,
  X = cbind(rep(1, N), decomp),
  tau0 = tau0(y, 5, ncol(decomp), N, fam = "gaussian"),
  slab_scl = 2,
  slab_df = 6,
  N_new = 0,
  X_new = matrix(nrow = 0, ncol = ncol(decomp) + 1)
)

mfit <- rstan::sampling(
  hs_reg, data = datlist,
  cores = 4, control = list(adapt_delta = 0.99)
)

# get shrunk coefficients
w <- colMeans(
  rstan::extract(mfit, pars = "beta")$beta
)
w_raw <- lm(ts ~ datlist$X - 1)$coefficients

theme <- theme_classic() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  )

large <- colnames(datlist$X)[which(round(abs(w), 1) > 0)]

large_ids <- which(colnames(datlist$X) %in% large)


raw <- ggplot(df, aes(x = time)) +
  geom_line(aes(y = y)) +
  theme +
  xlab("time")

# convert names
df <- cbind(
  df, decomp
)
df <- rename_with(df, ~ gsub("-", "_", .x, fixed = TRUE))
large <- gsub("-", "_", large)
small <- gsub("-", "_", small)
colnames(datlist$X) <- gsub("-", "_", colnames(datlist$X))

#### Waves panel ####

theme_trigs <- theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    axis.text.y = element_text(size = 12)
  )


df_long_waves <- df %>% select(c(time:C4_120, S10_120, C10_120, S50_120, C50_120)) %>%
  pivot_longer(
    cols = S1_120:C50_120,
    names_to = "wave",
    values_to = "y"
  )

# add colors
colors_key <- data.frame(
  wave = names(df)[c(3:10, 21:22, 101:102)],
  color = rep(PNWColors::pnw_palette("Bay", 6), each = 2)
)

df_long_waves <- left_join(
  df_long_waves,
  colors_key,
  by = "wave"
)
# make sure order is as I want
df_long_waves <- df_long_waves %>% mutate(
  wave = factor(wave, levels = names(df)[c(3:10, 21:22, 101:102)]),
  color = factor(color, levels = PNWColors::pnw_palette("Bay", 6))
)


# make the plot
waves_plot <- ggplot(df_long_waves, aes(x = time, y = y, color = color)) +
  facet_wrap(~wave, ncol = 1) +
  geom_line(linewidth = 0.75) +
  theme_trigs +
  xlab("") +
  ylab("") +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  theme(axis.line.x = element_blank()) +
  scale_color_manual(values = PNWColors::pnw_palette("Bay", 6)) +
  scale_y_continuous(breaks = c(-1, 1))




#### Wave combos ####

df_combo1_waves <- df_long_waves %>% filter(
  wave == "S1_120" | wave == "S10_120"
)

combo1_waves_plot <- ggplot(df_combo1_waves, aes(x = time, y = y, color = color)) +
  facet_wrap(~wave, ncol = 1, labeller = labeller(wave = c("", ""))) +
  geom_line(linewidth = 0.75) +
  theme_trigs +
  xlab("") +
  ylab("") +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  theme(axis.line.x = element_blank()) +
  scale_color_manual(values = unique(df_combo1_waves$color)) +
  scale_y_continuous(breaks = c(-1, 1))

df_combo2_waves <- df_long_waves %>% filter(
  wave == "C1_120" | wave == "S2_120"
)

combo2_waves_plot <- ggplot(df_combo2_waves, aes(x = time, y = y, color = color)) +
  facet_wrap(~wave, ncol = 1) +
  geom_line(linewidth = 0.75) +
  theme_trigs +
  xlab("") +
  ylab("") +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  theme(axis.line.x = element_blank()) +
  scale_color_manual(values = unique(df_combo2_waves$color)) +
  scale_y_continuous(breaks = c(-1, 1))

df_fit_waves <- df_long_waves %>% filter(
  wave %in% large
)

fit_waves_plot <- ggplot(df_fit_waves, aes(x = time, y = y, color = color)) +
  facet_wrap(~wave, ncol = 1) +
  geom_line(linewidth = 1) +
  theme_trigs +
  xlab("") +
  ylab("") +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  theme(axis.line.x = element_blank()) +
  scale_color_manual(values = unique(df_fit_waves$color)) +
  scale_y_continuous(breaks = c(-1, 1))

#### Linear combinations ####

df_combos <- df %>% mutate(
  combo_1 = w_raw[1] +  0.8 * S1_120 + 0.5 * S10_120,
  combo_2 = w_raw[1] + 0.5 * C1_120 + 0.2 * S2_120,
  y_hat = colMeans(
    rstan::extract(mfit, pars = "eta")$eta
  )
)

# find the color of the mixtures
w_c1 <- c(0.8, 0.5) / (0.8 + 0.5)
col_comb1 <- col2rgb(unique(df_combo1_waves$color)) %*% diag(w_c1) |>
  rowMeans() |>
  as.list() |> c(list(maxColorValue = 255)) |>
  do.call(rgb, args = _)

w_c2 <- c(0.5, 0.2) / (0.5 + 0.2)
col_comb2 <- col2rgb(unique(df_combo2_waves$color)) %*% diag(w_c2) |>
  rowMeans() |>
  as.list() |> c(list(maxColorValue = 255)) |>
  do.call(rgb, args = _)

w_fit <- abs(w[large_ids]) / sum(abs(w[large_ids]))
col_comb3 <- col2rgb(rep(unique(df_fit_waves$color), each = 2)) %*% diag(w_fit) |>
  rowMeans() |>
  as.list() |> c(list(maxColorValue = 255)) |>
  do.call(rgb, args = _)

comb1 <- ggplot(df_combos, aes(x = time)) +
  geom_line(data = df, aes(y = y), col = "grey") +
  geom_line(aes(y = combo_1), linewidth = 1.2, col = col_comb1) +
  theme

comb2 <- ggplot(df_combos, aes(x = time)) +
  geom_line(data = df, aes(y = y), col = "grey") +
  geom_line(aes(y = combo_2), linewidth = 1.2, col = col_comb2) +
  theme

fitted <- ggplot(df_combos, aes(x = time)) +
  geom_line(data = df, aes(y = y), col = "grey") +
  geom_line(aes(y = y_hat), color = col_comb3, linewidth = 1.2) +
  theme

#### Save sub-panels for manual assembly ####

if(!dir.exists(here::here("Figures/fourier_concept_subplots"))){
  dir.create(here::here("Figures/fourier_concept_subplots"))
}


ggsave(
  here::here("Figures/fourier_concept_subplots/data.png"),
  plot = raw, device = "png",
  width = 4, height = 2, units = "in",
  dpi = 300
)
ggsave(
  here::here("Figures/fourier_concept_subplots/fitted.png"),
  plot = fitted, device = "png",
  width = 4, height = 2, units = "in",
  dpi = 300
)
ggsave(
  here::here("Figures/fourier_concept_subplots/waves.png"),
  plot = waves_plot, device = "png",
  width = 4, height = 6, units = "in",
  dpi = 300
)
ggsave(
  here::here("Figures/fourier_concept_subplots/waves_comb1.png"),
  plot = combo1_waves_plot, device = "png",
  width = 4, height = 1, units = "in",
  dpi = 300
)
ggsave(
  here::here("Figures/fourier_concept_subplots/waves_comb2.png"),
  plot = combo2_waves_plot, device = "png",
  width = 4, height = 1, units = "in",
  dpi = 300
)
ggsave(
  here::here("Figures/fourier_concept_subplots/fitted_waves.png"),
  plot = fit_waves_plot, device = "png",
  width = 4, height = 2, units = "in",
  dpi = 300
)
ggsave(
  here::here("Figures/fourier_concept_subplots/comb1.png"),
  plot = comb1, device = "png",
  width = 4, height = 2, units = "in",
  dpi = 300
)
ggsave(
  here::here("Figures/fourier_concept_subplots/comb2.png"),
  plot = comb2, device = "png",
  width = 4, height = 2, units = "in",
  dpi = 300
)









