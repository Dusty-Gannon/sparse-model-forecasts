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
  slab_df = 6
)

mfit <- rstan::sampling(
  hs_reg, data = datlist,
  cores = 4, control = list(adapt_delta = 0.99)
)

# get unshrunk coefficients
w <- colMeans(
  rstan::extract(mfit, pars = "beta")$beta
)
w_raw <- lm(ts ~ datlist$X - 1)$coefficients

theme <- theme_classic() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  )

large <- colnames(datlist$X)[(order(abs(w[-1]), decreasing = T) + 1)[1:2]]
small <- colnames(datlist$X)[
  order(abs(w_raw), decreasing = T) %in% which(round(w, 2) == 0)
][1:2]
large_ids <- which(colnames(datlist$X) %in% large)
small_ids <- which(colnames(datlist$X) %in% small)

clrs <- c(PNWColors::pnw_palette("Bay", 2), "#974e62")

raw <- ggplot(df, aes(x = time)) +
  geom_line(aes(y = y)) +
  theme +
  xlab("") +
  ggtitle("raw data")

# convert names
df <- cbind(
  df, decomp
)
df <- rename_with(df, ~ gsub("-", "_", .x, fixed = TRUE))
large <- gsub("-", "_", large)
small <- gsub("-", "_", small)
colnames(datlist$X) <- gsub("-", "_", colnames(datlist$X))

l1 <- ggplot(df, aes_string(x = "time", large[1])) +
  geom_line(color = clrs[1]) +
  theme_void() +
  xlab("") +
  ylab(expression(w[1])) +
  ggtitle("example 1")

l2 <- ggplot(df, aes_string(x = "time", y = large[2])) +
  geom_line(color = clrs[1]) +
  theme_void() +
  ylab(expression(w[2])) +
  xlab("") +
  ggtitle("example 2")


s1 <- ggplot(df, aes_string(x = "time", y = small[1])) +
  geom_line(color = clrs[2]) +
  theme_void()

s2 <- ggplot(df, aes_string(x = "time", y = small[2])) +
  geom_line(color = clrs[2]) +
  theme_void()

# s32 <- ggplot(df, aes(x = time, y = `S32-100`)) +
#   geom_line(color = clrs[5]) +
#   theme +
#   ylab(expression(w[64])) +
#   xlab("")
#
# c33 <- ggplot(df, aes(x = time, y = `C33-100`)) +
#   geom_line(color = clrs[6]) +
#   theme +
#   ylab(expression(w[67]))


# start second panel with linear combos and fitted values

df <- df %>% mutate(
  combo_1 = w_raw[1] + w_raw[large_ids[1]] * datlist$X[,large[1]] + 0.5 * datlist$X[, small[1]],
  combo_2 = w_raw[1] + w_raw[large_ids[2]] * datlist$X[, large[2]] + 0.5 * datlist$X[, small[2]],
  y_hat = colMeans(
    rstan::extract(mfit, pars = "eta")$eta
  )
)

comb1 <- ggplot(df, aes(x = time)) +
  geom_line(data = df, aes(y = y), col = "grey") +
  geom_line(aes(y = combo_1), color = clrs[3], linewidth = 1) +
  theme_void()

comb2 <- ggplot(df, aes(x = time)) +
  geom_line(data = df, aes(y = y), col = "grey") +
  geom_line(aes(y = combo_2), color = clrs[3], linewidth = 1) +
  theme_void()

fitted <- ggplot(df, aes(x = time)) +
  geom_line(data = df, aes(y = y), col = "grey") +
  geom_line(aes(y = y_hat), color = "black", linewidth = 1) +
  theme +
  ggtitle("fitted model")

png(
  filename = here::here("Figures/concept_fig_fourier.png"),
  width = 6, height = 8, res = 300,
  units = "in"
)
(raw) / ((l1 / s1) | comb1) /
  ((l2 / s2) | comb2) /
  fitted
dev.off()



  # ggsave(
  #   here::here("Figures/concept_fig_fourier.png"),
  #   plot = .,
  #   device = "png",
  #   width = 6,
  #   height = 8,
  #   units = "in"
  # )

# panel1_dml <- rvg::dml(ggobj = panel1)
#
# officer::read_pptx() %>%
#   # add slide ----
# officer::add_slide() %>%
#   # specify object and location of object ----
# officer::ph_with(
#   panel1,
#   officer::ph_location(),
#   use_loc_size = F
# ) %>%
#   # export slide -----
# base::print(
#   target = here::here(
#     "Figures/concept_fig_editable.pptx"
#   )
# )






