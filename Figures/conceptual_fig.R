# constructing conceptual figures

library(tidyverse)
library(patchwork)
devtools::load_all()
set.seed(3400)

freq <- 100
n <- 5

N <- freq * n

ts <- sin(0.05 * 1:N) + 0.8 * sin(0.2 * 1:N) + 0.5*cos(0.5 * 1:N) +
  0.2 * sin(0.7 * (1:N)) + arima.sim(model = list(ar = c(0.6, 0.2)), n = N, sd = 0.3)

plot(ts, type = "l")

# make y a time series object
y <- stats::ts(ts, frequency = 100)

# do the fourier decomposition
decomp <- forecast::fourier(
  y, K = 50
)

# combine into dataframe
df <- data.frame(
  y = ts,
  time = 1:N
)

# get weights using ols
w <- lm(df$y ~ decomp)$coefficients

df <- cbind(df, as.data.frame(decomp))

theme <- theme_classic() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  )

large <- colnames(decomp)[order(w)[1:3] - 1]
small <- colnames(decomp)[order(w)[length(w):(length(w) - 2)] - 1]

clrs <- PNWColors::pnw_palette("Bay", n = 6)

raw <- ggplot(df, aes(x = time)) +
  geom_line(aes(y = y)) +
  theme +
  xlab("")

# convert names
df <- rename_with(df, ~ gsub("-", "_", .x, fixed = TRUE))
large <- gsub("-", "_", large)
small <- gsub("-", "_", small)

l1 <- ggplot(df, aes_string(x = "time", large[1])) +
  geom_line(color = clrs[1]) +
  theme_void() +
  xlab("") +
  ylab(expression(w[1]))

l2 <- ggplot(df, aes_string(x = "time", y = large[2])) +
  geom_line(color = clrs[2]) +
  theme_void() +
  ylab(expression(w[2])) +
  xlab("")


s11 <- ggplot(df, aes(x = time, y = `S11-100`)) +
  geom_line(color = clrs[3]) +
  theme +
  ylab(expression(w[22])) +
  xlab("")

c11 <- ggplot(df, aes(x = time, y = `S11-100`)) +
  geom_line(color = clrs[4]) +
  theme +
  ylab(expression(w[23])) +
  xlab("")

s32 <- ggplot(df, aes(x = time, y = `S32-100`)) +
  geom_line(color = clrs[5]) +
  theme +
  ylab(expression(w[64])) +
  xlab("")

c33 <- ggplot(df, aes(x = time, y = `C33-100`)) +
  geom_line(color = clrs[6]) +
  theme +
  ylab(expression(w[67]))


panel1 <- (raw / s1 / c1 / s11 / c11 / s32 / c33)

# start second panel with linear combos
all_ids <- c(2, 3, large_ids, small_ids)


df <- df %>% mutate(
  combo_1 = w[1] + w[2] * `S1-100` + w[large_ids[1]] * `S11-100`,
  combo_2 = w[1] + w[2] * `S1-100` + w[small_ids[1]] * `C11-100`,
  combo_3 = w[1] + w[large_ids[1]] * `S11-100` + w[large_ids[2]] * `S32-100`
)

comb1 <- ggplot(df, aes(x = time)) +
  geom_line(aes(y = y), col = "grey") +
  geom_line(aes(y = combo_1), color = clrs[1], linewidth = 1) +
  geom_line(aes(y = combo_1), color = clrs[3], linetype = "dotted", linewidth = 1) +
  theme

comb2 <- ggplot(df, aes(x = time)) +
  geom_line(aes(y = y), col = "grey") +
  geom_line(aes(y = combo_2), color = clrs[1], linewidth = 1) +
  geom_line(aes(y = combo_2), color = clrs[4], linetype = "dotted", linewidth = 1) +
  theme

comb3 <- ggplot(df, aes(x = time)) +
  geom_line(aes(y = y), color = "grey") +
  geom_line(aes(y = combo_3), color = clrs[3], linewidth = 1) +
  geom_line(aes(y = combo_3), color = clrs[5], linetype = "dotted", linewidth = 1) +
  theme

panel2 <- (comb1 / comb2 / comb3)

(panel1 | panel2) %>%
  ggsave(
    here::here("Figures/concept_fig_fourier.png"),
    plot = .,
    device = "png",
    width = 6,
    height = 8,
    units = "in"
  )

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






