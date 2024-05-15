# constructing conceptual figures

library(tidyverse)
library(patchwork)
devtools::load_all()
set.seed(3400)

freq <- 100
n <- 1

N <- freq * n

ts <- basic_timeseries(
  K = 20,
  num_strong = 2,
  n = n,
  freq = freq,
  trend_fraction = 0,
  prob_cycle = 1
)
plot(ts$y, type = "l")

# make y a time series object
y <- stats::ts(ts$y, frequency = 100)


decomp <- forecast::fourier(
  y, K = 50
)

# combine into dataframe

df <- data.frame(
  y = ts$y,
  time = 1:N
)

# get weights using ols
w <- lm(df$y ~ decomp)$coefficients

df <- cbind(df, as.data.frame(decomp))

theme <- theme_classic() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  )

raw <- ggplot(df, aes(x = time)) +
  geom_line(aes(y = y)) +
  theme +
  xlab("")


s1 <- ggplot(df, aes(x = time, y = `S1-100`)) +
  geom_line() +
  theme +
  xlab("") +
  ylab(expression(w[1]))

c1 <- ggplot(df, aes(x = time, y = `C1-100`)) +
  geom_line() +
  theme +
  ylab(expression(w[2])) +
  xlab("")

# now figure out which signals have largest weights
large <- names(sort(abs(w[-1]), decreasing = T)[1:2])
small <- names(sort(abs(w[-1]))[1:2])

large_ids <- which(names(w) %in% large)
small_ids <- which(names(w) %in% small)

s3 <- ggplot(df, aes(x = time, y = `S3-100`)) +
  geom_line() +
  theme +
  ylab(expression(w[6])) +
  xlab("")

s11 <- ggplot(df, aes(x = time, y = `S11-100`)) +
  geom_line() +
  theme +
  ylab(expression(w[22])) +
  xlab("")

s32 <- ggplot(df, aes(x = time, y = `S32-100`)) +
  geom_line() +
  theme +
  ylab(expression(w[64])) +
  xlab("")

c41 <- ggplot(df, aes(x = time, y = `C41-100`)) +
  geom_line() +
  theme +
  ylab(expression(w[83]))


(panel1 <- (raw / s1 / c1 / s3 / s11 / s32 / c41))

# start second panel with linear combos

wids <- which(
  names(w) %in%
    c("decompS1-100", "decompC1-100", "decompS3-100",)
)











