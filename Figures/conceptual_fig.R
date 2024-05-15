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

df <- cbind(df, as.data.frame(decomp))

raw <- ggplot(df, aes(x = time)) +
  geom_line(aes(y = y)) +
  theme_classic() +
  xlab("")

s1 <- ggplot(df, aes(x = time, y = `S1-100`)) +
  geom_line() +
  theme_classic() +
  xlab("") +
  ylab(expression(w[1]))

c1 <- ggplot(df, aes(x = time, y = `C1-100`)) +
  geom_line() +
  theme_classic() +
  xlab("")

s2 <- ggplot(df, aes(x = time, y = `S2-100`)) +
  geom_line() +
  theme_classic() +
  xlab("")

c2 <- ggplot(df, aes(x = time, y = `C2-100`)) +
  geom_line() +
  theme_classic() +
  xlab("")

s20 <- ggplot(df, aes(x = time, y = `S20-100`)) +
  geom_line() +
  theme_classic() +
  xlab("")

c20 <- ggplot(df, aes(x = time, y = `C20-100`)) +
  geom_line() +
  theme_classic()

(panel1 <- (raw / s1 / c1 / s2 / c2 / s20 / c20))

