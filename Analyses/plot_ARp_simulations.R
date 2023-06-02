# Plot summary data from AR-p beta-p simulation runs
# These simulations iterate through different time series lengths
# and different standard deviations for the random innovations
# where there are sparse covariates, lags in one covariate, and sparse AR terms
# comparing the fits of regularized and not regularized models

library(tidyverse)

mod_cols <- c("#a52a2aff", "#33406fff")

dd <- read_csv('Data/aquatic_sim_data/ARp_sims_5_31_condensed.csv')

dd <- dd %>%
  mutate(siglab = factor(paste0('sigma = ', sigma),
                         levels = c('sigma = 0.5', 'sigma = 2', 'sigma = 5')),
         prior = case_when(model == 'not_reg' ~ 'Gaussian',
                           model == 'reg' ~ 'Horseshoe'))

divergent_trans_cap <- 20

# Plot divergent transitions:
total_runs <- dd %>%
  group_by(n, sigma, prior) %>%
  summarize(total_runs = n())

con_runs <- filter(dd, divergent_trans <= divergent_trans_cap) %>%
  group_by(n, sigma, prior) %>%
  summarize(runs = n())

png('Manuscript/Figures/ARp_model_convergence.png',
    width = 8, height = 5, units = 'in', res = 300)
  left_join(con_runs, total_runs) %>%
    mutate(converged_runs = runs/total_runs,
           Sigma = factor(sigma)) %>%
    ggplot(aes(n, converged_runs, col = prior, lty = Sigma)) +
    geom_line(size = 0.9) +
    scale_color_manual('Prior', values = mod_cols) +
    theme_classic() +
    theme(legend.position = c(0.75, 0.25),
          legend.box = 'horizontal')+
    xlim(60,150)+
    ylab('Percent convergence') +
    xlab('Time series length')

dev.off()

# Subsample model fits so that the same number is in each category.
min_size = min(c(con_runs$runs, 50))
dat <- data.frame()
for(i in 1:nrow(total_runs)){
  tmp <- dd %>%
    filter(divergent_trans <= divergent_trans_cap,
           n == total_runs$n[i],
           sigma == total_runs$sigma[i],
           prior == total_runs$prior[i])
  rows <- sample(1:nrow(tmp), min_size)

  dat <- bind_rows(dat, tmp[rows,])
}

# Plot model results:
datr <- filter(dat, prior == 'Horseshoe')
datnr <- filter(dat, prior == 'Gaussian')

frmse <- dat %>%
  ggplot(aes(n, rmse_forecast, col = prior)) +
  geom_point(alpha = 0.3, size = 0.85) +
  geom_smooth(se = FALSE) +
  scale_color_manual('Prior', values = mod_cols) +
  scale_y_log10() +
  facet_grid(.~siglab)+
  ylab('Forecast RMSE')+
  xlab('Time Series Length')+
  theme_classic()+
  theme(panel.border = element_rect(fill = NA),
        panel.spacing = unit(0, 'line'),
        legend.position = c(0.92, 0.83),
        legend.title = element_text(size=8),
        legend.text = element_text(size=7),
        legend.spacing.y = unit(0.05, 'cm')
        # legend.key.size = unit(1, 'line')
        )

png('Manuscript/Figures/ARp_forecast_rmses.png',
    width = 6.5, height = 3.2, units = 'in', res = 300)
    frmse
dev.off()

ggplot(datr, aes(n, rmse_forecast)) +
  geom_violin(aes(group = cut_width(n, 5)), alpha = 0.5,
              fill = mod_cols[1]) +
  geom_violin(data = datnr, aes(group = cut_width(n, 5)), alpha = 0.5,
              fill = mod_cols[2]) +
  geom_smooth(aes(lty = siglab), col = mod_cols[1]) +
  geom_smooth(data = datnr, aes(lty = siglab), col = mod_cols[2]) +
  scale_color_manual(values = mod_cols) +
  scale_y_log10() +
  facet_wrap(siglab~., ncol = 1)+
  ylab('Model Forecast RMSE')+
  xlab('Time series length')+
  theme_classic()

# make a dummy plot to get a legend
dummy_df <- data.frame(Model = c('Regularized', 'Not Regularized'),
                       y = c(1,2))
dummy_plot <- ggplot(dummy_df, aes(Model, y, fill = Model)) +
  geom_bar(stat = 'identity', width = 1, alpha = 0.5) +
  scale_fill_manual(values = mod_cols,
                    name = "") +
  theme(legend.direction = 'horizontal')

lgnd <- ggpubr::get_legend(dummy_plot)
pp <- ggpubr::ggarrange(p, legend.grob = lgnd)

png(width = 10.5, height = 9, units = 'in', type = 'cairo', res = 300,
    filename = 'Figures/ARp_sim_fits2.png')
    ggpubr::annotate_figure(pp, top = ggpubr::text_grob('AR-p simulation model fits'))
dev.off()

dat %>%
  pivot_longer(cols = starts_with(c('beta', 'phi')),
               values_to = 'value',
               names_to = c('parameter', 'rate'),
               names_pattern = '(beta|phi)_([a-z_]+$)') %>%
  pivot_wider(values_from = 'value', names_from = 'rate') %>%
ggplot(aes(n, true_pos, col = model)) +
  geom_point() +
  facet_grid(siglab~parameter) +
  scale_color_manual(values = mod_cols) +
  theme_bw()

dd2 <- dat %>%
  pivot_longer(cols = starts_with(c('beta', 'phi')),
               values_to = 'value',
               names_to = c('parameter', 'rate'),
               names_pattern = '(beta|phi)_([a-z_]+$)') %>%
  pivot_wider(values_from = 'value', names_from = 'rate')

rates <-
  dd2 %>%
  pivot_longer(cols = ends_with('pos'),
               values_to = 'rate',
               names_to = 'category') %>%
  mutate(category = case_when(category == 'true_pos' ~ 'TPR',
                              category == 'false_pos' ~ 'FPR'),
         category = factor(category, levels = c('TPR', 'FPR'))) %>%
  group_by(prior, category, siglab, parameter, n) %>%
  summarize(rate_high = quantile(rate, 0.75, na.rm = T),
            rate_low = quantile(rate, 0.25, na.rm = T),
            rate = median(rate, na.rm = T)) %>%
  ungroup() %>%
  arrange(prior, category, siglab, parameter, n)

r <- rates %>%
  select(-rate_high, -rate_low) %>%
  tidyr::nest(data = c(rate, n)) %>%
  dplyr::mutate(
    m = purrr::map(data, loess,
                   formula = rate ~ n, span = 0.5),
    fitted = purrr::map(m, `[[`, 'fitted')
  )

r <- r %>%
  select(-m) %>%
  tidyr::unnest(cols = c(data, fitted)) %>%
  rename(rate_fit = fitted)

rh <- rates %>%
  select(-rate, -rate_low) %>%
  tidyr::nest(data = c(rate_high, n)) %>%
  dplyr::mutate(
    m = purrr::map(data, loess,
                   formula = rate_high ~ n, span = 0.5),
    fitted = purrr::map(m, `[[`, 'fitted')
  )

rh <- rh %>%
  select(-m) %>%
  tidyr::unnest(cols = c(data,fitted)) %>%
  rename(rate_high_fit = fitted)

rl <- rates %>%
  select(-rate, -rate_high) %>%
  tidyr::nest(data = c(rate_low, n)) %>%
  dplyr::mutate(
    m = purrr::map(data, loess,
                   formula = rate_low ~ n, span = 0.5),
    fitted = purrr::map(m, `[[`, 'fitted')
  )

rl <- rl %>%
  select(-m) %>%
  tidyr::unnest(cols = c(data, fitted)) %>%
  rename(rate_low_fit = fitted)

rr <- inner_join(r, rh) %>% inner_join(rl) %>%
  mutate(across(starts_with('rate'),
                ~case_when(. > 1 ~ 1,
                           . < 0 ~ 0,
                           TRUE ~ .)))

tpr <- rr %>%
  rename(Category = category) %>%
  ggplot(aes(n, rate_fit, color = prior, lty = Category))+
  geom_ribbon(aes(ymin = rate_low_fit, ymax = rate_high_fit, fill = prior),
              alpha = 0.25, color = NA)+
  geom_line(size = 0.75)+
  facet_grid(parameter~siglab) +
  scale_color_manual('Prior', values = mod_cols) +
  scale_fill_manual('Prior', values = mod_cols) +
  ylab('Rate')+
  xlab('Time series length') +
  theme_classic()+
  theme(panel.border = element_rect(fill = NA),
        panel.spacing = unit(0, 'line'),
        legend.position = 'top')
png('Manuscript/Figures/ARp_TPR_FPR.png',
    width = 6.5, height = 5, units = 'in', res = 300)
tpr
dev.off()

ggpubr::ggarrange(tpr, frmse, common.legend = TRUE,
                  nrow = 2, align = 'v', heights = c(2,1))
