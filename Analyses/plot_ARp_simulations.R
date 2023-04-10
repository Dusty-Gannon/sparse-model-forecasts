# Plot summary data from AR-p beta-p simulation runs
# These simulations iterate through different time series lengths
# and different standard deviations for the random innovations
# where there are sparse covariates, lags in one covariate, and sparse AR terms
# comparing the fits of regularized and not regularized models

library(tidyverse)
dd <- read_csv('Data/aquatic_sim_data/AR_p_sims_model_output_condensed.csv')

dat <- dd %>%
  # filter(divergent_trans < 20) %>%
  mutate(siglab = factor(paste0('sigma = ', sigma),
                         levels = c('sigma = 10', 'sigma = 1', 'sigma = 0.1'))) %>%
  select(-rmse_sigma, -unconverged_pars, -divergent_trans) %>%
  pivot_longer(cols = starts_with('rmse'),
               names_to = 'parameter', values_to = 'rmse',
               names_pattern = 'rmse_([a-z]+)') %>%
  mutate(parameter = factor(parameter, levels = c('forecast', 'phi', 'beta')))
datr <- filter(dat, model == 'reg')
datnr <- filter(dat, model == 'not_reg')
# not all simulations ran - that is something that needs to be figured out later
p <- datr %>%
  filter(rmse <50) %>%
ggplot(aes(n, rmse)) +
  geom_violin(aes(group = cut_width(n, 50)), alpha = 0.5, col = 'brown',
              fill =  'brown', position = 'dodge', scale = "width") +
  geom_violin(data = filter(datnr, rmse <50),
              aes(group = cut_width(n, 50)), alpha = 0.5, col = 'grey',
              fill =  'grey', position = 'dodge', scale = "width") +
  facet_grid(siglab~parameter)+
  scale_y_log10()+
  ylab('RMSE') + xlab('Time series length')+
  theme_bw()

col_pal <- c(alpha('grey', 0.5), alpha('brown', 0.5))

dummy_df <- data.frame(Model = c('Regularized', 'Not Regularized'),
                       y = c(1,2))
dummy_plot <- ggplot(dummy_df, aes(Model, y, fill = Model)) +
  geom_bar(stat = 'identity', width = 1) +
  scale_fill_manual(values = col_pal,
                    name = "") +
  theme(legend.direction = 'horizontal')

lgnd <- ggpubr::get_legend(dummy_plot)
pp <- ggpubr::ggarrange(p, legend.grob = lgnd)

png(width = 10.5, height = 9, units = 'in', type = 'cairo', res = 300,
    filename = 'Figures/ARp_sim_fits.png')
    ggpubr::annotate_figure(pp, top = ggpubr::text_grob('AR-p simulation model fits'))
dev.off()


length(which(dat$rmse_forecast >20))/nrow(dat)
# which models ran successfully?
plot(density(dd$rmse_forecast), xlim = c(0, 20))
plot(density(dd$n[dat$sigma == 0.1], na.rm = T))
lines(density(dd$n[dat$sigma == 1]), col = 2)
lines(density(dd$n[dat$sigma == 10]), col = 3)

library(RColorBrewer)
dd %>%
  mutate(sigma = factor(sigma))%>%
  group_by(n, sigma) %>%
  summarize(nsims = n()/2) %>%
  ggplot(aes(factor(n), nsims, fill = sigma))+
  geom_bar(position = 'dodge', stat = 'identity') +
  scale_fill_brewer(palette = 'Blues')+
  ylab('Number of simulations') + xlab('Time series length')+
  theme_minimal()


# quantify divergent transitions
