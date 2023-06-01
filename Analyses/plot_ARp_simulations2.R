# Plot summary data from AR-p beta-p simulation runs
# These simulations iterate through different time series lengths
# and different standard deviations for the random innovations
# where there are sparse covariates, lags in one covariate, and sparse AR terms
# comparing the fits of regularized and not regularized models

library(tidyverse)
dd <- read_csv('Data/aquatic_sim_data/ARp_sims_5_29b_condensed.csv')

dd <- dd %>%
  # filter(divergent_trans < 40) %>%
  mutate(siglab = factor(paste0('sigma = ', sigma),
                         levels = c('sigma = 5', 'sigma = 2', 'sigma = 0.5')))

dat <- dd %>%
  filter(divergent_trans <= 20) %>%
  select(-rmse_sigma, -unconverged_pars, -divergent_trans) %>%
  pivot_longer(cols = starts_with('rmse'),
               names_to = 'parameter', values_to = 'rmse',
               names_pattern = 'rmse_([a-z]+)') %>%
  mutate(parameter = factor(parameter, levels = c('forecast', 'phi', 'beta')))%>%
  filter(!is.na(siglab))
datr <- filter(dat, model == 'reg')
datnr <- filter(dat, model == 'not_reg')
# not all simulations ran - that is something that needs to be figured out later
p2 <- datr %>%
  filter(rmse <50) %>%
ggplot(aes(n, rmse)) +
  geom_violin(aes(group = cut_width(n, 60)), alpha = 0.5, col = 'brown',
              fill =  'brown', position = 'dodge', scale = "width") +
  geom_violin(data = filter(datnr, rmse <50),
              aes(group = cut_width(n, 60)), alpha = 0.5, col = 'grey',
              fill =  'grey', position = 'dodge', scale = "width") +
  facet_grid(siglab~parameter)+
  scale_y_log10()+
  ylab('RMSE') + xlab('Time series length')+
  theme_bw()
p <- datr %>%
  filter(rmse <50) %>%
ggplot(aes(n, rmse)) +
  geom_violin(aes(group = factor(n)), alpha = 0.5, col = 'brown',
              fill =  'brown', position = 'dodge', scale = "width") +
  geom_violin(data = filter(datnr, rmse <50),
              aes(group = factor(n)), alpha = 0.5, col = 'grey',
              fill =  'grey', position = 'dodge', scale = "width") +
  facet_grid(siglab~parameter)+
  scale_y_log10()+
  ylab('RMSE') + xlab('Time series length')+
  theme_bw()
datr %>%
  filter(rmse <50,
         parameter == 'forecast') %>%
ggplot(aes(n, rmse)) +
  geom_violin(aes(group = factor(n)), alpha = 0.5, col = 'brown',
              fill =  'brown', position = 'dodge', scale = "width") +
  geom_violin(data = filter(datnr, rmse <50, parameter == 'forecast'),
              aes(group = factor(n)), alpha = 0.5, col = 'grey',
              fill =  'grey', position = 'dodge', scale = "width") +
  facet_grid(siglab~.)+
  scale_y_log10()+
  ylab('RMSE') + xlab('Time series length')+
  theme_bw()
total_runs <- dd %>%
  group_by(n, sigma, model) %>%
  summarize(total_runs = n())
dd %>%
  filter(divergent_trans <= 20) %>%
  group_by(n, sigma, model) %>%
  summarize(runs = n()) %>%
  left_join(total_runs) %>%
  mutate(converged_runs = runs/total_runs) %>%
  ggplot(aes(factor(n), converged_runs, fill = factor(sigma))) +
  geom_bar(stat = 'identity', position = 'dodge') +
  facet_wrap(.~model) + theme_bw() + ylab('Fraction of runs without divergent transitions')

dd %>%
  ggplot(aes(factor(n), log10(divergent_trans), fill = siglab)) +
  geom_boxplot() + facet_wrap(.~model)


  filter(rmse <50,
         parameter == 'forecast') %>%
ggplot(aes(n, rmse)) +
  geom_violin(aes(group = factor(n)), alpha = 0.5, col = 'brown',
              fill =  'brown', position = 'dodge', scale = "width") +
  geom_violin(data = filter(datnr, rmse <50, parameter == 'forecast'),
              aes(group = factor(n)), alpha = 0.5, col = 'grey',
              fill =  'grey', position = 'dodge', scale = "width") +
  facet_grid(siglab~.)+
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
    filename = 'Figures/ARp_sim_fits2.png')
    ggpubr::annotate_figure(pp, top = ggpubr::text_grob('AR-p simulation model fits'))
dev.off()

dd %>%
  pivot_longer(cols = starts_with(c('beta', 'phi')),
               values_to = 'value',
               names_to = c('parameter', 'rate'),
               names_pattern = '(beta|phi)_([a-z_]+$)') %>%
  pivot_wider(values_from = 'value', names_from = 'rate') %>%
ggplot(aes(n, true_pos, col = model)) +
  geom_point() +
  facet_grid(siglab~parameter) +
  scale_color_manual(values = col_pal) +
  theme_bw()

dd2 <- dd %>%
  filter(divergent_trans <=20)%>%
  pivot_longer(cols = starts_with(c('beta', 'phi')),
               values_to = 'value',
               names_to = c('parameter', 'rate'),
               names_pattern = '(beta|phi)_([a-z_]+$)') %>%
  pivot_wider(values_from = 'value', names_from = 'rate')
ddsum <- dd2 %>% group_by(n, sigma, model, parameter) %>%
  summarize(across(starts_with(c('true', 'false')), mean))
dd2r <- filter(dd2, model == 'reg')
dd2nr <- filter(dd2, model == 'not_reg')
tp <- ggplot(dd2r, aes(n, true_pos)) +
  geom_violin(aes(group = factor(n)), col = 'brown', fill = 'brown', alpha = 0.5) +
  geom_violin(data = dd2nr, aes(group =  factor(n)), fill = 'grey', alpha = 0.5) +
  # geom_point(data = ddsum, aes(col = model))+
  # scale_color_manual(values = c('grey', 'brown')) +
  facet_grid(siglab~parameter) +
  theme_bw()
tp <- ggplot(dd2, aes(factor(n), true_pos)) +
  geom_boxplot(aes( fill = model), alpha = 0.5, outlier.color = NA) +
  facet_grid(siglab~parameter) +
  scale_fill_manual(values = c('brown', 'grey'))+
  ylab('True Positive Rate')+
  theme_bw()
fp <- ggplot(dd2, aes(factor(n), false_pos)) +
  geom_boxplot(aes( fill = model), alpha = 0.5, outlier.colour = NA) +
  facet_grid(siglab~parameter) +
  scale_fill_manual(values = c('brown', 'grey'))+
  ylab('False Positive Rate')+
  theme_bw()

ggpubr::ggarrange(tp, fp, common.legend = TRUE)
ggplot(dd2, aes(factor(n), false_neg)) +
  geom_boxplot(aes( fill = model), alpha = 0.5) +
  facet_grid(siglab~parameter) +
  scale_fill_manual(values = c('brown', 'grey'))+
  ylab('False Negative Rate')+
  theme_bw()
tn <- ggplot(dd2r, aes(n, true_neg)) +
  geom_violin(aes(group = cut_width(n, 50)), col = 'brown', fill = 'brown', alpha = 0.5) +
  geom_violin(data = dd2nr, aes(group = cut_width(n, 50)), fill = 'grey', alpha = 0.5) +
  facet_grid(siglab~parameter) +
  scale_color_manual(values = col_pal) +
  theme_bw()
fn <- ggplot(dd2r, aes(n, false_neg)) +
  geom_violin(aes(group =  factor(n)), col = 'brown', fill = 'brown', alpha = 0.5) +
  geom_violin(data = dd2nr, aes(group =  factor(n)), fill = 'grey', alpha = 0.5) +
  facet_grid(siglab~parameter) +
  scale_color_manual(values = col_pal) +
  theme_bw()
fp <- ggplot(dd2r, aes(n, false_pos)) +
  geom_violin(aes(group =  factor(n)), col = 'brown', fill = 'brown', alpha = 0.5) +
  geom_violin(data = dd2nr, aes(group =  factor(n)), fill = 'grey', alpha = 0.5) +
  facet_grid(siglab~parameter) +
  scale_color_manual(values = col_pal) +
  theme_bw()

ggpubr::ggarrange(tp, fp, fn)

# which models ran successfully?
library(RColorBrewer)
dd %>%
  filter(divergent_trans < 20) %>%
  mutate(sigma = factor(sigma))%>%
  group_by(n, sigma, model) %>%
  summarize(nsims = n()) %>%
  ggplot(aes(factor(n), nsims/300, fill = sigma))+
  geom_bar(position = 'dodge', stat = 'identity') +
  scale_fill_brewer(palette = 'Blues')+
  facet_wrap(.~model)+
  ylab('Number of simulations') + xlab('Time series length')+
  theme_minimal()


# quantify divergent transitions
dd %>% group_by(n, sigma, model) %>%
  filter(divergent_trans < 20) %>%
  summarize(min_div = min(divergent_trans),
    median_div = median(divergent_trans),
            count = n()) %>%
  arrange(count)
dd %>%
  ggplot(aes(n, divergent_trans)) +
  geom_violin(aes(group = c(cut_width(n, 50))), scale = 'width') +
  # geom_violin(aes(group = cut_width(n, 50)), fill = 'grey', alpha = 0.5) +
  facet_grid(siglab~model) +
  scale_color_manual(values = col_pal) +
  theme_bw() +ylim(0,500)

plot(density(dd$divergent_trans[dd$n == 50 & dd$sigma == 0.5]))
