################################################
###### Figures for AIC vs STAN comparison ######
################################################

library(here)
library(tidyr)
library(dplyr)
library(patchwork)
library(ggplot2)

################################################
# Figure 1: (Deprecated)
################################################

# source(here("Simulations/AICvsStanRMSE.R"))
#
# set.seed(3782309)
# n=100
# K=50
#
# timeSeries1=getTS(numTrials = 1, n = n, K = K, num_strong = 5,prob_cycle = 0.5, trend_fraction = 0.5,freq=1,sigma=0.5,probWeakCorr=0.8,numStrongCorr=3,strongSelf=T,corrLevel=0.5,corrChange=F, propChange=0.75, changeSize=0, changeTimeVar=0)
# cleanSeries1=cleanTS(timeSeries1[[1]])
# testSeries1=splitTS(cleanSeries1,set="test",n=n,nfit=round(0.6*n))
# trainSeries1=splitTS(cleanSeries1,set="train",n=n,nfit=round(0.6*n))
#
# AICseries1=AICselect(trainSeries1)
# STANseries1=STANselect(trainSeries1,testSeries1,nfit=60,n=100,K=50)
# AICpred=predict(AICseries1,testSeries1,se.fit = T)
#
#
# par(mfrow=c(3,1))
# par(mar=c(4,4,1,1))
#
# #95% confidence intervals for AIC
# plot(0,0,type="n",ylim=c(-6,4),xlim=c(0,100),xlab="Time",ylab="y (stepwise AIC)")
# polygon(c(61:100,100:61),c(AICpred$fit+1.96*AICpred$se.fit,rev(AICpred$fit-1.96*AICpred$se.fit)),col="lightblue",border=NA)
# lines(timeSeries1[[1]]$y,type="l",xlim=c(0,100))
# lines(61:100,AICpred$fit,col="blue",xlim=c(0,100),type="l")
#
# #95% confidence intervals for STAN
# STANfit=apply(STANgetpredict(STANseries1),2,median)
# STAN975=apply(STANgetpredict(STANseries1),2,quantile,0.975)
# STAN025=apply(STANgetpredict(STANseries1),2,quantile,0.025)
#
# plot(0,0,type="n",ylim=c(-6,4),xlim=c(0,100),xlab="Time",ylab="y (Horseshoe)")
# polygon(c(61:100,100:61),c(STAN975,rev(STAN025)),col="deeppink",border=NA)
# lines(timeSeries1[[1]]$y,type="l",xlim=c(0,100))
# lines(61:100,STANfit,col="deeppink4",xlim=c(0,100),type="l")
#
# # coefficient recovery plots
# plot(timeSeries1[[1]]$beta,xlab="parameter",ylab="value",pch=0,cex=2)
# abline(h=0)
#
# #get AIC coefficients in order
# coefAIC=numeric(51)
# seAIC=numeric(51)
#
# coefAIC[1]=summary(AICseries1)[4]$coefficients[1,1]
# seAIC[1]=summary(AICseries1)[4]$coefficients[1,2]
# for(i in 1:50){
#
#   finder=which(rownames(summary(AICseries1)[4]$coefficients)==paste0("driver_",i))
#   if(length(finder)==0){
#     coefAIC[i]=0
#     seAIC[i]=0
#   } else {
#     coefAIC[i]=summary(AICseries1)[4]$coefficients[finder,1]
#     seAIC[i]=summary(AICseries1)[4]$coefficients[finder,2]
#   }
#
# }
#
# # plot AIC coefficients
# points(x=1:51-0.2,coefAIC,pch=16,col="blue")
# segments(x0=1:51-0.2,y0=coefAIC-1.96*seAIC,y1=coefAIC+1.96*seAIC,col="blue")
#
# # get STAN coefficients in shape
#
#
# betapost=extract(STANseries1, pars = "beta")$beta
# means = apply(betapost, 2, mean)
#
# points(x=2:51+0.2,means,col="deeppink",pch=16)
# segments(x0=2:51+0.2,y0=apply(betapost, 2, quantile, probs = 0.025),y1=apply(betapost, 2, quantile, probs = 0.975),col="deeppink")


simData <- read.csv(here("Simulations/AICvsStanResults.csv"))

colors <- PNWColors::pnw_palette("Cascades", n = 4)
fills <- adjustcolor(colors, alpha.f = 0.6)

#
# ---- confusion plots ----
#

## ---- Option 1: Scatter of TPR vs TNR ----

conf_res <- simData %>%
  filter(corrLevel == 0.5) %>%
  dplyr::select(
    contains("_TPR") | contains("_TNR")
  )

conf_sub <- conf_res %>% select(
  contains("AIC") | contains("STAN")
) %>%
  rename(
    "RHS_TNR" = "STAN_beta_TNR",
    "RHS_TPR" = "STAN_beta_TPR"
  )

# now stack the columns
conf_long <- tibble(
  method = rep(c("stepwise AIC", "RHS"), each = nrow(conf_sub)),
  tpr = c(conf_sub$AIC_TPR, conf_sub$RHS_TPR),
  tnr = c(conf_sub$AIC_TNR, conf_sub$RHS_TNR)
)

ggplot(conf_long, aes(x = tnr, y = tpr, color = method)) +
  geom_jitter(shape = 1, width = 0.01, height = 0.03, alpha = 0.2) +
  theme_bw() +
  scale_color_manual(values = colors[c(2,4)]) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  theme(legend.title = element_blank()) +
  xlab("Rate of excluding weak variables") +
  ylab("Rate of including strong variables") +
  coord_fixed()


ggsave(
  here("Figures/TPRvsTNR_tradeoff.pdf"),
  height = 4, width = 5,
  device = "pdf",
  dpi = 300
)


## ---- Option 2: violins ----

# pivot long
conf_long2 <- conf_res %>%
  pivot_longer(
    cols = everything(),
    values_to = "rate",
    names_to = "group"
  ) %>%
  mutate(
    group = case_when(
      group == "STAN_beta_TNR" ~ "RHS_TNR",
      group == "STAN_beta_TPR" ~ "RHS_TPR",
      .default = group
    )
  ) %>%
  separate_wider_delim(
    cols = group,
    delim = "_",
    names = c("method", "metric")
  ) %>%
  filter(
    method != "GLM" & method != "modelAvg"
  ) %>%
  mutate(
    method = factor(
      method,
      levels = c("Full model", "AIC", "AIC model avg", "RHS")
    )
  )


pdf(file=here("Figures/confusion_rates_AICvRHS.pdf"), width = 5, height = 6)

ggplot(conf_long2, aes(x = method, y = rate)) +
  facet_wrap(vars(metric), ncol = 1) +
  geom_violin(aes(fill = method, color = method), adjust = 2) +
  stat_summary(
    color = "grey3",
    geom = "errorbar",
    fun.min = \(x){quantile(x, probs = 0.025)},
    fun.max = \(x){quantile(x, probs = 0.975)},
    width = 0,
    linewidth = 0.5
  ) +
  stat_summary(
    color = "grey3",
    geom = "errorbar",
    fun = mean,
    fun.min = \(x){quantile(x, probs = 0.1)},
    fun.max = \(x){quantile(x, probs = 0.9)},
    linewidth = 1,
    width = 0
  ) +
  stat_summary(
    geom = "point",
    aes(color = method),
    fun = mean,
    size = 1.5
  ) +
  stat_summary(
    geom = "point",
    fun = mean,
    color = "white",
    size = 0.7
  ) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = fills[c(1, 4)]) +
  scale_color_manual(values = colors[c(1, 4)]) +
  theme(
    legend.position = "none"
  )

dev.off()



#par(mfrow=c(3,1))

# ## ---- TPR plot ----
#
# vioplot::vioplot(cbind(
#   AIC1=corData6$TPRaic[which(corData6$corrLevel==0.1)],
#   STAN1=corData6$TPRstan[which(corData6$corrLevel==0.1)],
#   AIC5=corData6$TPRaic[which(corData6$corrLevel==0.5)],
#   STAN5=corData6$TPRstan[which(corData6$corrLevel==0.5)],
#   AIC9=corData6$TPRaic[which(corData6$corrLevel==0.9)],
#   STAN9=corData6$TPRstan[which(corData6$corrLevel==0.9)]
# ),
#   xlab="",
#   horizontal=T,
#   las=1,
#   names=c(
#     "Stepwise AIC\n correlation 0.1",
#     "RHS\n correlation 0.1",
#     "Stepwise AIC\n correlation 0.5",
#     "RHS\n correlation 0.5",
#     "Stepwise AIC\n correlation 0.9",
#     "RHS\n correlation 0.9"
#   ),
#   col=colors,
#   pchMed=20,
#   border=rep(c("tan4","royalblue4"), 3),
#   rectCol=rep(c("tan4","royalblue4"), 3),
#   lineCol=rep(c("tan4","royalblue4"), 3),
#   colMed=rep(c("tan4","royalblue4"), 3),
#   ylab=""
# )
# title(ylab="Density of TPR",line=7,cex.lab=1.2)
# title(xlab="True Positive Rate (TPR)",line=3,cex.lab=1.2)
# par(xpd=T)
# text(-0.16,7,"a)",cex=1.5)
# par(xpd=F)

# ## ---- TNR plot ----
#
# vioplot::vioplot(cbind(
#   AIC1=corData6$TNRaic[which(corData6$corrLevel==0.1)],
#   STAN1=corData6$TNRstan[which(corData6$corrLevel==0.1)],
#   AIC5=corData6$TNRaic[which(corData6$corrLevel==0.5)],
#   STAN5=corData6$TNRstan[which(corData6$corrLevel==0.5)],
#   AIC9=corData6$TNRaic[which(corData6$corrLevel==0.9)],
#   STAN9=corData6$TNRstan[which(corData6$corrLevel==0.9)]
# ),
#   xlab="",
#   horizontal=T,
#   las=1,
#   names=c(
#     "Stepwise AIC\n correlation 0.1",
#     "RHS\n correlation 0.1",
#     "Stepwise AIC\n correlation 0.5",
#     "RHS\n correlation 0.5",
#     "Stepwise AIC\n correlation 0.9",
#     "RHS\n correlation 0.9"
#   ),
#   col=c("wheat","skyblue","tan","steelblue1","tan3","steelblue4"),
#   pchMed=20,
#   border=rep(c("tan4","royalblue4"), 3),
#   rectCol=rep(c("tan4","royalblue4"), 3),
#   lineCol=rep(c("tan4","royalblue4"), 3),
#   colMed=rep(c("tan4","royalblue4"), 3),
#   ylab=""
# )
# title(ylab="Density of TNR",line=7,cex.lab=1.2)
# title(xlab="True Negative Rate (TNR)",line=3,cex.lab=1.2)
# par(xpd=T)
# text(-0.275,7,"b)",cex=1.5)
# par(xpd=F)

# ---- RMSE plot ----

rmse_res <- simData %>% filter(
  corrLevel == 0.5
) %>% dplyr::select(
  contains("RMSE")
)

# pivot long
rmse_long <- rmse_res %>% pivot_longer(
  cols = everything(),
  values_to = "RMSE",
  names_to = "method"
) %>% mutate(
  method = case_when(
    method == "RMSE_AIC" ~ "AIC",
    method == "RMSE_STAN" ~ "RHS",
    method == "RMSE_GLM" ~ "Full model",
    method == "RMSE_modelAvg" ~ "AIC model avg",
    .default = NA
  )
) %>% mutate(
  method = factor(
    method,
    levels = c("Full model", "AIC", "AIC model avg", "RHS")
  )
)

pdf(file = here("Figures/pred_rmse_corr_p5.pdf"), width = 5, height = 4)

ggplot(rmse_long, aes(x = RMSE, y = method)) +
  geom_violin(aes(color = method, fill = method), adjust = 1.5) +
  stat_summary(
    color = "grey3",
    geom = "errorbar",
    fun.min = \(x){quantile(x, probs = 0.025)},
    fun.max = \(x){quantile(x, probs = 0.975)},
    width = 0,
    linewidth = 0.5
  ) +
  stat_summary(
    color = "grey3",
    geom = "errorbar",
    fun = mean,
    fun.min = \(x){quantile(x, probs = 0.1)},
    fun.max = \(x){quantile(x, probs = 0.9)},
    linewidth = 1,
    width = 0
  ) +
  stat_summary(
    geom = "point",
    aes(color = method),
    fun = mean,
    size = 1.5
  ) +
  stat_summary(
    geom = "point",
    fun = mean,
    color = "white",
    size = 0.7
  ) +
  theme_bw() +
  scale_fill_manual(values = fills) +
  scale_color_manual(values = colors) +
  theme(
    legend.position = "none"
  ) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  xlab("Prediction RMSE")


dev.off()

# {
# plot(0, 0, type="n",
#      xlim=c(0, 4), ylim=c(0.5, 9.5),
#      xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
#
# rect(xleft=par("usr")[1], ybottom=0.5, xright=par("usr")[2], ytop=3.5,
#      col="grey92", border=NA)
# rect(xleft=par("usr")[1], ybottom=6.5, xright=par("usr")[2], ytop=9.5,
#      col="grey92", border=NA)
#
# vioplot::vioplot(cbind(
#   GLM1=corData6$RMSEglm[which(corData6$corrLevel==0.1)],
#   AIC1=corData6$RMSEaic[which(corData6$corrLevel==0.1)],
#   STAN1=corData6$RMSEstan[which(corData6$corrLevel==0.1)],
#   GLM5=corData6$RMSEglm[which(corData6$corrLevel==0.5)],
#   AIC5=corData6$RMSEaic[which(corData6$corrLevel==0.5)],
#   STAN5=corData6$RMSEstan[which(corData6$corrLevel==0.5)],
#   GLM9=corData6$RMSEglm[which(corData6$corrLevel==0.9)],
#   AIC9=corData6$RMSEaic[which(corData6$corrLevel==0.9)],
#   STAN9=corData6$RMSEstan[which(corData6$corrLevel==0.9)]
# ),
#   xlab="",
#   horizontal=T,
#   names=rep("", 9),
#   las=1,
#   col=colors,
#   pchMed=20,
#   add=T,
#   border=rep(c("darkgreen","tan4","royalblue4"), 3),
#   rectCol=rep(c("darkgreen","tan4","royalblue4"), 3),
#   lineCol=rep(c("darkgreen","tan4","royalblue4"), 3),
#   colMed=rep(c("darkgreen","tan4","royalblue4"), 3),
#   ylab=""
# )
# title(ylab="Correlation level",line=5,cex.lab=1.2)
# title(xlab="Prediction Root Mean Square Error (RMSE)",line=3,cex.lab=1.2)
#
# # correcting the axis labels
# axis(1,at=seq(0,4,by=1),labels=seq(0,4,by=1))
#
# axis(2, at = c(2, 5, 8),
#      labels = c(
#        expression(rho == 0.1),
#        expression(rho == 0.5),
#        expression(rho == 0.9)
#      ),
#      las=2)
#
# legend("topright",
#        legend = c("Full model", "Stepwise AIC", "RHS"),
#        fill = c("green4", "tan3", "steelblue4"),
#        border = c("darkgreen", "tan4", "royalblue4"),
#        bty = "n",
#        inset = c(0.01, 0.01))
#
# dev.off()
# }

# ---- Coverage results ----

## ---- Prediction interval coverage and width ----
predint_long <- simData %>%
  filter(corrLevel == 0.5) %>%
  select(contains("PI", ignore.case = F)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "group",
    values_to = "value"
  ) %>%
  separate_wider_delim(
    cols = "group",
    delim = "_",
    names = c("method", "metric"),
    too_few = "align_start"
  ) %>%
  pivot_wider(
    names_from = metric,
    values_from = value,
    values_fn = list
  ) %>%
  unnest_longer(
    col = c(PIcoverage, PIwidth)
  ) %>%
  mutate(
    method = case_when(
      method == "GLM" ~ "Full model",
      method == "MAvg" ~ "AIC model avg",
      method == "Stan" ~ "RHS",
      .default = method
    )
  ) %>%
  mutate(
    method = factor(method, levels = c("Full model", "AIC", "AIC model avg", "RHS"))
  )

predcov_plot <- ggplot(predint_long, aes(x = PIcoverage, y = method, color = method)) +
  stat_summary(
    geom = "pointrange",
    fun = mean,
    fun.min = \(x) quantile(x, probs = 0.025),
    fun.max = \(x) quantile(x, probs = 0.975)
  ) +
  geom_vline(xintercept = 0.95, linetype = "dashed") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    plot.title.position = "panel",
    plot.title = element_text(hjust = -0.2),
    plot.subtitle = element_text(hjust = 0)
  ) +
  scale_color_manual(values = colors) +
  xlab("Rate") +
  ylab("") +
  ggtitle("a)", subtitle = "Prediction interval coverage")

predwidth_plot <- ggplot(predint_long, aes(x = PIwidth/0.5, y = method, color = method)) +
  stat_summary(
    geom = "pointrange",
    fun = mean,
    fun.min = \(x) quantile(x, probs = 0.025),
    fun.max = \(x) quantile(x, probs = 0.975)
  ) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    plot.title.position = "panel",
    plot.title = element_text(hjust = -0.2),
    plot.subtitle = element_text(hjust = 0)
  ) +
  scale_color_manual(values = colors) +
  xlab("Units of irreducible error") +
  ylab("") +
  ggtitle("b)", subtitle = "Prediction interval width")

predint_combined_plot <- predcov_plot + predwidth_plot + plot_layout(guides = "collect")

ggsave(
  filename = here("Figures/prediction_interval_res_sim1.pdf"),
  plot = predint_combined_plot,
  device = "pdf",
  height = 4,
  width = 8,
  dpi = 300,
  units = "in"
)




## ---- coefficient CI coverage and width (appendix) ----

ci_cov_long <- simData %>%
  filter(corrLevel == 0.5) %>%
  select(matches("^(AIC|GLM|Stan|MAvg)_(coverage|width)(_escape|_noescape)?$")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "group",
    values_to = "value"
  ) %>%
  separate_wider_delim(
    cols = "group",
    delim = "_",
    names = c("method", "metric", "type"),
    too_few = "align_start"
  ) %>%
  pivot_wider(
    names_from = metric,
    values_from = value,
    values_fn = list
  ) %>%
  unnest_longer(
    col = c(coverage, width)
  ) %>%
  mutate(
    method = case_when(
      method == "GLM" ~ "Full model",
      method == "MAvg" ~ "AIC model avg",
      method == "Stan" ~ "RHS",
      .default = method
    )
  ) %>%
  mutate(
    method = factor(method, levels = c("Full model", "AIC", "AIC model avg", "RHS")),
    type2 = case_when(
      is.na(type) ~ "",
      type == "escape" ~ "(active)",
      type == "noescape" ~ "(not active)",
      .default = ""
    ),
    type = case_when(
      is.na(type) ~ "escape",
      .default = type
    ),
    grp = paste(method, type2) %>%
      factor(levels = c(
        "Full model ", "AIC ",
        "AIC model avg (active)", "AIC model avg (not active)",
        "RHS (active)", "RHS (not active)"
      ))
  )


ci_cov_plot <- ggplot(ci_cov_long, aes(x = coverage, y = grp)) +
  stat_summary(
    aes(color = grp, shape = type),
    geom = "pointrange",
    fun = mean,
    fun.min = \(x){quantile(x, probs = 0.025)},
    fun.max = \(x){quantile(x, probs = 0.975)}
  ) +
  geom_vline(xintercept = 0.95, linetype = "dashed") +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Rate") +
  ylab("") +
  scale_color_manual(values = c(colors[1:2], colors[c(3,3)], colors[c(4,4)])) +
  ggtitle("a)", subtitle = "CI coverage") +
  theme(
    plot.title.position = "panel",
    plot.title = element_text(hjust = -0.2),
    plot.subtitle = element_text(hjust = 0)
  )

ci_width_plot <- ggplot(ci_cov_long, aes(x = width, y = grp)) +
  stat_summary(
    aes(color = grp, shape = type),
    geom = "pointrange",
    fun = mean,
    fun.min = \(x){quantile(x, probs = 0.025)},
    fun.max = \(x){quantile(x, probs = 0.975)}
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Width") +
  ylab("") +
  scale_color_manual(values = c(colors[1:2], colors[c(3,3)], colors[c(4,4)])) +
  ggtitle("b)", subtitle = "CI width") +
  theme(
    plot.title.position = "panel",
    plot.title = element_text(hjust = -0.2),
    plot.subtitle = element_text(hjust = 0)
  )

ci_combined_plot <- ci_cov_plot + ci_width_plot

ggsave(
  filename = here("Figures/coefficient_CI_res_sim1.pdf"),
  plot = ci_combined_plot,
  device = "pdf",
  height = 4,
  width = 8,
  dpi = 300,
  units = "in"
)


# ---- Covariate shift results ----

decorData <- read.csv(here("Simulations/DecorrResults.csv"))

required_decor_cols <- c("RMSE_AIC", "RMSE_GLM", "RMSE_STAN", "RMSE_modelAvg")
if (!all(required_decor_cols %in% names(decorData))) {
  stop(
    "Simulations/DecorrResults.csv is missing ",
    paste(setdiff(required_decor_cols, names(decorData)), collapse = ", "),
    " — it predates the standardization/model-averaging changes to STANselect()/AICselect(). ",
    "Rerun Simulations/run_Decorr.sh and Simulations/run_combine_Decorr.sh on Beartooth to regenerate it."
  )
}

decor_rmse_long <- decorData %>%
  dplyr::select(corrLevel, corrChange, contains("RMSE")) %>%
  pivot_longer(
    cols = contains("RMSE"),
    values_to = "RMSE",
    names_to = "method"
  ) %>%
  mutate(
    method = case_when(
      method == "RMSE_AIC" ~ "AIC",
      method == "RMSE_STAN" ~ "RHS",
      method == "RMSE_GLM" ~ "Full model",
      method == "RMSE_modelAvg" ~ "AIC model avg",
      .default = NA
    ),
    method = factor(
      method,
      levels = c("Full model", "AIC", "AIC model avg", "RHS")
    ),
    shift = if_else(corrChange, "Covariate shift", "No shift"),
    shift = factor(shift, levels = c("Covariate shift", "No shift"))
  )

pdf(file = here("Figures/decorrelation_comparison.pdf"), width = 6, height = 5)

decor_rmse_long %>%
  filter(corrLevel == 0.5) %>%
  ggplot(., aes(x = RMSE, y = method)) +
    facet_wrap(
      vars(shift),
      ncol = 1
      #labeller = labeller(corrLevel = as_labeller(function(x) paste0("rho == ", x), label_parsed))
    ) +
    geom_violin(aes(color = method, fill = method), adjust = 1.5) +
    stat_summary(
      color = "grey3",
      geom = "errorbar",
      fun.min = \(x){quantile(x, probs = 0.025)},
      fun.max = \(x){quantile(x, probs = 0.975)},
      width = 0,
      linewidth = 0.5
    ) +
    stat_summary(
      color = "grey3",
      geom = "errorbar",
      fun = mean,
      fun.min = \(x){quantile(x, probs = 0.1)},
      fun.max = \(x){quantile(x, probs = 0.9)},
      linewidth = 1,
      width = 0
    ) +
    stat_summary(
      geom = "point",
      aes(color = method),
      fun = mean,
      size = 1.5
    ) +
    stat_summary(
      geom = "point",
      fun = mean,
      color = "white",
      size = 0.7
    ) +
    geom_vline(xintercept = 0.5, linetype = "dashed") +
    theme_bw() +
    scale_fill_manual(values = fills) +
    scale_color_manual(values = colors) +
    theme(
      legend.position = "none"
    ) +
    xlab("Prediction RMSE") +
    ylab("")

dev.off()


# ---- Covariate shift: coverage results ----

decor_shift <- decorData %>%
  filter(corrLevel == 0.5) %>%
  mutate(
    shift = if_else(corrChange, "Covariate shift", "No shift"),
    shift = factor(shift, levels = c("Covariate shift", "No shift"))
  )

## ---- Prediction interval coverage and width ----
predint_decor_long <- decor_shift %>%
  select(shift, contains("PI", ignore.case = F)) %>%
  pivot_longer(
    cols = -shift,
    names_to = "group",
    values_to = "value"
  ) %>%
  separate_wider_delim(
    cols = "group",
    delim = "_",
    names = c("method", "metric"),
    too_few = "align_start"
  ) %>%
  pivot_wider(
    names_from = metric,
    values_from = value,
    values_fn = list
  ) %>%
  unnest_longer(
    col = c(PIcoverage, PIwidth)
  ) %>%
  mutate(
    method = case_when(
      method == "GLM" ~ "Full model",
      method == "MAvg" ~ "AIC model avg",
      method == "Stan" ~ "RHS",
      .default = method
    )
  ) %>%
  mutate(
    method = factor(method, levels = c("Full model", "AIC", "AIC model avg", "RHS"))
  )

predcov_decor_plot <- ggplot(predint_decor_long, aes(x = PIcoverage, y = method, color = method)) +
  stat_summary(
    geom = "pointrange",
    fun = mean,
    fun.min = \(x) quantile(x, probs = 0.025),
    fun.max = \(x) quantile(x, probs = 0.975)
  ) +
  facet_wrap(vars(shift), ncol = 1) +
  geom_vline(xintercept = 0.95, linetype = "dashed") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    plot.title.position = "panel",
    plot.title = element_text(hjust = -0.2),
    plot.subtitle = element_text(hjust = 0)
  ) +
  scale_color_manual(values = colors) +
  xlab("Rate") +
  ylab("") +
  ggtitle("a)", subtitle = "Prediction interval coverage")

predwidth_decor_plot <- ggplot(predint_decor_long, aes(x = PIwidth/0.5, y = method, color = method)) +
  stat_summary(
    geom = "pointrange",
    fun = mean,
    fun.min = \(x) quantile(x, probs = 0.025),
    fun.max = \(x) quantile(x, probs = 0.975)
  ) +
  facet_wrap(vars(shift), ncol = 1) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    plot.title.position = "panel",
    plot.title = element_text(hjust = -0.2),
    plot.subtitle = element_text(hjust = 0)
  ) +
  scale_color_manual(values = colors) +
  xlab("Units of irreducible error") +
  ylab("") +
  ggtitle("b)", subtitle = "Prediction interval width")

predint_decor_combined_plot <- predcov_decor_plot + predwidth_decor_plot + plot_layout(guides = "collect")

ggsave(
  filename = here("Figures/prediction_interval_res_decorr.pdf"),
  plot = predint_decor_combined_plot,
  device = "pdf",
  height = 6,
  width = 8,
  dpi = 300,
  units = "in"
)


## ---- coefficient CI coverage and width (appendix) ----

# this gives identical results since the effect sizes don't change,
# just the correlations between them

ci_cov_decor_long <- decor_shift %>%
  select(shift, matches("^(AIC|GLM|Stan|MAvg)_(coverage|width)(_escape|_noescape)?$")) %>%
  pivot_longer(
    cols = -shift,
    names_to = "group",
    values_to = "value"
  ) %>%
  separate_wider_delim(
    cols = "group",
    delim = "_",
    names = c("method", "metric", "type"),
    too_few = "align_start"
  ) %>%
  pivot_wider(
    names_from = metric,
    values_from = value,
    values_fn = list
  ) %>%
  unnest_longer(
    col = c(coverage, width)
  ) %>%
  mutate(
    method = case_when(
      method == "GLM" ~ "Full model",
      method == "MAvg" ~ "AIC model avg",
      method == "Stan" ~ "RHS",
      .default = method
    )
  ) %>%
  mutate(
    method = factor(method, levels = c("Full model", "AIC", "AIC model avg", "RHS")),
    type2 = case_when(
      is.na(type) ~ "",
      type == "escape" ~ "(active)",
      type == "noescape" ~ "(not active)",
      .default = ""
    ),
    type = case_when(
      is.na(type) ~ "escape",
      .default = type
    ),
    grp = paste(method, type2) %>%
      factor(levels = c(
        "Full model ", "AIC ",
        "AIC model avg (active)", "AIC model avg (not active)",
        "RHS (active)", "RHS (not active)"
      ))
  )

ci_cov_decor_plot <- ggplot(ci_cov_decor_long, aes(x = coverage, y = grp)) +
  stat_summary(
    aes(color = grp, shape = type),
    geom = "pointrange",
    fun = mean,
    fun.min = \(x){quantile(x, probs = 0.025)},
    fun.max = \(x){quantile(x, probs = 0.975)}
  ) +
  facet_wrap(vars(shift), ncol = 1) +
  geom_vline(xintercept = 0.95, linetype = "dashed") +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Rate") +
  ylab("") +
  scale_color_manual(values = c(colors[1:2], colors[c(3,3)], colors[c(4,4)])) +
  ggtitle("a)", subtitle = "CI coverage") +
  theme(
    plot.title.position = "panel",
    plot.title = element_text(hjust = -0.2),
    plot.subtitle = element_text(hjust = 0)
  )

ci_width_decor_plot <- ggplot(ci_cov_decor_long, aes(x = width, y = grp)) +
  stat_summary(
    aes(color = grp, shape = type),
    geom = "pointrange",
    fun = mean,
    fun.min = \(x){quantile(x, probs = 0.025)},
    fun.max = \(x){quantile(x, probs = 0.975)}
  ) +
  facet_wrap(vars(shift), ncol = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Width") +
  ylab("") +
  scale_color_manual(values = c(colors[1:2], colors[c(3,3)], colors[c(4,4)])) +
  ggtitle("b)", subtitle = "CI width") +
  theme(
    plot.title.position = "panel",
    plot.title = element_text(hjust = -0.2),
    plot.subtitle = element_text(hjust = 0)
  )

ci_decor_combined_plot <- ci_cov_decor_plot + ci_width_decor_plot

ggsave(
  filename = here("Figures/coefficient_CI_res_decorr.pdf"),
  plot = ci_decor_combined_plot,
  device = "pdf",
  height = 6,
  width = 8,
  dpi = 300,
  units = "in"
)
