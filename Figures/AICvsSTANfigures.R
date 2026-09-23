################################################
###### Figures for AIC vs STAN comparison ######
################################################

library(here)
library(tidyverse)

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


################################################
# Figure 2:
################################################
corData2=read.csv(here("Simulations/AICvsSTANResults.csv"))

# Remove unwanted variables
# fix strongSelf=F
corData3=corData2[-which(corData2$strongSelf),]
# fix numStrongCorr=3
corData4=corData3[-which(corData3$numStrongCorr==1),]
# fix percent of weak covariates to 0.5
corData5=corData4[which(corData4$probWeakCorr==0.5),]
# fix correlations to 0.5 and 0.9
#corData6=corData5[-which(corData5$corrLevel==0.7),]
corData6=corData5

colors <- PNWColors::pnw_palette("Cascades", n = 4)
fills <- adjustcolor(colors, alpha.f = 0.6)

# ---- confusion plots ----

conf_res <- corData6 %>% dplyr::select(
  contains("_TPR") | contains("_TNR")
)

# pivot long
conf_long <- conf_res %>% pivot_longer(
  cols = everything(),
  values_to = "rate",
  names_to = "group"
) %>% mutate(
  group = case_when(
    group == "STAN_beta_TNR" ~ "RHS_TNR",
    group == "STAN_beta_TPR" ~ "RHS_TPR",
    .default = group
  )
) %>% separate_wider_delim(
  cols = group,
  delim = "_",
  names = c("method", "metric")
) %>% filter(
  method != "GLM" & method != "modelAvg"
) %>% mutate(
  method = factor(
    method,
    levels = c("Full model", "AIC", "AIC model avg", "RHS")
  )
)

### ---- load pdf device ----
pdf(file=here("Figures/confusion_rates_AICvRHS.pdf"), width = 5, height = 6)

ggplot(conf_long, aes(x = method, y = rate)) +
  facet_wrap(vars(metric), ncol = 1) +
  geom_violin(aes(fill = method, color = method), adjust = 1.5) +
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

# ---- confusion plots ----

rmse_res <- corData6 %>% filter(
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

coverage <- corData6 %>%
  select(c(corrLevel, contains("coverage"))) %>%
  pivot_longer(
    cols = AIC_coverage:MAvg_coverage_noescape,
    names_to = "method",
    values_to = "rate"
  ) %>%
  separate_wider_delim(
    cols = method,
    delim = "_",
    names = c("method", "what", "type"),
    too_few = "align_start"
  ) %>%
  filter(
    !(method == "MAvg" & type == "noescape")
  ) %>%
  mutate(
    type = case_when(is.na(type) ~ "escape", .default = type),
    method = case_when(
      method == "GLM" ~ "Full model",
      method == "Stan" & type == "escape" ~ "RHS active",
      method == "Stan" & type == "noescape" ~ "RHS inactive",
      method == "MAvg" ~ "AIC model avg",
      .default = method
    )
  ) %>%
  select(!what) %>%
  mutate(
    method = factor(
      method,
      levels = c("Full model", "AIC", "AIC model avg", "RHS active", "RHS inactive")
    )
  )


## ---- load pdf device ----

pdf(file = here("Figures/coverage_results_sim1.pdf"), width = 5, height = 7)

ggplot(coverage, aes(x = rate, y = method)) +
  facet_wrap(
    vars(corrLevel), ncol = 1,
    labeller = label_bquote(rows = rho == .(corrLevel))
  ) +
  stat_summary(
    aes(color = method, shape = type),
    geom = "pointrange",
    fun = mean,
    fun.min = \(x){quantile(x, probs = 0.025)},
    fun.max = \(x){quantile(x, probs = 0.975)}
  ) +
  geom_vline(xintercept = 0.95, linetype = "dashed") +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Coverage") +
  ylab("") +
  scale_color_manual(values = c(colors, colors[4]))

dev.off()


# ---- Covariate shift results ----

decorData=read.csv(file="Simulations/AICvsSTANdecorrelation.csv")
corr0.9_change=which(decorData$corrLevel==0.9&decorData$corrChange)
corr0.9_nochange=which(decorData$corrLevel==0.9&!decorData$corrChange)
corr0.5_change=which(decorData$corrLevel==0.5&decorData$corrChange)
corr0.5_nochange=which(decorData$corrLevel==0.5&!decorData$corrChange)
corr0.1_change=which(decorData$corrLevel==0.1&decorData$corrChange)
corr0.1_nochange=which(decorData$corrLevel==0.1&!decorData$corrChange)

{
pdf(file=here("Figures/decorrelation_comparison.pdf"), width = 7, height = 8)


par(mar=c(3,9,3,1))
par(mfrow=c(3,1))

# plot(1, type = "n", xlim = c(0, 12), ylim = c(0.5, 12.5),
#      xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")

## ---- correlation=0.9 ----
plot(0,0,type="n",ylim=c(0.5,6.5),xlim=c(0.5,40),xlab="",ylab="",xaxt="n",yaxt="n",bty = "n")

rect(xleft=-1.5,ybottom=-0.5,xright=45,ytop=3.5,col="gray",border=F)

vioplot::vioplot(cbind(
  GLMdecorr9=decorData$RMSEglm[corr0.9_change],
  AICdecorr9=decorData$RMSEaic[corr0.9_change],
  STANdecorr9=decorData$RMSEstan[corr0.9_change],
  GLMnorm9=decorData$RMSEglm[corr0.9_nochange],
  AICnorm9=decorData$RMSEaic[corr0.9_nochange],
  STANnorm9=decorData$RMSEstan[corr0.9_nochange]
),
  las=1,
  cex.axis=0.9,
  names=c(
    "GLM,\n decorrelated",
    "Stepwise AIC,\n decorrelated",
    "RHS,\ndecorrelated",
    "GLM,\nno shift",
    "Stepwise AIC,\nno shift",
    "RHS,\nno shift"
  ),
  xlab="",
  horizontal=T,
  ylab="",
  add=T,
  col=rep(colors[7:9], 2),
  border=rep(c("darkgreen","tan4","royalblue4"), 2),
  rectCol=rep(c("darkgreen","tan4","royalblue4"), 2),
  lineCol=rep(c("darkgreen","tan4","royalblue4"), 2),
  colMed=rep(c("darkgreen","tan4","royalblue4"), 2)
)

axis(1,at=seq(0,40,by=5),labels=seq(0,40,by=5))

axis(2, at = 1:6,
     labels = rep(c("Full model", "Stepwise AIC", "RHS"), 2),
     las=2)

text(38, 2, "Decorrelated", srt=90, cex=0.9, adj=0.5)
text(38, 5, "No shift",     srt=90, cex=0.9, adj=0.5)

title(ylab=expression(rho == 0.9),line=7,cex.lab=1.2)
par(xpd=T)
title(main="a)",cex=1.5, adj=0)
par(xpd=F)

## ---- correlation=0.5 ----
plot(0,0,type="n",ylim=c(0.5,6.5),xlim=c(0.5,40),xlab="",ylab="",xaxt="n",yaxt="n",bty = "n")

rect(xleft=-1.5,ybottom=-0.5,xright=45,ytop=3.5,col="gray",border=F)

vioplot::vioplot(cbind(
  GLMdecorr5=decorData$RMSEglm[corr0.5_change],
  AICdecorr5=decorData$RMSEaic[corr0.5_change],
  STANdecorr5=decorData$RMSEstan[corr0.5_change],
  GLMnorm5=decorData$RMSEglm[corr0.5_nochange],
  AICnorm5=decorData$RMSEaic[corr0.5_nochange],
  STANnorm5=decorData$RMSEstan[corr0.5_nochange]
),
  las=1,
  cex.axis=0.9,
  names=c(
    "GLM,\n decorrelated",
    "Stepwise AIC,\n decorrelated",
    "RHS,\ndecorrelated",
    "GLM,\nno shift",
    "Stepwise AIC,\nno shift",
    "RHS,\nno shift"
  ),
  xlab="",
  horizontal=T,
  ylab="",
  add=T,
  col=rep(colors[4:6], 2),
  border=rep(c("darkgreen","tan4","royalblue4"), 2),
  rectCol=rep(c("darkgreen","tan4","royalblue4"), 2),
  lineCol=rep(c("darkgreen","tan4","royalblue4"), 2),
  colMed=rep(c("darkgreen","tan4","royalblue4"), 2)
)

axis(1,at=seq(0,40,by=5),labels=seq(0,40,by=5))

axis(2, at = 1:6,
     labels = rep(c("Full model", "Stepwise AIC", "RHS"), 2),
     las=2)

text(38, 2, "Decorrelated", srt=90, cex=0.9, adj=0.5)
text(38, 5, "No shift",     srt=90, cex=0.9, adj=0.5)

title(ylab=expression(rho == 0.5),line=7,cex.lab=1.2)
par(xpd=T)
title(main = "b)",cex=1.5, adj = 0)
par(xpd=F)


## ---- correlation=0.1 ----
par(mar = c(5, 9, 3, 2))
plot(0,0,type="n",ylim=c(0.5,6.5),xlim=c(0.5,40),xlab="",ylab="",xaxt="n",yaxt="n",bty = "n")

rect(xleft=-1.5,ybottom=-0.5,xright=45,ytop=3.5,col="gray",border=F)

vioplot::vioplot(cbind(
  GLMdecorr1=decorData$RMSEglm[corr0.1_change],
  AICdecorr1=decorData$RMSEaic[corr0.1_change],
  STANdecorr1=decorData$RMSEstan[corr0.1_change],
  GLMnorm1=decorData$RMSEglm[corr0.1_nochange],
  AICnorm1=decorData$RMSEaic[corr0.1_nochange],
  STANnorm1=decorData$RMSEstan[corr0.1_nochange]
),
  las=1,
  cex.axis=0.9,
  names=c(
    "Stepwise AIC,\n decorrelated",
    "RHS,\ndecorrelated",
    "Stepwise AIC,\nno shift",
    "RHS,\nno shift"
  ),
  xlab="",
  horizontal=T,
  ylab="",
  add=T,
  col=rep(colors[1:3], 2),
  border=rep(c("darkgreen","tan4","royalblue4"), 2),
  rectCol=rep(c("darkgreen","tan4","royalblue4"), 2),
  lineCol=rep(c("darkgreen","tan4","royalblue4"), 2),
  colMed=rep(c("darkgreen","tan4","royalblue4"), 2)
)

axis(1,at=seq(0,40,by=5),labels=seq(0,40,by=5))

axis(2, at = 1:6,
     labels = rep(c("Full model", "Stepwise AIC", "RHS"), 2),
     las=2)

text(38, 2, "Decorrelated", srt=90, cex=0.9, adj=0.5)
text(38, 5, "No shift",     srt=90, cex=0.9, adj=0.5)

title(ylab=expression(rho == 0.1), line=7, cex.lab=1.2)
title(xlab = "Prediction RMSE", cex.lab = 1.2)
title(main = "c)", adj=0)
par(xpd=T)
text(-8,7,"c)",cex=1.5)
par(xpd=F)


dev.off()
}
