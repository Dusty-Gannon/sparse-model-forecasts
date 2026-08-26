################################################
###### Figures for AIC vs STAN comparison ######
################################################

library(here)

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
corData2=read.csv(here("Simulations/AICvsSTANcorrResults.csv"))

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

pdf(file=here("Figures/correlated_variable_effects_combined.pdf"), width = 5, height = 9)

par(mar=c(5,9,2,2))
par(mfrow=c(3,1))
vioplot::vioplot(cbind(AIC1=corData6$TPRaic[which(corData6$corrLevel==0.1)],STAN1=corData6$TPRstan[which(corData6$corrLevel==0.1)],AIC5=corData6$TPRaic[which(corData6$corrLevel==0.5)],STAN5=corData6$TPRstan[which(corData6$corrLevel==0.5)],AIC9=corData6$TPRaic[which(corData6$corrLevel==0.9)],STAN9=corData6$TPRstan[which(corData6$corrLevel==0.9)]), xlab="",
                 horizontal=T,las=1,names=c("Stepwise AIC\n correlation 0.1","RHS\n correlation 0.1","Stepwise AIC\n correlation 0.5","RHS\n correlation 0.5","Stepwise AIC\n correlation 0.9","RHS\n correlation 0.9"),col=c("wheat","skyblue","tan","steelblue1","tan3","steelblue4"),
                 pchMed=20,border=c("tan4","royalblue4","tan4","royalblue4","tan4","royalblue4"),rectCol=c("tan4","royalblue4","tan4","royalblue4","tan4","royalblue4"), lineCol=c("tan4","royalblue4","tan4","royalblue4","tan4","royalblue4"),colMed=c("tan4","royalblue4","tan4","royalblue4","tan4","royalblue4"),ylab="")
title(ylab="Density of TPR",line=7,cex.lab=1.2)
title(xlab="True Positive Rate (TPR)",line=3,cex.lab=1.2)
par(xpd=T)
text(-0.16,7,"a)",cex=1.5)
par(xpd=F)
vioplot::vioplot(cbind(AIC1=corData6$TNRaic[which(corData6$corrLevel==0.1)],STAN1=corData6$TNRstan[which(corData6$corrLevel==0.1)],AIC5=corData6$TNRaic[which(corData6$corrLevel==0.5)],STAN5=corData6$TNRstan[which(corData6$corrLevel==0.5)],AIC9=corData6$TNRaic[which(corData6$corrLevel==0.9)],STAN9=corData6$TNRstan[which(corData6$corrLevel==0.9)]), xlab="",
                 horizontal=T,las=1,names=c("Stepwise AIC\n correlation 0.1","RHS\n correlation 0.1","Stepwise AIC\n correlation 0.5","RHS\n correlation 0.5","Stepwise AIC\n correlation 0.9","RHS\n correlation 0.9"),col=c("wheat","skyblue","tan","steelblue1","tan3","steelblue4"),
                 pchMed=20,border=c("tan4","royalblue4","tan4","royalblue4","tan4","royalblue4"),rectCol=c("tan4","royalblue4","tan4","royalblue4","tan4","royalblue4"), lineCol=c("tan4","royalblue4","tan4","royalblue4","tan4","royalblue4"),colMed=c("tan4","royalblue4","tan4","royalblue4","tan4","royalblue4"),ylab="")
title(ylab="Density of TNR",line=7,cex.lab=1.2)
title(xlab="True Negative Rate (TNR)",line=3,cex.lab=1.2)
par(xpd=T)
text(-0.275,7,"b)",cex=1.5)
par(xpd=F)
vioplot::vioplot(cbind(GLM1=corData6$RMSEglm[which(corData6$corrLevel==0.1)],AIC1=corData6$RMSEaic[which(corData6$corrLevel==0.1)],STAN1=corData6$RMSEstan[which(corData6$corrLevel==0.1)],
                       GLM5=corData6$RMSEglm[which(corData6$corrLevel==0.5)],AIC5=corData6$RMSEaic[which(corData6$corrLevel==0.5)],STAN5=corData6$RMSEstan[which(corData6$corrLevel==0.5)],
                       GLM9=corData6$RMSEglm[which(corData6$corrLevel==0.9)],AIC9=corData6$RMSEaic[which(corData6$corrLevel==0.9)],STAN9=corData6$RMSEstan[which(corData6$corrLevel==0.9)]), xlab="",
                 horizontal=T,las=1,names=c("GLM\n correlation 0.1","Stepwise AIC\n correlation 0.1","RHS\n correlation 0.1",
                                            "GLM\n correlation 0.5","Stepwise AIC\n correlation 0.5","RHS\n correlation 0.5",
                                            "GLM\n correlation 0.9","Stepwise AIC\n correlation 0.9","RHS\n correlation 0.9"),
                 col=c("lightgreen","wheat","skyblue","green3","tan","steelblue1","forestgreen","tan3","steelblue4"),
                 pchMed=20,border=c("darkgreen","tan4","royalblue4","darkgreen","tan4","royalblue4","darkgreen","tan4","royalblue4"),
                 rectCol=c("darkgreen","tan4","royalblue4","darkgreen","tan4","royalblue4","darkgreen","tan4","royalblue4"),
                 lineCol=c("darkgreen","tan4","royalblue4","darkgreen","tan4","royalblue4","darkgreen","tan4","royalblue4"),
                 colMed=c("darkgreen","tan4","royalblue4","darkgreen","tan4","royalblue4","darkgreen","tan4","royalblue4"),ylab="")
title(ylab="Density of Prediction RMSE",line=7,cex.lab=1.2)
title(xlab="Prediction Root Mean Square Error (RMSE)",line=3,cex.lab=1.2)
par(xpd=T)
text(-0.16,10.5,"c)",cex=1.5)
par(xpd=F)

dev.off()


################################################
# Figure 4
################################################
decorData=read.csv(file="Simulations/AICvsSTANdecorrelation.csv")


pdf(file=here("Figures/decorrelation_comparison.pdf"), width = 7, height = 8)


par(mar=c(5,9,2,2))
par(mfrow=c(3,1))

corr0.9_change=which(decorData$corrLevel==0.9&decorData$corrChange)
corr0.9_nochange=which(decorData$corrLevel==0.9&!decorData$corrChange)
corr0.5_change=which(decorData$corrLevel==0.5&decorData$corrChange)
corr0.5_nochange=which(decorData$corrLevel==0.5&!decorData$corrChange)
corr0.1_change=which(decorData$corrLevel==0.1&decorData$corrChange)
corr0.1_nochange=which(decorData$corrLevel==0.1&!decorData$corrChange)

# plot(1, type = "n", xlim = c(0, 12), ylim = c(0.5, 12.5),
#      xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")

# correlation=0.9
plot(0,0,type="n",ylim=c(0.5,6.5),xlim=c(0.5,40),xlab="",ylab="",xaxt="n",yaxt="n",bty = "n")

rect(xleft=-1.5,ybottom=-0.5,xright=45,ytop=3.5,col="gray",border=F)

vioplot::vioplot(cbind(GLMdecorr9=decorData$RMSEglm[corr0.9_change],AICdecorr9=decorData$RMSEaic[corr0.9_change],STANdecorr9=decorData$RMSEstan[corr0.9_change],
                       GLMnorm9=decorData$RMSEglm[corr0.9_nochange],AICnorm9=decorData$RMSEaic[corr0.9_nochange],STANnorm9=decorData$RMSEstan[corr0.9_nochange]),las=1,cex.axis=0.9,
                 names=c("GLM,\n decorrelated","Stepwise AIC,\n decorrelated","RHS,\ndecorrelated","GLM,\nno shift","Stepwise AIC,\nno shift","RHS,\nno shift"),
                 xlab="Root Mean Square Error (RMSE)",horizontal=T,ylab="",add=T,
                 col=c("forestgreen","tan3","steelblue4","forestgreen","tan3","steelblue4"),
                 border=c("darkgreen","tan4","royalblue4","darkgreen","tan4","royalblue4"),
                 rectCol=c("darkgreen","tan4","royalblue4","darkgreen","tan4","royalblue4"),
                 lineCol=c("darkgreen","tan4","royalblue4","darkgreen","tan4","royalblue4"),
                 colMed=c("darkgreen","tan4","royalblue4","darkgreen","tan4","royalblue4"))

axis(1,at=seq(0,40,by=5),labels=seq(0,40,by=5))

axis(2, at = 1:6, labels = c("GLM,\n decorrelated","Stepwise AIC,\n decorrelated","RHS,\ndecorrelated","GLM,\nno shift","Stepwise AIC,\nno shift","RHS,\nno shift"),las=2)

title(ylab="Correlation 0.9",line=7,cex.lab=1.2)
par(xpd=T)
text(-8,7,"a)",cex=1.5)
par(xpd=F)
# correlation=0.5
plot(0,0,type="n",ylim=c(0.5,6.5),xlim=c(0.5,40),xlab="",ylab="",xaxt="n",yaxt="n",bty = "n")

rect(xleft=-1.5,ybottom=-0.5,xright=45,ytop=3.5,col="gray",border=F)

vioplot::vioplot(cbind(GLMdecorr5=decorData$RMSEglm[corr0.5_change],AICdecorr5=decorData$RMSEaic[corr0.5_change],STANdecorr5=decorData$RMSEstan[corr0.5_change],
                       GLMnorm5=decorData$RMSEglm[corr0.5_nochange],AICnorm5=decorData$RMSEaic[corr0.5_nochange],STANnorm5=decorData$RMSEstan[corr0.5_nochange]),las=1,cex.axis=0.9,
                 names=c("GLM,\n decorrelated","Stepwise AIC,\n decorrelated","RHS,\ndecorrelated","GLM,\nno shift","Stepwise AIC,\nno shift","RHS,\nno shift"),
                 xlab="Root Mean Square Error (RMSE)",horizontal=T,ylab="",add=T,
                 col=c("green3","tan","steelblue1","green3","tan","steelblue1"),
                 border=c("darkgreen","tan4","royalblue4","darkgreen","tan4","royalblue4"),
                 rectCol=c("darkgreen","tan4","royalblue4","darkgreen","tan4","royalblue4"),
                 lineCol=c("darkgreen","tan4","royalblue4","darkgreen","tan4","royalblue4"),
                 colMed=c("darkgreen","tan4","royalblue4","darkgreen","tan4","royalblue4"))

axis(1,at=seq(0,40,by=5),labels=seq(0,40,by=5))

axis(2, at = 1:6, labels = c("GLM,\n decorrelated","Stepwise AIC,\n decorrelated","RHS,\ndecorrelated","GLM,\nno shift","Stepwise AIC,\nno shift","RHS,\nno shift"),las=2)

title(ylab="Correlation 0.5",line=7,cex.lab=1.2)
par(xpd=T)
text(-8,7,"b)",cex=1.5)
par(xpd=F)
# correlation=0.1
plot(0,0,type="n",ylim=c(0.5,6.5),xlim=c(0.5,40),xlab="",ylab="",xaxt="n",yaxt="n",bty = "n")

rect(xleft=-1.5,ybottom=-0.5,xright=45,ytop=3.5,col="gray",border=F)

vioplot::vioplot(cbind(GLMdecorr1=decorData$RMSEglm[corr0.1_change],AICdecorr1=decorData$RMSEaic[corr0.1_change],STANdecorr1=decorData$RMSEstan[corr0.1_change],GLMnorm1=decorData$RMSEglm[corr0.1_nochange],AICnorm1=decorData$RMSEaic[corr0.1_nochange],STANnorm1=decorData$RMSEstan[corr0.1_nochange]),las=1,cex.axis=0.9,
                 names=c("Stepwise AIC,\n decorrelated","RHS,\ndecorrelated","Stepwise AIC,\nno shift","RHS,\nno shift"),
                 xlab="Root Mean Square Error (RMSE)",horizontal=T,ylab="",add=T,
                 col=c("lightgreen","wheat","skyblue","lightgreen","wheat","skyblue"),
                 border=c("darkgreen","tan4","royalblue4","darkgreen","tan4","royalblue4"),
                 rectCol=c("darkgreen","tan4","royalblue4","darkgreen","tan4","royalblue4"),
                 lineCol=c("darkgreen","tan4","royalblue4","darkgreen","tan4","royalblue4"),
                 colMed=c("darkgreen","tan4","royalblue4","darkgreen","tan4","royalblue4"))

axis(1,at=seq(0,40,by=5),labels=seq(0,40,by=5))

axis(2, at = 1:6, labels = c("GLM,\n decorrelated","Stepwise AIC,\n decorrelated","RHS,\ndecorrelated","GLM,\nno shift","Stepwise AIC,\nno shift","RHS,\nno shift"),las=2)

title(ylab="Correlation 0.1",line=7,cex.lab=1.2)
par(xpd=T)
text(-8,7,"c)",cex=1.5)
par(xpd=F)


dev.off()
