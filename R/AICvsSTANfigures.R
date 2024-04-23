################################################
###### Figures for AIC vs STAN comparison ######
################################################

library(here)

################################################
# Figure 1:
################################################




################################################
# Figure 2:
################################################
corData2=read.csv(here("Simulations/AICvsSTANcorrResults.csv"))

# Remove unwanted variables
corData3=corData2[-which(corData2$strongSelf),]
# fix numStrongCorr=3
corData4=corData3[-which(corData3$numStrongCorr==1),]
# fix percent of weak covariates to 0.5
corData5=corData4[which(corData4$probWeakCorr==0.5),]
# fix correlations to 0.5 and 0.9
corData6=corData5[-which(corData5$corrLevel==0.7),]


pdf(file=here("Figures/correlated_variable_effects_combined.pdf"), width = 5, height = 7.5)

par(mar=c(5,9,2,2))
par(mfrow=c(3,1))
vioplot::vioplot(cbind(AIC5=corData6$TPRaic[which(corData6$corrLevel==0.5)],AIC9=corData6$TPRaic[which(corData6$corrLevel==0.9)],STAN5=corData6$TPRstan[which(corData6$corrLevel==0.5)],STAN9=corData6$TPRstan[which(corData6$corrLevel==0.9)]), xlab="",
                 horizontal=T,las=1,names=c("Stepwise AIC\n correlation 0.5","Stepwise AIC\n correlation 0.9","Horseshoe\n correlation 0.5","Horseshoe\n correlation 0.9"),col=c("steelblue1","steelblue4","tan","tan3"),
                 pchMed=20,border=c("royalblue4","royalblue4","tan4","tan4"),rectCol=c("royalblue4","royalblue4","tan4","tan4"), lineCol=c("royalblue4","royalblue4","tan4","tan4"),colMed=c("royalblue4","royalblue4","tan4","tan4"),ylab="")
title(ylab="Density of TPR",line=7,cex.lab=1.2)
title(xlab="TPR",line=3,cex.lab=1.2)
vioplot::vioplot(cbind(AIC5=corData6$TNRaic[which(corData6$corrLevel==0.5)],AIC9=corData6$TNRaic[which(corData6$corrLevel==0.9)],STAN5=corData6$TNRstan[which(corData6$corrLevel==0.5)],STAN9=corData6$TNRstan[which(corData6$corrLevel==0.9)]), xlab="",
                 horizontal=T,las=1,names=c("Stepwise AIC\n correlation 0.5","Stepwise AIC\n correlation 0.9","Horseshoe\n correlation 0.5","Horseshoe\n correlation 0.9"),col=c("steelblue1","steelblue4","tan","tan3"),
                 pchMed=20,border=c("royalblue4","royalblue4","tan4","tan4"),rectCol=c("royalblue4","royalblue4","tan4","tan4"), lineCol=c("royalblue4","royalblue4","tan4","tan4"),colMed=c("royalblue4","royalblue4","tan4","tan4"),ylab="")
title(ylab="Density of TNR",line=7,cex.lab=1.2)
title(xlab="TNR",line=3,cex.lab=1.2)
vioplot::vioplot(cbind(AIC5=corData6$RMSEaic[which(corData6$corrLevel==0.5)],AIC9=corData6$RMSEaic[which(corData6$corrLevel==0.9)],STAN5=corData6$RMSEstan[which(corData6$corrLevel==0.5)],STAN9=corData6$RMSEstan[which(corData6$corrLevel==0.9)]), xlab="",
                 horizontal=T,las=1,names=c("Stepwise AIC\n correlation 0.5","Stepwise AIC\n correlation 0.9","Horseshoe\n correlation 0.5","Horseshoe\n correlation 0.9"),col=c("steelblue1","steelblue4","tan","tan3"),
                 pchMed=20,border=c("royalblue4","royalblue4","tan4","tan4"),rectCol=c("royalblue4","royalblue4","tan4","tan4"), lineCol=c("royalblue4","royalblue4","tan4","tan4"),colMed=c("royalblue4","royalblue4","tan4","tan4"),ylab="")
title(ylab="Density of Prediction RMSE",line=7,cex.lab=1.2)
title(xlab="Prediction RMSE",line=3,cex.lab=1.2)

dev.off()


################################################
# Figure 3
################################################
decorData=read.csv(file="Simulations/AICvsSTANdecorrelation.csv")


pdf(file=here("Figures/decorrelation_comparison.pdf"), width = 7, height = 4)


par(mar=c(5,9,2,2))
par(mfrow=c(1,1))
vioplot::vioplot(cbind(AICdecorr=decorData$RMSEaic[1:1800],STANdecorr=decorData$RMSEstan[1:1800],AICnorm=decorData$RMSEaic[1801:3600],STANnorm=decorData$RMSEstan[1801:3600]),las=1,cex.axis=0.9,names=c("Stepwise AIC\n decorrelated","Horseshoe\n decorrelated","Stepwise AIC\n normal","Horseshoe\n normal"),
                 xlab="",horizontal=T,col=c("steelblue1","tan","steelblue4","tan3"),pchMed=20,border=c("royalblue4","tan4","royalblue4","tan4"),rectCol=c("royalblue4","tan4","royalblue4","tan4"),
                 lineCol=c("royalblue4","tan4","royalblue4","tan4"),colMed=c("royalblue4","tan4","royalblue4","tan4"),ylab="")

dev.off()
