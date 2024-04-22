################################################
###### Figures for AIC vs STAN comparison ######
################################################

################################################
# Figure 1:
################################################




################################################
# Figure 2:
################################################
corData2=read.csv(file="Simulations/AICvsSTANcorrResults.csv")

# Remove unwanted variables
corData3=corData2[-which(corData2$strongSelf),]
# fix numStrongCorr=3
corData4=corData3[-which(corData3$numStrongCorr==1),]
# fix percent of weak covariates to 0.5
corData5=corData4[which(corData4$probWeakCorr==0.5),]
# fix correlations to 0.5 and 0.9
corData6=corData5[-which(corData5$corrLevel==0.7),]

par(mar=c(5,9,2,2))
par(mfrow=c(3,1))
vioplot::vioplot(cbind(AIC5=corData6$TPRaic[which(corData6$corrLevel==0.5)],AIC9=corData6$TPRaic[which(corData6$corrLevel==0.9)],STAN5=corData6$TPRstan[which(corData6$corrLevel==0.5)],STAN9=corData6$TPRstan[which(corData6$corrLevel==0.9)]), xlab="",
                 horizontal=T,las=1,names=c("Step-wise AIC\n correlation 0.5","Step-wise AIC\n correlation 0.9","Horseshoe\n correlation 0.5","Horseshoe\n correlation 0.9"),col=c("steelblue1","steelblue4","tan","tan3"),
                 pchMed=20,border=c("royalblue4","royalblue4","tan4","tan4"),rectCol=c("royalblue4","royalblue4","tan4","tan4"), lineCol=c("royalblue4","royalblue4","tan4","tan4"),colMed=c("royalblue4","royalblue4","tan4","tan4"),ylab="")
title(ylab="Density of TPR",line=7,cex.lab=1.2)
title(xlab="TPR",line=3,cex.lab=1.2)
vioplot::vioplot(cbind(AIC5=corData6$TNRaic[which(corData6$corrLevel==0.5)],AIC9=corData6$TNRaic[which(corData6$corrLevel==0.9)],STAN5=corData6$TNRstan[which(corData6$corrLevel==0.5)],STAN9=corData6$TNRstan[which(corData6$corrLevel==0.9)]), xlab="",
                 horizontal=T,las=1,names=c("Step-wise AIC\n correlation 0.5","Step-wise AIC\n correlation 0.9","Horseshoe\n correlation 0.5","Horseshoe\n correlation 0.9"),col=c("steelblue1","steelblue4","tan","tan3"),
                 pchMed=20,border=c("royalblue4","royalblue4","tan4","tan4"),rectCol=c("royalblue4","royalblue4","tan4","tan4"), lineCol=c("royalblue4","royalblue4","tan4","tan4"),colMed=c("royalblue4","royalblue4","tan4","tan4"),ylab="")
title(ylab="Density of TNR",line=7,cex.lab=1.2)
title(xlab="TNR",line=3,cex.lab=1.2)
vioplot::vioplot(cbind(AIC5=corData6$RMSEaic[which(corData6$corrLevel==0.5)],AIC9=corData6$RMSEaic[which(corData6$corrLevel==0.9)],STAN5=corData6$RMSEstan[which(corData6$corrLevel==0.5)],STAN9=corData6$RMSEstan[which(corData6$corrLevel==0.9)]), xlab="",
                 horizontal=T,las=1,names=c("Step-wise AIC\n correlation 0.5","Step-wise AIC\n correlation 0.9","Horseshoe\n correlation 0.5","Horseshoe\n correlation 0.9"),col=c("steelblue1","steelblue4","tan","tan3"),
                 pchMed=20,border=c("royalblue4","royalblue4","tan4","tan4"),rectCol=c("royalblue4","royalblue4","tan4","tan4"), lineCol=c("royalblue4","royalblue4","tan4","tan4"),colMed=c("royalblue4","royalblue4","tan4","tan4"),ylab="")
title(ylab="Density of Prediction RMSE",line=7,cex.lab=1.2)
title(xlab="Prediction RMSE",line=3,cex.lab=1.2)

################################################
# Figure 3
################################################
corData2=read.csv(file="Simulations/AICvsSTANdecorrelation.csv")


par(mar=c(5,9,2,2))
par(mfrow=c(1,1))
vioplot::vioplot(cbind(AICdecorr=corData2$RMSEaic[1:1800],AICnorm=corData2$RMSEaic[1801:3600],STANdecorr=corData2$RMSEstan[1:1800],STANnorm=corData2$RMSEstan[1801:3600]),las=1,cex.axis=0.9,names=c("Step-wise AIC\n decorrelated","Step-wise AIC\n normal","Horseshoe\n decorrelated","Horseshoe\n normal"),
                 xlab="",horizontal=T,col=c("steelblue1","steelblue4","tan","tan3"),pchMed=20,border=c("royalblue4","royalblue4","tan4","tan4"),rectCol=c("royalblue4","royalblue4","tan4","tan4"),
                 lineCol=c("royalblue4","royalblue4","tan4","tan4"),colMed=c("royalblue4","royalblue4","tan4","tan4"),ylab="")
