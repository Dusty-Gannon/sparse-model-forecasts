##################################################
## sequoia_climate_pca_plot.R
##################################################

## JPJ 27 i 25
## PURPOSE: to make a PCA showing yearly variation in climate data in the sequoia data
## USAGE: Rscript sequoia_climate_pca_plot.R


## load data and merge climate variables
ppt <- read.csv("BMG_ppt.csv", header=TRUE)
tdmean <- read.csv("BMG_td_mean.csv", header=TRUE)
tmax <- read.csv("BMG_tmax.csv", header=TRUE)
tmean <- read.csv("BMG_tmean.csv", header=TRUE)
tmin <- read.csv("BMG_tmin.csv", header=TRUE)
vpdmax <- read.csv("BMG_vpdmax.csv", header=TRUE)
vpdmin <- read.csv("BMG_vpdmin.csv", header=TRUE)
merged_data <- cbind(ppt[,6], tdmean[,6], tmax[,6], tmin[,6], vpdmax[,6], vpdmin[,6])


## Original try was based on all data, but patterns ruled by differences across months.
## Retry with a separate PCA for each month (LOOPS!!!)


colramp <- colorRampPalette(c("purple4", "white"))
colramp_118 <- colramp(118)

month_list <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")

pdf("sequoia_climate_pca_plot.pdf", height=9, width=12)
#quartz(height=9, width=12)
par(mar=c(5,5,2,2), mfrow=c(3,4))
for (i in 1:12){
	month_sub <- subset(merged_data, ppt[,5]==i)
	month_pc <- prcomp(month_sub, center=TRUE, scale.=TRUE)
	plot(month_pc$x[,1], month_pc$x[,2], type="n", main=month_list[i], cex.main=1.5, cex.lab=1.5, cex.axis=1.25, las=1,
		xlab=paste0("PC1 (", round(summary(month_pc)$importance[2,1]*100,1), "%)"),
		ylab=paste0("PC2 (", round(summary(month_pc)$importance[2,2]*100,1), "%)"))
	for (j in 1:dim(month_sub)[1]) {
		points(month_pc$x[j,1], month_pc$x[j,2], pch=21, bg=colramp_118[j], cex=2)
	}
}
dev.off()





