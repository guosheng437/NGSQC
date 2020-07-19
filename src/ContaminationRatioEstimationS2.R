######################################################
# Sheng Guo, PhD (guosheng@crownbio.com)
# July 2020
# Copyright (C) 2020 Sheng Guo, Crown Bioscience Inc.
######################################################
rm(list=ls())
library(MASS)
library(mclust)
library(EnvStats)

#--change directory to .\Program\testedcelllines_Table2 or .\Program\testedcelllines_Table3
file_list <- list.files(pattern="*.SNPratio")

sample2mc = read.table("Sample_MajorComponent.txt", head=TRUE,colClasses = "character")

results<-c('sample','ContaminationRatio(%)')
pdf(file='SNPratio.pdf')
for (i in 1:length(file_list)){
	maintitle = file_list[i]
	maintitle = gsub('\\.SNPratio', '', maintitle)
	maintitle = paste(maintitle, "(", sample2mc[i,2], ")", sep="") 
	
	X <- read.table(file_list[i])
	X <- X[,1]
	
	#--outlier detection--#
	rt<-rosnerTest(X, k = 10, warn = F)
	print(rt$n.outliers)
	if(rt$n.outliers>0){
		X<-X[1:(length(X)-rt$n.outliers)]
	}
		
	fit = Mclust(X,G=c(2,3),model="V", verbose=FALSE)
	if(is.null(fit)){
		fit = Mclust(X,G=1,model="V", verbose=FALSE)
	}
	
	if(is.null(fit)){
		ratio = median(X)
	}else{
		
		summary(fit)
	
		plot(fit, what="density", main=maintitle, xlab="Contamination ratio", xlim=c(0,0.2), cex.lab=1.667, cex.axis=1.667, cex.main=1.667, cex.sub=1.667)
		legend('topright', legend=maintitle,cex=1.667)
		rug(X)

		median(X[fit$classification!=1])
		median(X[fit$classification==1])		
		
		if(length(X[fit$classification==1])> 0.6*length(X)){
			ratio = median(X)
		}else{
			ratio = median(X[fit$classification!=1])
		}
	}
	
	ratio = ratio*100
	results <- rbind(results, c(file_list[i], ratio))
}
dev.off()
results<-as.data.frame(results[-1,])
names(results)<-c('sample','ContaminationRatio(%)')
write.csv(results, file="ContaminationRatio.csv")