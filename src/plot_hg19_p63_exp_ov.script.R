require('ggplot2')
require('reshape')
resultsFolder<-"results/tables/p63_exp_vs_ovs_overlap"
resultsFileNames<-list.files(resultsFolder,pattern="Encode")
resultsFiles<-paste(resultsFolder,resultsFileNames,sep="/")
resultsNames<-print(gsub(".-Y.bed","",gsub("hg19_p63_obs_vs_exp_wgEncodeBroadChipSeqPeaks","",resultsFileNames)))
myDataList<-lapply(resultsFiles,function(x) read.table(x,header=F,sep=" ",
                                                   col.names=c("Group","Iteration","Overlap")))
names(myDataList)<-resultsNames
myDataFrame<-melt.list(myDataList,id=1:3)
names(myDataFrame)<-c(names(myDataFrame)[-length(names(myDataFrame))],"Histone Modification")
myDataFrame[,"Histone Modification"]<-as.factor(myDataFrame[,"Histone Modification"])
str(myDataFrame)
p <- ggplot(myDataFrame,aes(get('Histone Modification'),log2(Overlap)))
p <- p + geom_boxplot() + geom_point(aes(color = Group)) 
p <- p + labs(list(title = "Observed vs. Expercted p63 BS - Histone BS overlap",
                   x = "Histone Modifications", 
                   y = "log2(Overlap)", colour = ""))
svg("results/plots/hg19_p63_hist_obs_vs_exp_all.svg",width = 14 ,height = 7)
p
dev.off()
png("results/plots/hg19_p63_hist_obs_vs_exp_all.png",width = 800, height = 400)
p
dev.off()