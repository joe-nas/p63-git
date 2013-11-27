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
ggplot(myDataFrame,aes(get('Histone Modification'),Overlap)) + geom_violin()