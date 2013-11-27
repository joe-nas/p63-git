require("Vennerable")
histoneFolder<- "../results/tables/venn_diagrams/h3k4me1_h3k4me3_h3k27ac/"
histoneFiles<-list.files("../results/tables/venn_diagrams/h3k4me1_h3k4me3_h3k27ac/",pattern="intersect")
histoneNames<-sapply(X=histoneFiles,FUN=function(x){gsub(x,pattern="intersect_union_wgEncodeBroadChipSeqPeaks",replacement="")})
print(histoneNames)[2,]
histoneFiles<- paste(histoneFolder,histoneFiles,sep="")
myData<-lapply(histoneFiles,function(x) read.table(x,header=F))