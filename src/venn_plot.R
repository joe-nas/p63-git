require("Vennerable")
setwd(".")

# shortcut to folder and filenames
histoneFolder<- "../results/tables/venn_diagrams/h3k4me1_h3k4me3_h3k27ac/"
histoneFiles<-list.files("../results/tables/venn_diagrams/h3k4me1_h3k4me3_h3k27ac/",pattern="intersect")
histoneFiles<- paste(histoneFolder,histoneFiles,sep="")
# loading data
myData<-lapply(histoneFiles,function(x) read.table(x,header=F))
# changing names
m<-regexpr("*NhekH3k{1,2}[0-9]{1,2}[a-z]{1,2}[a-z]{1,2}[0-9]*",histoneFiles)
histoneNames<-regmatches(histoneFiles,m)
names(myData)<-histoneNames
# compute and plot venn diagram
myVenn<-Venn(list(myData$NhekH3k27ac[,1],
                  myData$NhekH3k4me1[,1],
                  myData$NhekH3k4me3[,1]),
             SetNames=names(myData))
# save plot as svg then as png
svg("../results/plots/venn_h3k4me1_h3k4me3_h3k27ac.svg")
plot(myVenn)
dev.off()
png("../results/plots/venn_h3k4me1_h3k4me3_h3k27ac.png")
plot(myVenn)
dev.off()