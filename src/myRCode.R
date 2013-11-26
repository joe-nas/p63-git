 
 
collapse2gene_id <- function(txdb,which){
  object <- txdb
  gr <- biovizBase:::fetch(object, which, type = "all", 
    truncate.gaps = FALSE, truncate.fun = NULL, 
    ratio = 0.0025)

  gr.ids <- unique(levels(gr$gene_id))
  gr.reduce <- GRanges()
  mcols(gr.reduce)  <-  data.frame(gene_id=factor(),type=factor())
  for (i in gr.ids){
    gr.tmp <- gr[gr$gene_id %in% i]
    for ( j in c("cds","utr")){
      if (j %in% gr.tmp$type) {
        t <- reduce(gr.tmp[gr.tmp$type %in% j])
        t$gene_id <- i
        t$type <- j
        gr.reduce <- c(gr.reduce,t)
      }
    }
  }
  return(gr.reduce)
}

ens2name <- function(gr){
  require("biomaRt")
  ens_gene_ids <- unique(gr$gene_id)
  mart <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
  results <- getBM(attributes = c("external_gene_id","ensembl_gene_id"),
                   filters = "ensembl_gene_id", 
                   values = ens_gene_ids ,mart = mart)

  gr.ret <- GRanges()
  for (i in results$ensembl_gene_id){
    gr.tmp <- gr[gr$gene_id %in% i]
    mcols(gr.tmp)$external_gene_id <- results$external_gene_id[results$ensembl_gene_id %in% i]
    gr.ret <- c(gr.ret,gr.tmp)
  }
  return(gr.ret)
}
#"http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneNhekH3k4me1StdSig.bigWig"


getHistoneBW <- function(myConnection,myCoordinates,asRangedData=FALSE,name){
  library(rtracklayer)
  library(Rsamtools)
  asRangedData=asRangedData
  cat("downloading bigWig data from UCSC\n")
  df <- import(con=myConnection,which=myCoordinates,asRangedData=asRangedData)
  df$series <- name
  return(df)
}


get1KG <- function(snp,population,cleftGenes=NULL,txdb,titlegene,p63motifs,mystat,maf,rstart=250000,rend=250000){
  require("snpStats")
  require("biomaRt")
  require("GenomicRanges")
  require("ggplot2")
  require("AnnotationDbi")
  require("ggbio")
  
  mart <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
  snpmart <- useMart("snp",dataset="hsapiens_snp")

  makeGenomicRegion <- function(snp,snpmart){
    snpInfo <- snpInfo  <-  getBM(attributes = c('refsnp_id', 'chrom_start', 
                                                 'chr_name'),
                                  filters = 'snp_filter', values = snp, 
                                  mart = snpmart)
    genomicRegion <- GRanges(Rle(paste("chr",snpInfo$chr_name,sep="")),
                     IRanges(start=snpInfo$chrom_start-rstart,
                             end=snpInfo$chrom_start+rend),
                             name=snpInfo$refsnp_id,
                             snpLoc=snpInfo$chrom_start)
    return(genomicRegion)
  }

  cat("calculating Linkage Disequilibrium and LLR\n")
  getSnpStats <- function(ld.data,genomicRegion=NULL){
    if (is.null(genomicRegion)) {
      genomicRegion <- ld.data$genomicRegion
    }
    for (i in population){
      ld.data[[i]]$snp <- ld.data[[i]]$map[ld.data[[i]]$map$snp.names %in% snp,]
      ld.data[[i]]$rest <- ld.data[[i]]$map[!ld.data[[i]]$map$snp.names %in% snp,]
      if (nrow(ld.data[[i]]$snp) == "1") {
        cat(paste(paste("calculating Linkage Disequilibrium and LLR for target SNP in", i),".\n",sep=""))
        cat("calculating statistics\n")
        ld.data[[i]]$STAT <- ld(x=ld.data[[i]]$genotypes[,as.numeric(row.names(ld.data[[i]]$snp))],
                                y=ld.data[[i]]$genotypes[,as.numeric(row.names(ld.data[[i]]$rest))], 
                                stats=c("D.prime","LLR","R.squared"))

        cat("points to be plotted\n")
        switch(mystat,
        LD = ld.data[[i]]$MAF05 <- which(col.summary(ld.data[[i]]$genotype)$MAF[-as.numeric(row.names(ld.data[[i]]$snp))] >= maf & !is.na(ld.data[[i]]$STAT$D.prime | ld.data[[i]]$STAT$LLR)),
        RS = ld.data[[i]]$MAF05 <- which(col.summary(ld.data[[i]]$genotype)$MAF[-as.numeric(row.names(ld.data[[i]]$snp))] >= maf & !is.na(ld.data[[i]]$STAT$R.squared))
        )
        cat('FIRST DEBUG\n')
        print(ld.data[[i]]$MAF05)
        cat("points data frame\n")
        ld.data[[i]]$df.point <- data.frame(
          position=ld.data[[i]]$map$V2[ld.data[[i]]$MAF05],  
          D.prime=round(ld.data[[i]]$STAT$D.prime[ld.data[[i]]$MAF05],digits=10),
          R.squared=round(ld.data[[i]]$STAT$R.squared[ld.data[[i]]$MAF05],digits=10))
#         print(head(ld.data[[i]]$df.point))
        cat("smooth data frame\n")
        ld.data[[i]]$df.smooth <- data.frame(
          position=ld.data[[i]]$map$V2[ld.data[[i]]$MAF05],
          D.prime=c(smooth(round(ld.data[[i]]$STAT$D.prime[ld.data[[i]]$MAF05],digits=10))),
          R.squared=c(smooth(round(ld.data[[i]]$STAT$R.squared[ld.data[[i]]$MAF05],digits=10))))
#         print(head(ld.data[[i]]$df.smooth))  
        cat("p63 motifs\n")
        p63MotifHits <- findOverlaps(genomicRegion,p63motifs)
        p63motifs <- p63motifs[subjectHits(p63MotifHits)]
        
        cat("generate plot\n")
        ld.data$popOrder <- c(ld.data$popOrder,i)
        cat("pop order\n")
        switch(mystat,
          LD = ld.data$plots[[i]] <- (ggplot() + 
                                      geom_line(data=ld.data[[i]]$df.smooth,
                                                aes(x=position,y=D.prime)) + 
                                      geom_point(data=ld.data[[i]]$df.point,
                                                 aes(x=position,y=D.prime,color=D.prime),size=1.25)  + 
                                      scale_colour_gradient(low="yellow",high="red", limits=c(0, 1), breaks=c(0,0.5,1), guide = guide_colorbar(barheight= 1.8)) + 
                                      geom_rect(p63motifs,color="blue",rect.height=0.1) + 
                                      xlim(c(start(ld.data$genomicRegion),end(ld.data$genomicRegion)))),
          RS = ld.data$plots[[i]] <- (ggplot() + 
                                      geom_line(data=ld.data[[i]]$df.smooth,
                                                aes(x=position,y=R.squared)) + 
                                      geom_point(data=ld.data[[i]]$df.point,
                                                aes(x=position,y=R.squared,color=R.squared),size=1.25) + 
                                      scale_colour_gradient(low="blue", high="green",limits = c(0,1), breaks=c(0,0.5,1), guide = guide_colorbar(barheight= 1.8)) +
                                      geom_rect(p63motifs,color="blue",rect.height=0.1) + 
                                      xlim(c(start(ld.data$genomicRegion),end(ld.data$genomicRegion))) +
                                      geom_hline(yintercept=0.8, alpha=0.5))
        )
        cat("generate GRanges from SNPs for promoter overlap.\n")
        switch(mystat,
          LD = {ld.data[[i]]$h.LD <- na.omit(ld.data[[i]]$STAT$"D.prime"[,(round(ld.data[[i]]$STAT$D.prime,digits=10) >= 1 & ld.data[[i]]$STAT$LLR >= 2)])
                cat("HERE IS LD\n")
#                 print(ld.data[[i]]$h.LD)
               },
          RS = {print(ld.data[[i]]$STAT$R.squared >= 0.8)
                ld.data[[i]]$h.RS <- na.omit(ld.data[[i]]$STAT$"R.squared"[,(ld.data[[i]]$STAT$R.squared >= 0.8)])
                cat("HERE IS RS\n")
#                 print(ld.data[[i]]$h.RS)
               }
        )
        cat(paste(i,"\n"))
        
        #creates GRanges objects from h.LD data.frames

        switch(mystat,
          LD = if (!is.null(names(ld.data[[i]]$h.LD))) {
                try(ld.data[[i]]$h.LD <- GRanges(Rle(rep(seqnames(seqinfo(genomicRegion)),
                                                     length(ld.data[[i]]$h.LD))),
                                                 IRanges(
                start=ld.data[[i]]$map$V2[ld.data[[i]]$map$snp.names %in% names(ld.data[[i]]$h.LD)],
                end=ld.data[[i]]$map$V2[ld.data[[i]]$map$snp.names %in% names(ld.data[[i]]$h.LD)]),
                                                   names=rep(i,length(ld.data[[i]]$h.LD)),
                                                   snp=names(ld.data[[i]]$h.LD), 
                                                   D.prime=ld.data[[i]]$STAT$D.prime[,names(ld.data[[i]]$h.LD)]
                                           )
                )
               } else {
                cat("ld.data[[i]]$h.LD <- NULL\n")
                ld.data[[i]]$h.LD <- NULL
               },
        
        
        
          RS = if (!is.null(names(ld.data[[i]]$h.RS))) {
            cat('BEFORE SWITH WITHIN SWITCH RS \n')
            try(ld.data[[i]]$h.RS <- GRanges(
              Rle(rep(seqnames(seqinfo(genomicRegion)),length(ld.data[[i]]$h.RS))),
              IRanges(
                start=ld.data[[i]]$map$V2[ld.data[[i]]$map$snp.names %in% names(ld.data[[i]]$h.RS)],
                end=ld.data[[i]]$map$V2[ld.data[[i]]$map$snp.names %in% names(ld.data[[i]]$h.RS)]
                ),
              names=rep(i,length(ld.data[[i]]$h.RS)),
              snp=names(ld.data[[i]]$h.RS), 
              D.prime=ld.data[[i]]$STAT$R.squared[,names(ld.data[[i]]$h.RS)],
            ))
            cat("here is RS results as gr object\n")
#             print(ld.data[[i]]$h.RS)
          } else {
            cat("ld.data[[i]]$h.RS <- NULL\n")
            ld.data[[i]]$h.RS <- NULL
          }
        )

      } else {
        cat("Target SNP ins not known in ",i,", no calculation of Linkage Disequilibrium and LLR.\n")
      }
    }
    cat("return\n")
    print(summary(ld.data))
    return(ld.data)
  }

#gets 5' UTRs from biomaRt; overlaps them with h.LD SNPs; returns Arches plot
  getArches <- function(ld.data,cleftGenes,mart,genomicRegion=NULL){
    if (is.null(genomicRegion)) {
      genomicRegion <- ld.data$genomicRegion
    }
    
    cat("find h.LD SNP promoter overlap\n")
    prmtrs <- na.omit(getBM(attributes=c("chromosome_name","5_utr_start", 
          "5_utr_end","external_gene_id"), 
    filters=c("chromosome_name","start","end"), 
    values=list(chr_name=sub("chr","",seqnames(seqinfo(genomicRegion))),
        chrom_start=start(genomicRegion), 
        chrom_end=end(genomicRegion)
    ),mart=mart))
    prmtrs$Pstart <- prmtrs$"5_utr_start"-2000
    prmtrs$Pend <- prmtrs$"5_utr_end"+2000
    prmtrs <- try(GRanges(paste("chr",prmtrs$chromosome_name,sep=""),
          IRanges(prmtrs$Pstart,prmtrs$Pend),
          seqInfo=prmtrs$external_gene_id),
          silent=F)    
    if (!is.null(cleftGenes)) {
      prmtrs <- prmtrs[prmtrs$seqInfo %in% cleftGenes$V1]
    }
    ld.data$prmtrs <- prmtrs    
    tmp <- names(ld.data[lapply(ld.data,length) == 10])
    print(lapply(ld.data,length))
    print(summary(ld.data$CEU))
    print(head(ld.data$CEU$h.RS))
    print(tmp)
    tmp <- tmp[sapply(tmp,nchar) == 3]
    ld.data$promoterSNPs <- GRanges()
    cat(tmp)
    for ( i in tmp){
      cat(paste(paste("finding overlap between h.LD SNPs with 5'UTR +/- 2000bp in",i),"\n"))
      switch(mystat,
        LD = {cat("debug LD\n")
              cat(class(ld.data[[i]]$h.LD),"ld.data",class(prmtrs),"prmtrs\n")
              hits <- findOverlaps(ld.data[[i]]$h.LD,prmtrs)
              print(hits)
              if (length(hits) != 0) {
                ld.data$promoterSNPs <- append(ld.data$promoterSNPs,ld.data[[i]]$h.LD[unique(queryHits(hits))])
              }
        },
        
        RS = {cat("debug RS\n")
              cat(class(ld.data[[i]]$h.RS),"ld.data",class(prmtrs),"prmtrs\n")
              hits <- findOverlaps(ld.data[[i]]$h.RS,prmtrs)
              print(hits)
              cat('post problematic overlap\n')
#               print(d.data[[i]]$h.RS)
              if (length(hits) != 0) {
#               print(hits)
                ld.data$promoterSNPs <- append(ld.data$promoterSNPs,ld.data[[i]]$h.RS[unique(queryHits(hits))])
#                 print(ld.data$promoterSNPs)
                cat("after append\n")
              }
        }
        
      )
    }
  cat("promoter SNPs\n")

    #ld.data$prmtrs <- reduce(ld.data$prmtrs)
    snpInfo <- snpInfo  <-  getBM(attributes = c('refsnp_id', 'chrom_start', 'chr_name'),filters = 'snp_filter', values = snp, mart = snpmart)
    if (!is.null(cleftGenes)) {
      if (!is.null(ld.data$promoterSNPs) & length(ld.data$promoterSNPs) != "0") {
        arches.df <- data.frame(start=start(ld.data$promoterSNPs),
                                end=rep(rep(snpInfo$chrom_start,length(ld.data$promoterSNPs))),
                                names=ld.data$promoterSNPs$names)
        cat('DEBUG DOWNSTREAM\n')

        arches.df <- unique(t(apply(arches.df,1,sort.int,method='quick')))
        arches.df <- GRanges(
                  Rle(paste("chr",snpInfo$chr_name,sep="")),
                  IRanges(
                    start=as.numeric(arches.df[,1]),
                    end=as.numeric(arches.df[,2])
                  ),
                  names=arches.df[,3]
                )
        print(arches.df)
        cat("plotssssss\n")
        require("ggbio")
        ld.data$plots[["arches"]] <- ggplot(arches.df) + geom_arch(aes(color=names,alpha=0.25),facets=names ~ seqnames) + theme(legend.position="none")
      }
    }
    return(ld.data)
  }

  getMyreads <- function(ld.data,genomicRegion=NULL){
    cat("aggreating reads\n")
    if (is.null(genomicRegion)) {
      genomicRegion <- ld.data$genomicRegion
    }
    myreads <- list()
    h3k4me1 <- getHistoneBW(myCoordinates=genomicRegion,myConnection="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneNhekH3k4me1StdSig.bigWig",name="h3k4me1")
    #h3k4me3 <- getHistoneBW(myCoordinates=genomicRegion,myConnection="http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneNhekH3k4me3StdSig.bigWig",name="h3k4me3")
    h3k27ac <- getHistoneBW(myCoordinates=genomicRegion,myConnection="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneNhekH3k27acStdSig.bigWig",name="h3k27ac")
    #pol2 <- getHistoneBW(myCoordinates=genomicRegion,myConnection="http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneNhekPol2bStdSig.bigWig",name="pol2")
    p63reads <- getHistoneBW(myCoordinates=genomicRegion,myConnection="hg19p63.bigWig",name="p63")
    h3k4me1$score <- sqrt(h3k4me1$score)
    h3k27ac$score <- sqrt(h3k27ac$score)
    #pol2$score <- sqrt(pol2$score)
    #h3k4me3$score <- sqrt(h3k4me3$score)
    #myreads <- c(h3k4me1,h3k27ac,p63reads,pol2,h3k4me3)
    myreads <- c(h3k4me1,h3k27ac,p63reads)
    cat("start at myReads\n")
    ld.data$myReads <- myreads
    ld.data$plots$marks <- ggplot(ld.data$myReads) + stat_identity(geom='line',aes(x=start,y=score,color=series),facets=series ~ seqnames,size=0.4) + theme(legend.position="none") 
    cat("myReads are done\n")
    return(ld.data)
  }
  
  getGenes <- function(ld.data,txdb=NULL,genomicRegion=NULL){
    if (is.null(genomicRegion)) {
      genomicRegion <- ld.data$genomicRegion
    }
    if (is.null(txdb)) {
      txdb <- loadDb("TxDb.Hsapiens.BioMart.ensembl.GRCh37.p8.sqlite")
    }
    myTranscriptRanges <- GRanges(Rle(sub("chr","",seqnames(seqinfo(genomicRegion)))),IRanges(start=start(genomicRegion),end=end(genomicRegion)))
    transcripts <- collapse2gene_id(txdb,myTranscriptRanges)
    transcripts <- ens2name(transcripts)
    transcripts <- ggplot(transcripts)+geom_alignment(aes(group=external_gene_id,color=strand))
    ld.data$plots$Genes <- transcripts
    return(ld.data)
  }  

#loading or downloading and loading .ped and .info files  
  
  getFiles <- function(snp,population,genomicRegion=NULL){
    ld.data <- list()
    if (is.null(genomicRegion)) {
      genomicRegion <- makeGenomicRegion(snp,snpmart)
      ld.data$genomicRegion <- genomicRegion
    }
    ld.data$genomicRegion <- genomicRegion
    gRchrom <- sub("chr","",seqnames(seqinfo(genomicRegion)))
    gRstart <- start(genomicRegion)
    gRend <- end(genomicRegion)
    
    for (j in population){    
      tryError <- try(ld.data[[paste(j)]] <- read.pedfile(
#           file=paste("InfoPed/",gRchrom,".",gRstart,".",gRend,".",i,".ped",sep=""),
#           snps=paste("InfoPed/",gRchrom,".",gRstart,".",gRend,".",i,".info",sep=""))
          file=paste("files/",gRchrom,".",gRstart,".",gRend,".",j,".ped",sep=""),
          snps=paste("files/",gRchrom,".",gRstart,".",gRend,".",j,".info",sep=""))
          )
      cat("loading .ped and .info files\n")
      if (class(tryError) ==  "try-error") {
        cat("downloading .ped and .info files\n")      
#         system(command=paste("perl vcf_to_ped_converter.pl -vcf /scratch/jfalck/1000G/ALL.chr",gRchrom,".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz -sample_panel_file /scratch/jfalck/1000G/phase1_integrated_calls.20101123.ALL.panel -population ",i," -region ",gRchrom,":",gRstart,"-",gRend," -output_ped InfoPed/",gRchrom,".",gRstart,".",gRend,".",i,".ped -output_info InfoPed/",gRchrom,".",gRstart,".",gRend,".",i,".info",sep=""),wait=T)    
          system(command=paste("perl /scratch/jfalck/vcf_to_ped_converter.pl -vcf ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/ALL.chr",gRchrom,".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz -sample_panel_file /scratch/jfalck/1000G/phase1_integrated_calls.20101123.ALL.panel -population ",j," -region ",gRchrom,":",gRstart,"-",gRend," -output_ped /scratch/jfalck/files/",gRchrom,".",gRstart,".",gRend,".",j,".ped -output_info /scratch/jfalck/files/",gRchrom,".",gRstart,".",gRend,".",j,".info",sep=""),wait=T)
       cat("loading .ped and .info files\n")
        cat("InfoPed/",gRchrom,".",gRstart,".",gRend,".",j,".ped\n")
        ld.data[[paste(j)]] <- read.pedfile(
          file=paste("files/",gRchrom,".",gRstart,".",gRend,".",j,".ped",sep=""),
          snps=paste("files/",gRchrom,".",gRstart,".",gRend,".",j,".info",sep=""))
      }
    }
    return(ld.data)
  }
  
  
  population <- sort(population)
  ld.data <- getFiles(snp,population)
  cat("getSnpStats\n")
  ld.data <- getSnpStats(ld.data)
  cat("get genes\n")
  ld.data <- getGenes(ld.data)
  cat("getArches\n")
  ld.data <- getArches(ld.data,cleftGenes,mart)
  cat("getMyreads\n")
  ld.data <- getMyreads(ld.data)

  if (length(ld.data$plots$arches) == "0") {
    cat("if \n")
    ld.data$plots <- ld.data$plots[c("Genes","marks",ld.data$popOrder)]
  } else {
    cat("else\n")
    ld.data$plots <- ld.data$plots[c("Genes","marks","arches",ld.data$popOrder)]
  }
  cat("nach\n")
  #ld.data <- promoterOverlap(ld.data,cleftGenes,mart)
  #ld.data <- getGenes(ld.data,txdb,genomicRegion)
  #return(ld.data)
  p <- tracks(main=paste(snp,titlegene,sep=" "),ld.data$plots) + xlim(ranges(ld.data$genomicRegion)) + geom_vline(xintercept = c(ld.data$genomicRegion$snpLoc), alpha=0.5) + theme(axis.text.y=element_text(size=4.5),strip.text.y = element_text(size = 6,angle=90)) + labs(size=2)
#   
#   print(ld.data$popOrder)
#   print(summary(ld.data$plots))
#   print(IRanges(start=start(ld.data$genomicRegion),end=end(ld.data$genomicRegion)))
#   cat("ld.data$plots$YRI\n")
#   print(summary(ld.data$plots$YRI))
#   cat("ld.data$plots$marks\n")
#   print(summary(ld.data$plots$marks))
  return(p)
#   return(tracks(ld.data$plots))
  
#   print(summary(ld.data$plots))
#   print(ld.data$genomicRegion)
#   print(start(ld.data$genomicRegion))
#   p <- tracks(ld.data$plots$YRI) + xlim(IRanges(start=start(ld.data$genomicRegion),end=end(ld.data$genomicRegion)))
#   
#   return(p)
}

#get1KG("rs3809857",c("CEU","YRI","JPT","CHB"),cleftgenes,txdb,"test",p63Ranges)


# printfun <- function(){
#   require("ggplot2")
#   require("doMC")
#   require("foreach")
#   require("ggbio")
#   registerDoMC(cores=10)
#   require("AnnotationDbi")
#   for ( i in c(1:60)){
#     cat(paste(i,"\n",sep=""))
# #     cat(paste(paste(snplist$V1[i],snplist$V2[i],sep="_"),"\n",sep=""))
#     cat(paste(paste(schnips$lead_snp[i],schnips$name[i],sep="_"),"\n",sep=""))
#     #png(file=paste(paste(snplist$V1[i],snplist$V2[i],sep="_"),"png",sep="."))
# #     get1KG(schnips$lead_snp[i],c("CEU","CHB","JPT","YRI"),cleftGenes,txdb,schnips$name[i],p63Ranges)
#     ggsave(plot=get1KG(schnips$lead_snp[i],c("CEU","CHB","JPT","YRI"),cleftGenes,txdb,schnips$name[i],p63Ranges),filename=paste(paste(schnips$lead_snp[i],schnips$name[i],sep="_"),"pdf",sep="."),width=11, height=8.5)
#     cat("done\n")
#   }
# }

printfun.rs <- function(){
  require("ggplot2")
  require("doMC")
  require("foreach")
  require("ggbio")
  registerDoMC(cores=10)
  require("AnnotationDbi")
  cat("befor e for each\n")
  
  foreach ( i=c(1:nrow(unique(resAll$rsResults[,c(10,16)]))), .combine="c", .errorhandling="remove", .verbose=T) %dopar% {
#   for (i in c(1:60)){
    cat(paste(i,"\n",sep=""))
    cat(unique(resAll$rsResults[,c(10,16)])[i,1],"\n")  
    #png(file=paste(paste(snplist$V1[i],snplist$V2[i],sep="_"),"png",sep="."))
    ggsave(plot=get1KG(unique(resAll$rsResults[,c(10,16)])[i,1],c("CEU","CHB","JPT","YRI"),cleftGenes,txdb,paste(unique(resAll$rsResults[,c(10,16)])[i,2],sep=" "),p63Ranges,"RS",maf=0.00),filename=paste(unique(resAll$rsResults[,c(10,16)])[i,1],unique(resAll$rsResults[,c(10,16)])[i,2],"final.rs.maf0.00.pdf",sep="."),width=11, height=8.5)
    cat("done\n")
  }
}

printfun.ld <- function(){
  require("ggplot2")
  require("doMC")
  require("foreach")
  require("ggbio")
  registerDoMC(cores=10)
  require("AnnotationDbi")
  cat("befor e for each\n")
  
  foreach ( i=c(1:nrow(unique(resAll$ldResults[,c(10,16)]))), .combine="c", .errorhandling="remove", .verbose=T) %dopar% {
#   for (i in c(1:60)){
    cat(paste(i,"\n",sep=""))
    cat(unique(resAll$ldResults[,c(10,16)])[i,1],"\n")  
    #png(file=paste(paste(snplist$V1[i],snplist$V2[i],sep="_"),"png",sep="."))
    ggsave(plot=get1KG(unique(resAll$ldResults[,c(10,16)])[i,1],c("CEU","CHB","JPT","YRI"),cleftGenes,txdb,paste(unique(resAll$ldResults[,c(10,16)])[i,2],sep=" "),p63Ranges,"LD",maf=0.05),filename=paste(unique(resAll$ldResults[,c(10,16)])[i,1],unique(resAll$ldResults[,c(10,16)])[i,2],"final.ld.maf0.05.pdf",sep="."),width=11, height=8.5)
    cat("done\n")
  }
}






snpStufffun <- function(snps,populations){
  require("biomaRt")
  require("snpStats")
  require("doMC")
  require("GenomicRanges")
  require("foreach")

  registerDoMC(cores=10)

  snp.mart <- useMart("snp",dataset="hsapiens_snp")
  snps.gr <- makeGenomicRegions(snps=snps)
  print(snps.gr)
  downloadGenomicRegions <- function(snps.df,populations){
    rv <- 1:nrow(snps.df)
    foreach (i=rv , .combine='c' , .errorhandling='remove' ) %:% 
        foreach (j=populations, .combine='c' , .errorhandling='remove') %dopar% {

        gRchrom <- sub("chr","",snps.df$seqnames[i])
        gRstart <- snps.df$start[i]
        if (gRstart < 0) {
          gRstart <- 0
        }
        gRend <- snps.df$end[i]
        if (!file.exists(paste("/scratch/jfalck/files/",gRchrom,".",gRstart,".",gRend,".",j,".info",sep="")) | file.exists(paste("/scratch/jfalck/files/",gRchrom,".",gRstart,".",gRend,".",j,".ped"))) {
          cat("1.file does not exist\n")
#           system(command=paste("perl /scratch/jfalck/vcf_to_ped_converter.pl -vcf /scratch/jfalck/1000G/ALL.chr",gRchrom,".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz -sample_panel_file /scratch/jfalck/1000G/phase1_integrated_calls.20101123.ALL.panel -population ",j," -region ",gRchrom,":",gRstart,"-",gRend," -output_ped /scratch/jfalck/files/",gRchrom,".",gRstart,".",gRend,".",j,".ped -output_info /scratch/jfalck/files/",gRchrom,".",gRstart,".",gRend,".",j,".info",sep=""),wait=T) 
          system(command=paste("perl /scratch/jfalck/vcf_to_ped_converter.pl -vcf ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/ALL.chr",gRchrom,".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz -sample_panel_file /scratch/jfalck/1000G/phase1_integrated_calls.20101123.ALL.panel -population ",j," -region ",gRchrom,":",gRstart,"-",gRend," -output_ped /scratch/jfalck/files/",gRchrom,".",gRstart,".",gRend,".",j,".ped -output_info /scratch/jfalck/files/",gRchrom,".",gRstart,".",gRend,".",j,".info",sep=""),wait=T)
        } else if (!file.exists(paste("/scratch/jfalck/files/",gRchrom,".",gRstart,".",gRend,".",j,".info",sep="")) | file.exists(paste("/scratch/jfalck/files/",gRchrom,".",gRstart,".",gRend,".",j,".ped"))) {
          cat("2.file does not exist\n")
#           system(command=paste("perl /scratch/jfalck/vcf_to_ped_converter.pl -vcf /scratch/jfalck/1000G/ALL.chr",gRchrom,".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz -sample_panel_file /scratch/jfalck/1000G/phase1_integrated_calls.20101123.ALL.panel -population ",j," -region ",gRchrom,":",gRstart,"-",gRend," -output_ped /scratch/jfalck/files/",gRchrom,".",gRstart,".",gRend,".",j,".ped -output_info /scratch/jfalck/files/",gRchrom,".",gRstart,".",gRend,".",j,".info",sep=""),wait=T) 
          system(command=paste("perl /scratch/jfalck/vcf_to_ped_converter.pl -vcf ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/ALL.chr",gRchrom,".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz -sample_panel_file /scratch/jfalck/1000G/phase1_integrated_calls.20101123.ALL.panel -population ",j," -region ",gRchrom,":",gRstart,"-",gRend," -output_ped /scratch/jfalck/files/",gRchrom,".",gRstart,".",gRend,".",j,".ped -output_info /scratch/jfalck/files/",gRchrom,".",gRstart,".",gRend,".",j,".info",sep=""),wait=T)
        } 
        else if (!file.exists(paste("/scratch/jfalck/files/",gRchrom,".",gRstart,".",gRend,".",j,".info",sep="")) | file.exists(paste("/scratch/jfalck/files/",gRchrom,".",gRstart,".",gRend,".",j,".ped"))) {
          cat("3.file does not exist\n")
#           system(command=paste("perl /scratch/jfalck/vcf_to_ped_converter.pl -vcf /scratch/jfalck/1000G/ALL.chr",gRchrom,".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz -sample_panel_file /scratch/jfalck/1000G/phase1_integrated_calls.20101123.ALL.panel -population ",j," -region ",gRchrom,":",gRstart,"-",gRend," -output_ped /scratch/jfalck/files/",gRchrom,".",gRstart,".",gRend,".",j,".ped -output_info /scratch/jfalck/files/",gRchrom,".",gRstart,".",gRend,".",j,".info",sep=""),wait=T)
          system(command=paste("perl /scratch/jfalck/vcf_to_ped_converter.pl -vcf ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/ALL.chr",gRchrom,".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz -sample_panel_file /scratch/jfalck/1000G/phase1_integrated_calls.20101123.ALL.panel -population ",j," -region ",gRchrom,":",gRstart,"-",gRend," -output_ped /scratch/jfalck/files/",gRchrom,".",gRstart,".",gRend,".",j,".ped -output_info /scratch/jfalck/files/",gRchrom,".",gRstart,".",gRend,".",j,".info",sep=""),wait=T)
        } else if (file.exists(paste("/scratch/jfalck/files/",gRchrom,".",gRstart,".",gRend,".",j,".info",sep="")) | file.exists(paste("/scratch/jfalck/files/",gRchrom,".",gRstart,".",gRend,".",j,".ped"))) {
          cat("file exists\n")
        }
#~           system(command=paste("perl vcf_to_ped_converter.pl -vcf /media/jonas/extern/1000G/ALL.chr",gRchrom,".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz -sample_panel_file /media/jonas/extern/1000G/phase1_integrated_calls.20101123.ALL.panel -population ",j," -region ",gRchrom,":",gRstart,"-",gRend," -output_ped batch2/",gRchrom,".",gRstart,".",gRend,".",j,".ped -output_info batch2/",gRchrom,".",gRstart,".",gRend,".",j,".info",sep=""),wait=T)
#          cat("downloading .ped and .info files jfalck\n")
#        }
        return(paste(i,j))
    }
  }

  downloadGenomicRegions(snps.df=snps.gr$regions.df,populations=populations)
}

makeGenomicRegions <- function(snps){
  require("biomaRt")
  cat("genomicRegions\n")
  snpmart <- useMart("snp",dataset="hsapiens_snp")
  genomicRegions <- list()
  cat("biomaRt: snsp 2 genomicRegions\n")
  genomicRegions$snps.df <-  getBM(attributes = c('refsnp_id', 'chrom_start', 'chr_name'),filters = 'snp_filter', values = snps, mart = snpmart)
  regions.gr <- GRanges(Rle(paste("chr",genomicRegions$snps.df$chr_name,sep="")),
          IRanges(start=genomicRegions$snps.df$chrom_start-250000,
              end=genomicRegions$snps.df$chrom_start+250000),
                name=genomicRegions$snps.df$refsnp_id,snpLoc=genomicRegions$snps.df$chrom_start)
  regions.df <- as.data.frame(regions.gr)
  
  return(list("regions.gr"=regions.gr,"regions.df"=regions.df))
}


getSNPs <- function(snps){
  require("biomaRt")
  require("GenomicRanges")
  
  snpmart <- useMart("snp",dataset="hsapiens_snp")
  
  mySNPs <- getBM(attributes = c('refsnp_id', 
                'chrom_start', 
                'chr_name'),
          filters = 'snp_filter', 
          values = snps, 
          mart = snpmart)
  
  mySNPs <- GRanges(Rle(paste("chr",mySNPs$chr_name,sep="")),
          IRanges(start=mySNPs$chrom_start,
              end=mySNPs$chrom_start),
          name=mySNPs$refsnp_id)
  
  return(mySNPs)
}

getGenesInRegions <- function(gr){
  require("biomaRt")
  
  mart <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
  
  prmtrs.df <- na.omit(getBM(attributes=c("chromosome_name",
                      "5_utr_start",
                      "5_utr_end",
                      "external_gene_id",
                      "strand"), 
                filters=c("chromosome_name",
                      "start",
                      "end"), 
                values=list(chr_name=sub("chr","",seqnames(seqinfo(gr))),
                      chrom_start=start(gr), 
                      chrom_end=end(gr)),
                mart=mart))
  
  
  myExt<-function(x){
  for(i in seq(nrow(x))){
      if(x[i,]$"strand" == -1){
        x[i,]$"5_utr_end" <- x[i,]$"5_utr_end"+1000
      }else if(x[i,]$strand == 1){
        x[i,]$"5_utr_start" <- x[i,]$"5_utr_start"-1000
      }
    }
   return(x)
  }
  
  prmtrs.df<-myExt(prmtrs.df)
  prmtrs.df$"5_utr_start_ext" <- prmtrs.df$"5_utr_start"
  prmtrs.df$"5_utr_end_ext" <- prmtrs.df$"5_utr_end"
  
  prmtrs.gr <- GRanges(Rle(paste("chr",prmtrs.df$chromosome_name,sep="")),
              IRanges(start=c(prmtrs.df$"5_utr_start"),
                  end=c(prmtrs.df$"5_utr_end")),
              strand=prmtrs.df$strand,
              name=prmtrs.df$external_gene_id)
  print(prmtrs.gr)
  return(prmtrs.gr)
}

getGenesFromExternelGeneID<-function(myHgnc_symbols){
  require("biomaRt")
  
  mart <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

  geneLocations.df <- na.omit(getBM(attributes=c("chromosome_name",
                      "5_utr_start",
                      "5_utr_end",
                      "external_gene_id",
                      "strand"), 
                filters="hgnc_symbol", 
                values=myHgnc_symbols,
                mart=mart))
                
  geneLocations.df$"5_utr_start_ext" <- geneLocations.df$"5_utr_start" - 250000
  geneLocations.df$"5_utr_end_ext" <- geneLocations.df$"5_utr_end" + 250000
  
  
  myExt<-function(x){
  for(i in seq(nrow(x))){
      if(x[i,]$"strand" == -1){
        x[i,]$"5_utr_end" <- x[i,]$"5_utr_end"+2000
      }else if(x[i,]$strand == 1){
        x[i,]$"5_utr_start" <- x[i,]$"5_utr_start"-2000
      }
    }
   return(x)
  }
  
  geneLocations.df<-myExt(geneLocations.df)
  
  
  utr5Locations.gr <- GRanges(Rle(paste("chr",geneLocations.df$chromosome_name,sep="")),
              IRanges(start=geneLocations.df$"5_utr_start",
                  end=geneLocations.df$"5_utr_end"),
              strand=geneLocations.df$strand,
              name=geneLocations.df$external_gene_id)
              
  windowLocations.gr <- GRanges(Rle(paste("chr",geneLocations.df$chromosome_name,sep="")),
              IRanges(start=geneLocations.df$"5_utr_start_ext",
                  end=geneLocations.df$"5_utr_end_ext"),
              strand=geneLocations.df$strand,
              name=geneLocations.df$external_gene_id)
  
  
  return(list("utr5"=utr5Locations.gr,"windows"=windowLocations.gr))
}

loadInfoPed <- function(regions,populations){
  require("doMC")
  require("foreach")
  require("snpStats")
  registerDoMC(cores=10)

  ld.data <- list()
  for (pop in populations){
#~     cat(pop,"\n")
    ld.data[[paste(pop)]] <- foreach (i=1:nrow(regions), .combine='c', .errorhandling='remove') %dopar% {
          if (regions$start[i] < 0) {
              regions$start[i] <- 0
          }
          structure(list(read.pedfile(file=paste("/scratch/jfalck/files/",sub("chr","",regions$seqnames[i]),".",regions$start[i],".",regions$end[i],".",pop,".ped",sep=""),
                        snps=paste("/scratch/jfalck/files/",sub("chr","",regions$seqnames[i]),".",regions$start[i],".",regions$end[i],".",pop,".info",sep="")
                        )
                  )
                ,names=as.character(paste(sub("chr","",regions$seqnames[i]),".",regions$start[i],".",regions$end[i],"_",regions$name[i],sep=""))
              )
      }
    }
    return(ld.data)
}

calcLD <- function(ld.data){
  require("doMC")
  require("foreach")
  require("snpStats")
  require("GenomicRanges")
  registerDoMC(cores=14)
  for (n in names(ld.data)) {
    ld.data[[n]] <- foreach (i=names(ld.data[[n]]), .combine='c', .errorhandling='remove') %dopar% {
#       cat("for loop",n,i,"\n")
      cat('.')
      str <- unlist(strsplit(i,"_"))
      loc <- str[1]
      chr <- paste("chr",unlist(strsplit(loc,"\\."))[1],sep="")
      snp <- str[2]
      snp <- ld.data[[n]][[i]]$map[ld.data[[n]][[i]]$map$snp.names %in% snp,]
      rest <- ld.data[[n]][[i]]$map[!ld.data[[n]][[i]]$map$snp.names %in% snp,]

      if (nrow(snp) == "1") {
#         cat(paste(paste("calculating Linkage Disequilibrium and LLR for target SNP in", i),".\n",sep=""))
#         cat("calculating statistics\n")
        LD <- ld(  x=ld.data[[n]][[i]]$genotypes[,as.numeric(row.names(snp))],
              y=ld.data[[n]][[i]]$genotypes[,as.numeric(row.names(rest))], 
              stats=c("D.prime","LLR","R.squared"))
        dp.df <- data.frame(t(LD$D.prime))
        names(dp.df) <- "D.prime"
        llr.df <- data.frame(t(LD$LLR))
        names(llr.df) <- "LLR"
        rs.df <- data.frame(t(LD$R.squared))
        names(rs.df) <- "R.squared"
        snps <- data.frame("chromosome"=rep(chr,nrow(rest)),
                "start"=rest$V2,
                "end"=rest$V2,
                "name"=rest$snp.names,
                "D.prime"=dp.df,
                "LLR"=llr.df,
                "R.squared"=rs.df,
                "population"=rep(n,nrow(rest)),
                "lead_snp"=snp$snp.names,
                row.names=NULL,
                stringsAsFactors=FALSE)
        snps.high.ld <- snps[snps["D.prime"] == 1 & snps["LLR"] >= 2,]
        snps.high.rs <- snps[snps["R.squared"] >= 0.8,]
        STAT <- list("LD"=LD,"snp"=snp,"rest"=rest,"snps"=snps,"snps.high.ld"=snps.high.ld,"snps.high.rs"=snps.high.rs)
        structure(list(STAT),names=as.character(paste(i)))
      }
    }
  }
  return(ld.data)
}


getHighLD <- function(pop,myld){
  require(data.table)
  require(plyr)
  snp.list.main.ld <- list()
  snp.list.main.rs <- list()
  for ( i in pop ){
      snp.df.tmp.ld <- sapply(myld[[i]],"[","snps.high.ld")
      snp.df.tmp.rs <- sapply(myld[[i]],"[","snps.high.rs")
      snp.list.main.ld[[i]] <- ldply(snp.df.tmp.ld,data.frame)
      snp.list.main.rs[[i]] <- ldply(snp.df.tmp.rs,data.frame)
  }
#~   snp.df.main <- rbindlist(snp.list.main)
#~   snp.df.main <- do.call("rbind",snp.list.main)
  snp.df.main.ld <- ldply(snp.list.main.ld,data.frame)
  snp.df.main.rs <- ldply(snp.list.main.rs,data.frame)
  print('return(snp.df.main)')
#~   print(nrow(snp.df.main))
#~   snp.df.main <- na.omit(snp.df.main)
#~   snp.df.main <- snp.df.main
  return(list("snp.df.main.ld"=snp.df.main.ld,"snp.df.main.rs"=snp.df.main.rs))
}
 

myScript <- function(numberOfSnps,numberOfBatches,snps){
require("batch")
require("GenomicRanges")
batches <- msplit(1:numberOfSnps,numberOfBatches)
mynumbers <- c()
print(summary(batches))
allhighld <- data.frame()
allhighrs <- data.frame()
  for ( i in batches){
    cat(i,"\n")
    cat("make genomic regions\n")
    gr <- makeGenomicRegions(unique(snps))
    cat('getGeneesInRegion\n')
    myPromoters <- getGenesInRegions(gr$regions.gr)
    cat('laodInfoPed\n')
    myINFO <- loadInfoPed(gr$regions.df,c("CEU","CHB","YRI","JPT"))
    cat('calcLD\n')
    myLD <- calcLD(myINFO)
    cat('getHighLD\n')
    highsnps <- getHighLD(pop=c("CEU","CHB","YRI","JPT"),myld=myLD)
    allhighld <- rbind(allhighld,highsnps$snp.df.main.ld)
    allhighrs <- rbind(allhighrs,highsnps$snp.df.main.rs)
  }
return(list("LD"=allhighld,"RS"=allhighrs))
}


newFun<-function(genes,snps){
  require("biomaRt")
  require("GenomicRanges")
  require("rtracklayer")
  require("batch")
  
  allhighld <- data.frame()
  allhighrs <- data.frame()

  mySNPS <- getSNPs(snps)
  myCleftWindows <- getGenesFromExternelGeneID(genes)
  
  
  print(myCleftWindows)
  ov <- findOverlaps(mySNPS,myCleftWindows$windows)
  
  print(mySNPS[queryHits(ov)])
  print(myCleftWindows$windows[subjectHits(ov)])
  
  snpWindowCoordinates<-makeGenomicRegions(as.vector(mySNPS[queryHits(ov)]$name))
  myInfo <- loadInfoPed(snpWindowCoordinates$regions.df,c("CEU","CHB","YRI","JPT"))
  myLD <- calcLD(myInfo)
  
  highsnps <- getHighLD(pop=c("CEU","CHB","YRI","JPT"),myld=myLD)
  allhighld <- highsnps$snp.df.main.ld[complete.cases(highsnps$snp.df.main.ld),]
  allhighrs <- highsnps$snp.df.main.rs[complete.cases(highsnps$snp.df.main.rs),]
  
      
    cat('convert high ld and rs data.frames to GRanges objects\n')
    cat('ld.gr\n')
#     print(allhighld[3089:3150,])
    ld.gr<-GRanges(Rle(allhighld$chromosome),IRanges(start=allhighld$start,end=allhighld$end),name=allhighld$name)
    
    cat('rs.gr\n')
    rs.gr<-GRanges(Rle(allhighrs$chromosome),IRanges(start=allhighrs$start,end=allhighrs$end),name=allhighrs$name)

#     print(class(targetGenes))
    cat('finding overlapping snps between CL/P promoters and high rs/ld SNPs.\n')
    print(ov.tg.ld<-findOverlaps(myCleftWindows$utr5,ld.gr))
    print(ov.tg.rs<-findOverlaps(myCleftWindows$utr5,rs.gr))
    
    cat('cbind ld\n')
    ld.results<-cbind(as.data.frame(allhighld[subjectHits(ov.tg.ld),]),as.data.frame(myCleftWindows$utr5[queryHits(ov.tg.ld),]))
    cat('cbind rs\n')
    rs.results<-cbind(as.data.frame(allhighrs[subjectHits(ov.tg.rs),]),as.data.frame(myCleftWindows$utr5[queryHits(ov.tg.rs),]))

    return(list("ldResults"=ld.results, "rsResults"=rs.results))
  
#   return(list("TargetGenes"=myCleftWindows[subjectHits(ov)],"TargetSNPs"=mySNPS[queryHits(ov)],"overlap"=ov,"allGenes"=myCleftWindows,"allSNPs"=mySNPS))
  
}


snpsNearInterestingGenes<-function(snps,genes){
  require("biomaRt")
  require("GenomicRanges")
  require("rtracklayer")
  require("batch")
  
  allhighld <- data.frame()
  allhighrs <- data.frame()
  
#   lsnps<-1:length(snps)
  
#   nbatches<-length(snps)/2500
  batches <- msplit(snps,5)
  print(summary(batches))
  for( snps in batches ){
    snpRegions<-makeGenomicRegions(as.vector(snps))
    genesInSnpRegions<-getGenesInRegions(snpRegions$regions.gr)
    genesInSnpRegions<-list("regions.df"=as.data.frame(genesInSnpRegions),"regions.gr"=genesInSnpRegions)
    targetGenes<-genesInSnpRegions$regions.df[genesInSnpRegions$regions.df$name %in% genes,]
    print(head(targetGenes,n=50))
    targetGenes<-list("regions.df"=targetGenes,
                      "regions.gr"=GRanges(Rle(targetGenes$seqnames),
                                           IRanges(start=targetGenes$start,
                                                  end=targetGenes$end),
                                           name=targetGenes$name))
    ov<-findOverlaps(snpRegions$regions.gr,targetGenes$regions.gr)
    print(ov)
    myInfo <- loadInfoPed(snpRegions$regions.df[queryHits(ov),],c("CEU","CHB","YRI","JPT"))
    cat("myInfo done\n")
    myLD <- calcLD(myInfo)
    
    highsnps <- getHighLD(pop=c("CEU","CHB","YRI","JPT"),myld=myLD)
    allhighld <- rbind(allhighld,highsnps$snp.df.main.ld[complete.cases(highsnps$snp.df.main.ld),])
    allhighrs <- rbind(allhighrs,highsnps$snp.df.main.rs[complete.cases(highsnps$snp.df.main.rs),])
  }
    cat('convert high ld and rs data.frames to GRanges objects\n')
    cat('ld.gr\n')
#     print(allhighld[3089:3150,])
    ld.gr<-GRanges(Rle(allhighld$chromosome),IRanges(start=allhighld$start,end=allhighld$end),name=allhighld$name)
    
    cat('rs.gr\n')
    rs.gr<-GRanges(Rle(allhighrs$chromosome),IRanges(start=allhighrs$start,end=allhighrs$end),name=allhighrs$name)

    print(class(targetGenes))
    cat('finding overlapping snps between CL/P promoters and high rs/ld SNPs.\n')
    print(ov.tg.ld<-findOverlaps(targetGenes$regions.gr,ld.gr))
    print(ov.tg.rs<-findOverlaps(targetGenes$regions.gr,rs.gr))
    
    
    ld.results<-cbind(allhighld[subjectHits(ov.tg.ld),],targetGenes$regions.df[queryHits(ov.tg.ld),])
    rs.results<-cbind(allhighrs[subjectHits(ov.tg.rs),],targetGenes$regions.df[queryHits(ov.tg.rs),])

    return(list("ldResults"=ld.results, "rsResults"=rs.results))
}
                      

#snpsNearInterestingGenes(snps=myOverlap$V7[1:100], genes=cleftGenes$V1)
