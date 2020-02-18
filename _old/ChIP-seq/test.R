library(rGADEM,BSgenome.Hsapiens.UCSC.hg19,IRanges)

path<- system.file("extdata/Test_100.bed",package="rGADEM")
BED<-read.table(path,header=FALSE,sep="\t")
BED<-data.frame(chr=as.factor(BED[,1]),start=as.numeric(BED[,2]),end=as.numeric(BED[,3]))
rgBED<-IRanges(start=BED[,2],end=BED[,3])
Sequences<-RangedData(rgBED,space=BED[,1])
# DHS= BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm9.masked, DHS_rng$V1, start=DHS_rng$start, end=DHS_rng$end)
gadem<-GADEM(Sequences,verbose=1,genome=Hsapiens)