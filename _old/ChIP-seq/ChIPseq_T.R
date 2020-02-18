# find the signal of T
verbose=F
test=F


args_mac=c(
  DHS_file="/Volumes/Nutcase_DataZhu/forFan/DHS_score/5k/forhit_l8b4/tmphide/mSpleen_intersection.bed",
  kmerTabFile="/Volumes/Nutcase_DataZhu/Analysis/DHS_prediction/BeiData/kcnt/cycle4_Spleen_sig_10mer_u_cmp_cycle1_spleen_sig_10mer_u.txt",

  outDir="/Volumes/Nutcase_DataZhu/Analysis/DHS_prediction/CENTER_PRED_FINAL_v1/test/center_pred_10mer_Spleen/"
)

# head-----------
scriptName<-"maxBiasPair_posSpecific"
isMac=grepl("apple",version[1])
RlibDir= ifelse(isMac, "/Volumes/Nutcase_DataZhu/lib/Rlib/", "/wrk/data/fangjie/lib/Rlib/")

#-----------
usage=paste(' ~/R/R320/bin/Rscript ', scriptName," '<DHS_file>' '<kmerTabFile>' <outputDir> \n",sep = "")
args=args_mac; if(!isMac) {args=commandArgs(trailingOnly=TRUE)}
if (!isMac && length(args)==0) {cat("usage:  ",usage); stop() } 
filePattern=args[1]
kmerTabFile= args[2]
outDir=ifelse(is.na(args[3]), paste("./",scriptName,"/",sep = ""), args[3])

#---------------------------------------------------
system(paste("mkdir -p ", outDir,sep = ""));
system(paste("mkdir -p ", outDir,"/plot/",sep = "")); # make plot dir
system(paste("mkdir -p ", outDir,"/data/",sep = "")); # make plot dir
cmd= paste("echo 'Rscript ",scriptName," ", paste(args,collapse = " "),"' >",outDir,"/Rcmd_",scriptName,".txt", sep = ""); system(cmd) # memo cmd
#---------------------end head----------------------#

targetDir=dirname(filePattern)
allfiles<- list.files(path=targetDir, pattern = basename(filePattern),full.names = F)#read in all txt and fastq files
tmpDir=paste(outDir,"/tmp/",sep = ""); system(paste("mkdir -p ",tmpDir,sep = ""))

#--------------------end fixed part------------------


# source("http://bioconductor.org/biocLite.R")  # for nutcase server
if (!require(pacman,quietly = T)){ install.packages(pacman) } 
pacman::p_load(readr,stringr,ggplot2,dplyr,reshape2) #dplyr

source(paste0(RlibDir,"/Mylib/ggplot2_common.R"))
source(paste0(RlibDir,"/Mylib/fileIO.R")) # for oneOffCalc
# sourceCpp(paste0(RlibDir,"/common/kcnt.cpp")) 



################
# functions ---------

k_col_gen <- function(allseq1,seqLen1,kmerLen1)
{  
  # construct df cols with kmerLen
  allkmer1=data.frame() 
  kmerNo1=as.integer(1+seqLen1-kmerLen1); if(test){kmerNo1=kNotest}
  for (i in 1:(kmerNo1)) 
  {
    curr=vector()
    for (j in 1:kmerLen1) {curr=paste(curr,allseq1[[i+j-1]],sep = "")}
    
    #create df if not exist, add to it for later cycs
    if(nrow(allkmer1)==0) {allkmer1=as.data.frame(as.character(curr))} else {allkmer1[[i]]= curr} 
    if(verbose) {cat(i)}
  }
  colnames(allkmer1)[1]="V1"; allkmer1[[1]]=as.character(allkmer1[[1]])
  allkmer1
}

calcMaxBias <- function(topNo)
{
  if(!isMac)
  {
    # rmdup_3point
    tmpDir=paste(outDir,"/tmp/",sep = ""); system(paste("mkdir -p ",tmpDir,sep = ""));
    cmd=(paste( "rmdup_3point_tmp.pl -f ",file," -k 20 -o ",tmpDir , sep = "")); system(cmd);
    # cat (cmd); cat("\n"); 
    tmpFile= paste(outDir,"/tmp/",filei, sep = ""); # for tmp fasta file
    file=tmpFile
  }
  
  # begin to do work
  allseq=read_fwf(file, fwf_widths(  rep(1,seqLen) , as.character(c(1:seqLen)) ), col_types = cols(.default = col_character())  ) # seperation, colname
  
  allkmer=k_col_gen(allseq,seqLen,kmerLen)
  pseudoCnt=nrow(allkmer)/4^(kmerLen*2)/20 
  
  
  allcomb= combn(c(1:kmerNo),2) 
  allcomb= allcomb[ ,(allcomb[2,]-allcomb[1,]>=kmerLen)] # only spaced >= kmerLen has no interference
  df <- data.frame(pos1=allcomb[1,], pos2=allcomb[2,])
  df[,c("maxfoldchnTop","maxfoldchnBottom")]= apply(df, 1, function(x) {calcMaxBiasSingle(x[1],x[2],topNo,allkmer,pseudoCnt)} ) %>% t
  
  if(!isMac){system(paste("rm ",tmpFile,sep = ""))}
  return(df)
}

calcMaxBiasSingle <- function(col1,col2,topNo,allkmer,pseudoCnt)
{
  actFreq= table(allkmer[[col1]],allkmer[[col2]])
  actFreq= (actFreq[!grepl("N",rownames(actFreq)),!grepl("N",colnames(actFreq))] + pseudoCnt) %>% prop.table 
  expFreq= margin.table(actFreq,1) %o% margin.table(actFreq,2) %>% as.matrix() %>% melt
  actFreq= melt(actFreq)
  actFreq= mutate(actFreq, foldchn= log2(actFreq$value/ expFreq$value)) %>% arrange(foldchn)
  bottomDevMean= head(actFreq,topNo)$foldchn %>% mean
  topDevMean= tail(actFreq,topNo)$foldchn %>% mean
  return(c(topDevMean,bottomDevMean))
}



# body ========