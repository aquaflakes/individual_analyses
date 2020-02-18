# highest kmer pair MI against distance (for TF motif length and cross-strand binding)
# 2016-10-24
kmerLen=3

test=F
calcOverWrite=F

maxGap=ifelse(test,10,0)
testFileAdd= ifelse(test,"_test","")
verbose=F


# parameters to use when running on Mac
args_mac=c(
         reads_file="/Volumes/Nutcase_DataZhu/seqFiles2/FJ4.4CAP2_PE/allreads/FJ4.4_PE_147/dual_trim/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-4-CAP2-J18-TF345-DLX1IIIc4_S618_L002_R1_001.peared_trimmed.fq.gz",
         # reads_file="/Volumes/Nutcase_DataZhu/seqFiles2/FJ4.4CAP2_PE/allreads/FJ4.4_PE_147/dual_trim/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-4-CAP2-A06-TF99-RFX5IIIc4_S390_L002_R1_001.peared_trimmed.fq.gz",
         outDir="/Volumes/Nutcase_DataZhu/seqFiles2/FJ4.4CAP2_PE/analysis/maxBiasPair_allpos/147c4/"
         )

# head-----------
scriptName<-"mutualInfo.R"
    isMac=grepl("apple",version[1])
    RlibDir= ifelse(isMac, "/Volumes/Nutcase_DataZhu/lib/Rlib/", "/wrk/data/fangjie/lib/Rlib/")
       
#-----------
usage=paste(' ~/R/R320/bin/Rscript ', scriptName," '<readsfile>' <outputDir> \n",sep = "")
    args=args_mac; if(!isMac) {args=commandArgs(trailingOnly=TRUE)}
    if (!isMac && length(args)==0) {cat("usage:  ",usage); stop() } 
filePattern=args[1]
outDir=ifelse(is.na(args[2]), paste("./",scriptName,"/",sep = ""), args[2])

#---------------------------------------------------
    system(paste("mkdir -p ", outDir,sep = ""));
    system(paste("mkdir -p ", outDir,"/plot/",sep = ""))
    system(paste("mkdir -p ", outDir,"/one_off_calc/",sep = ""))
    system(paste("mkdir -p ", outDir,"/data/",sep = ""))
    cmd= paste("echo 'Rscript ",scriptName," ", paste(args,collapse = " "),"' >",outDir,"/Rcmd_",scriptName,".txt", sep = ""); system(cmd) # memo cmd
#---------------------end head----------------------#

targetDir=dirname(filePattern)
allfiles<- list.files(path=targetDir, pattern = basename(filePattern),full.names = F) #read in all txt and fastq files
tmpDir=paste(outDir,"/tmp/",sep = ""); system(paste("mkdir -p ",tmpDir,sep = ""));
#--------------------end fixed part------------------


# source("http://bioconductor.org/biocLite.R")  # for nutcase server
if (!require(pacman,quietly = T)){ install.packages("pacman") } 
pacman::p_load(readr,stringr,ggplot2,dplyr,reshape2,data.table,Rcpp,entropy)
source(paste0(RlibDir,"/Mylib/ggplot2_common.R"))
source(paste0(RlibDir,"/Mylib/fileIO.R")) # for oneOffCalc
source(paste0(RlibDir,"/Mylib/seqFile.R")) # for dedup
sourceCpp(paste0(RlibDir,"/common/kcnt.cpp")) 


################
# functions ---------
gapkCount_c0 <-function(c0File)
{
  allreads_c0= read_table(c0File,col_names = F) %>% .[[1]]
  if (test && (allreads_c0 %>% length() > 10000)) {allreads_c0=allreads_c0[1:10000]}
  
  pseudoCnt= min((length(allreads_c0)/4^(kmerLen*2)/20) %>% as.integer(),10)
  
  gapped_kCnt_c0= kmerCnt_allgap(allreads_c0, k = kmerLen, maxGap = maxGap) #maxGap=0 if no limit
  gapped_kCnt_df_c0= gapped_kCnt_c0 %>% 
    dplyr::group_by(kmer1) %>% mutate(counts= counts+ pseudoCnt) %>%
    mutate(prop_act_c0= counts/sum(counts)) %>% 
    ungroup
  return(gapped_kCnt_df_c0)
}

gapkCount <-function(file_full)
{

  if(isMac||test) allreads= read_table(file_full,col_names = F) else allreads= dedup(file_full) # dedup if running on server
  allreads= allreads[[1]]

  if (test && (allreads %>% length() > 10000)) {allreads=allreads[1:10000]}
  readsLen= allreads[1] %>% nchar()
  
  pseudoCnt= min((length(allreads)/4^(kmerLen*2)/20) %>% as.integer(),10)
  
  # bk from current cycle
  allkCnt= kmerCnt(allreads,kmerLen,collapse = T) + pseudoCnt  #counts for half of the gapped kmer
  prop_allkCnt=allkCnt/sum(allkCnt)
  
  # takd minimun diviation from c0 and curr cyc, calc fold change
  gapped_kCnt= kmerCnt_allgap(allreads, k = kmerLen, maxGap = maxGap)
  
  # calc foldchn_VS_mindiff
  gapped_kCnt_df= gapped_kCnt %>% 
    dplyr::group_by(kmer1) %>% mutate(counts= counts+ pseudoCnt) %>%
    mutate(prop_act= counts/sum(counts)) %>% 
    ungroup %>% 
    mutate(prop_exp_curr= prop_allkCnt[kmer2], prop_exp_c0=gapped_kCnt_df_c0[match(.$kmer,gapped_kCnt_df_c0$kmer), ]$prop_act_c0) %>%
    mutate(min_diff_exp= ifelse(abs(prop_exp_curr-prop_act)> abs(prop_exp_c0-prop_act), prop_exp_c0, prop_exp_curr), foldchn_VS_mindiff= log2(prop_act/min_diff_exp)) %>% 
    arrange(foldchn_VS_mindiff)
  return(gapped_kCnt_df)
}

MIcalc <- function (gap_curr,punish,data) 
{
  df=data %>% dplyr::filter(gapLen==gap_curr)
  t= xtabs(counts~kmer1+kmer2, data=df) %>% prop.table()
  # t= t[!grepl("N",rownames(t)),!grepl("N",colnames(t))]
  mi.empirical(t,unit ="log2")+punish
}

MI_calc_wrap <- function(gapLen, data)
{
  MI_bits=rep(NA,length(gapLen))

  # punish= -(4^(kmerLen*2)-1)/(2*log(2)*nrow(allseq)) + 2* ((4^(kmerLen)-1)/(2*log(2)*nrow(allseq)))
  punish=0
  for (i in 1:length(gapLen)) 
  {
    gap_curr=gapLen[i]
    if(verbose) {cat(paste0("curr gapLen: ",gap_curr,"\n"))}
    MI_bits[i]= MIcalc(gap_curr,punish, data)
  } 
  return(MI_bits)
}

file_name=allfiles[1];
  cat(paste("processing ",file_name,"\n",sep = "")); 
  file_full= paste(targetDir, file_name,sep = "/") #full path
  
  saved_oneOff_filename= paste0("/data/",basename(file_full),testFileAdd,"_k",kmerLen,"_maxGap",maxGap,".csv")
  gapped_kCnt_df_CAP= oneOffCalc(saved_oneOff_filename, calcFun =gapkCount, calcFunParamList =list(file_full) )
  gapped_kCnt_df_TF= read_csv("/Volumes/Nutcase_DataZhu/seqFiles2/FJ4.4CAP2_PE/analysis/maxBiasPair_allpos/TFctrlc3/data/Trulig147v1IIIFJ4-4-TFctrl-J18-TF345-DLX1IIIc3_S234_L004_R1_001.peared_trimmed.fq.gz_k3_maxGap0.csv",col_names = T)
  gapped_kCnt_df_Nu= read_csv("/Volumes/Nutcase_DataZhu/seqFiles2/FJ6.1Methyl_PE/analysis/maxBiasPair_allpos/FJ6.1AllPosBiasC5/data/Trulig147v1IIIFJ6-1-OriMe1-A02-1xCAGOct64ngIIIc4-Bq4-TTGGTGGT_S770_L008_R1_001.peared_trimmed.fq.gz_k3_maxGap0.csv",col_names = T)
  
  
  # begin to do work
  seqLen= 101;  kmerNo=as.integer(1+seqLen-kmerLen)
  # allseq=read_csv(file_full,col_names = F)
  gapLen=0:(max(gapped_kCnt_df_CAP$gapLen))

  
  MI_df= data.frame(gapLen=gapLen, DLX1_and_Nucleosome=MI_calc_wrap(gapLen,gapped_kCnt_df_CAP), DLX1=MI_calc_wrap(gapLen,gapped_kCnt_df_TF), nucleosome=MI_calc_wrap(gapLen,gapped_kCnt_df_Nu))
  plot_df=melt(MI_df,id="gapLen")
  colnames(plot_df)[2:3]=c("exp","MI")
  
  
  
  
  
  
  
  
  
  
  p=ggplot() + trans_bk_xyaxis + geom_line(data=plot_df,aes(x=gapLen,y=log10(MI), color=exp)) +
    scale_y_continuous(expand = c(0,0))+ scale_x_continuous(expand = c(0,0))+
    geom_vline(xintercept = c(10,20,30,40,50,60,70,80,90)-3,size=0.4,color="grey")+theme(legend.position = c(.60, .80)) + labs(color="experiment")+
    scale_color_manual(values=c("red","blue", "black"))

  print(p)
  
  # system(paste("mkdir -p ", outDir,"/MI_gap_plot/",sep = ""))
  img="/Users/zhu/Desktop/Nucleosome_manuscript/plotting_Scripts/MI_gapLen/DLX.pdf"
  ggsave(filename =img ,plot = p,width = 4,height =3 )
  system( paste0("inkscape -f ",img," -b white -d 200 --export-png=",img,".png") )


# }
  
