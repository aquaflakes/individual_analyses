# run maxBiasPair_posSpecific.R first to generate loess of maxfoldchnTop

kmerLen=3


setwd("/Users/zhu/Desktop/Nucleosome_manuscript/plotting_Scripts/penetration_clustering/")

# source("http://bioconductor.org/biocLite.R")  # for nutcase server
if (!require(pacman,quietly = T)){ install.packages(pacman) } 
pacman::p_load(readr,stringr,ggplot2,dplyr,reshape2) #dplyr

isMac=grepl("apple",version[1])
RlibDir= ifelse(isMac, "/Volumes/Nutcase_DataZhu/lib/Rlib/", "/wrk/data/fangjie/lib/Rlib/")
source(paste0(RlibDir,"/Mylib/ggplot2_common.R"))
source(paste0(RlibDir,"/Mylib/fileIO.R")) # for oneOffCalc
# sourceCpp(paste0(RlibDir,"/common/kcnt.cpp")) 


guide200=read_csv("guide200.csv", col_names = T) %>% dplyr::filter(SigClustering==1)
guide200$mean_penetration= "NA"


for (i in 1:nrow(guide200))
{
  file=Sys.glob(paste0("/Volumes/Nutcase_DataZhu/seqFiles2/FJ4.4CAP2_PE/analysis/maxBiasPair_posSpecific/200c4/loess_enrich/",guide200[i,]$pos,"_*"))
  curr_loess=read_csv(file, col_names = T)
  left= curr_loess %>% dplyr::filter(pos1<=74)
  left_penetration= left %>% dplyr::filter(loess0.45>= (min(left$loess0.45)+max(left$loess0.45))/2) %>% .$pos1 %>% max + kmerLen-1+0.5
  
  right= curr_loess %>% dplyr::filter(pos1>=74)
  right_penetration= 154-(right %>% dplyr::filter(loess0.45>= (min(right$loess0.45)+max(right$loess0.45))/2) %>% .$pos1 %>% min + kmerLen-1+0.5)
  
  avg_penetration= mean(c(right_penetration,left_penetration))
  guide200[i,"mean_penetration"]=avg_penetration
  
}

guide200=guide200 %>% mutate(family=as.factor(family), mean_penetration=as.numeric(mean_penetration))
p=ggplot(guide200, aes(family, mean_penetration, color=class)) +geom_boxplot() + geom_point(aes(text=TF)) + theme_bw() + scale_y_continuous(expand = c(0,0)) + 
  labs(color=("TF class")) + xlab("TF family")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
# library(plotly)
# p1=ggplotly(p)
# print(p1)
img_p1=paste0(getwd(),"/TF_penetration.pdf")
ggsave(file=img_p1, plot = p, width = 10,height = 8)
system( paste0("inkscape -f ",img_p1," -b white -d 200 --export-png=",img_p1,".png") )

