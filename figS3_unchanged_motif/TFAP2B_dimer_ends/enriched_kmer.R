fjComm::clear_()

file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4-B20-TFAP2BoooNIIIc4-CGTTCTTT_S44_L003_R1_001.peared_trimmed.fq.gz"
TF_file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-TFctrl147c4-B20-TFAP2BoooNIIIc4-CGTTCTTT_S1580_L007_R1_001.peared_trimmed.fq.gz"
# file="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_147/Trulig147v1IIIFJ4-4-CAP2-E21-TF35-GSX2IIIc4_S501_L002_R1_001.peared_trimmed.fq_u.gz"


# seq=read_csv(file,col_names = F)
# seq=rmdup(seq,1)
# add_args=list(
#   seqs = seq[[1]],kmerLen = 3,filter_for_spacing =F ,spacing = 0:30,
#   verbose=F, pseudo=10, type = "maxBias", maxBias_dimer_Params=list(type="topMI",topNo=10L)
# )
# resultDf= oneOffCalc(paste0("_tmp/",basename(file)), calcFun = ic_related_calc, calcFunParamList =add_args, useScriptPath = T,asObj = T)
#
# pp=ggplot(data = resultDf)+ geom_point(aes(pos1-0.2,pos2,color=topk1),size=0.5)+ geom_point(aes(pos1+0.2,pos2,color=topk2),size=0.5)
# plotly::ggplotly(pp)
#
# gkcount=gkmerCntBit(strings = seq[[1]],gapNo = 1,k = 4L,gapMins = 0,gapMaxs = 30,pseudo = 10,diffLen = F,posInfo = F,all_possible_k = T) %>%  melt()
# gkcount= gkcount %>% dplyr::filter(Var2 %in% paste0(9:15,"n")) %>% arrange(desc(value))

xyp=xyplot(xfile = file,yfile = TF_file,test = F,kmerLen = 3,gapped = T,gapNo = 2,gapMins = c(0,0),gapMaxs = c(5,5),topKmersEach = 300)
xyp=xyp+xlab("Counts (nucleosomal DNA)")+ylab("Counts (free DNA)")
print(xyp)
gg_save_plotly(xyp)
xyp$layers[[1]]$mapping=NULL #delete the color
xyp$layers[[2]]=NULL
gg_save_pdf(xyp,width = 7,height = 6)

