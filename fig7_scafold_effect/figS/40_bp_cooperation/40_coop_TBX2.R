fjComm::clear_()

minSpacing=28; maxSpacing=50
TFname="TBX2"
use_nth_top_as_seed=1

# seq=read_csv("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4-H6-TBX2oooN0-5IIIc4-ACTGTCAC_S174_L003_R1_001.peared_trimmed.fq.gz")
# seq=rmdup(seq,1)[[1]]
# result=fjComm::ic_related_calc(seqs = seq,filter_for_spacing = F,kmerLen = 3,type = "maxBias",  maxBias_dimer_Params=list(type="topMI",topNo=10L))
seq=read_csv("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c5Nu-H6-TBX2oooN0-5IIIc5-ACTGTCAC_S558_L004_R1_001.peared_trimmed.fq.gz")
seq=rmdup(seq,1)[[1]]
all147=fjComm::nusel_get_147() %>% dplyr::filter(use==1 & order==1 & TF=="TBX2")
result=read_csv(all147$c5_EMI_Nu_file[1]) %>% as.data.frame()
p=fjComm::gg_heat2D_MI(result,grad_colors = gg_steelblue_red)+guides(fill=F)
gg_save_pdf(p,5.8,5.5,filename = "TBX2_c5Nu")
ggsave(file=img_p1, plot = p_enrich, width = 5.8/2.54,height = 5.5/2.54)
# print(plotly::ggplotly(p))

targetRng = result %>%  dplyr::filter(pos2-pos1<=maxSpacing & pos2-pos1>=minSpacing) %>% mutate(topk1=as.character(topk1),topk2=as.character(topk2))
topRec= targetRng%>% arrange(desc(topMIsum)) %>% .[use_nth_top_as_seed,]

gapLen = topRec[[2]][1]-topRec[[1]][1]-nchar(topRec[[5]][1])
# Rcpp::sourceCpp("posed_seed.cpp")
pfm=pfm_from_seed(seqs = seq,seed1 = topRec[1,"topk1"],gapLen = gapLen, seed2 = topRec[1,"topk2"],seed1_start = topRec$pos1,flankLen = 4)

seed=paste0( topRec[[1]][1],"_",topRec[[2]][1],"_",topRec[[5]][1],gapLen,"N",topRec[[6]][1])
pdf(paste0(TFname,"_",seed,"_motif.pdf"),height = 3)
fjComm::plotMotif_pfmMat(pfm)
dev.off()



minSpacing=70; maxSpacing=90
targetRng = result %>%  dplyr::filter(pos2-pos1<=maxSpacing & pos2-pos1>=minSpacing) %>% mutate(topk1=as.character(topk1),topk2=as.character(topk2))
topRec= targetRng%>% arrange(desc(topMIsum)) %>% head(1)#%>% dplyr::filter(topk1=="GTG"& topk2=="GTG") %>% head(1) #.[3,]

gapLen = topRec[[2]][1]-topRec[[1]][1]-nchar(topRec[[5]][1])
# Rcpp::sourceCpp("posed_seed.cpp")
pfm=pfm_from_seed(seqs = seq,seed1 = topRec[1,"topk1"],gapLen = gapLen, seed2 = topRec[1,"topk2"],seed1_start = topRec$pos1,flankLen = 4)

seed=paste0( topRec[[1]][1],"_",topRec[[2]][1],"_",topRec[[5]][1],gapLen,"N",topRec[[6]][1])
pdf(paste0(TFname,"_",seed,"_motif.pdf"),height = 3)
fjComm::plotMotif_pfmMat(pfm)
dev.off()

