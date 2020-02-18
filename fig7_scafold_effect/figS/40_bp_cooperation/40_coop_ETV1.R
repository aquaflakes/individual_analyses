fjComm::clear_()

minSpacing=28; maxSpacing=41
TFname="ETV1"
use_nth_top_as_seed=3

seq=read_csv("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4-H11-ETV1IIIc4-GCTTTGGA_S179_L003_R1_001.peared_trimmed.fq.gz")
seq=rmdup(seq,1)[[1]]
result=fjComm::ic_related_calc(seqs = seq,filter_for_spacing = F,kmerLen = 3,type = "maxBias",  maxBias_dimer_Params=list(type="foldchn",topNo=10L))
p=fjComm::gg_heat2D_MI(result)
# print(plotly::ggplotly(p))

targetRng = result %>%  dplyr::filter(pos2-pos1<=maxSpacing & pos2-pos1>=minSpacing) %>% mutate(topk1=as.character(topk1),topk2=as.character(topk2))
topRec= targetRng%>% arrange(desc(foldchnTop)) %>% .[use_nth_top_as_seed,]

gapLen = topRec[[2]][1]-topRec[[1]][1]-nchar(topRec[[5]][1])
Rcpp::sourceCpp("posed_seed.cpp")
pfm=pfm_from_seed(seqs = seq,seed1 = topRec[1,"topk1"],gapLen = gapLen, seed2 = topRec[1,"topk2"],seed1_start = topRec$pos1,flankLen = 4,all_start_with_specified_gap = F)

seed=paste0( topRec[[1]][1],"_",topRec[[2]][1],"_",topRec[[5]][1],gapLen,"N",topRec[[6]][1])
pdf(paste0(TFname,"_",seed,"_motif.pdf"),height = 3)
fjComm::plotMotif_pfmMat(pfm)
dev.off()

