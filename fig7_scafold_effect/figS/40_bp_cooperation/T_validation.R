fjComm::clear_()

fg=read_csv("~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/FJ4.4_PE_200/dual_trim/2_Trim_adaptor_kmer/Trulig200v1IIIFJ4-4-CAP2-P10-TF377-TIIIc4_S754_L002_R1_001.peared_trimmed.fq.gz")
bg=read_csv("~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/FJ4.4_PE_200/dual_trim/2_Trim_adaptor_kmer/Trulig200v1IIIFJ4-4-CAP2-P10-TF377-TIIIc3_S370_L001_R1_001.peared_trimmed.fq.gz")


# seq=read_csv("~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/FJ4.4_PE_200/dual_trim/2_Trim_adaptor_kmer/Trulig200v1IIIFJ4-4-CAP2-P10-TF377-TIIIc4_S754_L002_R1_001.peared_trimmed.fq.gz")
# seq=rmdup(seq,1)[[1]]

result=oneOffCalc("/tmp/T_maxBias.Robj", fjComm::ic_related_calc, list(seqs = seq,filter_for_spacing = F,kmerLen = 3,type = "maxBias",  maxBias_dimer_Params=list(type="foldchn",topNo=10L)) )

p=fjComm::gg_heat2D_MI(result)
# print(plotly::ggplotly(p))

targetRng = result %>%  dplyr::filter(pos2-pos1<=90 & pos2-pos1>=70) %>% mutate(topk1=as.character(topk1),topk2=as.character(topk2))
topRec= targetRng%>% arrange(desc(foldchnTop)) %>% .[1,]

pfm=pfm_from_seed(seqs_or_file = fg, seqs_or_file_bg = NA, seed1 = "CAC", gapLen = topRec[[2]][1]-topRec[[1]][1]-nchar(topRec[[5]][1])-1, seed2 = "GTG", seed1_start = topRec$pos1,flankLen = 4,all_start_with_specified_gap=T,two_strands = T)


# pfm=pfm_from_seed(seqs = seq,seed1 = "TAA", gapLen = 0, seed2 = "CACCT",flankLen = 0, all_start_with_specified_gap=T,two_strands = T)

fjComm::plotMotif_pfmMat(pfm )#+(pfm %>% fjComm::matRevComp()))

pdf(paste0("T_motif.pdf"),height = 3)
fjComm::plotMotif_pfmMat(pfm)
dev.off()
