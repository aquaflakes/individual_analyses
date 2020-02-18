fjComm::clear_()

max_pos_num=10

file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/D2_Trulig147v1IIIFJ4-5-CAP2r-147c4-D2-EOMESoooNIIIc4-TTAGAGCT_S74_L003_R1_001.peared_trimmed.fq.gz_3mer_foldchn_maxBias.csv"
seqfile= "~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4-D2-EOMESoooNIIIc4-TTAGAGCT_S74_L003_R1_001.peared_trimmed.fq.gz"
# c0file= "~/Nut_zhuData/seqFiles2/c0/subset_100000_Trulig147v1IIIPEx200IIIc0_S3_R1_001.peared_trimmed.fq.gz"

resultDf= read_csv(file)
  result_diag= resultDf %>% dplyr::filter(pos2-pos1==3) %>% dplyr::filter(pos1>5 & pos1<max(pos1)-4)
  max_pos1_starts= result_diag %>% top_n(10,topMIsum) %>% .$pos1

seq= SELEXFile(seqfile)
seq$getSeq(dup_rm = T)
seq$count_k(k = 6,collapse = F,diffLen = F,asDf = T,pseudo = 5,all_possible_k = T,cmp_c0 = F,count_c0 = T,rc = T)

kcnt_topMI_pos=seq$kmerCnt %>% slice(max_pos1_starts) %>% colSums()
top_kmer= kcnt_topMI_pos %>% sort() %>% tail(20) %>% names()
print(top_kmer)

rio::export(kcnt_topMI_pos %>% sort() %>% tail(20) %>% {data.frame(kmer=names(.),count=.)}, file = paste0(script_path_from_fun,"/EOMES_top_kmers.csv"),format = "csv")


# p=ggplot(data = resultDf)+ geom_point(aes(pos1-0.2,pos2,color=topk1),size=0.5)+ geom_point(aes(pos1+0.2,pos2,color=topk2),size=0.5)
# plotly::ggplotly(p)

# p=gg_heat2D_MI(resultDf)+ geom_text(data=resultDf,aes(pos1,pos2,label=paste0(topk1,"_",pos2-pos1-3,"N","_",topk2)),alpha=0)
# plotly::ggplotly(p)

file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/P14_Trulig147v1IIIFJ4-5-CAP2r-147c4-P14-PITX2IIIc4-GGGTCGTG_S374_L003_R1_001.peared_trimmed.fq.gz_3mer_foldchn_maxBias.csv"
seqfile= "~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4-P14-PITX2IIIc4-GGGTCGTG_S374_L003_R1_001.peared_trimmed.fq.gz"
# c0file= "~/Nut_zhuData/seqFiles2/c0/subset_100000_Trulig147v1IIIPEx200IIIc0_S3_R1_001.peared_trimmed.fq.gz"

resultDf= read_csv(file)
result_diag= resultDf %>% dplyr::filter(pos2-pos1==3) %>% dplyr::filter(pos1>5 & pos1<max(pos1)-4)
max_pos1_starts= result_diag %>% top_n(10,topMIsum) %>% .$pos1

seq= SELEXFile(seqfile)
seq$getSeq(dup_rm = T)
seq$count_k(k = 6,collapse = F,diffLen = F,asDf = T,pseudo = 5,all_possible_k = T,cmp_c0 = F,count_c0 = T,rc = T)

kcnt_topMI_pos=seq$kmerCnt %>% slice(max_pos1_starts) %>% colSums()
top_kmer= kcnt_topMI_pos %>% sort() %>% tail(20) %>% names()
print(top_kmer)

rio::export(kcnt_topMI_pos %>% sort() %>% tail(20) %>% {data.frame(kmer=names(.),count=.)}, file = paste0(script_path_from_fun,"/PITX_top_kmers.csv"),format = "csv")
