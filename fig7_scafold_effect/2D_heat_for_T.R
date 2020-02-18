fjComm::clear_()
# file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4-C18-RFX5oooNIIIc4-TATATCTG_S66_L003_R1_001.peared_trimmed.fq.gz"
file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4-E3-GSX1IIIc4-TAAGGTCG_S99_L003_R1_001.peared_trimmed.fq.gz"
# file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c5Nu-P8-ToooN0-5IIIc5-TCTTAGGT_S752_L004_R1_001.peared_trimmed.fq.gz"
# file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c5Nu-H11-ETV1IIIc5-GCTTTGGA_S563_L004_R1_001.peared_trimmed.fq.gz"

# T
file="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_200/Trulig200v1IIIFJ4-4-CAP2-P10-TF377-TIIIc4_S754_L002_R1_001.peared_trimmed.fq_u.gz"
seq= readr::read_csv(file,col_names = F)
seq=rmdup(seq,1)
result_TFMI=ic_related_calc(seqs = seq[[1]],kmerLen = 3,filter_for_spacing = F,spacing = 0:10,type = "maxBias",maxBias_dimer_Params=list(type="topMI",topNo=10L))
FJ4.4_lig200_TFMI_P10_T= gg_heat2D_MI(result_TFMI,grad_colors = gg_steelblue_red)+guides(fill=F)
print(FJ4.4_lig200_TFMI_P10_T)
gg_save_all(FJ4.4_lig200_TFMI_P10_T,width = 5.8,height = 5.5)

# T-TFctrl
file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig200v1IIIFJ4-5-CAP2r-TFctrl200c3-P10-TIIIc3-CATCACGC_S2674_L008_R1_001.peared_trimmed.fq.gz"
seq= readr::read_csv(file,col_names = F)
seq=rmdup(seq,1)
result_TFMI=ic_related_calc(seqs = seq[[1]],kmerLen = 3,filter_for_spacing = F,spacing = 0:10,type = "maxBias",maxBias_dimer_Params=list(type="topMI",topNo=10L))
FJ4.5_lig200_TFMI_P10_T_TFctrl= gg_heat2D_MI(result_TFMI,grad_colors = gg_steelblue_red)+guides(fill=F)
print(FJ4.5_lig200_TFMI_P10_T_TFctrl)
gg_save_all(FJ4.4_lig200_TFMI_P10_T,width = 5.8,height = 5.5)
