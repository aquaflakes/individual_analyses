fjComm::clear_()

# T
file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/N24_Trulig147v1IIIFJ4-5-CAP2r-147c5Nu-N24-ToooN2XIIIc5-GACGTATG_S720_L004_R1_001.peared_trimmed.fq.gz_3mer_foldchn_maxBias.csv"
# seq= readr::read_csv(file,col_names = F)
# seq=rmdup(seq,1)
# result_TFMI=ic_related_calc(seqs = seq[[1]],kmerLen = 3,filter_for_spacing = F,spacing = 0:10,type = "maxBias",maxBias_dimer_Params=list(type="topMI",topNo=10L))
result_TFMI=readr::read_csv(file,col_names = T)
FJ4.4_lig200_TFMI_P10_T= gg_heat2D_MI(result_TFMI,grad_colors = gg_steelblue_red)+guides(fill=F)
print(FJ4.4_lig200_TFMI_P10_T)
gg_save_all(FJ4.4_lig200_TFMI_P10_T,width = 5.8,height = 5.5,newName = "T_lig147_2DEMI")

# # T-TFctrl
# file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig200v1IIIFJ4-5-CAP2r-TFctrl200c3-P10-TIIIc3-CATCACGC_S2674_L008_R1_001.peared_trimmed.fq.gz"
# seq= readr::read_csv(file,col_names = F)
# seq=rmdup(seq,1)
# result_TFMI=ic_related_calc(seqs = seq[[1]],kmerLen = 3,filter_for_spacing = F,spacing = 0:10,type = "maxBias",maxBias_dimer_Params=list(type="topMI",topNo=10L))
# FJ4.5_lig200_TFMI_P10_T_TFctrl= gg_heat2D_MI(result_TFMI,grad_colors = gg_steelblue_red)+guides(fill=F)
# print(FJ4.5_lig200_TFMI_P10_T_TFctrl)
# gg_save_all(FJ4.4_lig200_TFMI_P10_T,width = 5.8,height = 5.5)

readsfile="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c5Nu-N24-ToooN2XIIIc5-GACGTATG_S720_L004_R1_001.peared_trimmed.fq.gz"
seq= readr::read_csv(readsfile,col_names = F)
# seq=rmdup(seq,1)
motif=pfm_from_seed(seq[[1]],seed1 = "CAC",gapLen = 76,two_strands = T,seed2 = "GTG",flankLen = 4,all_start_with_specified_gap = F,seed1_start = 11,rmdup_fg = F)
pdf("T_lig147_motif.pdf",height = 4)
fjComm::plotMotif_pfmMat(motif)
dev.off()
