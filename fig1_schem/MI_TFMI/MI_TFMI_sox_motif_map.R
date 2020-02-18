fjComm::clear_()
# file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4-C18-RFX5oooNIIIc4-TATATCTG_S66_L003_R1_001.peared_trimmed.fq.gz"
# file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4-E3-GSX1IIIc4-TAAGGTCG_S99_L003_R1_001.peared_trimmed.fq.gz"
# file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4-K3-VSX2IIIc4-CTGCATTT_S243_L003_R1_001.peared_trimmed.fq.gz"
# file="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_147/Trulig147v1IIIFJ4-4-CAP2-E21-TF35-GSX2IIIc4_S501_L002_R1_001.peared_trimmed.fq_u.gz"
# file="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_147/Trulig147v1IIIFJ4-4-CAP2-O19-TF94-EVX1IIIc4_S739_L002_R1_001.peared_trimmed.fq_u.gz"
file="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_147/Trulig147v1IIIFJ4-4-CAP2-C05-TF15-SOX11IIIc4_S437_L002_R1_001.peared_trimmed.fq_u.gz"


# pfmfile="SOX.pfm.txt"
# pfmfile="~/Nut_zhuData/SELEXphp/data/pfm/shortspacek/*_AR_GAACAATG*"
# pfmfile= Sys.glob(pfmfile)[1]
pfmfile= motif_getpfm_from_TFname("SOX11")$file[1]
SOX11_motif_map=fjComm::motif_plot_moods_singlePWM_hit(readsfile = file,pfmfile = pfmfile,dup_rm = T,combine_2strand_hits = T)
gg_save_diag(SOX11_motif_map)
print(SOX11_motif_map)

pdf("SOX11_motif.pdf",height = 3)
fjComm::plotMotif_pfmFile(pfmfile,ic.scale = T)
dev.off()


seq= readr::read_csv(file,col_names = F)
seq=rmdup(seq,1)
result_MI= oneOffCalc(paste0(basename(file),"_MI"),calcFun = ic_related_calc, calcFunParamList = list(seqs = seq[[1]],kmerLen = 3,filter_for_spacing = F,spacing = 0:10,type = "MI",maxBias_dimer_Params=list(type="topMI",topNo=10L)),asObj = T, useScriptPath = T   )
MI_plot= gg_heat2D_MI(result_MI,grad_colors = gg_steelblue_red)+guides(fill=F)
print(MI_plot)
gg_save_all(MI_plot,width = 5.2,height = 5)


result_TFMI= oneOffCalc(paste0(basename(file),"_TFMI"),calcFun = ic_related_calc, calcFunParamList = list(seqs = seq[[1]],kmerLen = 3,filter_for_spacing = F,spacing = 0:10,type = "maxBias",maxBias_dimer_Params=list(type="topMI",topNo=10L)), asObj = T, useScriptPath = T)
# result_TFMI=ic_related_calc(seqs = seq[[1]],kmerLen = 3,filter_for_spacing = F,spacing = 0:10,type = "maxBias",maxBias_dimer_Params=list(type="topMI",topNo=10L))
TFMI_plot= gg_heat2D_MI(result_TFMI,grad_colors = gg_steelblue_red)+guides(fill=F)+gg_theme_transparent
print(TFMI_plot)
gg_save_all(TFMI_plot,width = 5.2,height = 5.2)
TFMI_diag= gg_heat2D_diag(result_TFMI,grad_colors = gg_steelblue_red)+gg_theme_bordered_diag
gg_save_diag(TFMI_diag)

