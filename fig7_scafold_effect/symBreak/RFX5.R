library(fjComm)
setwd(get_scriptpath())

pfmfile= "RFX5_motif_FJ4.4.txt"
file="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_147/Trulig147v1IIIFJ4-4-CAP2-A06-TF99-RFX5IIIc4_S390_L002_R1_001.peared_trimmed.fq_u.gz"
RFX5_map_both_str= motif_plot_moods_singlePWM_hit(file,pfmfile,dup_rm = F,p_moods = "0.0001",colorGrad = c("#ceffce","red"),strand_filter = c("+","-"),combine_2strand_hits = F)
print(RFX5_map_both_str)
gg_save_all(RFX5_map_both_str,width = 4,height = 1.5)

# pdf("motif.pdf",height = 3)
# plotMotif_pfmFile(pfmfile,ic.scale = T)
# dev.off()

file="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_TFctrl/Trulig147v1IIIFJ4-4-TFctrl-A06-TF99-RFX5IIIc3_S6_L004_R1_001.peared_trimmed.fq_u.gz"
RFX5_TFctrl_map_both_str= motif_plot_moods_singlePWM_hit(file,pfmfile,dup_rm = F,p_moods = "0.0001",colorGrad = c("#ceffce","red"),strand_filter = c("+","-"),combine_2strand_hits = F)
print(RFX5_TFctrl_map_both_str)
gg_save_all(RFX5_TFctrl_map_both_str,width = 4,height = 1.5)
