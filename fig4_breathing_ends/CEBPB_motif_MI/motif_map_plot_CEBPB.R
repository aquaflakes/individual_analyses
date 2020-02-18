fjComm::clear_()


MIfile="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/O16_Trulig147v1IIIFJ4-5-CAP2r-147c4-O16-CEBPBIIIc4-ACGGTAAT_S352_L003_R1_001.peared_trimmed.fq.gz_3mer_foldchn_maxBias.csv"
df=read_csv(MIfile)
CEBPB_MI_diag= gg_heat2D_diag(df,grad_colors = gg_steelblue_red)+gg_theme_bordered_diag
gg_save_diag(CEBPB_MI_diag)

MIfile="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/O16_Trulig147v1IIIFJ4-5-CAP2r-TFctrl147c4-O16-CEBPBIIIc4-ACGGTAAT_S1888_L007_R1_001.peared_trimmed.fq.gz_3mer_foldchn_maxBias.csv"
df=read_csv(MIfile)
CEBPB_MI_diag_TFctrl= gg_heat2D_diag(df,grad_colors = gg_steelblue_red)+gg_theme_bordered_diag
gg_save_diag(CEBPB_MI_diag_TFctrl)



motif_idx_to_use=1
TFname= "CEBPB"
file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4-O16-CEBPBIIIc4-ACGGTAAT_S352_L003_R1_001.peared_trimmed.fq.gz"
dup_rm= T

library(fjComm)
setwd(get_scriptpath())


all_pfmfile= motif_getpfm_from_TFname(TFname)
pfmfile=all_pfmfile[motif_idx_to_use,]$file
pfm_title=all_pfmfile[motif_idx_to_use,]$pfmName
motif_map=fjComm::motif_plot_moods_singlePWM_hit(readsfile = file,pfmfile = pfmfile,dup_rm = dup_rm,combine_2strand_hits = T)
gg_save_diag(motif_map)
print(motif_map)

pdf(paste0(TFname,"_motif.pdf"),height = 3)
fjComm::plotMotif_pfmFile(pfmfile,ic.scale = T,title = pfm_title)
dev.off()




motif_idx_to_use=1
TFname= "CEBPB"
file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-TFctrl147c4-O16-CEBPBIIIc4-ACGGTAAT_S1888_L007_R1_001.peared_trimmed.fq.gz"
dup_rm= T

library(fjComm)
setwd(get_scriptpath())


all_pfmfile= motif_getpfm_from_TFname(TFname)
pfmfile=all_pfmfile[motif_idx_to_use,]$file
pfm_title=all_pfmfile[motif_idx_to_use,]$pfmName
motif_map=fjComm::motif_plot_moods_singlePWM_hit(readsfile = file,pfmfile = pfmfile,dup_rm = dup_rm,combine_2strand_hits = T)
gg_save_diag(motif_map,newName = "motif_map_TFctrl")
print(motif_map)

pdf(paste0(TFname,"_motif_TFctrl.pdf"),height = 3)
fjComm::plotMotif_pfmFile(pfmfile,ic.scale = T,title = pfm_title)
dev.off()
