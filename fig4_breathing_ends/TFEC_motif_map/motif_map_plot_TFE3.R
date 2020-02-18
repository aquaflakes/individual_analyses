fjComm::clear_()
library(fjComm)
setwd(get_scriptpath())

motif_idx_to_use=1
TFname= "TFE3"
file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4-N4-TFE3IIIc4-AAATTTTG_S316_L003_R1_001.peared_trimmed.fq.gz"
dup_rm= T

library(fjComm)
setwd(get_scriptpath())


all_pfmfile= motif_getpfm_from_TFname(TFname)
pfmfile=all_pfmfile[motif_idx_to_use,]$file
pfm_title=all_pfmfile[motif_idx_to_use,]$pfmName
motif_map=fjComm::motif_plot_moods_singlePWM_hit(readsfile = file,pfmfile = pfmfile,dup_rm = dup_rm,combine_2strand_hits = T)
gg_save_pdf(motif_map,width = 4,height = 1.5)
print(motif_map)

pdf(paste0(TFname,"_motif.pdf"),height = 3)
fjComm::plotMotif_pfmFile(pfmfile,ic.scale = T,title = pfm_title)
dev.off()




motif_idx_to_use=1
TFname= "TFE3"
file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-TFctrl147c4-N4-TFE3IIIc4-AAATTTTG_S1852_L007_R1_001.peared_trimmed.fq.gz"
dup_rm= T

library(fjComm)
setwd(get_scriptpath())


all_pfmfile= motif_getpfm_from_TFname(TFname)
pfmfile=all_pfmfile[motif_idx_to_use,]$file
pfm_title=all_pfmfile[motif_idx_to_use,]$pfmName
motif_map=fjComm::motif_plot_moods_singlePWM_hit(readsfile = file,pfmfile = pfmfile,dup_rm = dup_rm,combine_2strand_hits = T)
gg_save_pdf(motif_map,width = 4,height = 1.5,filename = "motif_map_TFctrl")
print(motif_map)

pdf(paste0(TFname,"_motif_TFctrl.pdf"),height = 3)
fjComm::plotMotif_pfmFile(pfmfile,ic.scale = T,title = pfm_title)
dev.off()
