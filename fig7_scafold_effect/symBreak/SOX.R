fjComm::clear_()

motif_idx_to_use=1
TFname= c("SOX10","SOX11","SOX12","SOX18","SOX7")
files=c("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4-L13-SOX10oooNIIIc4-GGCCAGCT_S277_L003_R1_001.peared_trimmed.fq.gz",
        "~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_147/Trulig147v1IIIFJ4-4-CAP2-C05-TF15-SOX11IIIc4_S437_L002_R1_001.peared_trimmed.fq_u.gz",
        "~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_147/Trulig147v1IIIFJ4-4-CAP2-C09-TF17-SOX12IIIc4_S441_L002_R1_001.peared_trimmed.fq_u.gz",
        "~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4-D8-SOX18IIIc4-AGGTGTAA_S80_L003_R1_001.peared_trimmed.fq.gz",
        "~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4-P5-SOX7IIIc4-AAAATTGC_S365_L003_R1_001.peared_trimmed.fq.gz"
        )
dup_rm= T


for (i in 1:5)
{
  all_pfmfile= motif_getpfm_from_TFname("SOX11")
  pfmfile=all_pfmfile[motif_idx_to_use,]$file
    pfmfile=("SOX_motif.txt")
  pfm_title=all_pfmfile[motif_idx_to_use,]$pfmName
  motif_map=fjComm::motif_plot_moods_singlePWM_hit(readsfile = files[i],pfmfile = pfmfile,dup_rm = dup_rm,combine_2strand_hits = F) + theme(axis.title.y = element_blank())
  gg_save_pdf(motif_map,width = 4,height = 1.5,filename = paste0(TFname[i],"_motif_map"))
  print(motif_map)
}  


# # pdf(paste0(TFname,"_motif.pdf"),height = 3)
# # fjComm::plotMotif_pfmFile(pfmfile,ic.scale = T,title = pfm_title)
# # dev.off()
# 
# 
# 
# 
# motif_idx_to_use=1
# TFname= "SOX11"
# file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-TFctrl147c4-L13-SOX10oooNIIIc4-GGCCAGCT_S1813_L007_R1_001.peared_trimmed.fq.gz"
# dup_rm= T
# 
# library(fjComm)
# setwd(get_scriptpath())
# 
# 
# all_pfmfile= motif_getpfm_from_TFname(TFname)
# pfmfile=all_pfmfile[motif_idx_to_use,]$file
# pfm_title=all_pfmfile[motif_idx_to_use,]$pfmName
# motif_map=fjComm::motif_plot_moods_singlePWM_hit(readsfile = file,pfmfile = pfmfile,dup_rm = dup_rm,combine_2strand_hits = F)
# gg_save_pdf(motif_map,width = 4,height = 1.5,filename = "motif_map_TFctrl")
# print(motif_map)
# 
# # pdf(paste0(TFname,"_motif_TFctrl.pdf"),height = 3)
# # fjComm::plotMotif_pfmFile(pfmfile,ic.scale = T,title = pfm_title)
# # dev.off()
