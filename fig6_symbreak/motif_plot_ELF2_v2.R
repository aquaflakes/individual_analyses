fjComm::clear_()


motif_idx_to_use=1
TFname= "ELF2"
file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig200v1IIIFJ4-5-CAP2r-200c4-M24-ELF2IIIc4-AGCGCCTC_S1464_L006_R1_001.peared_trimmed.fq.gz"
dup_rm= T


all_pfmfile= motif_getpfm_from_TFname(TFname)
pfmfile=all_pfmfile[motif_idx_to_use,]$file
pfm_title=all_pfmfile[motif_idx_to_use,]$pfmName
motif_map=fjComm::motif_plot_moods_singlePWM_hit(readsfile = file,pfmfile = pfmfile,dup_rm = T,combine_2strand_hits = F)#,colorGrad = c("white","red"))
gg_save_pdf(motif_map,width = 5.53,height = 2.33)
print(motif_map)

pdf(paste0(TFname,"_motif.pdf"),height = 3)
fjComm::plotMotif_pfmFile(pfmfile,ic.scale = T,title = pfm_title)
dev.off()



file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig200v1IIIFJ4-5-CAP2r-TFctrl200c4-M24-ELF2IIIc4-AGCGCCTC_S2232_L008_R1_001.peared_trimmed.fq.gz"
# file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-TFctrl147c4-M24-ELF2IIIc4-AGCGCCTC_S1848_L007_R1_001.peared_trimmed.fq.gz"
# file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig200v1IIIFJ4-5-CAP2r-TFctrl200c1-M24-ELF2IIIc2-AGCGCCTC_S2232_L007_R1_001.peared_trimmed.fq.gz"


# all_pfmfile= motif_getpfm_from_TFname(TFname)
pfmfile=all_pfmfile[motif_idx_to_use,]$file
pfm_title=all_pfmfile[motif_idx_to_use,]$pfmName
motif_map=fjComm::motif_plot_moods_singlePWM_hit(readsfile = file,pfmfile = pfmfile,dup_rm = T,combine_2strand_hits = F)
gg_save_pdf(motif_map,width = 5.53,height = 2.33,filename = "ELF2_TFctrl_motif_map")
print(motif_map)

pdf(paste0("ELF2_TFctrl_","_motif.pdf"),height = 3)
fjComm::plotMotif_pfmFile(pfmfile,ic.scale = T,title = pfm_title)
dev.off()
