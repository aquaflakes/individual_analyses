fjComm::clear_()

motif_idx_to_use=1
TFname= "CREB3L4"
# file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4-C4-CREB3L1IIIc4-TTCTTCCC_S52_L003_R1_001.peared_trimmed.fq.gz"
file="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_200/Trulig200v1IIIFJ4-4-CAP2-A24-TF108-CREB3L4IIIc4_S408_L002_R1_001.peared_trimmed.fq_u.gz"
# file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig200v1IIIFJ4-5-CAP2r-200c4-A24-CREB3L4IIIc4-GAAACGTG_S1176_L006_R1_001.peared_trimmed.fq.gz"


all_pfmfile= motif_getpfm_from_TFname(TFname)
pfmfile=all_pfmfile[motif_idx_to_use,]$file
  pfmfile="CREB3L4_motif.txt"
pfm_title=all_pfmfile[motif_idx_to_use,]$pfmName
motif_map=fjComm::motif_plot_moods_singlePWM_hit(readsfile = file,pfmfile = pfmfile,dup_rm = T,combine_2strand_hits = F)
gg_save_pdf(motif_map,width = 4,height = 1.5)
print(motif_map)

pdf(paste0(TFname,"_motif.pdf"),height = 3)
fjComm::plotMotif_pfmFile(pfmfile,ic.scale = T,title = pfm_title)
dev.off()



# file="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_TFctrl/Trulig147v1IIIFJ4-4-TFctrl-A24-TF108-CREB3L4IIIc3_S24_L004_R1_001.peared_trimmed.fq_u.gz"
file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig200v1IIIFJ4-5-CAP2r-TFctrl200c4-A24-CREB3L4IIIc4-GAAACGTG_S1944_L008_R1_001.peared_trimmed.fq.gz"

all_pfmfile= motif_getpfm_from_TFname(TFname)
pfmfile=all_pfmfile[motif_idx_to_use,]$file
  pfmfile="CREB3L4_motif.txt"
pfm_title=all_pfmfile[motif_idx_to_use,]$pfmName
motif_map=fjComm::motif_plot_moods_singlePWM_hit(readsfile = file,pfmfile = pfmfile,dup_rm = T,combine_2strand_hits = F)
gg_save_pdf(motif_map,width = 4,height = 1.5,filename = "CREB3L4_TFctrl_motif_map")
print(motif_map)

pdf(paste0("CREB3L4_TFctrl_","_motif.pdf"),height = 3)
fjComm::plotMotif_pfmFile(pfmfile,ic.scale = T,title = pfm_title)
dev.off()
