fjComm::clear_()

motif_idx_to_use=1
TFname= "CREB1"
file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig200v1IIIFJ4-5-CAP2r-200c4-J21-CREB1IIIc4-AGAATGAC_S1389_L006_R1_001.peared_trimmed.fq.gz"
# file="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_200/Trulig200v1IIIFJ4-4-CAP2-J21-TF251-CREB1IIIc4_S621_L002_R1_001.peared_trimmed.fq_u.gz"

all_pfmfile= motif_getpfm_from_TFname(TFname)
pfmfile=all_pfmfile[motif_idx_to_use,]$file
pfm_title=all_pfmfile[motif_idx_to_use,]$pfmName
motif_map=fjComm::motif_plot_moods_singlePWM_hit(readsfile = file,pfmfile = pfmfile,dup_rm = T,combine_2strand_hits = F)
gg_save_pdf(motif_map,width = 4,height = 2)
print(motif_map)

pdf(paste0(TFname,"_motif.pdf"),height = 3)
fjComm::plotMotif_pfmFile(pfmfile,ic.scale = T,title = pfm_title)
dev.off()



# file="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_TFctrl/Trulig147v1IIIFJ4-4-TFctrl-J21-TF251-CREB1IIIc3_S237_L004_R1_001.peared_trimmed.fq_u.gz"
file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig200v1IIIFJ4-5-CAP2r-TFctrl200c4-J21-CREB1IIIc4-AGAATGAC_S2157_L008_R1_001.peared_trimmed.fq.gz"
# file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-TFctrl147c4-J21-CREB1IIIc4-AGAATGAC_S1773_L007_R1_001.peared_trimmed.fq.gz"


all_pfmfile= motif_getpfm_from_TFname(TFname)
pfmfile=all_pfmfile[motif_idx_to_use,]$file
pfm_title=all_pfmfile[motif_idx_to_use,]$pfmName
motif_map=fjComm::motif_plot_moods_singlePWM_hit(readsfile = file,pfmfile = pfmfile,dup_rm = T,combine_2strand_hits = F)
gg_save_pdf(motif_map,width = 4,height = 2,filename = "CREB1_TFctrl_motif_map")
print(motif_map)

pdf(paste0("CREB1_TFctrl_","_motif.pdf"),height = 3)
fjComm::plotMotif_pfmFile(pfmfile,ic.scale = T,title = pfm_title)
dev.off()
