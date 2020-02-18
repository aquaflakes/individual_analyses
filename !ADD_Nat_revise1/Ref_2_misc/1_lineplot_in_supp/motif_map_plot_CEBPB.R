fjComm::clear_()

# lig147=nusel_get_147() %>% dplyr::filter(use==1&order==1)
# MIfile=lig147 %>% dplyr::filter(TF=="CEBPB") %>% c4_EMI_file

MIfile="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/O16_Trulig147v1IIIFJ4-5-CAP2r-147c4-O16-CEBPBIIIc4-ACGGTAAT_S352_L003_R1_001.peared_trimmed.fq.gz_3mer_foldchn_maxBias.csv"
df1=read_csv(MIfile) %>% mutate(TF="NCAP") %>% dplyr::filter(pos2-pos1==3)
CEBPB_MI_diag= gg_heat2D_diag(df1,grad_colors = gg_steelblue_red)+gg_theme_bordered_diag
# gg_save_diag(CEBPB_MI_diag)

MIfile="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/O16_Trulig147v1IIIFJ4-5-CAP2r-TFctrl147c4-O16-CEBPBIIIc4-ACGGTAAT_S1888_L007_R1_001.peared_trimmed.fq.gz_3mer_foldchn_maxBias.csv"
df2=read_csv(MIfile) %>% mutate(TF="HT") %>% dplyr::filter(pos2-pos1==3)
CEBPB_MI_diag_TFctrl= gg_heat2D_diag(df2,grad_colors = gg_steelblue_red)+gg_theme_bordered_diag
# gg_save_diag(CEBPB_MI_diag_TFctrl)

MIdf=rbind(df1,df2)
MI_plot_CEBPB=ggplot(MIdf)+geom_line(aes(pos1,topMIsum,color=TF))+xlab("(bp)")+ylab("E-MI (bit)")+scale_color_manual(values = c("steelblue","orange"),name="E-MI")+scale_x_continuous(expand = c(0,0),breaks = c(1,96))+scale_y_continuous(breaks = c(0,0.15)) 
gg_save_pdf(MI_plot_CEBPB,7,3)


lig147=nusel_get_147() %>% dplyr::filter(use==1&order==1)
temp=lig147 %>% dplyr::filter(TF=="CEBPB")

motif_idx_to_use=1
TFname= "CEBPB"
all_pfmfile= motif_getpfm_from_TFname(TFname)
pfmfile=all_pfmfile[motif_idx_to_use,]$file

file=temp$c4_file
motif_map=fjComm::motif_plot_moods_singlePWM_hit(readsfile = file,pfmfile = pfmfile,dup_rm = T,combine_2strand_hits = T)
readsNum=file %>% read_csv(col_names=F) %>% rmdup() %>% nrow()
df1=motif_map$data %>% mutate(TF="NCAP",count=count/readsNum)

file=temp$TFctrl_file
motif_map=fjComm::motif_plot_moods_singlePWM_hit(readsfile = file,pfmfile = pfmfile,dup_rm = T,combine_2strand_hits = T)
readsNum=file %>% read_csv(col_names=F) %>% rmdup() %>% nrow()
df2=motif_map$data %>% mutate(TF="HT",count=count/readsNum)
Motifdf=rbind(df1,df2)
Motif_plot_CEBPB=ggplot(Motifdf)+geom_line(aes(pos+1,count,color=TF))+xlab("(bp)")+ylab("Motif Freq.")+scale_color_manual(values = c("steelblue","orange"),name="")+scale_x_continuous(expand = c(0,0),breaks = c(1,92))+scale_y_continuous(limits = c(0,0.03),breaks = c(0,0.03)) 
gg_save_pdf(Motif_plot_CEBPB,7,3)

