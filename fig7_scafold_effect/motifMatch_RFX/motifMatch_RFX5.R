fjComm::clear_()
# match reads with pfm file, and plot
motif_plot_moods_singlePWM_hit <- function(readsfile, pfmfile, dup_rm=TRUE, p_moods="0.0001", colorGrad= c("#ceffce","red"), strand_filter=c("+","-"), combine_2strand_hits=FALSE)
{
  selex= SELEXFile(readsfile)
  selex$getSeq(dup_rm = dup_rm)
  selex$moodsMap(pfmfile, p=p_moods )
  selex$moodsResult= selex$moodsResult %>% group_by(pos,strand) %>% summarise(count=n()) %>% ungroup() # bin counts with the same strand and pos

  if (length(strand_filter)==2 & (!combine_2strand_hits))
  {
    selex$moodsResult= selex$moodsResult %>% group_by(strand) %>% mutate(ncount= prop.table(count)) %>% ungroup()
    plot=ggplot(selex$moodsResult)+ geom_tile(aes(pos+1,strand,fill=ncount))+ scale_fill_gradientn(colors =colorGrad)+ guides(fill=F)  +ylab("Strand")+
      gg_axis_x_noExp(breaks=gg_breaks(selex$moodsResult$pos+1))+ scale_y_discrete(expand = c(0,0)) + gg_theme_Publication() + theme(axis.line= element_blank(),axis.ticks.y = element_blank())
  }else{
    if (combine_2strand_hits) strand_filter=c("+","-")
    if (length(strand_filter)==1) {selex$moodsResult= selex$moodsResult %>% dplyr::filter(strand %in% strand_filter) %>% mutate(ncount= prop.table(count))}
    else if (length(strand_filter)==2){selex$moodsResult= selex$moodsResult %>% group_by(pos) %>% summarise(count=sum(count)) %>% ungroup() %>% mutate(ncount= prop.table(count)) }
    else {stop("strand filter length is not correct !!")}
    plot=ggplot(selex$moodsResult)+ geom_tile(aes(pos+1,1,fill=ncount))+ scale_fill_gradientn(colors =colorGrad)+ guides(fill=F) +
      gg_axis_x_noExp(breaks=gg_breaks(selex$moodsResult$pos+1,end_plus_middle = T))+gg_axis_y_noExp()+ gg_theme_Publication_diag
  }
  return(plot+xlab("(bp)")+scale_x_continuous(expand = c(0, 0),breaks =  gg_breaks(selex$moodsResult$pos+1,end_plus_middle = T)))
}

MIfile="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/analysis/2d_IC_related_pub/data/A6_Trulig147v1IIIFJ4-4-CAP2-A06-TF99-RFX5IIIc4_S390_L002_R1_001.peared_trimmed.fq_u.gz_3mer_foldchn_maxBias.csv"
df=read_csv(MIfile)
RFX5_MI_diag= gg_heat2D_diag(df,grad_colors = gg_steelblue_red)+gg_theme_bordered_diag
gg_save_diag(RFX5_MI_diag)

MIfile="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/analysis/2d_IC_related_pub/data/A6_Trulig147v1IIIFJ4-4-TFctrl-A06-TF99-RFX5IIIc3_S6_L004_R1_001.peared_trimmed.fq_u.gz_3mer_foldchn_maxBias.csv"
df=read_csv(MIfile)
RFX5_MI_diag_TFctrl= gg_heat2D_diag(df,grad_colors = gg_steelblue_red)+gg_theme_bordered_diag
gg_save_diag(RFX5_MI_diag_TFctrl)


pfmfile= "RFX5_motif_FJ4.4.txt"
file="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_147/Trulig147v1IIIFJ4-4-CAP2-A06-TF99-RFX5IIIc4_S390_L002_R1_001.peared_trimmed.fq_u.gz"
RFX5_map_both_str= motif_plot_moods_singlePWM_hit(file,pfmfile,dup_rm = F,p_moods = "0.0001",colorGrad = c("#ceffce","red"),strand_filter = c("+","-"),combine_2strand_hits = T)
print(RFX5_map_both_str)
gg_save_diag(RFX5_map_both_str)

pdf("motif.pdf",height = 3)
plotMotif_pfmFile(pfmfile,ic.scale = T)
dev.off()

file="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_TFctrl/Trulig147v1IIIFJ4-4-TFctrl-A06-TF99-RFX5IIIc3_S6_L004_R1_001.peared_trimmed.fq_u.gz"
RFX5_TFctrl_map_both_str= motif_plot_moods_singlePWM_hit(file,pfmfile,dup_rm = F,p_moods = "0.0001",colorGrad = c("#ceffce","red"),strand_filter = c("+","-"),combine_2strand_hits = T)
print(RFX5_TFctrl_map_both_str)
gg_save_diag(RFX5_TFctrl_map_both_str)

