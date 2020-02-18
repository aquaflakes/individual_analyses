fjComm::clear_()

# lig147=nusel_get_147() %>% dplyr::filter(use==1& tag=="Flag" &order==1&batch=="FJ45")

# well="N16"

calcRatio <- function(well)
{
# well="N24"
  boundfile= paste0( "~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/",well,"_Trulig147v1IIIFJ4-5-CAP2r-147c5Nu-*_foldchn_maxBias.csv") %>% Sys.glob()
  unboundfile= paste0( "~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/",well,"_Trulig147v1IIIFJ4-5-CAP2r-147c5Free*_foldchn_maxBias.csv") %>% Sys.glob()
  
  bound=boundfile %>% read_csv()
  unbound=unboundfile %>% read_csv()
  
  free_Nu_diff=bound %>% mutate(topMIsum=log2(unbound$topMIsum/bound$topMIsum))

  
  # c5_all_sig_free_Nu_diff_ratio= ggbasic + stat_summary_2d(data = free_Nu_diff,aes_string("pos1","pos2",z=names(free_Nu_diff)[3]),fun = mean,binwidth = 1) + scale_fill_gradientn(name= expression("log"[2]*"(E-MI ratio)"),colors = c("#7dbb78","white","red"), limits=quantile(free_Nu_diff[[3]], c(0.01, 0.99)), oob=scales::squish,values = gg_gen_scale_with_quan(free_Nu_diff[[3]], 0.01, 0.99) )+ gg_bkcolor("white") + xlab("(bp)") + ylab("(bp)")+
  #   gg_theme_Publication()+theme(panel.grid.major = element_line(size = 0.4,color="grey",linetype = "dashed"),panel.grid.minor = element_line(size = 0.4,color="grey",linetype = "dashed"))+
  #   theme(legend.justification=c(1,0), legend.position=c(1,0))+scale_x_continuous(breaks = gg_breaks(free_Nu_diff$pos1,end_plus_middle = T),expand = c(0,0))+scale_y_continuous(breaks = gg_breaks(free_Nu_diff$pos2,end_plus_middle = T),expand = c(0,0))
  # print(c5_all_sig_free_Nu_diff_ratio)
  
  TFcurr= extract_TF(boundfile,sufix_rm = T)
  c5_all_sig_free_Nu_diff_ratio= ggbasic + stat_summary_2d(data = free_Nu_diff,aes_string("pos1","pos2",z=names(free_Nu_diff)[3]),fun = mean,binwidth = 1) + scale_fill_gradientn(name= expression("Log"[2]*"(E-MI ratio)"),colors = c("green","white","red"),guide = F, limits=c(-0.8,0.8), oob=scales::squish,values = c(0,0.5,1),breaks=c(-0.8,0,0.8) )+ gg_bkcolor("white") + xlab("(bp)") + ylab("(bp)")+
    gg_theme_Publication()+theme(panel.grid.major = element_line(size = 0.4,color="grey",linetype = "dashed"),panel.grid.minor = element_line(size = 0.4,color="grey",linetype = "dashed"))+
    theme(legend.justification=c(1,0), legend.position=c(1,0))+scale_x_continuous(breaks = gg_breaks(free_Nu_diff$pos1,end_plus_middle = T),expand = c(0,0))+scale_y_continuous(breaks = gg_breaks(free_Nu_diff$pos2,end_plus_middle = T),expand = c(0,0))+
    ggtitle(TFcurr)
  # print(c5_all_sig_free_Nu_diff_ratio)
  
  return(c5_all_sig_free_Nu_diff_ratio)
}


# wells=lig147$pos %>% sample(20)
wells=qw("M24 O23 P3 O21 G19 I16 P14 D22 P19 J2 G18 A16 B2 C19 C4 M7 A19 I8 B21 O6")
resultCtrl=lapply(wells,calcRatio); resultCtrl$cols=7
do.call(fjComm::gg_multiplot, resultCtrl)

resultT=lapply(c("N16","N24","P8","P12"),calcRatio); resultT$cols=7
do.call(fjComm::gg_multiplot, resultT)


