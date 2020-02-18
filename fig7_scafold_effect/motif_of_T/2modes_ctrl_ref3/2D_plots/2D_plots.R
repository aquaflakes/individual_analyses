fjComm::clear_()

# lig147=nusel_get_147() %>% dplyr::filter(use==1& tag=="Flag" &order==1&batch=="FJ45")

# well="N16"

calcRatio <- function(well)
{
# well="N24"
  boundfile= paste0( "~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/",well,"_Trulig147v1IIIFJ4-5-CAP2r-147c5Nu-*_foldchn_maxBias.csv") %>% Sys.glob()
  # unboundfile= paste0( "~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/",well,"_Trulig147v1IIIFJ4-5-CAP2r-147c5Free*_foldchn_maxBias.csv") %>% Sys.glob()
  
  bound=boundfile %>% read_csv()
  TFname=extract_TF(boundfile,sufix_rm = T)
  plot=gg_heat2D_MI(bound,grad_colors = gg_steelblue_red)+theme(legend.position = "none",plot.title = element_text(size = 9,face = "bold"))+ggtitle(TFname)
  
  return(plot)
}


# wells=lig147$pos %>% sample(20)
wells=qw("M24 O23 P3 O21 G19 I16 P14 D22 P19 J2 G18 A16 B2 C19 C4 M7 A19 I8 B21 O6")
resultCtrl=lapply(wells,calcRatio); resultCtrl$cols=5
pdf("ctrlTF_2d_EMI.pdf",width = 20/2.54,height = 17.3/2.54)
do.call(fjComm::gg_multiplot, resultCtrl)
dev.off()

# resultT=lapply(c("N16","N24","P8","P12"),calcRatio); resultT$cols=7
# do.call(fjComm::gg_multiplot, resultT)


