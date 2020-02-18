fjComm::clear_()

# lig147=nusel_get_147() %>% dplyr::filter(use==1& tag=="Flag" &order==1&batch=="FJ45")

# well="N16"

calcRatio <- function(well)
{
  boundfile= paste0( "~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/",well,"_Trulig147v1IIIFJ4-5-CAP2r-147c5Nu-*_foldchn_maxBias.csv") %>% Sys.glob()
  unboundfile= paste0( "~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/",well,"_Trulig147v1IIIFJ4-5-CAP2r-147c5Free*_foldchn_maxBias.csv") %>% Sys.glob()

  bound_mode2= rio::import(boundfile) %>% dplyr::filter(pos2-pos1 >=77 & pos2-pos1 <=83 ) %>% {mean(.$topMIsum)}
  unbound_mode2= rio::import(unboundfile) %>% dplyr::filter(pos2-pos1 >=77 & pos2-pos1 <=83 ) %>% {mean(.$topMIsum)}

  bound_mode1= rio::import(boundfile) %>% dplyr::filter(pos2-pos1 >=50 & pos2-pos1 <=70 ) %>% {mean(.$topMIsum)}
  unbound_mode1= rio::import(unboundfile) %>% dplyr::filter(pos2-pos1 >=50 & pos2-pos1 <=70 ) %>% {mean(.$topMIsum)}

  mode2_ratio= (bound_mode2 / unbound_mode2)
  mode1_ratio= (bound_mode1 / unbound_mode1)

  return(c(mode2_ratio=mode2_ratio, mode1_ratio=mode1_ratio,TF=extract_TF(boundfile,sufix_rm = T)))
}


# wells=lig147$pos %>% sample(20)
wells=qw("M24 O23 P3 O21 G19 I16 P14 D22 P19 J2 G18 A16 B2 C19 C4 M7 A19 I8 B21 O6")
resultCtrl=lapply(wells,calcRatio)  %>% {do.call(rbind,.)} %>% as_data_frame() %>% mutate_at(vars(ends_with("ratio")), funs(as.numeric)) %>% mutate(lable="ctrlTFs")

resultT=lapply(c("N16","N24","P8","P12"),calcRatio) %>% {do.call(rbind,.)} %>% as_data_frame() %>% mutate_at(vars(ends_with("ratio")), funs(as.numeric))%>% mutate(lable="T")

TBX2Nu="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/*_Trulig147v1IIIFJ4-5-CAP2r-147c5Nu-*TBX2[(III)|(ooo)]*_foldchn_maxBias.csv" %>% Sys.glob() %>% .[1:3]
TBXwells=TBX2Nu %>% lapply(extract_well) %>% unlist()
resultTBX=lapply(TBXwells,calcRatio) %>% {do.call(rbind,.)} %>% as_data_frame() %>% mutate_at(vars(ends_with("ratio")), funs(as.numeric))%>% mutate(lable="TBX2")


plotdf=rbind(resultCtrl,resultT,resultTBX) %>% mutate(diff=mode2_ratio-mode1_ratio)
bk_dist=plotdf$diff[1:20]
  alltests=numeric(27)
  for (i in 1:27)
  {
    t_test=t.test(bk_dist,mu=plotdf$diff[i],alternative = "less")
    alltests[i]=t_test$p.value[1]
    
  }
conf_interval=t_test$conf.int %>% paste0(collapse=" ~ ")
plotdf$p_value=alltests %>% as.numeric() %>% prettyNum(digits=2)

plotdf %<>% mutate(mode2_ratio=round(mode2_ratio,2),mode1_ratio=round(mode1_ratio,2),diff=round(diff,4))
plotdf1=plotdf %>% select(c(3,1,2,6)) %>% set_colnames(c("TF","E-MI ratio (Gs)","E-MI ratio (bk)","p_value"))
plotdf1 %>% writeClipboard()
