fjComm::clear_()

well="N16"

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

  return(c(mode2_ratio=mode2_ratio, mode1_ratio=mode1_ratio))
}

result=lapply(c("N16","N24","P8","P12"),calcRatio) %>% {do.call(rbind,.)} %>% as_data_frame()
  t_test=with(result,t.test(mode2_ratio,mode1_ratio,paired = T))
pval=t_test$p.value

result=as.matrix(result) %>% melt %>% set_colnames(c("var1","mode","log2_b_ub") ) %>% mutate(log2_b_ub= log2(log2_b_ub))
# ggplot(data = result,aes(mode,log2_b_ub))+ geom_point()

mode2_enrich_in_bound= ggplot(data = result,aes(mode,log2_b_ub))+ geom_point(size=0.8,alpha=0.99)+ stat_summary(color="red",alpha=0.5) + gg_axis_x_labels(c("Mode 2", "Background")) +
  ylab(expression("Log"[2]*"(E-MI"[b]*"/E-MI"[ub]*")"))+
  geom_path(data=data.frame(x=c(1,1,2,2),y=c(0.27,0.28,0.28,0.27)),aes(x,y)) + annotate("text",x=1.5,y=0.33,size=9/3,label= paste0("p=",signif(pval,2)))+
  gg_theme_Publication()+scale_y_continuous(breaks = gg_breaks(result$log2_b_ub))+ theme(axis.title.x = element_blank(),axis.text.x = element_text(size=9))+ gg_axis_x_label_angle()
gg_save_all(mode2_enrich_in_bound,width = 4,height = 4.5)
