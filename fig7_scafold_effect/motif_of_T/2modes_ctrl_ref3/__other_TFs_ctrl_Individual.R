fjComm::clear_()

# lig147=nusel_get_147() %>% dplyr::filter(use==1& tag=="Flag" &order==1&batch=="FJ45")

# well="N16"

TBX2Nu="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/*_Trulig147v1IIIFJ4-5-CAP2r-147c5Nu-*TBX2[(III)|(ooo)]*_foldchn_maxBias.csv" %>% Sys.glob()
TBX2Nu[4] %>% read_csv() %>% gg_heat2D_MI()

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
  t_testCtrl=with(resultCtrl,t.test(mode2_ratio,mode1_ratio,paired = T))
  # resultCtrl %<>% as.matrix() %>% melt %>% set_colnames(c("var1","mode","log2_b_ub") ) %>% mutate(log2_b_ub= log2(log2_b_ub)) %>% mutate(TF="Ctrls")

resultT=lapply(c("N16","N24","P8","P12"),calcRatio) %>% {do.call(rbind,.)} %>% as_data_frame() %>% mutate_at(vars(ends_with("ratio")), funs(as.numeric))%>% mutate(lable="T")
  t_testT=with(resultT,t.test(mode2_ratio,mode1_ratio,paired = T))
  # resultT %<>% as.matrix() %>% melt %>% set_colnames(c("var1","mode","log2_b_ub") ) %>% mutate(log2_b_ub= log2(log2_b_ub)) %>% mutate(TF="T")

  pval=t_testT$p.value
  pvalCtrl=t_testCtrl$p.value

plotdf=rbind(resultCtrl,resultT) %>% rename(label="lable")


#%>% mutate(label=paste0(mode,"_",TF)) 
p=ggplot(data = plotdf,aes(factor(label,levels = c("mode2_ratio_Ctrls","mode1_ratio_Ctrls","mode2_ratio_T","mode1_ratio_T")),log2_b_ub))+ geom_point(size=0.8,alpha=0.99)
p=p+ stat_summary(color="red",alpha=0.5) + gg_axis_x_labels(c("Mode 2 (Rnd)", "Bk (Rnd)", "Mode 2 (T)", "Bk (T)")) +
  ylab(expression("Log"[2]*"(E-MI"[b]*"/E-MI"[ub]*")"))+
  geom_path(data=data.frame(x=c(1,1,2,2),y=c(0.27,0.28,0.28,0.27)),aes(x,y)) + annotate("text",x=1.5,y=0.33,size=9/3,label= paste0("p=",signif(pval,2)))+
  gg_theme_Publication()+scale_y_continuous()+ theme(axis.title.x = element_blank(),axis.text.x = element_text(size=9))+ gg_axis_x_label_angle()+
  geom_path(data=data.frame(x=c(3,3,4,4),y=c(0.27,0.28,0.28,0.27)),aes(x,y)) + annotate("text",x=3.5,y=0.33,size=9/3,label= paste0("p=",signif(pvalCtrl,2)))

gg_save_pdf(p,width = 7,height = 6,filename = "ctrl_2_modes")
