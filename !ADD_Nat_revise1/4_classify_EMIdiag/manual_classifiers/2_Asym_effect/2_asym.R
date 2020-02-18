fjComm::clear_()
dataDir="~/Nut_zhuData/Analysis/Nusel/FJ4.4/DNA_rotation_symbreak/result/201802_allTFs/"

# for use
library(plotly,ggplot2)

df=read_tsv(paste0(dataDir,"/output_200cmpTFctrl_mac.txt")) %>% mutate(family=lapply(TF,TF_classify) %>% unlist())
# c4_tTest_mean_diff is mean(abs(ddG))
p=ggplot()+geom_point(data = df,aes(`log10(t_energy$p.value)`,c4_tTest_mean_diff,  comment=TF,well=pos,size=I(max_foldchn_sig*0.5) ), alpha=ifelse(df$family=="ETS",1,0.3),color=ifelse(df$family=="ETS","red","black"))
p=p+ theme_bw()+
  ylab("Orientation asymmetry")+guides(size=F,color=F)+
  scale_y_continuous(expand = c(0.002,0.002))+
  # geom_text(data = df,aes(`log10(t_energy$p.value)`,t_energy_mean_diff+0.02, label=TF),alpha=1/3)+
  gg_theme_Publication()+scale_y_continuous(breaks = gg_breaks(df$c4_tTest_mean_diff,fromMin = F),limits = c(0,df$c4_tTest_mean_diff %>% max),expand = c(0,0.004))+
  scale_x_continuous(breaks = c(1,-10,-20),trans = "reverse")
p1=ggplotly(p)
print(p1)
p=p+xlab(expression("log"*""["10"]*"(p value)"))
# gg_save_pdf(p,width=6.2,height=6,path =dataDir, filename = "c4_asym_corr_ctrl_p")

output= df %>% select(TF,c4_tTest_mean_diff,`log10(t_energy$p.value)`) %>% mutate(TF=TF %>% str_replace("ooo.*","")) %>% set_colnames(c("TF","Orientational_asym","p_value_log10")) %>%
  mutate(Asym=ifelse(Orientational_asym>=0.1 & p_value_log10<=-5.7,"Yes","No"))
saveRDS(output,"2_asym")
