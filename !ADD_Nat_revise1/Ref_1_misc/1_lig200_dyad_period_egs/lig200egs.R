fjComm::clear_()

lig200=fjComm::nusel_get_200() %>% dplyr::filter(use==1&order==1&(TF=="RFX5"|TF=="SHOX"))
lig147=fjComm::nusel_get_147() %>% dplyr::filter(use==1&order==1&(TF=="RFX5"|TF=="SHOX"))

TFs=c(lig200$TF,lig147$TF)
lig=qw("Lig200 Lig200 Lig147 Lig147")
EMIs=c(lig200$c4_EMI_file,lig147$c4_EMI_file)
plotdf=vector("list",4)

for (i in 1:4)
{
  EMIcurr=EMIs[i] %>% read_csv() %>% dplyr::filter(pos2-pos1==3) %>% mutate(TF=TFs[i],lig=lig[i])
  plotdf[[i]]=EMIcurr
}

plotdf=do.call(rbind,plotdf)

p=ggplot(plotdf)+geom_line(aes(pos1,topMIsum))+facet_grid(TF~lig,scales = "free")+scale_x_continuous(breaks = c(25,75,125)) +scale_y_continuous(breaks = c(0.01,0.03))+
  xlab("(bp)")+ylab("")+theme(strip.text = element_text(size = 8, family = "Helvetica"))
gg_save_pdf(p,8,6,filename = "147_200_corr")

