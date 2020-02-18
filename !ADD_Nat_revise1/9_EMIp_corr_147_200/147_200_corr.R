fjComm::clear_()

EMIp147="147EMIp" %>% readRDS()
EMIp200="200EMIp" %>% readRDS() %>% dplyr::rename(EMIp200="EMIp147")

plotdf=dplyr::inner_join(EMIp147,EMIp200,by=c("TF","class"))

p=ggplot(plotdf)+geom_point(aes(EMIp147,EMIp200))+xlab("Lig147 E-MI penetration (bp)")+ylab("Lig200 E-MI penetration (bp)")+scale_y_continuous(breaks = c(10,30,50))+gg_theme_Publication()
gg_save_pdf(p,6.7,6.7,filename = "EMI_corr_147_200")

with(plotdf,cor.test(EMIp147,EMIp200))
