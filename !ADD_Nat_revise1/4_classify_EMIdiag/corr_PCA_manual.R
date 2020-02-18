fjComm::clear_()
manual=readRDS("manual_classifiers/all_classes")
automatic= readRDS("PCA/PCA_EMI_cor_dim")

automatic=left_join(automatic,manual,by="TF")

cor(automatic$cor_dim1,automatic$EMI.penetration..lig147.)
cor(automatic$cor_dim5,automatic$Dyad.EMI.intensity)



my_grob=grobTree(textGrob("r = " %>% paste0(cor(automatic$cor_dim5,automatic$Dyad.EMI.intensity) %>% prettyNum(digits=3)),x=0.1, y=0.95,hjust=0, gp=gpar(fontsize=9,fontface="italic")))
Dyad_with_dim5=ggplot(automatic)+geom_point(aes(cor_dim5,`Dyad.EMI.intensity`,color=`Dyad.preference`))+gg_theme_Publication()+xlab("Corr(EMI_diag, Dim_5)")+ylab("Dyad EMI intensity") +
  scale_x_continuous(breaks = c(-0.5,0,0.5))+scale_y_continuous(breaks = c(0,0.6))+scale_color_manual(values = c("black","red"),guide=F)+annotation_custom(my_grob)
gg_save_pdf(Dyad_with_dim5,5.5,5,path ="plots/")


my_grob=grobTree(textGrob("r = " %>% paste0(cor(automatic$cor_dim1,automatic$EMI.penetration..lig147.) %>% prettyNum(digits=3)),x=0.1, y=0.95,hjust=0, gp=gpar(fontsize=9,fontface="italic")))
End_with_dim1= ggplot(automatic)+geom_point(aes(cor_dim1,EMI.penetration..lig147.,color=End.preference))+gg_theme_Publication()+xlab("Corr(EMI_diag, Dim_1)")+ylab("EMI penetration (bp)") +
  scale_x_continuous(breaks = c(-0.5,0,0.5))+scale_y_continuous(breaks = c(10,40))+scale_color_manual(values = c("red","black"),guide=F)+annotation_custom(my_grob)
gg_save_pdf(End_with_dim1,5.5,5,path ="plots/")


# my_grob=grobTree(textGrob("r = " %>% paste0(cor(automatic$cor_dim2,automatic$EMI.penetration..lig147.) %>% prettyNum(digits=3)),x=0.1, y=0.95,hjust=0, gp=gpar(fontsize=9,fontface="italic")))
# End_with_dim2= ggplot(automatic)+geom_point(aes(cor_dim2,EMI.penetration..lig147.,color=End.preference))+gg_theme_Publication()+xlab("Corr(EMI_diag, Dim_1)")+ylab("EMI penetration (bp)") +
#   scale_x_continuous()+scale_y_continuous()+scale_color_manual(values = c("red","black"),guide=F)+annotation_custom(my_grob)
# gg_save_pdf(End_with_dim2,5.5,5,path ="plots/")

my_grob=grobTree(textGrob("r = " %>% paste0(cor(automatic$cor_dim3,automatic$Periodicity.strength) %>% prettyNum(digits=3)),x=0.6, y=0.95,hjust=0, gp=gpar(fontsize=9,fontface="italic")))
Period_with_dim3= ggplot(automatic)+geom_point(aes(cor_dim3,Periodicity.strength,color=Periodic.preference))+gg_theme_Publication()+xlab("Corr(EMI_diag, Dim_3)")+ylab("Periodicity strength\n(FFT-AUC)") +
  scale_x_continuous(breaks = c(-0.5,0,0.5))+scale_y_continuous(breaks = c(0,0.005))+scale_color_manual(values = c("black","red"),guide=F)+annotation_custom(my_grob)
gg_save_pdf(Period_with_dim3,6.5,5,path ="plots/")

my_grob=grobTree(textGrob("r = " %>% paste0(cor(automatic$cor_dim4,automatic$Periodicity.strength) %>% prettyNum(digits=3)),x=0.6, y=0.95,hjust=0, gp=gpar(fontsize=9,fontface="italic")))
Period_with_dim4= ggplot(automatic)+geom_point(aes(cor_dim4,Periodicity.strength,color=Periodic.preference))+gg_theme_Publication()+xlab("Corr(EMI_diag, Dim_4)")+ylab("Periodicity strength\n(FFT-AUC)") +
  scale_x_continuous(breaks = c(-0.5,0,0.5))+scale_y_continuous(breaks = c(0,0.005))+scale_color_manual(values = c("black","red"),guide=F)+annotation_custom(my_grob)
gg_save_pdf(Period_with_dim4,6.5,5,path ="plots/")





