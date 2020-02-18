fjComm::clear_()

cyc5b="~/Nut_zhuData/seqFiles2/FJ6.1Methyl_PE/allreads/FJ6.1_OriMe1_PE/dual_trim/2_Trim_adaptor_kmer/Trulig147v1IIIFJ6-1-OriMe1-B02-1xCAGOct32ngIIIc5-Bq1-GATCCACG_S14_L003_R1_001.peared_trimmed.fq.gz"
cyc5ub="~/Nut_zhuData/seqFiles2/FJ6.1Methyl_PE/allreads/FJ6.1_OriMe1_PE/dual_trim/2_Trim_adaptor_kmer/Trulig147v1IIIFJ6-1-OriMe1-A06-MeStdOct64ngIIIc5-Bq1-AAAAAACG_S6_L003_R1_001.peared_trimmed.fq.gz"
# cyc5ub="~/Nut_zhuData/seqFiles2/FJ6.1Methyl_PE/allreads/FJ6.1_OriMe1_PE/dual_trim/2_Trim_adaptor_kmer/Trulig147v1IIIFJ6-1-OriMe1-C05-StdUbOct64ngIIIc5-Bq1-CTGCCATG_S29_L003_R1_001.peared_trimmed.fq.gz"


kcounts=kmerCntBit_file(file = cyc5b,k = 9L,diffLen = F,collapse = T,asDf = T,all_possible_k = T,pseudo = 5,rmdup = F,bkfile = cyc5ub)
kcounts %<>% mutate(C=str_count(kmer,"C"),G=str_count(kmer,"G"),CG=str_count(kmer,"CG")) %>% .[sample(1:nrow(kcounts)),]

# CG
p=ggplot(kcounts)+geom_point(aes(counts,bkcounts,color=CG),alpha=1/3,size=0.1)+scale_color_gradientn(colours = c("grey","red"))+gg_theme_Publication()+xlab("Bound")+ylab("Unbound")+scale_x_continuous(breaks = c(0,2000))+scale_y_continuous(breaks = c(0,1000))
gg_save_png(p,6,4,filename = "CG_xyplot")

p=ggplot(kcounts)+geom_boxplot(aes(CG %>% factor(),counts/bkcounts,color=CG),outlier.size = 0.1)+scale_color_gradientn(colours = c("grey","red"),guide = F)+gg_theme_Publication()+xlab("CG content")+scale_y_continuous(breaks = c(0,8))+theme(axis.title.y = element_blank())#ylab("Count ratio (bound/unbound)")
gg_save_pdf(p,4.5,4,filename = "CG_box_ratio")
# p=ggplot(kcounts)+geom_boxplot(aes(CG %>% factor(),bkcounts,color=CG),outlier.size = 0.1)+scale_color_gradientn(colours = c("grey","red"),guide = F)+gg_theme_Publication()+ylab("Bound")+xlab("CG number")+scale_y_continuous(breaks = c(0,1000))+theme(axis.title.y = element_blank())
# gg_save_pdf(p,4.5,4,filename = "CG_box_unbound")
#

# C
p=ggplot(kcounts)+geom_point(aes(counts,bkcounts,color=C),alpha=1/3,size=0.1)+scale_color_gradientn(colours = c("grey","red"))+gg_theme_Publication()+xlab("Bound")+ylab("Unbound")+scale_x_continuous(breaks = c(0,2000))+scale_y_continuous(breaks = c(0,1000))
gg_save_png(p,6,4,filename = "C_xyplot")
p=ggplot(kcounts)+geom_boxplot(aes(C %>% factor(),counts/bkcounts,color=C),outlier.size = 0.1)+scale_color_gradientn(colours = c("grey","red"),guide = F)+gg_theme_Publication()+xlab("C content")+scale_y_continuous(breaks = c(0,8))+theme(axis.title.y = element_blank())#ylab("Count ratio (bound/unbound)")
gg_save_pdf(p,4.5,4,filename = "C_box_ratio")


# p=ggplot(kcounts)+geom_boxplot(aes(C %>% factor(),counts,color=C),outlier.size = 0.1)+scale_color_gradientn(colours = c("grey","red"),guide = F)+gg_theme_Publication()+xlab("C number")+scale_y_continuous(breaks = c(0,2000))+theme(axis.title.y = element_blank())
# gg_save_pdf(p,4.5,4,filename = "C_box_bound")
# p=ggplot(kcounts)+geom_boxplot(aes(C %>% factor(),bkcounts,color=C),outlier.size = 0.1)+scale_color_gradientn(colours = c("grey","red"),guide = F)+gg_theme_Publication()+ylab("Bound")+xlab("C number")+scale_y_continuous(breaks = c(0,1000))+theme(axis.title.y = element_blank())
# gg_save_pdf(p,4.5,4,filename = "C_box_unbound")


# G
p=ggplot(kcounts)+geom_point(aes(counts,bkcounts,color=G),alpha=1/3,size=0.1)+scale_color_gradientn(colours = c("grey","red"))+gg_theme_Publication()+xlab("Bound")+ylab("Unbound")+scale_x_continuous(breaks = c(0,2000))+scale_y_continuous(breaks = c(0,1000))
gg_save_png(p,6,4,filename = "G_xyplot")
p=ggplot(kcounts)+geom_boxplot(aes(G %>% factor(),counts/bkcounts,color=G),outlier.size = 0.1)+scale_color_gradientn(colours = c("grey","red"),guide = F)+gg_theme_Publication()+xlab("G content")+scale_y_continuous(breaks = c(0,8))+theme(axis.title.y = element_blank())#ylab("Count ratio (bound/unbound)")
gg_save_pdf(p,4.5,4,filename = "G_box_ratio")


# p=ggplot(kcounts)+geom_boxplot(aes(G %>% factor(),counts,color=G),outlier.size = 0.1)+scale_color_gradientn(colours = c("grey","red"),guide = F)+gg_theme_Publication()+xlab("G number")+scale_y_continuous(breaks = c(0,2000))+theme(axis.title.y = element_blank())
# gg_save_pdf(p,4.5,4,filename = "G_box_bound")
# p=ggplot(kcounts)+geom_boxplot(aes(G %>% factor(),bkcounts,color=G),outlier.size = 0.1)+scale_color_gradientn(colours = c("grey","red"),guide = F)+gg_theme_Publication()+ylab("Bound")+xlab("G number")+scale_y_continuous(breaks = c(0,1000))+theme(axis.title.y = element_blank())
# gg_save_pdf(p,4.5,4,filename = "G_box_unbound")
