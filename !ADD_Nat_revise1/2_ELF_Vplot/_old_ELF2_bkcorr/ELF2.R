fjComm::clear_()

# pallete=RColorBrewer::brewer.pal(n = 11, name = "Spectral") %>% rev()

# ggplot(sigdf)+geom_tile(aes(X11,X4,color=count))

library(MASS)

evalsig<-function()
{
  sigdf=read_tsv("signal.txt",col_names = F) %>% dplyr::select(X11,X4) %>% dplyr::filter(X4<=200)
  sigdf %<>% dplyr::group_by(X11,X4) %>% summarise(sigcount=n()); sigdf
}
sigdf=fjComm::oneOffCalc("oneoff/sigdensity",evalsig) %>% mutate(sigcount=sigcount/mean(sigcount))

evalbk<-function()
{
  bkdf=read_tsv("bk.txt",col_names = F) %>% dplyr::select(X11,X4) %>% dplyr::filter(X4<=200)
  bkdf %<>% dplyr::group_by(X11,X4) %>% summarise(bkcount=n()); bkdf
}
bkdf=fjComm::oneOffCalc("oneoff/bkdensity",evalbk) %>% mutate(bkcount=bkcount/mean(bkcount))
ggplot(bkdf)+geom_tile(aes(X11,X4,fill=bkcount))+scale_x_continuous(limits = c(-200,200),expand = c(0,0),breaks = c(-200,0,200))+scale_fill_gradientn(colours =c("white","#ead100","red","red"))
ggplot(sigdf)+geom_tile(aes(X11,X4,fill=sigcount))+scale_x_continuous(limits = c(-200,200),expand = c(0,0),breaks = c(-200,0,200))+scale_fill_gradientn(colours =c("white","#ead100","red","red"))

bkcorrdf=dplyr::inner_join(sigdf,bkdf,by=c("X4","X11")) %>% mutate(count=sigcount-bkcount)
ggplot(bkcorrdf)+geom_tile(aes(X11,X4,fill=count))+scale_x_continuous(limits = c(-200,200),expand = c(0,0),breaks = c(-200,0,200))+scale_fill_gradientn(colours =c("white","#ead100","red","red"),na.value = "grey")

bkcorrdf=dplyr::inner_join(sigdf,bkdf,by=c("X4","X11")) %>% mutate(count=bkcount-sigcount)
ggplot(bkcorrdf)+geom_tile(aes(X11,X4,fill=count))+scale_x_continuous(limits = c(-200,200),expand = c(0,0),breaks = c(-200,0,200))+scale_fill_gradientn(colours =c("white","#ead100","red","red"))



# p=ggplot(sigdf)+stat_bin_2d(aes(X11,X4),binwidth = 2)+scale_fill_gradientn(colours =c("white","#ead100","red","red"),limits=c(0,30),name="Count",breaks=c(0,30))+scale_y_continuous(limits = c(0,200),expand = c(0,0),breaks = c(0,100,200))+scale_x_continuous(limits = c(-200,200),expand = c(0,0),breaks = c(-200,0,200))+
#     xlab("Distance from motif (bp)")+ ylab("Fragment length (bp)")+gg_theme_Publication()+theme(legend.direction = "horizontal",legend.position = c(0.8,0.2))
# print(p)
# # gg_save_pdf(p,6,5.5,filename = "ELF2_Vplot")
# 
# ggplot()+geom_density(aes(df %>% dplyr::filter(X4>100 & X4<180) %>% .$X11))+scale_x_continuous(limits = c(-200,200))
# 
# 
# 
# # range(sigdf$X4)
# image(k,useRaster = F)
# 
