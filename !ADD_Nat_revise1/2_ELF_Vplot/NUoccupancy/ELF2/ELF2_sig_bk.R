fjComm::clear_()

nu_cut_low=140
nu_cut_high=170

evalsig<-function()
{
  sigdf=read_tsv("signal.txt",col_names = F) %>% dplyr::select(X11,X4) %>% dplyr::filter(X4<=200)
  sigdf %<>% dplyr::group_by(X11,X4) %>% summarise(count=n()) %>% ungroup(); sigdf
}
sigdf=fjComm::oneOffCalc("oneoff/sigdensity",evalsig) %>% mutate(freq=count/sum(count))

evalbk<-function()
{
  bkdf=read_tsv("bk.txt",col_names = F) %>% dplyr::select(X11,X4) %>% dplyr::filter(X4<=200)
  bkdf %<>% dplyr::group_by(X11,X4) %>% summarise(count=n()) %>% ungroup(); bkdf
}
bkdf=fjComm::oneOffCalc("oneoff/bkdensity",evalbk) %>% mutate(freq=count/sum(count))

evalsalt<-function()
{
  bkdf=read_tsv("500mM37degELF2.txt",col_names = F) %>% dplyr::select(X11,X4) %>% dplyr::filter(X4<=200)
  bkdf %<>% dplyr::group_by(X11,X4) %>% summarise(count=n()) %>% ungroup(); bkdf
}
saltdf=fjComm::oneOffCalc("oneoff/saltdensity",evalsalt) %>% mutate(freq=count/sum(count))


sigdf %<>% dplyr::filter(X4>nu_cut_low & X4<nu_cut_high) %>% group_by(X11) %>% summarise(freq=sum(freq)) %>% mutate(exp="Sites")
# ggplot(sigdf)+geom_line(aes(X11,freq))

bkdf %<>% dplyr::filter(X4>nu_cut_low & X4<nu_cut_high) %>% group_by(X11) %>% summarise(freq=sum(freq)) %>% mutate(exp="Non-sites")
# ggplot(bkdf)+geom_line(aes(X11,freq))

saltdf %<>% dplyr::filter(X4>nu_cut_low & X4<nu_cut_high) %>% group_by(X11) %>% summarise(freq=sum(freq)) %>% mutate(exp="Sites (salt treatment)")
# ggplot(saltdf)+geom_line(aes(X11,freq))
plotdf=rbind(sigdf,bkdf,saltdf) %>% mutate(exp=factor(exp,levels=c("Sites","Non-sites","Sites (salt treatment)")))


p=ggplot(plotdf)+geom_smooth(aes(X11,freq*10000,color=exp),span=0.05,se=T,size=0.7)+facet_wrap(~exp,ncol = 1)+scale_y_continuous(expand = c(0,0),breaks = c(0.0002,0.0004)*10000)+
  scale_x_continuous(expand = c(0,0))+scale_color_manual(values = c("#F8766D", "grey40", "grey40"),guide=F)+xlab("Dist. from motif (bp)")+ylab("Nucleosome occupancy       ")+gg_theme_Publication()+
  theme(strip.background = element_blank())
gg_save_pdf(p,5.15,6.4,filename = "ELF2_asym")



