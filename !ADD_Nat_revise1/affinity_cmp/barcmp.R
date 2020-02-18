fjComm::clear_()
pseudo=5
alldiss1=Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.6CAPFL/allreads/barSep147c5/*bTAAT*147c5after2*")
dir.create("Reads/")

  testfile=alldiss1[alldiss1 %>% str_match("-H2-") %>% {!is.na(.)[,1]}] #RFX5
freefile=testfile
nufile=freefile %>% str_replace("bTAAT","bATTA") %>% Sys.glob()
freebeforefile=str_replace(testfile,"after2(-.{2,3}-).*","before\\1*") %>% Sys.glob()
nubeforefile=freebeforefile %>% str_replace("bTAAT","bATTA") %>% Sys.glob()

  # # cp files
  # file.copy(c(freefile,nufile,freebeforefile,nubeforefile),"Reads")
  # allcycles= "/Users/zhu/Nut_zhuData/seqFiles2/FJ4.6CAPFL/allreads/barSep147c[1-4]/Trulig147*-CAPFL-147c*-H2-*" %>% Sys.glob()
  # dir.create("Reads/allcycles/")
  # file.copy(allcycles, "Reads/allcycles/")


  pwmMap<-function(file)
  {
    nu=SELEXFile(file); nu$getSeq(dup_rm = F); nu$moodsMap(pwmFile = "pwm/RFX5_motif_FJ4.4.txt", p = "0.0001",batch = T)
    attr(nu,"posCnt")=nu$moodsResult %>% group_by(pos) %>% summarise(count=n()) %>% arrange(desc(pos))
    nu
  }

nu=pwmMap(nufile); nu=attr(nu,"posCnt") %>% mutate(time="af",type="nu")
free=pwmMap(freefile); free=attr(free,"posCnt") %>% mutate(time="af",type="free")
nubf=pwmMap(nubeforefile); nubf=attr(nubf,"posCnt") %>% mutate(time="bf",type="nu")
freebf=pwmMap(freebeforefile); freebf=attr(freebf,"posCnt") %>% mutate(time="bf",type="free")
  ratio_bf=nubf$count/freebf$count
  ratio_af=nu$count/free$count


data=rbind(nu,free,nubf,freebf)
p_bf=ggplot()+geom_line(data=data %>% dplyr::filter(time=="bf"),aes(pos,count,color=type))+scale_color_manual(values = c("steelblue","orange"));
#p_bf=p_bf+geom_tile(data=NULL,aes(nubf$pos,240,fill=ratio_bf %>% log2),height=20)+scale_fill_distiller(palette = "YlOrRd",direction = 1,limits=c(-5,0.5))#scale_fill_gradientn(colours = c("green","white","red"),values =c(0,0.8,1),limits=c(-4,1))
p_bf=p_bf+scale_y_continuous(expand = c(0,0),breaks = c(0,500))+coord_cartesian(ylim= c(0,650))+scale_x_continuous(expand = c(0,0),breaks = c(0,80))+gg_theme_Publication()+xlab("(bp)")+ylab("Motif hits")
print(p_bf)
gg_save_pdf(p_bf,6,6,filename = "for_scale")
gg_save_pdf(p_bf+theme(legend.position = "none"),4.5,3.7,filename = "before_hits")

p_af=ggplot()+geom_line(data=data %>% dplyr::filter(time=="af"),aes(pos,count,color=type))+ scale_color_manual(values = c("steelblue","orange"));
#p_af=p_af+geom_tile(data=NULL,aes(nubf$pos,240,fill=ratio_af %>% log2),height=20)+ scale_fill_distiller(palette = "YlOrRd",direction = 1,limits=c(-5,0.5))#  scale_fill_gradientn(colours = c("green","white","red"),values =c(0,0.8,1),limits=c(-4,1))
p_af=p_af+scale_y_continuous(expand = c(0,0),breaks = c(0,2000))+coord_cartesian(ylim= c(0,2400))+scale_x_continuous(expand = c(0,0),breaks = c(0,80))+gg_theme_Publication()+xlab("(bp)")+ylab("Motif hits")
print(p_af)
gg_save_pdf(p_af+theme(legend.position = "none"),4.7,3.7,filename = "after_hits")

affinity_data=tibble(pos=nubf$pos, log_ratio=log(ratio_af/ratio_bf)) %>% arrange(pos)
aff_diff=fjComm::gg_heat2D_diag_vect(c(0,0,0,0,affinity_data$log_ratio))+scale_fill_gradientn(colours = c("blue","white","red"),limits=c(-3,3),values = c(0,0.5,1))
gg_save_diag(aff_diff)
aff_diff_scale=ggplot()+geom_raster(aes(1:92,1,fill=c(0,0,0,0,affinity_data$log_ratio)))+scale_fill_gradientn(colours = c("blue","white","red"),limits=c(-3,3),values = c(0,0.5,1))
gg_save_diag(aff_diff_scale)
