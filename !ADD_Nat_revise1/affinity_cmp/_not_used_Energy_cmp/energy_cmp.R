fjComm::clear_()

  binders=fjComm::nusel_binders()
Periodic_binders=binders[binders=="Periodic/dyad"] %>% names()
FJ46guide=read_tsv("~/Nut_zhuData/SELEXphp/guide/FJ4.6CAPFL_lig147TAAT_Fangj201709.txt") %>% dplyr::filter(qPCR_reliable!="0"&successful!=0)
FJ46guide %<>% dplyr::filter(TF %in% c("RFX5")) #"PAX7",,"TGIF2","MYOG","LMX1B", ,"HOXC10"
FJ46TFs=FJ46guide$TF %>% str_replace("ooo.*$","")
FJ46guide=FJ46guide %>% dplyr::filter(FJ46TFs %in% Periodic_binders)




alldiss1=Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.6CAPFL/allreads/barSep147c5/*bTAAT*147c5after2*")

single_well <- function(rowInGuide)
{
  TFname=FJ46guide$TF[rowInGuide]; well=FJ46guide$pos[rowInGuide]
    testfile=alldiss1[alldiss1 %>% str_match(glue::glue("-{well}-")) %>% {!is.na(.)[,1]}] #RFX5
  freefile=testfile
  nufile=freefile %>% str_replace("bTAAT","bATTA") %>% Sys.glob()
  freebeforefile=str_replace(testfile,"after2(-.{2,3}-).*","before\\1*") %>% Sys.glob()
  nubeforefile=freebeforefile %>% str_replace("bTAAT","bATTA") %>% Sys.glob()

  bk_count=tibble(pos=0:87,dummy=0); pseudo=0

    pwmMap<-function(file)
    {
      nu=SELEXFile(file); nu$getSeq(dup_rm = F); nu$moodsMap(pwmFile = "../pwm/RFX5_motif_FJ4.4.txt", p = "0.0001",batch = T)
      attr(nu,"posCnt")=nu$moodsResult %>% group_by(pos) %>% summarise(count=n()) %>% arrange(pos) %>% dplyr::right_join(bk_count,by="pos") %>% mutate(count= ifelse(is.na(count),pseudo, count)) %>% select(pos,count)
      nu
    }

  nu=pwmMap(nufile); nu=attr(nu,"posCnt") %>% mutate(time="af",type="nu")
  free=pwmMap(freefile); free=attr(free,"posCnt") %>% mutate(time="af",type="free")
  nubf=pwmMap(nubeforefile); nubf=attr(nubf,"posCnt") %>% mutate(time="bf",type="nu")
  freebf=pwmMap(freebeforefile); freebf=attr(freebf,"posCnt") %>% mutate(time="bf",type="free")
    ratio_bf=nubf$count/freebf$count
    ratio_af=nu$count/free$count
     # browser()
  tibble(pos=nu$pos, log_ratio=log(ratio_af/ratio_bf))
}

data=FJ46guide %>% select("TF") %>% mutate(result=purrr::map(1:nrow(FJ46guide),single_well))
data=data %>% unnest(result)
ggplot(data)+geom_tile(aes(pos,1,fill=log_ratio))+facet_grid(TF~.)+scale_fill_gradientn(colours = c("green","white","red"),limits=c(-4,4),values = c(0,0.5,1))
# p_e=ggplot()+geom_tile(data=NULL,aes(nubf$pos,240,fill=ratio_bf %>% log2),height=20)+scale_fill_distiller(palette = "YlOrRd",direction = 1,limits=c(-4,1))#scale_fill_gradientn(colours = c("green","white","red"),values =c(0,0.8,1),limits=c(-4,1))
# fjComm::gg_heat2D_diag_vect(log(ratio_af/ratio_bf))+scale_fill_gradientn(colours = c("green","white","red"),limits=c(-1.5,1.5),values = c(0,0.5,1))
# (ratio_af/ratio_bf) %>% plot(type="l")

