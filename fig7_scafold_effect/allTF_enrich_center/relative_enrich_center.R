fjComm::clear_()

label384= fjComm::gen_384_label()

get_TFMI_diag <- function(well)
{
  file=paste0("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/",well,"_*-147c4-*foldchn_maxBias.csv")
  file=Sys.glob(file)
  well= extract_well(file)
  if(length(file)==0) return(NA)
  TFname= extract_TF(file)
  tfmi=rio::import(file) %>% dplyr::filter(pos2-pos1==3)
  return(tfmi$topMIsum)
}

result= oneOffCalc("all_TFMI_diag",calcFun =function(x){
  fjComm::parallel(label384,get_TFMI_diag,workers = 25) %>% {do.call(rbind,.)}
} ,asObj = T,useScriptPath = T)

TFMIall=data.frame(pos=20:76,TFMI=colMeans(result,na.rm = T)[20:76])
# TFMIall=data.frame(pos=1:ncol(result),TFMI=colMeans(result,na.rm = T))
allTFMI_center= TFMIall %>% ggplot(aes(pos,TFMI))+geom_line()+ xlab("position (bp)") + ylab("TF-MI (bit)")
gg_save_pdf(allTFMI_center,width = 5,height = 3.5)




# hwm=15; TFMIall %>% ggplot(aes(pos,TFMI))+geom_line()+geom_line(aes(pos,TFMIall$TFMI %>% math_get_baseline(hwm)),color="red")+geom_line(aes(pos,TFMIall$TFMI %>% math_sub_baseline(hwm)),color="red")
