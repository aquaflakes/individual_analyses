fjComm::clear_()

label384= fjComm::gen_384_label()

get_TFMI_diag <- function(well,add="Free")
{
  file=paste0("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/",well,"_*-147c5",add,"-*foldchn_maxBias.csv")
  file=Sys.glob(file)
  well= extract_well(file)
  if(length(file)==0) return(NA)
  TFname= extract_TF(file)
  tfmi=rio::import(file) %>% dplyr::filter(pos2-pos1==3)
  return(tfmi$topMIsum)
}

resultFree= oneOffCalc("all_TFMI_diag_Free",calcFun =function(x){
  fjComm::parallel(label384,get_TFMI_diag,workers = 25,add="Free") %>% {do.call(rbind,.)}
} ,asObj = T,useScriptPath = T)

TFMIall=data.frame(pos=20:76,TFMI=colMeans(resultFree,na.rm = T)[20:76])
# TFMIall=data.frame(pos=1:ncol(result),TFMI=colMeans(result,na.rm = T))
allTFMI_center_ub= TFMIall %>% ggplot(aes(pos,TFMI))+geom_line()+ xlab("position (bp)") + ylab("TF-MI (bit)")
gg_save_pdf(allTFMI_center_ub,width = 5,height = 3.5)



resultNu= oneOffCalc("all_TFMI_diag_Nu",calcFun =function(x){
  fjComm::parallel(label384,get_TFMI_diag,workers = 25,add="Nu") %>% {do.call(rbind,.)}
} ,asObj = T,useScriptPath = T)

TFMIall=data.frame(pos=20:76,TFMI=colMeans(resultNu,na.rm = T)[20:76])
# TFMIall=data.frame(pos=1:ncol(result),TFMI=colMeans(result,na.rm = T))
allTFMI_center_b= TFMIall %>% ggplot(aes(pos,TFMI))+geom_line()+ xlab("position (bp)") + ylab("TF-MI (bit)")
gg_save_pdf(allTFMI_center_b,width = 5,height = 3.5)



TFMIall_ub=data.frame(pos=1:ncol(resultFree),TFMI=colMeans(resultFree,na.rm = T))
TFMIall_b=data.frame(pos=1:ncol(resultNu),TFMI=colMeans(resultNu,na.rm = T))

TFMI_b_vs_ub= TFMIall_b %>% mutate(TFMI=TFMI-TFMIall_ub$TFMI)
TFMI_b_vs_ub_diffplot= TFMI_b_vs_ub$TFMI %>% gg_heat2D_diag_vect()+ xlab("position (bp)") + ylab("TF-MI (bit)")
gg_save_pdf(TFMI_b_vs_ub_diffplot,4,1.5)


