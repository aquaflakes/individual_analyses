fjComm::clear_()

label384=fjComm::gen_384_label()

posCnt <- function(well,isTFctrl=F,dint="TA")
{
# check "TA" phase
readsfile= paste0("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147*-CAP2r-147c5Nu-",well,"-*.fq.gz")
# if(isTFctrl) readsfile= paste0("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147*-CAP2r-TFctrl147c4-",well,"-*.fq.gz")
readsfile= Sys.glob(readsfile)
seq= read_csv(readsfile,col_names = F) %>% rmdup(1)
seq_lines= nrow(seq)
sig_kcnt= kmerCntBit(strings = seq[[1]],k = 2L,diffLen = F,collapse = F)

c0_file= "~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/FJ4.4_PE_c0/dual_trim/2_Trim_adaptor_kmer/subset_100000_Trulig147v1IIIPEx200IIIc0_S3_R1_001.peared_trimmed.fq.gz"
seq_c0= read_csv(c0_file,col_names = F)
c0_lines= nrow(seq_c0)
c0_kcnt= oneOffCalc("c0_count/c0_count",calcFun = function(x){kmerCntBit(strings = seq_c0[[1]],k = 2L,diffLen = F,collapse = F)},asObj = T,overWrite = F,useScriptPath = T)

# corr for bk
sub_mat= sweep(c0_kcnt,2,colMeans(c0_kcnt),"-") / c0_lines * seq_lines # (each di-nt  -  mean) / c0_lines * seq_lines
sig_kcnt_corr= sig_kcnt- sub_mat
return(sig_kcnt_corr[,"TA"])
}

result= parallel(label384,fun = function(well) oneOffCalc(
  saved_oneOff_filename = paste0("indiv_TA_posCnt/",well,"_.robj"),
  calcFun = posCnt,
  calcFunParamList = list(well,isTFctrl=F, dint="TA"),asObj = T,useScriptPath = T
  ),workers = 25)
result= do.call(rbind,result) %>% colMeans()
TA_posCnt_allavg= gg_heat2D_diag(data_frame(pos1=seq_along(result),pos2=seq_along(result),result))+ theme(axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())
gg_save_all(TA_posCnt_allavg,width = 4.157,height = 0.6)

forScale=data_frame(pos1=seq_along(result),pos2=seq_along(result),result) %>% gg_heat2D_MI()
gg_save_pdf(forScale,width = 6,height = 6)
