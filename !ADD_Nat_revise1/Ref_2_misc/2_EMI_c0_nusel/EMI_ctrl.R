fjComm::clear_()
# file="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_147/Trulig147v1IIIFJ4-4-CAP2-C05-TF15-SOX11IIIc4_S437_L002_R1_001.peared_trimmed.fq_u.gz"
# c4file="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/FJ4.4_PE_TFctrl/Trulig147v1IIIFJ4-4-TFctrl-C05-TF15-SOX11IIIc3_S53_L004_R1_001.peared_trimmed.fq.gz"
# file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4-P16-DMBX1IIIc4-AACAAGCA_S376_L003_R1_001.peared_trimmed.fq.gz"
# c4file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-TFctrl147c4-P16-DMBX1IIIc4-AACAAGCA_S1912_L007_R1_001.peared_trimmed.fq.gz"
file="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_147/Trulig147v1IIIFJ4-4-CAP2-G14-TF139-HSF1IIIc4_S542_L002_R1_001.peared_trimmed.fq_u.gz"
c4file="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/FJ4.4_PE_TFctrl/Trulig147v1IIIFJ4-4-TFctrl-G14-TF139-HSF1IIIc3_S158_L004_R1_001.peared_trimmed.fq.gz"


c0file="~/Nut_zhuData/seqFiles2/c0/subset_500000_Trulig147v1IIIPEx200IIIc0_S3_R1_001.peared_trimmed.fq.gz"
nufile="~/Nut_zhuData/seqFiles2/FJ6.1Methyl_PE/allreads/FJ6.1_OriMe1_PE/dual_trim/2_Trim_adaptor_kmer/Trulig147v1IIIFJ6-1-OriMe1-B02-1xCAGOct32ngIIIc4-Bq4-ATCACTGA_S782_L008_R1_001.peared_trimmed.fq.gz"


plotEMI <-function(file,TFname="HSF1")
{
  seq= readr::read_csv(file,col_names = F)
  seq=rmdup(seq,1)
  result_TFMI= oneOffCalc(paste0(basename(file),"_TFMI"),calcFun = ic_related_calc, calcFunParamList = list(seqs = seq[[1]],kmerLen = 3,filter_for_spacing = F,spacing = 0:10,type = "maxBias",maxBias_dimer_Params=list(type="topMI",topNo=10L)), asObj = T, useScriptPath = T)
  # result_TFMI=ic_related_calc(seqs = seq[[1]],kmerLen = 3,filter_for_spacing = F,spacing = 0:10,type = "maxBias",maxBias_dimer_Params=list(type="topMI",topNo=10L))
  TFMI_plot= gg_heat2D_MI(result_TFMI,grad_colors = gg_steelblue_red)#+guides(fill=F)
  print(TFMI_plot)
  gg_save_pdf(TFMI_plot,width = 7.2,height = 5.2,filename = TFname)
  result_TFMI
}

# applist=list(c(file,"SOX12"))
result= purrr::map2(c(file,c0file,nufile,c4file),c("HSF1","cyc0","Nucyc4","HSF1_TFctrl"),plotEMI)
result=purrr::map2(result,c("HSF1","cyc0","Nucyc4","HSF1_TFctrl"), function(x,y)x %>% mutate(TF=y)) 

plotdf=do.call(rbind,result) #%>% dplyr::filter(TF!="HSF1")
TF_wrap=factor(plotdf$TF,levels =c("cyc0","Nucyc4","HSF1","HSF1_TFctrl"),labels =c("Cycle 0","Nucl. SELEX","NCAP","HT")  )
plotdf$TF=TF_wrap

p=ggplot(plotdf)+geom_tile(aes(pos1,pos2,fill=topMIsum))+facet_wrap(~TF)+scale_fill_gradientn(colours =gg_steelblue_red,limits=c(0,0.05),oob=scales::squish,breaks=c(0,0.05),name="E-MI (bits)")+expand_limits(fill=0)+scale_x_continuous(expand = c(0,0),breaks = c(1,96))+scale_y_continuous(expand = c(0,0),breaks = c(4,99))+xlab("(bp)")+ylab("(bp)")+gg_theme_Publication()
gg_save_pdf(p,12.5,10,filename = "EMI_nucl_c0_ctrl")
