fjComm::clear_()
TFfiles="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147*-147c5Free*"
Nufiles="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147*-147c5Nu*"
outDir="1_FJ4.5_2gap_c5/"

kmerLen=3L #if gapped, kmerLen becomes the segment length
gapped=T;  gapNo = 2; gapMins = c(0,0); gapMaxs = c(5,5)

test=F

library(fjComm)
# kmerLen=9L; gapped=T
topPercent=0.001
  topRows= (4^kmerLen*topPercent) %>% as.integer()
  if (gapped) topRows= (4^(kmerLen*(gapNo+1))*topPercent) %>% as.integer()

allwells= expand.grid(letters[1:16] %>% toupper(),1:24)
allwells= paste0(allwells[[1]],allwells[[2]])

well="D23"
  well_alt= well %>% str_replace("(\\w)(\\d)","\\10\\2")
  TFctrl_file= c( Sys.glob(paste0(TFfiles,well,"-*") ), Sys.glob(paste0(TFfiles,well_alt,"-*") ))
  Nu_file= c(Sys.glob(paste0(Nufiles,well,"-*") ), Sys.glob(paste0(Nufiles,well_alt,"-*") ))
  if((length(TFctrl_file)+length(Nu_file))!=2) return(NA)

  xyplot=xyplot(xfile = Nu_file, yfile=TFctrl_file, test = test, kmerLen=kmerLen, gapped=gapped, gapNo = gapNo, gapMins = gapMins, gapMaxs = gapMaxs, topKmersEach = topRows)
  # print(xyplot)
  xyplot$layers[[1]]$mapping=NULL
  xyplot=xyplot+xlab("Counts (bound)")+ylab("Counts (unbound)")+scale_y_continuous(breaks = c(1000,3000,5000))+gg_theme_Publication()

  TF=extract_TF(TFctrl_file); well=extract_well(TFctrl_file)
  gg_save_pdf(xyplot,width = 7,height = 6,filename = "ETV1_free_vs_Nu")
  # gg_save_png(xyplot, ifelse(gapped,11,9),8, filename = paste0("plot/",well,"_",TF,"_kmerCmp_",basename(Nu_file)),path = outDir)
  # gg_save_plotly(xyplot,filename = paste0("plot/",well,"_",TF,"_kmerCmp_",basename(Nu_file)))

# lapply("H11",plotSingle)
# if(!test) fjComm::parallel(allwells,plotSingle,workers = 15) else lapply(allwells,plotSingle)


