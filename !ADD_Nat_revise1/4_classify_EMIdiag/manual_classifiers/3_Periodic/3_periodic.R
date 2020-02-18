# > median value of sig_bk_ratio are periodic


# minus mean to reduce signal around 0
# FFT power normalized by dividing (mean * length of the vector)
# fit line segment
# take 0.08-0.12 as signal
# take the range 0.14-0.3 /bp-1 as background to normalize

# The diagonal of TF-MI for each TF was windowed by Welch’s function, subtracted with the mean, and then subjected to FFT. The obtained power spectrum was further normalized with the mean of TF-MI and the length of the diagonal. The strength of the 10-bp periodicity was estimated with the FFT-AUC of 0.08–0.12 bp-1 (corrected for the background level estimated from 0.14–0.3 bp-1). The phase of FFT was examined at 0.102 bp-1. The same process was applied to the TA dinucleotide counts across all positions of each TF’s library.


fjComm::clear_()
# lig147= fjComm::nusel_get_147() %>% dplyr::filter(use==1&order==1)
c4EMI= oneOffCalc("../../PCA/oneOff/c4EMI",calcFun = function(){ lig147$c4_EMI_file %>% purrr::map(read_csv) %>% set_names(lig147$TF) })
c4EMICtrl= oneOffCalc("../../PCA/oneOff/c4EMICtrl",calcFun = function(){ lig147$c4_EMI_TFctrl_file %>% purrr::map(read_csv) %>% set_names(lig147$TF) })

TFnames=names(c4EMI)
c4EMI_diag= oneOffCalc("../../PCA/oneOff/c4EMI_d",calcFun = function(){ c4EMI %>% purrr::map(~dplyr::filter(.,pos2-pos1==3)) %>% set_names(TFnames) })
c4EMICtrl_diag= oneOffCalc("../../PCA/oneOff/c4EMICtrl_d",calcFun = function(){ c4EMICtrl %>% purrr::map(~dplyr::filter(.,pos2-pos1==3)) %>% set_names(TFnames) })

NCAPdiags= vector("list",length(c4EMI_diag))
for(i in seq_along(NCAPdiags))
{
  NCAPdiags[[i]]= c4EMI_diag[[i]][,3]%>% set_colnames(TFnames[i])
}
NCAPdiags= do.call(cbind,NCAPdiags)
NCAPdiags %<>% mutate_all(function(x)squishQuan_normalize(x))
NCAPdiags %<>% t %>% set_colnames(1:ncol(.))



fft_result= apply(NCAPdiags,1,function(tfmi)fjComm::fft_snRatio_phase(c(0,0,tfmi,0,0),periodicity = 10.2))
fft_result= fft_result %>% t %>% .[,1:2] %>% as.data.frame()
  median_TA_phase=-0.04850917
fft_result$phase_rel_TA_in_pi= (fft_result[,2]-median_TA_phase) %>% {ifelse(abs(.-2)<abs(.),.-2,.)}
fft_result$TF=rownames(fft_result)

fft_result %<>% select(TF, sig_bk_ratio, phase_rel_TA_in_pi) %>% as_tibble()
  periodic_thr=fft_result$sig_bk_ratio %>% quantile(0.5)
fft_result %<>% mutate(isPeriodic=ifelse(sig_bk_ratio>=periodic_thr,"Yes","No"))
fft_result %<>% mutate(groovePref=ifelse(abs(phase_rel_TA_in_pi)>=0.5,"Minor","Major")) %>% mutate(groovePref=ifelse(isPeriodic=="Yes",groovePref,NA))
  fft_result %<>% .[,c(1,2,4,3,5)]
  # periodic_TFs=fft_result %>% dplyr::filter(isPeriodic=="Yes") %>% .$TF
  # ggheat(NCAPdiags[periodic_TFs,],clustering = "row",rescaling = "row")+scale_fill_gradientn(colours = gg_steelblue_red)
  # fft_raw(NCAPdiags["PITX2",]) %>% .$norm_pow %>% plot(type="l")
output=fft_result

saveRDS(output,"3_periodicity")
