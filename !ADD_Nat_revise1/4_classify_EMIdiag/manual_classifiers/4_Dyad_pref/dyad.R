# intensity of 44-56 bp / total intensity

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


middle_str=NCAPdiags[,44:56] %>% rowSums()
allsig=NCAPdiags %>% rowSums()

output= (middle_str/allsig) %>% fjComm::namedVect_to_df("TF","dyad_sig_frac") #%>% arrange(desc(.[[2]]))
  # topTFs=output$TF[1:20]
  # ggheat(NCAPdiags[topTFs,])
output %<>% mutate(dyadPref=ifelse(.[[2]]>=0.217,"Yes","No"))

saveRDS(output,"4_dyad")
