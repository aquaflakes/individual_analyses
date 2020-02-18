# intensity of pos1-pos2=75-85 bp / intensity of pos1-pos2=50-70 bp

fjComm::clear_()
# lig147= fjComm::nusel_get_147() %>% dplyr::filter(use==1&order==1)
c4EMI= oneOffCalc("../../PCA/oneOff/c4EMI",calcFun = function(){ lig147$c4_EMI_file %>% purrr::map(read_csv) %>% set_names(lig147$TF) })
c4EMICtrl= oneOffCalc("../../PCA/oneOff/c4EMICtrl",calcFun = function(){ lig147$c4_EMI_TFctrl_file %>% purrr::map(read_csv) %>% set_names(lig147$TF) })

# TFnames=names(c4EMI)
# c4EMI_diag= oneOffCalc("../../PCA/oneOff/c4EMI_d",calcFun = function(){ c4EMI %>% purrr::map(~dplyr::filter(.,pos2-pos1==3)) %>% set_names(TFnames) })
# c4EMICtrl_diag= oneOffCalc("../../PCA/oneOff/c4EMICtrl_d",calcFun = function(){ c4EMICtrl %>% purrr::map(~dplyr::filter(.,pos2-pos1==3)) %>% set_names(TFnames) })
#
#     NCAPdiags= vector("list",length(c4EMI_diag))
#       for(i in seq_along(NCAPdiags))
#       {
#         NCAPdiags[[i]]= c4EMI_diag[[i]][,3]%>% set_colnames(TFnames[i])
#       }
# NCAPdiags= do.call(cbind,NCAPdiags)
# NCAPdiags %<>% mutate_all(function(x)squishQuan_normalize(x))
# NCAPdiags %<>% t %>% set_colnames(1:ncol(.))

gyre_mod_str= purrr::map(c4EMI,  function(x)(x%>% dplyr::filter((pos2-pos1) %in% 75:85 ) %>% {sum(.$topMIsum)})/(x %>% dplyr::filter((pos2-pos1) %in% 50:70 ) %>% {sum(.$topMIsum)}))

output=gyre_mod_str %>% unlist() %>% namedVect_to_df("TF","gyre_span_mod_str") %>% mutate(Gyre_span_mode=ifelse(gyre_span_mod_str>0.31,"Yes","No"))
# c4EMI[["ZGLP1"]] %>% fjComm::gg_heat2D_MI()

saveRDS(output,"5_gyres")
