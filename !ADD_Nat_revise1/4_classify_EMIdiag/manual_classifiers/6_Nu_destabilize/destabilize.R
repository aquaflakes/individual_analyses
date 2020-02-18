# c5EMI Free > Nu

fjComm::clear_()
# lig147= fjComm::nusel_get_147() %>% dplyr::filter(use==1&order==1)
c5EMI= oneOffCalc("../../PCA/oneOff/c5EMI_Nu",calcFun = function(){ lig147$c5_EMI_Nu_file %>% purrr::map( function(x) if(is.na(x)) NA else read_csv(x) ) %>% set_names(lig147$TF) })
c5EMIFree= oneOffCalc("../../PCA/oneOff/c5EMI_free",calcFun = function(){ lig147$c5_EMI_Free_file %>% purrr::map(function(x) if(is.na(x)) NA else read_csv(x)) %>% set_names(lig147$TF) })

result= vector("list",length(c5EMI))
for (i in 1:length(c5EMI))
{
  if ((!is.data.frame(c5EMI[[i]]))) {result[[i]]=NA}
  else
  {
    EMI_Free= c5EMIFree[[i]] %>% {sum(.$topMIsum)}
    EMI_Nu= c5EMI[[i]] %>% {sum(.$topMIsum)}
    result[[i]]=log2(EMI_Free/EMI_Nu)
  }
}

output=result %>% set_names(names(c5EMI)) %>% unlist() %>% namedVect_to_df("TF","log2_Free_Nu_ratio") %>% mutate(nucl_diss=ifelse(.[[2]]>0,"Destabilizer","Stabilizer"))

saveRDS(output,"6_diss")
