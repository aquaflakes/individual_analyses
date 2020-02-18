fjComm::clear_()
lig200= fjComm::nusel_get_200() %>% dplyr::filter(use==1&order==1)
c4EMI= oneOffCalc("PCA/oneOff/c4EMI",calcFun = function(){ lig200$c4_EMI_file %>% purrr::map(read_csv) %>% set_names(lig200$TF) })
# c4EMICtrl= oneOffCalc("PCA/oneOff/c4EMICtrl",calcFun = function(){ lig200$c4_EMI_TFctrl_file %>% purrr::map(read_csv) %>% set_names(lig200$TF) })

TFnames=names(c4EMI)
c4EMI_diag= oneOffCalc("PCA/oneOff/c4EMI_d",calcFun = function(){ c4EMI %>% purrr::map(~dplyr::filter(.,pos2-pos1==3)) %>% set_names(TFnames) })
# c4EMICtrl_diag= oneOffCalc("PCA/oneOff/c4EMICtrl_d",calcFun = function(){ c4EMICtrl %>% purrr::map(~dplyr::filter(.,pos2-pos1==3)) %>% set_names(TFnames) })

    NCAPdiags= vector("list",length(c4EMI_diag))
      for(i in seq_along(NCAPdiags))
      {
        NCAPdiags[[i]]= c4EMI_diag[[i]][,3]%>% set_colnames(TFnames[i])
      }
NCAPdiags= do.call(cbind,NCAPdiags)
# NCAPdiags %<>% mutate_all(function(x)squishQuan_normalize(x))
NCAPdiags %<>% t %>% set_colnames(1:ncol(.))



loessSingle <- function(pos,curr_col)
{
  # if (applyEnv$applyCnt==21) browser()
  # loess and penetration
  loessCurr= loess(curr_col~ pos,span=0.45)
  leftCutOff= (length(curr_col)/2-1) %>% floor()
  leftLoess= predict(loessCurr,seq_len(leftCutOff))
  rightLoess= predict(loessCurr,pos[-(seq_len(leftCutOff))])
  left_penetration= which(leftLoess>= (min(leftLoess)+max(leftLoess))/2) %>% max
  right_penetration= which(rightLoess>= (min(rightLoess)+max(rightLoess))/2) %>% length()
  penetration= mean(c(left_penetration,right_penetration))
  # applyEnv$applyCnt= 1+ applyEnv$applyCnt
  return(penetration)
}

EMIpVect=apply(NCAPdiags,1,function(x)loessSingle(seq_along(x),x))
EMIp = namedVect_to_df(EMIpVect, "TF","EMIp")
EMIp %<>% mutate(class=lapply(TF,TF_classify) %>% unlist()) %>% mutate(classbz=ifelse(class=="bZIP","bZIP","0"))
EMIp %<>% arrange(desc(EMIp)) %>% mutate(endPref=cut(EMIp,breaks = c(-Inf,20,Inf),labels = c("Yes","No")))
#TFs with a EMI penetration (lig200) less than 20 are defined as having End preference
EMIp %<>% mutate(binder= fjComm::nusel_binders()[TF])
# EMIp %<>% arrange(desc(EMIp)) %>% mutate(endPref=cut(EMIp,breaks = c(-Inf,10,25,Inf),labels = c("End binder","Intermediate","Central binder")))
# According to the EMI penetration (lig200), TFs are classified into "Central binder" (>= 25), "Intermediate" (>=10), and "End binder" (<10)
output=EMIp %>% select(TF,class,EMIp,endPref) %>% dplyr::rename(EMIp147=EMIp)

saveRDS(output,"1_EMI_p")
