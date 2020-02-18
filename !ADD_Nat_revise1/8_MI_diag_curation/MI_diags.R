fjComm::clear_()

lig147=fjComm::nusel_get_147() %>% mutate(upload=ifelse((use==1 & order==1),1,0))
  allrep147 <- function(TFname) lig147 <<- lig147 %>%  mutate_cond(TF==TFname & use==1, upload=1) %>% mutate_cond(TF==TFname & use==1, TF=paste0(TF,"_rep",order))
for (i in c("T","TBX2","EOMES")) allrep147(i)
lig147 %<>% dplyr::filter(upload==1)

c4EMI= oneOffCalc("oneOff/c4EMI",calcFun = function(){ lig147$c4_EMI_file %>% purrr::map(function(x)read_csv(x) %>% dplyr::filter(pos2-pos1==3)) %>% set_names(lig147$TF) })
lig147df=purrr::map(c4EMI,function(x)x$topMIsum) %>% {do.call(rbind,.)}
lig147df=tibble::rownames_to_column(lig147df %>% as.data.frame() %>% set_colnames(1:96),"TF") %>% arrange(TF)
lig147df %>% writeClipboard()


lig200=fjComm::nusel_get_200() %>% mutate(upload=ifelse((use==1 & order==1),1,0))
lig200 %<>% dplyr::filter(upload==1)

c4EMI_200= oneOffCalc("oneOff/c4EMI_200",calcFun = function(){ lig200$c4_EMI_file %>% purrr::map(function(x)read_csv(x) %>% dplyr::filter(pos2-pos1==3)) %>% set_names(lig200$TF) })
lig200df=purrr::map(c4EMI_200,function(x)x$topMIsum) %>% {do.call(rbind,.)}
lig200df=tibble::rownames_to_column(lig200df %>% as.data.frame() %>% set_colnames(1:149),"TF") %>% arrange(TF)
lig200df %>% writeClipboard()


success_TFctrls_147="~/Nut_zhuData/SELEXphp/guide/FJ_HT_lig147_all_success.txt" %>% read_tsv() %>% dplyr::filter(is.na(comment))
  batch_=success_TFctrls_147$Batch %>% str_extract("FJ4\\..") %>% str_replace("\\.","")
  pos_=success_TFctrls_147$pos
lig147 %<>% dplyr::filter(paste0(batch,pos) %in% paste0(batch_,pos_) ) %>% mutate(TF=TF %>% str_replace("_rep.*",""))

c4EMICtrl= oneOffCalc("oneOff/c4EMICtrl147",calcFun = function(){ lig147$c4_EMI_TFctrl_file %>% purrr::map(function(x)read_csv(x) %>% dplyr::filter(pos2-pos1==3)) %>% set_names(lig147$TF) })
lig147dfctrl=purrr::map(c4EMICtrl,function(x)x$topMIsum) %>% {do.call(rbind,.)}
lig147dfctrl=tibble::rownames_to_column(lig147dfctrl %>% as.data.frame() %>% set_colnames(1:96),"TF") %>% arrange(TF)
lig147dfctrl %>% writeClipboard()



