fjComm::clear_()

End_pref= readRDS("1_End_pref/1_EMI_p")
Asym= readRDS("2_Asym_effect/2_asym")
Period= readRDS("3_Periodic/3_periodicity")
Dyad= readRDS("4_Dyad_pref/4_dyad")
Gyres= readRDS("5_Gyre_span/5_gyres")
Diss= readRDS("6_Nu_destabilize/6_diss")


classify= End_pref %>% dplyr::full_join(Period,by="TF") %>%
  dplyr::full_join(Dyad,by="TF") %>%
  dplyr::full_join(Gyres,by="TF") %>%
  dplyr::full_join(Diss,by="TF") %>%
  dplyr::full_join(Asym,by="TF")
classify$class= classify$TF %>% lapply(TF_classify) %>% unlist()
classify %>% writeClipboard()

saveRDS(classify,"all_classes")
