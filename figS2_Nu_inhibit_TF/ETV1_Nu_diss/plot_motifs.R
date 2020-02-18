fjComm::clear_()
Nu_seed_motif=read_lines("~/Nut_zhuData/Analysis/Analysis2/Nusel/2017_0817_ETV1_diffmotif_Nu_Free/Freefile_Nuseed.txt")[19:22]
Nu_seed_motif %<>% str_split("\t") %>% as.data.frame() %>% t %>% set_rownames(.[,2]) %>% .[,-c(1,2)] %>% as.data.frame() %>%
  mutate_all(.funs = funs(as.character(.))) %>% mutate_all(.funs = funs(as.numeric(.)))
pdf(file = "Nu_seed_motif.pdf",height = 2.5)
fjComm::plotMotif_pfmMat(Nu_seed_motif)
dev.off()

Free_seed_motif=read_lines("~/Nut_zhuData/Analysis/Analysis2/Nusel/2017_0817_ETV1_diffmotif_Nu_Free/Freefile_Freeseed.txt")[19:22]
Free_seed_motif %<>% str_split("\t") %>% as.data.frame() %>% t %>% set_rownames(.[,2]) %>% .[,-c(1,2)] %>% as.data.frame() %>%
  mutate_all(.funs = funs(as.character(.))) %>% mutate_all(.funs = funs(as.numeric(.)))
pdf(file = "Free_seed_motif.pdf",height = 2.5)
fjComm::plotMotif_pfmMat(Free_seed_motif)
dev.off()
