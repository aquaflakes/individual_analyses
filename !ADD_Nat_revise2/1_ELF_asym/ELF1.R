fjComm::clear_()

ELF5="~/Nut_zhuData/seqFiles2/FJ4.6CAPFLLin/allreads/pear_and_trimmed/*1nu50L*B21*ELF5IIIc4*" %>% Sys.glob() %>% .[1]
  motif_lin="~/Nut_zhuData/Analysis/Analysis2/Nusel_revise/1_motif_and_HTctrl/manual_curations/3_curated_motifs_and_make_table/NCAP_motifs/*ELF5*" %>% Sys.glob()
    # motif %>% read_tsv(col_names = F) %>% as.matrix() %>% fjComm::ggseqlogo_lab()
    pwm=motif_lin %>% read_tsv(col_names = F) %>% as.matrix() 
result=fjComm::motif_plot_moods_singlePWM_hit(ELF5,motif_lin)
result_mut=result
result_mut$data %<>% mutate(pos=ifelse(pos>60,pos+117,pos))
result_mut=result_mut+scale_x_continuous(breaks = c(1,141+117),expand = c(0,0))+geom_vline(xintercept = c(45,258-45),color="grey")
ELF5Lin=result_mut %>% gg_save_diag(width = 6.5, height = 2,newName = "ELF5Lin")

ELF="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/*-147c4-*ELF*IIIc4*" %>% Sys.glob()# %>% .[1]
motif="~/Nut_zhuData/Analysis/Analysis2/Nusel_revise/1_motif_and_HTctrl/manual_curations/3_curated_motifs_and_make_table/NCAP_motifs/*ELF*" %>% Sys.glob()
  order1=ELF %>% purrr::map(fjComm::extract_TF) %>% unlist()
  order2=motif %>% purrr::map(~str_extract(.,"(?<=\\.)ELF.(?=\\.c)")) %>% unlist()
results=purrr::map(1:4,~fjComm::motif_plot_moods_singlePWM_hit(ELF[.],motif[.]))

ELF1=results[[2]] %>% gg_save_diag(width = 3.5, height = 2,newName = "ELF1")
ELF2=results[[4]] %>% gg_save_diag(width = 3.5, height = 2,newName = "ELF2")
