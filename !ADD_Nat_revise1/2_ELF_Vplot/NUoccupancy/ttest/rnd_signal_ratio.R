fjComm::clear_()

nu_cut_low=140
nu_cut_high=170

# center point between 83 to 183
cut_low=83
cut_high=183

calc_nu_occ_ratio<-function(file,TF="rnd")
{
    sigdf=read_tsv(file,col_names = F) %>% dplyr::select(X11,X4) %>% dplyr::filter(X4<=200)
    sigdf %<>% dplyr::group_by(X11,X4) %>% summarise(count=n()) %>% ungroup() %>% mutate(freq=count/sum(count))
    
  sigdf %<>% dplyr::filter(X4>nu_cut_low & X4<nu_cut_high) %>% group_by(X11) %>% summarise(freq=sum(freq)) #%>% mutate(exp="ELF1")
  leftoccu=sigdf %>% dplyr::filter(X11>=-cut_high&X11<=-cut_low) %>% {sum(.$freq)}
  rightoccu=sigdf %>% dplyr::filter(X11>=cut_low & X11<=cut_high) %>% {sum(.$freq)}
  return(list(c(leftoccu=leftoccu,rightoccu=rightoccu,ratio=leftoccu/rightoccu,TF=TF)))
}

ratio_ELF1=calc_nu_occ_ratio("../ELF1/LoVoFixed.fragment_length.closest_dist_to_center_within_500.ELF1_Control_TAAGCA40NGTC_KW_NACCCGGAAGTN_m1_c3_short_overlap_LoVo_ChIP.no_blacklist.txt","ELF1")
ratio_ELF2=calc_nu_occ_ratio("../ELF2/signal.txt",TF = "ELF2")

randomFiles=Sys.glob("rnd_result_link/*.txt")
ratio_rnd=oneOffCalc("oneoff/rnd_results",function(){lapply(randomFiles,calc_nu_occ_ratio)})

ratio_rnd=purrr::map(ratio_rnd,~.[[1]]) %>% rbind_all()
p_ELF1=t.test(ratio_rnd$ratio %>% as.numeric,mu=ratio_ELF1[[1]]["ratio"] %>% as.numeric())
p_ELF2=t.test(ratio_rnd$ratio %>% as.numeric,mu=ratio_ELF2[[1]]["ratio"] %>% as.numeric())

# p=ggplot(sigdf)+geom_line(aes(X11,freq))
# plotly::ggplotly(p)
