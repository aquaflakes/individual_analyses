fjComm::clear_()

lig200=fjComm::nusel_get_200() %>% dplyr::filter(use==1 & order==1)
SOXs=lig200 %>% dplyr::filter(stringr::str_detect(TF,"SOX"))
SOXs=SOXs[1:5,] %>% dplyr::filter(TF %in% c("SOX9","SOX10")) #"SOX7"

TFnames=SOXs$TF
NCAP=SOXs$c4_EMI_file %>%  purrr::map2(TFnames, function(x,y) read_csv(x,col_names = T) %>% mutate(TF=y) %>% dplyr::filter(pos2-pos1==3)) %>% {do.call(rbind,.)} %>% mutate(type="NCAP")
HT=SOXs$c4_EMI_TFctrl_file %>%  purrr::map2(TFnames, function(x,y) read_csv(x,col_names = T) %>% mutate(TF=y) %>% dplyr::filter(pos2-pos1==3)) %>% {do.call(rbind,.)} %>% mutate(type="HT")
results=rbind(NCAP,HT)

p=ggplot(results)+geom_line(aes(pos1,topMIsum,color=type))+facet_wrap(~TF)+ylab("E-MI (bit)")+xlab("(bp)") +scale_x_continuous(breaks=c(1,149))+ theme(panel.spacing = unit(1, "lines"))
gg_save_pdf(p,8,4)
