fjComm::clear_()

wells50L= "P5 E5 O11 I20 K4" %>% str_split(" ") %>% unlist
  files50L="~/Nut_zhuData/seqFiles2/FJ4.6CAPFLLin/analysis/2d_IC_related_pub_1Nu50L/data/*maxBias.csv" %>% Sys.glob()
  filter=(files50L %>% str_extract("(?<=1Nu50Lc4-)[^-]*")) %in% wells50L
files50L=files50L[filter]

data=tibble(TF=purrr::map(files50L,extract_TF) %>% unlist, result=purrr::map(files50L,read_csv)) %>% unnest(result) %>% dplyr::filter(pos2-pos1==3 & (pos1<=65 | pos1>=81)) %>% mutate(pos1=ifelse(pos1>=81,pos1+117,pos1))
data %<>% group_by(TF) %>% mutate(topMIsum=squishQuan_normalize(topMIsum)) %>% ungroup()
p=ggplot(data)+geom_tile(aes(pos1,TF,fill=topMIsum),height=0.5)+scale_fill_gradientn(colours = gg_steelblue_red,guide = F)+xlab("(bp)")+gg_theme_Publication()+scale_x_continuous(breaks = c(1,262))+
  theme(axis.title.y = element_blank(),axis.line.x = element_blank())+
  geom_vline(xintercept = c(56,198+9),color="grey")
p %T>% print()%>% gg_save_pdf(9,5,filename = "NuLin50")

