fjComm::clear_()

# df=read_tsv("data/K562_rep1_20.6U_79.2U_304U_MNase.fragment_length.closest_dist_to_center_within_500.YBX1_Control_TTGCTC40NGAA_KR_YGTWCCAYC_m1_c4_short_overlap_K562_ChIP.no_blacklist.txt",col_names = F)
# df=read_tsv("data/K562_rep1_20.6U_79.2U_304U_MNase.fragment_length.closest_dist_to_center_within_500.CTCF_AJ_TAGCGA20NGCT_NGCGCCMYCTAGYGGTN_m2_c4_Cell2013_overlap_K562_ChIP.no_blacklist.txt",col_names = F)

calc_cvg <- function(df,TFname)
{
  df %<>% mutate(left_edge=X11-X4/2, right_edge=X11+X4/2) %>% dplyr::filter(left_edge >-1000 & right_edge<1000) %>% dplyr::filter(X4>140)
  # browser()
  rngs=IRanges::IRanges(start = df$left_edge,end = df$right_edge)
  cvgs=IRanges::coverage(rngs,shift = 1000); cvgLen=cvgs %>% length()
  cvgdf=tibble(pos=1:cvgLen-1000,cvg= cvgs %>% as.integer()) %>% mutate(TF=TFname)
}


df1=read_tsv("data/*YBX1*" %>% Sys.glob() %>% .[1],col_names = F) %>% calc_cvg("YBX1")
df2=read_tsv("data/*HMBOX1*" %>% Sys.glob() %>% .[1],col_names = F)%>% calc_cvg("HMBOX1")

df3=read_tsv("data/*PKNOX1*" %>% Sys.glob() %>% .[1],col_names = F)%>% calc_cvg("PKNOX1")
# df4=read_tsv("data/*YY1*" %>% Sys.glob() %>% .[1],col_names = F)%>% calc_cvg("YY1")
df4=read_tsv("data/*ATF2*" %>% Sys.glob() %>% .[1],col_names = F)%>% calc_cvg("ATF2")
df5=read_tsv("data/*ELF1*" %>% Sys.glob() %>% .[1],col_names = F)%>% calc_cvg("ELF1")
df6=read_tsv("data/*ATF3*" %>% Sys.glob() %>% .[1],col_names = F)%>% calc_cvg("ATF3")
df=rbind(df1,df2,df3,df4,df5,df6)


    


p=ggplot(df %>% dplyr::filter(pos>-350&pos<350))+geom_line(aes(pos,cvg))+
  facet_wrap(~TF,2,3,scales = "free_y")+
  scale_x_continuous(limits = c(-350,350),expand = c(0,0))+xlab("(bp)")+ylab("MNase-seq coverage")+gg_theme_Publication()
print(p)

gg_save_pdf(p,14,6,filename = "corr_examples")


