fjComm::clear_()

# pallete=RColorBrewer::brewer.pal(n = 11, name = "Spectral") %>% rev()
df=read_tsv("K562_rep1_20.6U_79.2U_304U_MNase.fragment_length.closest_dist_to_center_within_500.YBX1_Control_TTGCTC40NGAA_KR_YGTWCCAYC_m1_c4_short_overlap_K562_ChIP.no_blacklist.txt",col_names = F)
df %<>% mutate(left_edge=X11-X4/2, right_edge=X11+X4/2) %>% dplyr::filter(left_edge >-1000 & right_edge<1000) %>% dplyr::filter(X4>140)
p=ggplot(df)+stat_bin_2d(aes(X11,X4),binwidth = 2)+scale_fill_gradientn(colours =c("white","#ead100","red"),limits=c(0,200),oob=scales::squish, name="Count",breaks=c(0,25))+scale_y_continuous(limits = c(0,200),expand = c(0,0),breaks = c(0,100,200))+scale_x_continuous(limits = c(-600,600),expand = c(0,0),breaks = c(-200,0,200))+
    xlab("Distance from motif (bp)")+ ylab("Fragment length (bp)")+gg_theme_Publication()+theme(legend.direction = "horizontal",legend.position = c(0.8,0.2))
print(p)
gg_save_pdf(p,6,5.5,filename = "YBX1_Vplot")

data=rbind( tibble(side="5'",edge=df$left_edge),tibble(side="3'",edge=df$right_edge) )
p=ggplot(data)+geom_density(aes(edge,color=side),bw=5)+geom_vline(xintercept = 0)+scale_x_continuous(limits = c(-350,350))
gg_save_pdf(p,12,5.5,filename = "YBX1_Vplot_cutpos")

rngs=IRanges::IRanges(start = df$left_edge,end = df$right_edge)
cvgs=IRanges::coverage(rngs,shift = 1000); cvgLen=cvgs %>% length()
cvgdf=tibble(pos=1:cvgLen-1000,cvg= cvgs %>% as.integer())
p=ggplot(cvgdf)+geom_line(aes(pos,cvg))+scale_x_continuous(limits = c(-350,350))+scale_y_continuous(limits = range(cvgdf %>% dplyr::filter(pos>-300&pos<300) %>% .$cvg ))+geom_vline(xintercept = 0)
gg_save_pdf(p,12,5.5,filename = "YBX1_Vplot_cvg")

p=ggplot(cvgdf)+geom_line(aes(pos,cvg))+scale_x_continuous(limits = c(-350,350))+scale_y_continuous(limits = range(cvgdf %>% dplyr::filter(pos>-300&pos<300) %>% .$cvg ))+xlab("(bp)")+theme(axis.line.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())
gg_save_pdf(p,3,2,filename = "YBX1_Vplot_cvg")
