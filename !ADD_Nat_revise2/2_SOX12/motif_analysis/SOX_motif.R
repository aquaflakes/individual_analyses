fjComm::clear_()

lig200=fjComm::nusel_get_200() %>% dplyr::filter(use==1 & order==1)
SOXs=lig200 %>% dplyr::filter(stringr::str_detect(TF,"SOX"))
# SOXs=SOXs[1:5,] %>% dplyr::filter(TF %in% c("SOX9","SOX10")) #"SOX7"
SOXs=SOXs[1:5,] %>% dplyr::filter(TF %in% c("SOX11","SOX12")) #"SOX7"




nu11=fjComm::motif_plot_moods_singlePWM_hit(SOXs[1,]$c4_file,"SOX11.pfm",combine_2strand_hits = T)
TF11=fjComm::motif_plot_moods_singlePWM_hit(SOXs[1,]$TFctrl_file,"SOX11.pfm",combine_2strand_hits = T)

nu12=fjComm::motif_plot_moods_singlePWM_hit(SOXs[2,]$c4_file,"SOX12.pfm",combine_2strand_hits = T)
TF12=fjComm::motif_plot_moods_singlePWM_hit(SOXs[2,]$TFctrl_file,"SOX12.pfm",combine_2strand_hits = T)

EMI11=SOXs[1,]$c4_EMI_TFctrl_file %>% read_csv() %>% fjComm::gg_heat2D_diag(grad_colors = gg_steelblue_red)
EMI12=SOXs[2,]$c4_EMI_TFctrl_file %>% read_csv() %>% fjComm::gg_heat2D_diag(grad_colors = gg_steelblue_red)

fjComm::gg_multiplot(nu11,nu12,TF11,TF12,EMI11, EMI12,cols = 3)
fjComm::gg_multiplot(nu11,nu12,TF11,TF12,cols = 2)
fjComm::gg_multiplot(EMI11, EMI12,cols = 1)

kcnt11=kmerCntBit_file(SOXs[1,]$TFctrl_file) %>% mutate(pos=1:nrow(.)) %>% melt(id="pos") %>% set_colnames(c("pos","nt","count"))
kcnt12=kmerCntBit_file(SOXs[2,]$TFctrl_file) %>% mutate(pos=1:nrow(.)) %>% melt(id="pos") %>% set_colnames(c("pos","nt","count"))
ggplot(kcnt12)+geom_line(aes(pos,count,color=nt))
p=kcnt12 %>% dplyr::filter(nt=="GG") %>% .$count %>% fjComm::gg_heat2D_diag_vect(grad_colors = gg_steelblue_red)
gg_save_diag(p,newName = "SOX12_GGcount")
