fjComm::clear_()

SOX1245="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/*200*SOX12*IIIc4*" %>% Sys.glob()
SOX1145="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/*200*SOX11*IIIc4*" %>% Sys.glob()



fjComm::gg_heat2D_diag (SOX1145[7] %>% read_csv(),  grad_colors =gg_steelblue_red )
fjComm::gg_heat2D_diag (SOX1145[8] %>% read_csv(),  grad_colors =gg_steelblue_red )


plotdf=rbind(SOX1245[1] %>% read_csv() %>% mutate(type="NCAP"), SOX1245[2] %>% read_csv() %>% mutate(type="HT")) %>% dplyr::filter(pos2-pos1==3)

ggplot(plotdf) + geom_line(aes(pos1,topMIsum,color=type))+ylab("E-MI (bit)")+xlab("(bp)")
# gg_save_pdf()
