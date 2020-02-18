fjComm::clear_()

MIfile="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/D4_Trulig147v1IIIFJ4-5-CAP2r-147c4-D4-EOMESoooN2XIIIc4-ATCACTGA_S76_L003_R1_001.peared_trimmed.fq.gz_3mer_foldchn_maxBias.csv"
df=read_csv(MIfile)
EOMES_MI_diag= gg_heat2D_diag(df,grad_colors = gg_steelblue_red)
gg_save_diag(EOMES_MI_diag)

MIfile="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/P14_Trulig147v1IIIFJ4-5-CAP2r-147c4-P14-PITX2IIIc4-GGGTCGTG_S374_L003_R1_001.peared_trimmed.fq.gz_3mer_foldchn_maxBias.csv"
df=read_csv(MIfile)
PITX2_MI_diag= gg_heat2D_diag(df,grad_colors = gg_steelblue_red)
gg_save_diag(PITX2_MI_diag)
