fjComm::clear_()

MIfile="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/analysis/2d_IC_related_pub/data/C5_Trulig147v1IIIFJ4-4-CAP2-C05-TF15-SOX11IIIc4_S437_L002_R1_001.peared_trimmed.fq_u.gz_3mer_foldchn_maxBias.csv"
df=read_csv(MIfile)
SOX11_MI_diag= gg_heat2D_diag(df,grad_colors = gg_steelblue_red)+gg_theme_bordered_diag
gg_save_diag(SOX11_MI_diag)

MIfile="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/analysis/2d_IC_related_pub/data/C9_Trulig147v1IIIFJ4-4-CAP2-C09-TF17-SOX12IIIc4_S441_L002_R1_001.peared_trimmed.fq_u.gz_3mer_foldchn_maxBias.csv"
df=read_csv(MIfile)
SOX12_MI_diag= gg_heat2D_diag(df,grad_colors = gg_steelblue_red)+gg_theme_bordered_diag
gg_save_diag(SOX12_MI_diag)

MIfile="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/D8_Trulig147v1IIIFJ4-5-CAP2r-147c4-D8-SOX18IIIc4-AGGTGTAA_S80_L003_R1_001.peared_trimmed.fq.gz_3mer_foldchn_maxBias.csv"
df=read_csv(MIfile)
SOX18_MI_diag= gg_heat2D_diag(df,grad_colors = gg_steelblue_red)+gg_theme_bordered_diag
gg_save_diag(SOX18_MI_diag)

MIfile="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/L13_Trulig147v1IIIFJ4-5-CAP2r-147c4-L13-SOX10oooNIIIc4-GGCCAGCT_S277_L003_R1_001.peared_trimmed.fq.gz_3mer_foldchn_maxBias.csv"
df=read_csv(MIfile)
SOX10_MI_diag= gg_heat2D_diag(df,grad_colors = gg_steelblue_red)+gg_theme_bordered_diag
gg_save_diag(SOX10_MI_diag)

MIfile="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/P5_Trulig147v1IIIFJ4-5-CAP2r-147c4-P5-SOX7IIIc4-AAAATTGC_S365_L003_R1_001.peared_trimmed.fq.gz_3mer_foldchn_maxBias.csv"
df=read_csv(MIfile)
SOX7_MI_diag= gg_heat2D_diag(df,grad_colors = gg_steelblue_red)+gg_theme_bordered_diag
gg_save_diag(SOX7_MI_diag)
