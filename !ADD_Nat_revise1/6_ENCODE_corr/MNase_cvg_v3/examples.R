fjComm::clear_()

df=read_tsv("K562_merged_rep1_20.6U_79.2U_304U_MNase.sorted.rmdup.autosomes.over140fragments.fragment_coverage.All_motifs_overlap_K562_ChIP.autosomes_top500.v2.200bp_flank_from_center.average_per_TF.txt")
df %<>% dplyr::filter(TF %in% qw("YBX1 HMBOX1 PKNOX1 ATF2 ELF1 ATF3")) %>% set_colnames(qw("TF pos cvg sites"))
p=ggplot(df)+geom_line(aes(pos-199,cvg))+facet_wrap(~TF,2,3,scales = "free_y")+
  scale_x_continuous(limits = c(-200,200),expand = c(0,0))+xlab("(bp)")+ylab("MNase-seq coverage")+gg_theme_Publication()
print(p)

gg_save_pdf(p,14,6,filename = "corr_examples_not_used")


