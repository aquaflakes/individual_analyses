fjComm::clear_()

seq="CTGGAGAATCCCGGTCTGCAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGGTATTGTTTATTTTGTTCCTCCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGT"
Biostrings::writeXStringSet(seq %>% Biostrings::DNAStringSet(),"SOXlig.fasta")

pfm=read_tsv("sox11_motif.pfm",col_names = F) %>% as.matrix() 
ggseqlogo_lab(pfm)

system("moods_dna.py -m sox11_motif.pfm -s SOXlig.fasta -p 1 >match_results.csv",intern = F)
result=read_csv("match_results.csv",col_names = F)
result %<>% mutate(X7=revComp(X6))
p=ggplot(result)+geom_line(aes(X3+1,X5,color=X4)) + scale_x_continuous(expand=c(0,0),breaks = c(1,140)) + scale_colour_hc(name="Strand")+xlab("(bp)")+ylab("Score")
gg_save_pdf(p,10,4,filename = "SOX_score_along_lig")

pdf("SOX11_motif.pdf",height = 4)
fjComm::plotMotif_pfmMat(pfm)
dev.off()
