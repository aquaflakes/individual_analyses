library(fjComm)
setwd(fjComm::get_scriptpath())

monomer= read_tsv("monomer.txt",col_names = F)
pdf("monomer.pdf",height = 4)
fjComm::plotMotif_pfmMat(monomer)
dev.off()

dimer= read_tsv("dimer.txt",col_names = F)
pdf("dimer.pdf",height = 4)
fjComm::plotMotif_pfmMat(dimer)
dev.off()

dual_helices= read_tsv("CAC76NGTGmul_n.txt",col_names = F)
pdf("dual_helices.pdf",height = 4)
fjComm::plotMotif_pfmMat(dual_helices)
dev.off()
