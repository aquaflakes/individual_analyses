fjComm::clear_()
pacman::p_load(QuasR)

R1_file= "~/Nut_zhuData/forEevi/FJ7.1-1_MNase_test/Fangjie_FJ7.1-1_MNaseT_PE1/170222_K00110_0117_BHGLWJBBXX_Fangjie_MNase/FJ7.1-1-MNaseTest/FJ7-1-1-MNaseTest-E7-LoVo-x16MNase-ATATTCCG_S2_L007_R1_001.fastq.gz"
R2_file= "~/Nut_zhuData/forEevi/FJ7.1-1_MNase_test/Fangjie_FJ7.1-1_MNaseT_PE1/170222_K00110_0117_BHGLWJBBXX_Fangjie_MNase/FJ7.1-1-MNaseTest/FJ7-1-1-MNaseTest-E7-LoVo-x16MNase-ATATTCCG_S2_L007_R2_001.fastq.gz"
outDir=script_path_from_fun; outDir=paste0(outDir,"/trimmed/"); system(paste0("mkdir -p ",outDir))
fastqc=F


## trim_garole
# cmd=paste0("trim_galore -q 30 --paired --stringency 5 --illumina ",if(fastqc)"--fastqc"else""," --gzip -o ",outDir," ",R1_file," ",R2_file)
# system(cmd)

library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
sampleDf=data.frame(FileName1="trimmed/FJ7-1-1-MNaseTest-E7-LoVo-x16MNase-ATATTCCG_S2_L007_R1_001_val_1.fq.gz",	FileName2="trimmed/FJ7-1-1-MNaseTest-E7-LoVo-x16MNase-ATATTCCG_S2_L007_R2_001_val_2.fq.gz",	SampleName="MNase")
write_tsv(sampleDf,paste0(script_path_from_fun,"/samples_tmp.tsv"),col_names = T)

# sampleFile <- "extdata/samples_chip_single.txt"
# auxFile <- "extdata/auxiliaries.txt"
# genomeFile <- "extdata/hg19sub.fa"
proj1 <- qAlign("samples_tmp.tsv", genome= "BSgenome.Hsapiens.UCSC.hg19.masked")

library(GenomicAlignments)
library(Rsamtools)
bamFile= BamFile(file = "~/Nut_zhuData/fromOthers/Eevi/LoVo_DNase/LoVo_DNase_rmdup.bam",yieldSize = 10000)
# bamFile= BamFile(file = "trimmed/FJ7-1-1-MNaseTest-E7-LoVo-x16MNase-ATATTCCG_S2_L007_R1_001_val_1_3d850add768.bam",yieldSize = 1000000,asMates = T)
open(bamFile)
chunkPE <- readGAlignmentPairs(bamFile,use.names = )
close(bamFile)


library(GenomicAlignments)
library(Rsamtools)
bamFile= BamFile(file = "~/Nut_zhuData/seqFiles2/FJ7.1MNase/analysis/align/testout/FJ7-1-1-MNaseTest-B5-GP5dFixed-x4MNase-GCCGAATT_S1_L004.sorted.rmdup.bam",yieldSize = 10000)
# bamFile= BamFile(file = "trimmed/FJ7-1-1-MNaseTest-E7-LoVo-x16MNase-ATATTCCG_S2_L007_R1_001_val_1_3d850add768.bam",yieldSize = 1000000,asMates = T)
open(bamFile)
chunkPE <- readGAlignmentPairs(bamFile,use.names = )
close(bamFile)




# file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)
# sampleFile <- "extdata/samples_chip_single.txt"
# auxFile <- "extdata/auxiliaries.txt"
# genomeFile <- "extdata/hg19sub.fa"
# proj1 <- qAlign(sampleFile, genome=genomeFile, auxiliaryFile=auxFile)
