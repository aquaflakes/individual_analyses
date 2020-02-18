fjComm::clear_()
file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)
library(QuasR)
library(BSgenome)
library(Rsamtools)
library(rtracklayer)
library(GenomicFeatures)
library(Gviz)

sampleFile <- "extdata/samples_chip_single.txt"
auxFile <- "extdata/auxiliaries.txt"
genomeFile <- "extdata/hg19sub.fa"
proj1 <- qAlign(sampleFile, genome=genomeFile, auxiliaryFile=auxFile)

qQCReport(proj1, pdfFilename="extdata/qc_report.pdf")
alignmentStats(proj1)
qExportWig(proj1, binsize=100L, scaling=TRUE, collapseBySample=TRUE)


# qCount
pacman::p_load(GenomicFeatures,Rsamtools)
chrLen <- scanFaIndex(genomeFile)
chrominfo <- data.frame(chrom=as.character(seqnames(chrLen)),
                          length=width(chrLen),
                          is_circular=rep(FALSE, length(chrLen)))
annotFile <- "extdata/hg19sub_annotation.gtf"
txdb <- makeTxDbFromGFF(file=annotFile, format="gtf",
                          chrominfo=chrominfo,
                          dataSource="Ensembl",
                          organism="Homo sapiens")
promReg <- promoters(txdb, upstream=1000, downstream=500,
                       columns=c("gene_id","tx_id"))
gnId <- sapply(mcols(promReg)$gene_id, paste, collapse=",")
promRegSel <- promReg[ match(unique(gnId), gnId) ]
names(promRegSel) <- unique(gnId)

cnt <- qCount(proj1, promRegSel)


gr1 <- import("Sample1.wig.gz")
gr2 <- import("Sample2.wig.gz")
library(Gviz)
axisTrack <- GenomeAxisTrack()
dTrack1 <- DataTrack(range=gr1, name="Sample 1", type="h")
dTrack2 <- DataTrack(range=gr2, name="Sample 2", type="h")
txTrack <- GeneRegionTrack(txdb, name="Transcripts", showId=TRUE)
plotTracks(list(axisTrack, dTrack1, dTrack2, txTrack),
           chromosome="chr3", extend.left=1000)

# useful
library(GenomicAlignments)
bam <- readGAlignments("extdata/chip_1_1_d6c1510d82a.bam")
distance(MNaseGrange, TFmotif_pos) #distance(query, subject)
tt=as(bam,"GRanges")
