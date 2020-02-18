# pacman::p_load(Rsamtools)
# which <- RangesList(seq1=IRanges(1000, 2000), seq2=IRanges(c(100, 1000), c(1000, 2000)))
# what <- c("rname", "strand", "pos", "qwidth", "seq")
# param <- ScanBamParam(which=which, what=what)
#
# bamFile <- system.file("extdata", "ex1.bam", package="Rsamtools")
# bam <- scanBam(bamFile, param=param)


pacman::p_load(nucleoSim)

val.num       <- 50     ### Number of well-positioned nucleosomes
val.del       <- 10     ### Number of well-positioned nucleosomes to delete
val.var       <- 30     ### variance associated to well-positioned nucleosomes
val.fuz       <- 10     ### Number of fuzzy nucleosomes
val.fuz.var   <- 50     ### variance associated to fuzzy nucleosomes
val.max.cover <- 70     ### Maximum coverage for one nucleosome
val.nuc.len   <- 147    ### Distance between nucleosomes
val.len.var   <- 10     ### Variance associated to the length of the reads
val.lin.len   <- 20     ### The length of the DNA linker
val.rnd.seed  <- 100    ### Set seed when result needs to be reproducible
val.offset    <- 10000  ### The number of bases used to offset
### all nucleosomes and reads

## Create sample using a Normal distribution
sample <- nucleoSim::syntheticNucReadsFromDist(wp.num=val.num,
                                               wp.del=val.del, wp.var=val.var,
                                               fuz.num=val.del, fuz.var=val.fuz.var,
                                               max.cover=val.max.cover,
                                               nuc.len=val.nuc.len,
                                               len.var=val.len.var,
                                               lin.len=val.lin.len,
                                               rnd.seed=val.rnd.seed,
                                               distr="Normal", offset=val.offset)

## Create visual representation of the synthetic nucleosome sample
# plot(sample)

library(GenomicRanges)

## Transform sample dataset into GRanges object
sampleGRanges <- GRanges(seqnames = sample$dataIP$chr,
                         ranges = IRanges(start = sample$dataIP$start,
                                          end = sample$dataIP$end),
                         strand = sample$dataIP$strand)

library(devtools)
install_github("Bioconductor-mirror/RJMCMCNucleosomes")
## Segment sample into candidate regions
sampleSegmented <- segmentation(reads = sampleGRanges, zeta = 147,
                                delta = 40, maxLength = 1000)

## Number of segments created
length(sampleSegmented)
