fjComm::clear_()
allTFs=fjComm::nusel_get_147()
PITX=allTFs %>% dplyr::filter(TF=="PITX2")
PITX_147_file=PITX[1,]$c5Nu_file; c0file="~/Nut_zhuData/seqFiles2/c0/subset_100000_Trulig147v1IIIPEx200IIIc0_S3_R1_001.peared_trimmed.fq.gz"

SELEX=SELEXFile(PITX_147_file,c0_SELEXFile_obj = SELEXFile(c0file))
SELEX$getSeq(dup_rm = T)
SELEX$count_k(k = 6,collapse = T,diffLen = F,asDf = T,pseudo = 10,all_possible_k = T,cmp_c0 = T,count_c0 = T)
top_kmers=SELEX$kmerCnt %>% arrange(desc(foldchn)) %>% head(6) %>% .$kmer

kmer_curr=top_kmers[1]
motif_hits=str_locate_all(SELEX$seq,kmer_curr) %>% {do.call(rbind,.)} %>% as.data.frame()
p=ggplot(motif_hits)+geom_histogram(aes(start),binwidth = 1)+ggtitle(kmer_curr)
gg_save_plotly(p,filename = "PITX2" %>% paste0(kmer_curr))

# kmer_curr="TA"
# motif_hits=str_locate_all(SELEX$seq,kmer_curr) %>% {do.call(rbind,.)} %>% as.data.frame()
# p=ggplot(motif_hits)+geom_histogram(aes(start),binwidth = 1)+ggtitle(kmer_curr)
# gg_save_plotly(p,filename = "PITX2" %>% paste0(kmer_curr))


EOMES=allTFs %>% dplyr::filter(TF=="EOMES")
EOMES_147_file=EOMES[2,]$c5Nu_file; c0file="~/Nut_zhuData/seqFiles2/c0/subset_100000_Trulig147v1IIIPEx200IIIc0_S3_R1_001.peared_trimmed.fq.gz"

SELEX=SELEXFile(EOMES_147_file,c0_SELEXFile_obj = SELEXFile(c0file))
SELEX$getSeq(dup_rm = T)
SELEX$count_k(k = 6,collapse = T,diffLen = F,asDf = T,pseudo = 10,all_possible_k = T,cmp_c0 = T,count_c0 = T)
top_kmers=SELEX$kmerCnt %>% arrange(desc(foldchn)) %>% head(6) %>% .$kmer

kmer_curr=top_kmers[1]
motif_hits=str_locate_all(SELEX$seq,kmer_curr) %>% {do.call(rbind,.)} %>% as.data.frame()
p=ggplot(motif_hits)+geom_histogram(aes(start),binwidth = 1)+ggtitle(kmer_curr)+scale_x_continuous(limits = c(0,90))
gg_save_plotly(p,filename = "EOMES" %>% paste0(kmer_curr))
