overWrite=F

library(fjComm)
# reads_file="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_147/Trulig147v1IIIFJ4-4-CAP2-P14-TF379-PITX2IIIc4_S758_L002_R1_001.peared_trimmed.fq_u.gz"
# reads_file="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_147/Trulig147v1IIIFJ4-4-CAP2-C03-TF14-SCXBIIIc4_S435_L002_R1_001.peared_trimmed.fq_u.gz"

TF_reads_file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4-P24-SHOXIIIc4-GGTATGAA_S384_L003_R1_001.peared_trimmed.fq.gz"
# TF_reads_file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4-P14-PITX2IIIc4-GGGTCGTG_S374_L003_R1_001.peared_trimmed.fq.gz"
TFctrl_reads_file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-TFctrl147c1-P24-SHOXIIIc1-GGTATGAA_S384_L004_R1_001.peared_trimmed.fq.gz"

Nu_reads_file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4-F20-NOTFIIIc4-GCTCGATA_S140_L003_R1_001.peared_trimmed.fq.gz"

c0File="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/FJ4.4_PE_c0/dual_trim/2_Trim_adaptor_kmer/subset_100000_Trulig147v1IIIPEx200IIIc0_S3_R1_001.peared_trimmed.fq.gz"


gapLen_enrich_plot <- function(reads_file,c0file,dup_rm=T,topNo=10L)
{
  selex=SELEXFile(reads_file,SELEXFile(c0file)); selex$getSeq(dup_rm = dup_rm)
  selex$count_gpk(k = 3L,gapNo = 1,gapMins = 0,gapMaxs = selex$seqLen-6,diffLen = F,posInfo = F,all_possible_k = T,pseudo = 10,cmp_c0 = T,count_c0 = T)

          # # foldchn from curr cycle
          # selex$count_k(k = 3L, collapse = T,asDf = T,diffLen = F,pseudo = 10,all_possible_k = T)
          # selex$kmerCnt= selex$kmerCnt %>% mutate(freq=prop.table(counts))
          # selex$kmerCntGap= selex$kmerCntGap %>% mutate(freq_act=prop.table(counts))
          # selex$kmerCntGap$freq_exp=selex$kmerCnt[match(selex$kmerCntGap$kmer1,selex$kmerCnt$kmer),]$freq * selex$kmerCnt[match(selex$kmerCntGap$kmer2,selex$kmerCnt$kmer),]$freq
          # selex$kmerCntGap= selex$kmerCntGap %>% mutate(foldchn=log2(freq_act/freq_exp))

  plotdf=selex$kmerCntGap %>% group_by(gapLen) %>% top_n(topNo,foldchn) %>% summarise(foldchn= mean(foldchn))


  return(list(plotdf=plotdf,resultdf=selex$kmerCntGap))
}

TFresult=oneOffCalc(saved_oneOff_filename = paste0("tmp/",basename(TF_reads_file)),calcFun = gapLen_enrich_plot, calcFunParamList = list(TF_reads_file,c0File), asObj = T, useScriptPath = T,overWrite = overWrite)
TFctrlresult=oneOffCalc(saved_oneOff_filename = paste0("tmp/",basename(TFctrl_reads_file)),calcFun = gapLen_enrich_plot, calcFunParamList = list(TFctrl_reads_file,c0File), asObj = T, useScriptPath = T,overWrite = overWrite)
Nuresult=oneOffCalc(saved_oneOff_filename = paste0("tmp/",basename(Nu_reads_file)),calcFun = gapLen_enrich_plot, calcFunParamList = list(Nu_reads_file,c0File), asObj = T, useScriptPath = T,overWrite = overWrite)

gap_plot= ggplot()+ geom_line(data = TFresult$plotdf,aes(gapLen,foldchn,color="CAP"))+geom_line(data = TFctrlresult$plotdf,aes(gapLen,foldchn,color="TFctrl"))+ #geom_line(data = Nuresult$plotdf,aes(gapLen,foldchn,color="Nuctrl"))+
  # scale_colour_manual(name="",values=c(CAP="orange",TFctrl="steelblue", Nuctrl="grey"),label=c("nucleosomal DNA","free_DNA","noTF")) +
  scale_colour_manual(name="",values=c(CAP="orange",TFctrl="steelblue"),label=c("nucleosomal DNA","free_DNA")) + xlab("Gap length of kmers (bp)") + ylab(expression("log"[2]*"(fold-change) of top kmers"))+
  gg_theme_Publication()+
 scale_x_continuous(limits = c(0,50)) + scale_y_continuous(limits = c(0.2,1.1))

print(gap_plot)
gg_save_all(gap_plot,width = 8,height = 5,newName = "enrichment_gaplen_plot")



