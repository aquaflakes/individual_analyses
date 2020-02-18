fjComm::clear_()
# pacman::p_load(ShortRead)

ncapfile="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/FJ4.4_PE_147/dual_trim/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-4-CAP2-C05-TF15-SOX11IIIc4_S437_L002_R1_001.peared_trimmed.fq.gz"
htfile="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/FJ4.4_PE_TFctrl/Trulig147v1IIIFJ4-4-TFctrl-C05-TF15-SOX11IIIc3_S53_L004_R1_001.peared_trimmed.fq.gz"
c0file="~/Nut_zhuData/seqFiles2/c0/subset_1000000_Trulig147v1IIIPEx200IIIc0_S3_R1_001.peared_trimmed.fq.gz"

motifhits<- function(ncapfile,padlenTotal=200)
{
  ncap=SELEXFile(ncapfile)
  ncap$getSeq(dup_rm = T)
  ncap$moodsMap(pwmFile = "../sox11_motif.pfm")
  moods=ncap$moodsResult %>% dplyr::filter(pos>=36 & pos<=58) %>% mutate(seq=ncap$seq[rangeNo] %>% as.character())
      motiflen=moods$match[1] %>% nchar
      liglen=moods$seq[1] %>% nchar
        upper_index=liglen-motiflen
    moods %<>% mutate(pos=ifelse(strand=="+",pos,upper_index-pos)) %>% mutate(seq=ifelse(strand=="+",seq,revComp(seq)))
    
    moods %<>% mutate(X7=str_sub(seq,pos+1,pos+motiflen))
      
      Npad=padlenTotal-liglen
      motif_left=moods$pos; motif_right=93-moods$pos; lrdiff=motif_left-motif_right
      Npad_base=(Npad-abs(lrdiff))/2; Npad_add=Npad-Npad_base*2
      Npad_left=ifelse(lrdiff<=0,Npad_base+Npad_add,Npad_base); Npad_right=Npad-Npad_left
    moods %<>% mutate(seq=paste0(strrep("N",Npad_left),seq,strrep("N",Npad_right))) 
    moods %<>% mutate(posN=pos+Npad_left,match=str_sub(seq,posN+1,posN+motiflen)) 
  return(moods %>% select(pos,score,match,seq,posN))
}

saveRDS(motifhits(ncapfile),"NCAP.Robj")
saveRDS(motifhits(htfile),"HT.Robj")
saveRDS(motifhits(c0file),"cyc0.Robj")

readRDS("HT.Robj") %>% .$seq %>% DNAStringSet() %>% writeXStringSet('HT.fa')
readRDS("NCAP.Robj") %>% .$seq %>% DNAStringSet() %>% writeXStringSet('NCAP.fa')
readRDS("cyc0.Robj") %>% .$seq %>% DNAStringSet() %>% writeXStringSet('cyc0.fa')
