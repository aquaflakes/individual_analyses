fjComm::clear_()
pacman::p_load(ShortRead)

seqs=read_csv("~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/FJ4.4_PE_147/dual_trim/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-4-CAP2-C05-TF15-SOX11IIIc4_S437_L002_R1_001.peared_trimmed.fq.gz",col_names = F)
seqs=rmdup(seqs) %>% head(20000)
seqs[[1]] %>% DNAStringSet() %>% writeXStringSet('SOX11c4.fa')
# install_github("TsuPeiChiu/DNAshapeR", ref = "ext-features",force=T)

library(DNAshapeR)
fn <- "SOX11c4.fa"
# pred <- getShape(fn)
pred <- getShape(fn,parse = T)#shapeType=c("MGW")# shapeType = c("Stretch", "Tilt", "Buckle", "Shear","Opening", "Rise", "Shift", "Stagger", "Slide"), parse = T)
plotShape(pred$MGW)
plotShape(pred$HelT)
plotShape(pred$ProT)
plotShape(pred$Roll)
plotShape(pred$EP)




seqs=read_csv("~/Nut_zhuData/seqFiles2/FJ6.1Methyl_PE/allreads/FJ6.1_OriMe1_PE/dual_trim/2_Trim_adaptor_kmer/Trulig147v1IIIFJ6-1-OriMe1-B02-1xCAGOct32ngIIIc4-Bq4-ATCACTGA_S782_L008_R1_001.peared_trimmed.fq.gz",col_names = F)
seqs=rmdup(seqs) %>% head(20000)
seqs[[1]] %>% DNAStringSet() %>% writeXStringSet('CAGc4.fa')
# install_github("TsuPeiChiu/DNAshapeR", ref = "ext-features",force=T)

library(DNAshapeR)
fn <- "CAGc4.fa"
# pred <- getShape(fn)
pred <- getShape(fn,parse = T)#shapeType=c("MGW")# shapeType = c("Stretch", "Tilt", "Buckle", "Shear","Opening", "Rise", "Shift", "Stagger", "Slide"), parse = T)
plotShape(pred$MGW)
plotShape(pred$HelT)
plotShape(pred$ProT)
plotShape(pred$Roll)
plotShape(pred$EP)


MGW=readShape("CAGc4.fa.MGW") %>% colMeans(na.rm = T) %>% {tibble(type="Nu",feature="MGW",pos=1:length(.),meanVal=.)}
MGW=readShape("SOX11c4.fa.MGW") %>% colMeans(na.rm = T) %>% {tibble(type="SOX11",feature="MGW",pos=1:length(.),meanVal=.)} %>% rbind(MGW)
MGW %<>% group_by(type) %>% mutate(meanVal=scale(meanVal)) %>% ungroup()
# ggplot(MGW)+geom_line(aes(pos,meanVal,color=type))+geom_smooth(aes(pos,meanVal,color=type))
ggplot(MGW)+geom_smooth(aes(pos,meanVal,color=type),span = 0.15)

Roll=readShape("CAGc4.fa.Roll") %>% colMeans(na.rm = T) %>% {tibble(type="Nu",feature="Roll",pos=1:length(.),meanVal=.)}
Roll=readShape("SOX11c4.fa.Roll") %>% colMeans(na.rm = T) %>% {tibble(type="SOX11",feature="Roll",pos=1:length(.),meanVal=.)} %>% rbind(Roll)
Roll %<>% group_by(type) %>% mutate(meanVal=scale(meanVal)) %>% ungroup()
# ggplot(Roll)+geom_line(aes(pos,meanVal,color=type))+geom_smooth(aes(pos,meanVal,color=type))
ggplot(Roll)+geom_smooth(aes(pos,meanVal,color=type),span = 0.15)





