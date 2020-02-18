fjComm::clear_()

# library(DNAshapeR)
# c0= "../1_reads_collection/cyc0.fa"
# ncap="../1_reads_collection/NCAP.fa"
# ht="../1_reads_collection/HT.fa"
# 
# shapeType=c("MGW","HelT","ProT","Roll","EP","Stretch", "Tilt", "Buckle", "Shear","Opening", "Rise", "Shift", "Stagger", "Slide")
# ht <- getShape(ht,parse = T,shapeType = shapeType)#shapeType=c("MGW")# shapeType = c("Stretch", "Tilt", "Buckle", "Shear","Opening", "Rise", "Shift", "Stagger", "Slide"), parse = T)
# c0 <- getShape(c0,parse = T,shapeType = shapeType)
# ncap <- getShape(ncap,parse = T,shapeType = shapeType)
# 
# plotShape(ht$MGW)
# plotShape(c0$MGW)
# plotShape(ncap$MGW)
# 
# 
# collectAll <-function(feature="MGW")
# {
#   MGW=ht[[feature]] %>% colMeans(na.rm = T) %>% {tibble(type="HT",feature=feature,pos=1:length(.),meanVal=.)}
#   MGW=ncap[[feature]]%>% colMeans(na.rm = T) %>% {tibble(type="NCAP",feature=feature,pos=1:length(.),meanVal=.)} %>% rbind(MGW)
#   MGW=c0[[feature]]%>% colMeans(na.rm = T) %>% {tibble(type="c0",feature=feature,pos=1:length(.),meanVal=.)} %>% rbind(MGW)
#   MGW %<>% mutate(pos=pos-97) #%>% dplyr::filter(pos<0 | pos>7)
#   p=ggplot(MGW)+geom_line(aes(pos,meanVal,color=type))+scale_x_continuous(expand = c(0,0),limits = c(-50,50))+geom_vline(xintercept = c(0,7),alpha=0.3)
#   # print(plotly::ggplotly(p))
#   print(p)
#   # ggplot(MGW)+geom_smooth(aes(pos,meanVal,color=type),span = 0.15)
#   return(MGW)
# }
# 
# plotdata=do.call(rbind,lapply(shapeType,collectAll))
# saveRDS(plotdata,"plotdata.Robj")
plotdata=readRDS("plotdata.Robj")

# MGW=collectAll("MGW")
# HelT=collectAll("HelT")
# ProT=collectAll("ProT")
# Roll=collectAll("Roll")
# EP=collectAll("EP")
# 
# plotdata=rbind(MGW,HelT,ProT,Roll,EP)
p=ggplot(plotdata)+geom_line(aes(pos,meanVal,color=type),size=0.3)+scale_x_continuous(expand = c(0,0),limits = c(-50,50),breaks = c(-50,0,50))+geom_vline(xintercept = c(0,7),alpha=0.35,linetype=8,color="black")+
  facet_wrap(~feature,2,7,scales = "free_y")+scale_color_manual(values = c("grey","#5997d1","black"),name="",guide=F) + xlab("(bp)")+ylab("Mean value")
print(p)
# stop()

gg_save_pdf(p,20,7,filename = "allfeatures")
# tt=readRDS("../1_reads_collection/cyc0.Robj")

limitsCalc<- function(x){
  browser()
  limits <- c(min(x),max(x))/2
  limits
}


