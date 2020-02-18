fjComm::clear_()

### manually copy table S5
# df=readClipBoard(colNames = T)
# saveRDS(df,"rawdata")
df=readRDS("rawdata")
df200only=df[196:nrow(df),]
df=df[1:195,]
classify_df=df %>% select(EMI.penetration..lig147.,Periodicity.strength,Dyad.EMI.intensity,Gyre.spanning.mode.intensity,Orientational.asymmetry,Log2.EMIub.EMIb.)
TFnames=df$TF
classify_df=classify_df %>% as.matrix() %>% gtools::na.replace(0) %>% set_rownames(TFnames)

# normalize col
classify_df[,1:5]=apply(classify_df[,1:5],2,squishQuan_normalize)
# classify_df[,6]=apply(classify_df[,6],2,squishQuan_normalize)
lig147=fjComm::nusel_get_147() %>% dplyr::filter(use==1 & order==1 & !is.na(c5_EMI_Nu_file)& tag!="SBP")
# tt=classify_df[,6]
classify_df[,6][!(rownames(classify_df) %in% lig147$TF )]=0


# ylabel=df[hclust(dist(classify_df))$order, ]$TF %>% as.character()
p=ggheat(classify_df %>% as.matrix(),clustering = "row",labRow = T)+scale_x_continuous(expand = c(0,0),limits = c(0.5,6.5),breaks = 1:6,labels = qw("E P D Gs Asym Ds"))+
 scale_fill_gradientn(colours = c("green", "white", "red"), guide=F, limits= c(-1,1))+theme(axis.text.y =element_text(size=4))
gg_save_pdf(p,10,28,filename = "all_class")

# masAsym=df$Orientational.asymmetry %>% is.nan%>%  max#  %>% quantile(0.995)



classify_df=df200only %>% select(Orientational.asymmetry)
TFnames=df200only$TF %>% as.character()
  ordering=order(classify_df$Orientational.asymmetry)
classify_df=classify_df %>% dplyr::slice(ordering) %>% as.matrix() %>% gtools::na.replace(0) 
  rownames(classify_df)=TFnames[ordering]
# normalize col
classify_df=apply(classify_df,2,squishQuan_normalize) %>% cbind(0)

# ylabel=df[hclust(dist(classify_df))$order, ]$TF %>% as.character()
p=ggheat(classify_df %>% as.matrix(),labRow = T)+scale_x_continuous(expand = c(0,0),limits = c(0.5,6.5),breaks = 1:6,labels = qw("Asym E P D Gs Diss"))+
  scale_fill_gradientn(colours = c("white", "red"),guide=F)+theme(axis.text.y =element_text(size=4))
gg_save_pdf(p,10,28/nrow(df)*nrow(df200only)+3,filename = "all_class_200only")
