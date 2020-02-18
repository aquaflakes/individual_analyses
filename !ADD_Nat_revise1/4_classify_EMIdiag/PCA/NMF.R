fjComm::clear_()
# lig147= fjComm::nusel_get_147() %>% dplyr::filter(use==1&order==1)
c4EMI= oneOffCalc("oneOff/c4EMI",calcFun = function(){ lig147$c4_EMI_file %>% purrr::map(read_csv) %>% set_names(lig147$TF) })
c4EMICtrl= oneOffCalc("oneOff/c4EMICtrl",calcFun = function(){ lig147$c4_EMI_TFctrl_file %>% purrr::map(read_csv) %>% set_names(lig147$TF) })

TFnames=names(c4EMI)
c4EMI_diag= oneOffCalc("oneOff/c4EMI_d",calcFun = function(){ c4EMI %>% purrr::map(~dplyr::filter(.,pos2-pos1==3)) %>% set_names(TFnames) })
c4EMICtrl_diag= oneOffCalc("oneOff/c4EMICtrl_d",calcFun = function(){ c4EMICtrl %>% purrr::map(~dplyr::filter(.,pos2-pos1==3)) %>% set_names(TFnames) })

    NCAPdiags= vector("list",length(c4EMI_diag))
      for(i in seq_along(NCAPdiags))
      {
        NCAPdiags[[i]]= c4EMI_diag[[i]][,3]%>% set_colnames(TFnames[i])
      }
NCAPdiags= do.call(cbind,NCAPdiags)
# NCAPdiags %<>% mutate_all(function(x)squishQuan_normalize(x))
NCAPdiags %<>% t %>% set_colnames(1:ncol(.))

p_load(NMF)
set.seed(10)
nmfR=nmf(NCAPdiags, 5)

vect=nmfR@fit@H
ind=nmfR@fit@W %>% as.data.frame(); ind$TF=rownames(ind); ind$binder=fjComm::nusel_binders()[ind$TF]; ind %<>% mutate(sub=V1-V4)


#
rownames(vect)=1:5 %>% as.character()

# ggheat(vect %>% as.matrix(),rescaling = "row")+scale_fill_gradientn(colours = gg_steelblue_red)#+scale_y_continuous(breaks=0:9,labels=1:10)
NMF_result=ggheat(vect %>% as.matrix(),rescaling = "row")+scale_fill_gradientn(colours = gg_steelblue_red,guide=F)+scale_x_continuous(expand = c(0,0),breaks = c(1,96))+ylab("NMF components")+xlab(("(bp)"))
gg_save_pdf(NMF_result,5.2,5)


stop()

subset1=ind %>% arrange(desc(V5)) %>% head(20) %>% .$TF
subset2=ind %>% arrange(desc(V8)) %>% head(20) %>% .$TF
subset=c(subset1,subset2)
subsetData= NCAPdiags[subset,]
ggheat(subsetData,rescaling = "row",clustering = function(x)x)+scale_fill_gradientn(colours = gg_steelblue_red)

# EOMES= tibble(pos=1:96,TF="_CREB5",value=NCAPdiags["CREB5",] %>% squishQuan_normalize())
# ETV4= tibble(pos=1:96,TF="_MLX",value=NCAPdiags["MLX",] %>% squishQuan_normalize())
# pca3= pcs[,1:2] %>% apply(2,squishQuan_normalize) %>% melt() %>% setnames(c("pos","TF","value"))
# plotdata=rbind(EOMES,ETV4,pca3)
# ggplot(plotdata)+geom_tile(aes(pos,TF,fill=value))+scale_fill_gradientn(colours = gg_steelblue_red)

