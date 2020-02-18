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
NCAPdiags %<>% mutate_all(function(x)squishQuan_normalize(x))
NCAPdiags %<>% t %>% set_colnames(1:ncol(.))

p_load(fastICA)

icaR=fastICA(NCAPdiags, 5)

# num=5:10
#
# for (i in num)
# {
#   pca = PCA(NCAPdiags, graph = F,ncp = i)
#   pcs=pca$var$coord
#   p=ggheat(pcs %>% as.matrix() %>% t,rescaling = "row")+scale_fill_gradientn(colours = gg_steelblue_red)
#   gg_save_pdf(p,width = 10, height = i*2,filename = "PCA_n_" %>% paste0(i))
# }
ggheat(icaR$K %>% as.matrix() %>% t,rescaling = "row")+scale_fill_gradientn(colours = gg_steelblue_red)

ind=icaR$S %>% as.data.frame() %>% {.$TF=rownames(.);.}
#
# pca = PCA(NCAPdiags, graph = F,ncp = 6)
#
# p_load(Factoshiny)
#   # PCAshiny(Mydata)
# resshiny = PCAshiny(pca)
#
# pcs=pca$var$coord
# pcs=pca$var$cor
#
# ind=pca$ind$coord %>% as.data.frame();  ind$TF=rownames(ind)
# ggheat(pcs %>% as.matrix() %>% t,rescaling = "row")+scale_fill_gradientn(colours = gg_steelblue_red)
#
# EOMES= tibble(pos=1:96,TF="EOMES",value=NCAPdiags["EOMES",] %>% squishQuan_normalize())
# ETV4= tibble(pos=1:96,TF="ETV1",value=NCAPdiags["TBX4",] %>% squishQuan_normalize())
# pca3= pcs[,1:4] %>% apply(2,squishQuan_normalize) %>% melt() %>% setnames(c("pos","TF","value"))
# plotdata=rbind(EOMES,ETV4,pca3)
# ggplot(plotdata)+geom_tile(aes(pos,TF,fill=value))+scale_fill_gradientn(colours = gg_steelblue_red)

