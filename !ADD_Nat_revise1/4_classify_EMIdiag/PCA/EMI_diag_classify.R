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

p_load(FactoMineR)

# num=5:10
#
# for (i in num)
# {
#   pca = PCA(NCAPdiags, graph = F,ncp = i)
#   pcs=pca$var$coord
#   p=ggheat(pcs %>% as.matrix() %>% t,rescaling = "row")+scale_fill_gradientn(colours = gg_steelblue_red)
#   gg_save_pdf(p,width = 10, height = i*2,filename = "PCA_n_" %>% paste0(i))
# }


pca = PCA(NCAPdiags, graph = F,ncp = 5)


# p_load(Factoshiny)
  # PCAshiny(Mydata)
# resshiny = PCAshiny(pca)

# pcs=pca$var$coord
pcs=pca$var$cor
 # pcs %>% t %>% writeClipboard()

cor_dim1=apply(NCAPdiags,1,function(x)cor(pcs[,1],x))
cor_dim2=apply(NCAPdiags,1,function(x)cor(pcs[,2],x))
cor_dim3=apply(NCAPdiags,1,function(x)cor(pcs[,3],x))
cor_dim4=apply(NCAPdiags,1,function(x)cor(pcs[,4],x))
cor_dim5=apply(NCAPdiags,1,function(x)cor(pcs[,5],x))

cor_dim=tibble(TF=names(cor_dim1),cor_dim1,cor_dim2,cor_dim3,cor_dim4,cor_dim5)


cor_dim %<>% mutate(dim1s2=cor_dim1-cor_dim2) %>% mutate(binder=fjComm::nusel_binders()[TF])

saveRDS(cor_dim,"PCA_EMI_cor_dim")


PCA_result=ggheat(pcs %>% as.matrix() %>% t,rescaling = "row")+scale_fill_gradientn(colours = gg_steelblue_red,guide=F)+scale_x_continuous(expand = c(0,0),breaks = c(1,96))+ylab("PCA components")+xlab(("(bp)"))
gg_save_pdf(PCA_result,6,5)




# ind=pca$ind$coord %>% as.data.frame();
# tt=ind[c("CREB5","MLX"),]
# # ind$TF=rownames(ind); ind$binder=fjComm::nusel_binders()[ind$TF]; ind %<>% mutate(sub=Dim.1-Dim.2)
# ggheat(pcs %>% as.matrix() %>% t,rescaling = "row")+scale_fill_gradientn(colours = gg_steelblue_red)
#
# EOMES= tibble(pos=1:96,TF="_CREB5",value=NCAPdiags["CREB5",] %>% squishQuan_normalize())
# ETV4= tibble(pos=1:96,TF="_MLX",value=NCAPdiags["MLX",] %>% squishQuan_normalize())
# pca3= pcs[,1:2] %>% apply(2,squishQuan_normalize) %>% melt() %>% setnames(c("pos","TF","value"))
# plotdata=rbind(EOMES,ETV4,pca3)
# ggplot(plotdata)+geom_tile(aes(pos, TF, fill=value))+scale_fill_gradientn(colours = gg_steelblue_red)

