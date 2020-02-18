fjComm::clear_()

EMIdiff <- function(x,y)
{
  allTF_MI_diags_free= readr::read_csv(x) %>% dplyr::filter(pos2-pos1==3) %>% arrange(pos1)
  allTF_MI_diags_Nu= readr::read_csv(y) %>% dplyr::filter(pos2-pos1==3) %>% arrange(pos1)
  allTF_MI_diags_free[[3]]= log2(allTF_MI_diags_free[[3]]/allTF_MI_diags_Nu[[3]])
  allTF_MI_diags_free
}

lig147=fjComm::nusel_get_147() %>% dplyr::filter(batch=="FJ45"&use==1&tag!="SBP")

double_conc="B18 I14 L8 M6 B22 L13 D14" %>% str_split(" ") %>% unlist()
normal_conc="A23 M2 L23 D18 A19 G8 B20" %>% str_split(" ") %>% unlist()

double_conc_info= lig147 %>% dplyr::filter(pos %in% double_conc) %>% arrange(TF) %>% as_tibble() %>%  mutate(EMIdiff=purrr::map2(c5_EMI_Free_file, c5_EMI_Nu_file,EMIdiff), conc="double")
normal_conc_info= lig147 %>% dplyr::filter(pos %in% normal_conc) %>% arrange(TF) %>% as_tibble() %>%  mutate(EMIdiff=purrr::map2(c5_EMI_Free_file, c5_EMI_Nu_file,EMIdiff),conc="normal")

# # cp files
# targetDir="Reads/diss/"
# dir.create(targetDir)
# file.copy(c(double_conc_info$c5Nu_file,double_conc_info$c5Free_file,double_conc_info$c4_file %>% str_replace("147c4-([^-]*)-.*$","147c[1-4]-\\1-\\*") %>% Sys.glob()),targetDir)



data=rbind(double_conc_info,normal_conc_info) %>% dplyr::select(TF,EMIdiff,conc)
data=unnest(data) %>% dplyr::filter(pos2-pos1==3)
data1= data %>% group_by(TF,conc)  %>% ungroup() %>% mutate(conc=factor(conc,levels = c("normal","double"),labels = c("1X","2X"))) %>% mutate(topMIsum=fjComm::squishQuan_normalize(topMIsum,quan=c(0.05,0.95),normalize = F))

p=ggplot(data1 )+ geom_tile(aes(pos1,TF, fill=topMIsum), height=0.5) + facet_wrap(~conc) + gg_theme_Publication() + scale_fill_gradientn(colours = c("#7dbb78", "white", "red"),guide = F,values = (c(-0.55,0,0.99)+0.55)/(0.55+0.99) ) +
  theme(axis.title.x = element_blank(),axis.line.x = element_blank(),  axis.title.y = element_blank() , strip.text.x = element_text(size = 10),strip.background = element_blank())+scale_x_continuous(breaks = c(1,96))
print(p)
gg_save_pdf(p, 10, 5, filename = "conc_effect_diss")

data1$topMIsum %>% range() %>% print()
