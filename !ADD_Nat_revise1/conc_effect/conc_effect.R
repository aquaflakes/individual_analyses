fjComm::clear_()

lig200=fjComm::nusel_get_200() %>% dplyr::filter(batch=="FJ45")

double_conc="B18 I14 L8 M6 B22 L13 D14 P13" %>% str_split(" ") %>% unlist()
normal_conc="A23 M2 L23 D18 A19 G8 B20 L11" %>% str_split(" ") %>% unlist()

double_conc_info= lig200 %>% dplyr::filter(pos %in% double_conc) %>% arrange(TF) %>% as_tibble() %>%  mutate(EMI=purrr::map(c4_EMI_file,read_csv), conc="double")
normal_conc_info= lig200 %>% dplyr::filter(pos %in% normal_conc) %>% arrange(TF) %>% as_tibble() %>%  mutate(EMI=purrr::map(c4_EMI_file,read_csv),conc="normal")

# # cp files
# targetDir="Reads/lig200/"
# dir.create(targetDir)
# file.copy(c(double_conc_info$c4_file %>% str_replace("200c4-([^-]*)-.*$","200c[1-4]-\\1-\\*") %>% Sys.glob()),targetDir)


data=rbind(double_conc_info,normal_conc_info) %>% dplyr::select(TF,EMI,conc)
data=unnest(data) %>% dplyr::filter(pos2-pos1==3)
data1= data %>% group_by(TF,conc) %>% mutate(topMIsum= squishQuan_normalize(topMIsum))  %>% ungroup() %>% mutate(conc=factor(conc,levels = c("normal","double"),labels = c("1X","2X")))

p=ggplot(data1)+ geom_tile(aes(pos1,TF, fill=topMIsum), height=0.5) + facet_wrap(~conc) + gg_theme_Publication() + scale_fill_gradientn(colours = gg_steelblue_red,guide = F) +
  theme(axis.title.x = element_blank(),axis.line.x = element_blank(),  axis.title.y = element_blank() , strip.text.x = element_text(size = 10))+scale_x_continuous(breaks = c(1,149))
print(p)

# gg_save_pdf(p, 10, 8, filename = "conc_effect")
gg_save_pdf(p, 8.4, 6.7, filename = "conc_effect")
