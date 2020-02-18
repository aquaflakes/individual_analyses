fjComm::clear_()

allTFs=readRDS("allTFinfo_S5.RData")

TF147=fjComm::nusel_get_147() %>% dplyr::filter(use==1& order==1) %>% dplyr::slice(match(allTFs$V1 %>% as.character(),TF))

allTFs$Tag=TF147$tag
allTFs %>% dplyr::filter(Tag=="SBP")
