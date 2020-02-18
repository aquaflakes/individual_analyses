fjComm::clear_()

ori=readxl::read_excel("PITX_contacts_arr.xlsx",1,col_names = F)
mut=readxl::read_excel("PITX_contacts_arr.xlsx",2,col_names = F)

ori=ori %>% arrange(X__1,X__2)
oriout=tibble(AA=paste0(ori[[1]],"-",ori[[2]] %>% str_replace_all("[A-Z]","")), Base=ori[[5]] %>% str_replace("D","") %>% toupper(),Distance=ori[[9]])
oriout %<>% group_by(AA) %>% summarise(Base= paste0(Base," (",Distance,")") %>% paste0(collapse=", "))
# oriout %>% writeClipboard(sep = "\t\t")

mut=mut %>% arrange(X__1,X__2)
mutout=tibble(AA=paste0(mut[[1]],"-",mut[[2]] %>% str_replace_all("[A-Z]","")), Base=mut[[5]] %>% str_replace("D","") %>% toupper(),Distance=mut[[9]])
mutout %<>% group_by(AA) %>% summarise(Base= paste0(Base," (",Distance,")") %>% paste0(collapse=", "))
# mutout %>% writeClipboard(sep = "\t\t")

output=dplyr::full_join(oriout,mutout,"AA") %>% set_names(c("AA", "Contacts (free DNA)","Contacts (nucl. DNA)"))
output %>% writeClipboard("\t")
