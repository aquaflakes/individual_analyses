
library(fjComm)
setwd(get_scriptpath())

files=Sys.glob("*.txt")
for(file in files)
{
  pdf(paste0(file,"_motif.pdf"),height = 4)
  fjComm::plotMotif_pfmFile(file,ic.scale = T,title = file)
  dev.off()
}


