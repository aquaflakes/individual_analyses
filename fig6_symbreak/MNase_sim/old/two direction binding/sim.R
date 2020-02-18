library(Rcpp)
library(ggplot2)
library(dplyr)
sourceCpp("../Desktop/simCpp.cpp")

TFpos=0; TFsize=15; TFblock=500
Nusize=147; Nublockf=50; Nublockc=1000 # bloking factor for Nu, flank and center
# TFpenetration=20 # half decay point
linkerLen=10:30

TF_Nu_spacing_fw= -200:200
TF_Nu_spacing_rev=-200:200
# TF_Nu_spacing=0:50

seqNum=2000L; seqLenHalf=600

freeDNA_cut_prob=0.1; cut_times=1

sd_nu= 10
side_bias=0.8


# Nu win for multiply
Nu_exp=log(Nublockf/Nublockc)/73
d=-73:73 #dist to dyad
y=Nublockc*exp(Nu_exp*abs(d)) # block factor
Nu_win=1/y

# TF win for multiply
TF_win= rep(1/TFblock,TFsize)


# gen seq
seq=matrix(1,seqNum,seqLenHalf*2+1); colnames(seq)= -seqLenHalf:seqLenHalf

# gen Nu positions
Nupos_p1= ((TFsize-1)/2+(Nusize-1)/2+mean(TF_Nu_spacing_fw)+1) %>% {rnorm(seqNum,mean(.),sd_nu)} %>% as.integer()
Nupos_p2= Nupos_p1 + (rnorm(seqNum,mean(linkerLen),sd_nu) %>% as.integer()) + Nusize +1
  Nupos_m1_Nup1= Nupos_p1 - (rnorm(seqNum,mean(linkerLen),sd_nu) %>% as.integer())  -Nusize -1
  Nupos_m1_TF=((TFsize-1)/2+(Nusize-1)/2+TF_Nu_spacing_rev+1) %>% {rnorm(seqNum,mean(.),sd_nu)} %>% as.integer() %>% -.
Nupos_m1=  pmin(Nupos_m1_Nup1,Nupos_m1_TF)
Nupos_m2= Nupos_m1 - (rnorm(seqNum,mean(linkerLen),sd_nu) %>% as.integer()) -Nusize -1
  Nupos_m1[sample(c(T,F),length(Nupos_m1),replace = T,prob = c(side_bias,1-side_bias))]=NA
  Nupos_m2[sample(c(T,F),length(Nupos_m2),replace = T,prob = c(side_bias,1-side_bias))]=NA
  hist(c(Nupos_p1,Nupos_p2,Nupos_m1,Nupos_m2),breaks = 1000)
  
  # # gen seq
  # seq=matrix(1,seqNum,seqLenHalf*2+1); colnames(seq)= -seqLenHalf:seqLenHalf
  # 
  # # gen Nu positions
  # # Nupos_p1= ((TFsize-1)/2+(Nusize-1)/2+mean(TF_Nu_spacing_fw)+1) %>% rnorm(.,mean(.),sd_nu) %>% as.integer()
  # Nupos_p1= ((TFsize-1)/2+(Nusize-1)/2+TF_Nu_spacing_fw+1) %>% sample(seqNum,replace=T)
  # Nupos_p2= Nupos_p1 + sample(linkerLen, seqNum,replace=T) + Nusize +1
  # Nupos_m1_Nup1= Nupos_p1 - sample(linkerLen, seqNum,replace=T) -Nusize -1
  # Nupos_m1_TF=((TFsize-1)/2+(Nusize-1)/2+TF_Nu_spacing_rev+1) %>% sample(seqNum,replace=T) %>% -.
  # Nupos_m1=  pmin(Nupos_m1_Nup1,Nupos_m1_TF)
  # Nupos_m2= Nupos_m1 - sample(linkerLen, seqNum,replace=T) -Nusize -1
  # Nupos_m1[sample(c(T,F),length(Nupos_m1),replace = T,prob = c(side_bias,1-side_bias))]=NA
  # Nupos_m2[sample(c(T,F),length(Nupos_m2),replace = T,prob = c(side_bias,1-side_bias))]=NA
  # hist(c(Nupos_p1,Nupos_p2,Nupos_m1,Nupos_m2),breaks = 1000)

# # gen Nu positions
# Nupos_m1= ((TFsize-1)/2+(Nusize-1)/2+TF_Nu_spacing_rev+1) %>% sample(seqNum,replace=T) %>% -.
# Nupos_m2= Nupos_m1 - sample(linkerLen, seqNum,replace=T) -Nusize -1
#   Nupos_p1_Num1= Nupos_m1 + sample(linkerLen, seqNum,replace=T) +Nusize +1
# #   Nupos_p1_TF=((TFsize-1)/2+(Nusize-1)/2+TF_Nu_spacing_fw+1) %>% sample(seqNum,replace=T)
# # Nupos_p1=  pmax(Nupos_p1_Num1,Nupos_p1_TF)
#   Nupos_p1=  Nupos_p1_Num1
# Nupos_p2= Nupos_p1 + sample(linkerLen, seqNum,replace=T) + Nusize +1
# hist(c(Nupos_p1,Nupos_p2,Nupos_m1,Nupos_m2),breaks = 1000)



# gen Nu positions
# Nupos_p0= (0:100) %>% sample(seqNum,replace=T)
# seq=place(seq,Nupos_p0,Nu_win)
# 
# Nupos_p1= ((TFsize-1)/2+(Nusize-1)/2+TF_Nu_spacing+1) %>% sample(seqNum,replace=T)
#   Nupos_p1= pmax(Nupos_p1,Nupos_p0+ sample(linkerLen, seqNum,replace=T) + Nusize +1 )
# Nupos_p2= Nupos_p1 + sample(linkerLen, seqNum,replace=T) + Nusize +1
# Nupos_m1= ((TFsize-1)/2+(Nusize-1)/2+TF_Nu_spacing+1) %>% sample(seqNum,replace=T) %>% -.
#   Nupos_m1= pmin(Nupos_m1,Nupos_p0- sample(linkerLen, seqNum,replace=T) -Nusize -1 )
# Nupos_m2= Nupos_m1 - sample(linkerLen, seqNum,replace=T) -Nusize -1
#   hist(c(Nupos_p1,Nupos_p2,Nupos_m1,Nupos_m2,Nupos_p0),breaks = 1000)


# add digestion barrier
seq=place(seq,rep(0,seqNum),TF_win)
seq=place(seq,Nupos_p1,Nu_win)
seq=place(seq,Nupos_p2,Nu_win)
seq=place(seq,Nupos_m1,Nu_win)
seq=place(seq,Nupos_m2,Nu_win)


MNaseCuts=cutseq(seq,freeDNA_cut_prob,cut_times)
  # View(MNaseCuts[,350:450])
  # colSums(MNaseCuts) %>% plot
fragments= frag_mid_and_pos(MNaseCuts)
p=ggplot(fragments,aes(fragMid,fragLen)) + geom_point(alpha=0.3)+ theme_bw() + scale_x_continuous(expand = c(0,0))+ scale_y_continuous(expand = c(0,0))
print(p+scale_x_continuous(expand = c(0,0),limits = c(-200,200))+scale_y_continuous(expand = c(0,0),limits = c(0,200)))





