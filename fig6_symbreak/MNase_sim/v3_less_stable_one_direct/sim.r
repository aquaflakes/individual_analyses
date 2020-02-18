fjComm::clear_()
Rcpp::sourceCpp("simCpp.cpp")

TFpos=0; TFsize=15; TFblockf=50; TFblockc=100; TF_half_size= as.integer(TFsize/2)
Nusize=147; Nublockf=50; Nublockc=100;  Nu_half_size= as.integer(Nusize/2) # bloking factor for Nu, flank and center
# TFpenetration=20 # half decay point
# linkerLen=10:30

# TF_Nu_spacing_fw= -200:200
# TF_Nu_spacing_rev=-200:200
# TF_Nu_spacing=0:50

seqNum=2000L; seqLenHalf=600

freeDNA_cut_prob=0.02; cut_times=5

sd_nu= 10
side_bias=0.8


# Nu win for multiply
Nu_exp=log(Nublockf/Nublockc)/Nu_half_size
d=-Nu_half_size:Nu_half_size #dist to dyad
y=Nublockc*exp(Nu_exp*abs(d)) # block factor
Nu_win=1/y

# TF win for multiply
TF_exp=log(TFblockf/TFblockc)/TF_half_size
d=-TF_half_size:TF_half_size #dist to dyad
y=TFblockc*exp(TF_exp*abs(d)) # block factor
TF_win=1/y


# gen seq (cut prob at each pos)
seq=matrix(1,seqNum,seqLenHalf*2+1); colnames(seq)= -seqLenHalf:seqLenHalf

# gen Nu positions
Nupos_p1=
Nupos_p1= ((TFsize-1)/2+(Nusize-1)/2+TF_Nu_spacing_fw+1) %>% sample(seqNum,replace=T)
# Nupos_p1= ((TFsize-1)/2+(Nusize-1)/2+mean(TF_Nu_spacing_fw)+1) %>% {rnorm(seqNum,mean(.),sd_nu)} %>% as.integer()
Nupos_p2= Nupos_p1 + (rnorm(seqNum,mean(linkerLen),sd_nu) %>% as.integer()) + Nusize +1
  Nupos_m1_Nup1= Nupos_p1 - (rnorm(seqNum,mean(linkerLen),sd_nu) %>% as.integer())  -Nusize -1
Nupos_m1_TF=((TFsize-1)/2+(Nusize-1)/2+TF_Nu_spacing_rev+1) %>% {rnorm(seqNum,mean(.),sd_nu)} %>% as.integer() %>% -.
Nupos_m1=  pmin(Nupos_m1_Nup1,Nupos_m1_TF)
  # Nupos_m1=Nupos_m1_Nup1
Nupos_m2= Nupos_m1 - (rnorm(seqNum,mean(linkerLen),sd_nu) %>% as.integer()) -Nusize -1
  Nupos_m1[sample(c(T,F),length(Nupos_m1),replace = T,prob = c(side_bias,1-side_bias))]=NA
  Nupos_m2[sample(c(T,F),length(Nupos_m2),replace = T,prob = c(side_bias,1-side_bias))]=NA
  hist(c(Nupos_p1,Nupos_p2,Nupos_m1,Nupos_m2),breaks = 1000)

  # # gen seq
  # seq=matrix(1,seqNum,seqLenHalf*2+1); colnames(seq)= -seqLenHalf:seqLenHalf
  #
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
  overlap_pos_upper=(TFsize-1)/2+(Nusize-1)/2+1;  overlap_pos_upper= overlap_pos_upper+10;
  overlap_pos_lower=50
  filters=(Nupos_m1>overlap_pos_lower & Nupos_m1<overlap_pos_upper) | (Nupos_p1>overlap_pos_lower & Nupos_p1<overlap_pos_upper)| (Nupos_p2>overlap_pos_lower & Nupos_p2<overlap_pos_upper)

  TFpos= rep(0,seqNum); TFpos[filters]=NA
seq=place(seq,TFpos,TF_win)
seq=place(seq,Nupos_p1,Nu_win)
seq=place(seq,Nupos_p2,Nu_win)
seq=place(seq,Nupos_m1,Nu_win)
seq=place(seq,Nupos_m2,Nu_win)


MNaseCuts=cutseq(seq,freeDNA_cut_prob,cut_times)
fragments= frag_mid_and_pos(MNaseCuts)
MNaseCuts=cutseq(seq,freeDNA_cut_prob,10)
fragments= rbind(fragments, frag_mid_and_pos(MNaseCuts))
MNaseCuts=cutseq(seq,freeDNA_cut_prob,2)
fragments= rbind(fragments, frag_mid_and_pos(MNaseCuts))

# p=ggplot(fragments,aes(fragMid,fragLen)) + stat_bin_2d (binwidth = 1)+ theme_bw()+scale_fill_gradientn(colours = c("white","black"),limits=c(0,10)) + scale_x_continuous(expand = c(0,0))+ scale_y_continuous(expand = c(0,0))
p=ggplot(fragments,aes(fragMid,fragLen)) + geom_point(alpha=0.3)+ theme_bw() + scale_x_continuous(expand = c(0,0))+ scale_y_continuous(expand = c(0,0))
p=p+scale_x_continuous(expand = c(0,0),limits = c(-200,200))+scale_y_continuous(expand = c(0,0),limits = c(0,200))
print(p)

gg_save_png(p,width = 10,height = 10,filename = "sim_result")



# 2 direction opening
TFpos= rep(0,seqNum);
seq=place(seq,TFpos,TF_win)
MNaseCuts=cutseq(seq,freeDNA_cut_prob,cut_times)
fragments= frag_mid_and_pos(MNaseCuts)
MNaseCuts=cutseq(seq,freeDNA_cut_prob,10)
fragments= rbind(fragments, frag_mid_and_pos(MNaseCuts))
MNaseCuts=cutseq(seq,freeDNA_cut_prob,2)
fragments= rbind(fragments, frag_mid_and_pos(MNaseCuts))

# p=ggplot(fragments,aes(fragMid,fragLen)) + stat_bin_2d (binwidth = 1)+ theme_bw()+scale_fill_gradientn(colours = c("white","black"),limits=c(0,10)) + scale_x_continuous(expand = c(0,0))+ scale_y_continuous(expand = c(0,0))
p=ggplot(fragments,aes(fragMid,fragLen)) + geom_point(alpha=0.3)+ theme_bw() + scale_x_continuous(expand = c(0,0))+ scale_y_continuous(expand = c(0,0))
p=p+scale_x_continuous(expand = c(0,0),limits = c(-200,200))+scale_y_continuous(expand = c(0,0),limits = c(0,200))
print(p)

gg_save_png(p,width = 10,height = 10,filename = "sim_result_2direction")


