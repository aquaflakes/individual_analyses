fjComm::clear_()

df1= "ELF2_Vplot_bin_counts.tsv" %>% read_tsv()
plotdata=df1[,3:4] %>% melt()

testObj=t.test(df1$Count_left,df1$Count_right)
testObj$p.value





# cvg_plot=ggplot()+geom_line(aes(-200:200,cvg))+ scale_x_continuous(limits = c(-200,200),expand = c(0,0),breaks = c(-200,0,200))+ scale_y_continuous(breaks = c(12000,14000))+xlab("Distance from motif (bp)")+ylab("Coverage")
# gg_save_pdf(cvg_plot,6,3,filename = "ELF2_cut_coverage")

p= ggplot(plotdata,aes(variable,value))+geom_point(size=0.8,alpha=0.99)+ stat_summary(color="red",alpha=0.5) + gg_axis_x_labels(c("3' end", "5' end")) +
  ylab("Count")+ scale_y_continuous(limits = c(20,150), breaks = c(50, 100)) + xlab("")+
  geom_path(data=data.frame(x=c(1,1,2,2),y=c(120,140,140,130)),aes(x,y))+
  annotate("text",x=1.5,y=150,size=9/3,label= paste0("p=",signif(testObj$p.value,2)))+ gg_theme_Publication()+ gg_axis_x_label_angle()

gg_save_pdf(p,4,6,filename = "ELF1_t_test")
