B<-ggplot(All,aes(x = V1, y = V2,color=group,alpha=0.5)) +
  geom_line(size = 0.5) +
  #scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 300))+theme_classic()+
  geom_vline(xintercept=150,lty=2,color="grey")

A<-ggplot(All,aes(x = sampleInfo, y = V1, weight=Weight, fill = group)) +
  geom_violin(bw = 5) +
  #scale_y_continuous(breaks = seq(0, 800, 50)) +
  #scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  #scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 20) +
  #ggpubr::rotate_x_text(angle = 20) +
  ylim(0,800)+theme_classic()+
  ylab("Fragment Length") +
  xlab("")

pdf("Fragment_dist.pdf",width = 20,height = 10)
ggarrange(A,B,ncol=2)
dev.off()

redup_fragsize<-read.table("/alldata/zhouying/CUT&RUN/GSE136332_SALL4/05_align/SRR10022375.fragmentLen.txt")
for(i in colnames(redup_fragsize)){
  redup_fragsiz[,i]<-as.numeric(redup_fragsize[,i])
}
BB2<-mutate(redup_fragsize, Weight = as.numeric(V2)/sum(as.numeric(V2)),sampleInfo=rep("SRR10022375",17942))
BB2<-BB2[-1,]
BB=rbind(BB,BB2)

sampleList = c("SRR10022374","SRR10022375","SRR10022376","SRR10022377","SRR10022378","SRR10022379")
groupList = c(rep("Igg",3),rep("SALL4",3))
All<-data.frame(matrix(ncol=6,nrow=0))
colnames(All)<-c("V1","V2","Weight","group","sampleInfo")
for(j in 1:length(sampleList)){
  #histInfo = strsplit(hist, "_")[[1]]
  AA = read.table(paste0(sampleList[j],".fragmentLen.txt"), header = FALSE) 
  for(i in colnames(AA)){
    AA[,i]<-as.numeric(AA[,i])
  }
  BB<-mutate(AA, Weight = as.numeric(V2)/sum(as.numeric(V2)), group = groupList[j], sampleInfo = sampleList[j]) 
  BB<-BB[-1,]
  
  All=rbind(All,BB)
}
