#RNAseq DEG analysis
rootpath="/data/RNAseq/Fan_20230412/04_align_star/05_DGE/List/"
setwd(rootpath)
raw_counts_list<-dir()

#combined all counts results
counts_matrix=data.frame()
for (i in 1:length(raw_counts_list)){
  tmp<-read.table(raw_counts_list[i],sep='\t',row.names = 1)
  tmp[,"Gene"]<-rownames(tmp)
  tmp<-tmp[,-c(2:3)]
  tmp<-tmp[-c(1:4),]
  file_name<-raw_counts_list[i]
  s<-strsplit(file_name,split = "_",fixed = T)[[1]][1:3]
  sample_name<-paste(s[1],s[2],s[3],sep = "_")
  colnames(tmp)<-c(sample_name,"Gene")
  if (dim(counts_matrix)[1]==0){
    counts_matrix<-tmp
  }else{
    counts_matrix<-merge(counts_matrix,tmp,by="Gene",all.x = T)
  }
}
write.csv(counts_matrix,"raw_counts.csv",row.names = F)
#split count_mat
counts_mat<-read.csv("raw_counts.csv")

library(rtracklayer)
library(dplyr)
gtf<-import("/data/ref_data/anno/gencode.v36.annotation.gtf")%>%as.data.frame()
anno<-gtf[gtf$type=="gene",]
anno<-anno[anno$gene_type=="protein_coding",c(10,12,4)]
counts_mat<-merge(anno,counts_mat,by.x="gene_id",by.y="Gene")
counts_mat<-counts_mat[which(rowSums(counts_mat[,-c(1:2)])>3),]
dup_name<-counts_mat[duplicated(counts_mat$gene_name),]$gene_name
for (i in dup_name){
  tmp<-counts_mat[counts_mat$gene_name==i,]
  counts_mat<-counts_mat[counts_mat$gene_name!=i,]
  tmp1<-as.data.frame(c(tmp[1,1:2],colSums(tmp[,-c(1:2)])))
  counts_mat<-rbind(counts_mat,tmp1)
}

counts_mat1<-counts_mat
counts_mat1<-counts_mat[,c(1:3,4:27)]
counts_mat1<-counts_mat1[which(rowSums(counts_mat1[,-c(1:3)])>3),]
counts_mat2<-counts_mat[,c(1:3,28:51)]
counts_mat2<-counts_mat2[which(rowSums(counts_mat2[,-c(1:3)])>3),]

rownames(counts_mat1)<-counts_mat1$gene_name
counts_mat1<-counts_mat1[,-c(1:3)]
counts_mat1<-counts_mat1[,-c(1:3,23,24,27,16:18)]
rownames(counts_mat2)<-counts_mat2$gene_name
counts_mat2<-counts_mat2[,-c(1:3,16,19,20,22:24)]
#genering DElist data &normalization
library(edgeR)
group1<-c(rep("KASUMI_1LACZ_1",3),rep("KASUMI_1sg1",3),rep("KKASUMI_1sg4",3),
          rep("Kasumi_1_PWsh4_dox",3),rep("Kasumi_1_shEV",3),rep("Kasumi_1_shEV_dox",3))
Kasumi <- DGEList(counts=counts_mat1,group=group1)
group2<-c(rep("NB4_PWsh1",3),rep("NB4_PWsh1_dox",3),
          rep("NB4_PWsh4",3),rep("NB4_PWsh4_dox",3),
          rep("NB4_shEV",3),rep("NB4_shEV_dox",3))
NB4 <- DGEList(counts=counts_mat2,group=group2)

filtNorm<-function(dge){
keep <- filterByExpr(dge)
dge <- dge[keep,,keep.lib.sizes=FALSE] 
dge <- calcNormFactors(dge)
return(dge)
}

kasumi_norm<-filtNorm(Kasumi)
NB4_norm<-filtNorm(NB4)
#sample heatmap
library(pheatmap)
library("RColorBrewer")
sampleHeatmap<-function(logCPM){
  #sampledist<-dist(t(logCPM))
  samplecor <- cor(logCPM)
  samplecorMatrix <- as.matrix(samplecor)
  colnames(samplecorMatrix) <- NULL
  #colors <- colorRampPalette( brewer.pal(9, "Blues") )(255)
  p<-pheatmap(samplecorMatrix,
              display_numbers = T,silent = T)
  return(p)
}

#HepG2_logCPM<-cpm(HepG2_MHCC_norm[,1:9],log=TRUE,prior.count=2)
#MHCC_logCPM<-cpm(HepG2_MHCC_norm[,10:18],log=TRUE,prior.count=2)
Kasumi_logCPM<-cpm(kasumi_norm,log=TRUE,prior.count=2)
NB4_logCPM<-cpm(NB4_norm,log=TRUE,prior.count=2)

fpkm_cal<-function(count_mat,anno){
  tmp<-data.frame(gene_name=rownames(count_mat))
  tmp<-merge(tmp,anno,by="gene_name")
  tmp<-tmp[!duplicated(tmp$gene_name),]
  lengths<-tmp$width
  names(lengths)<-tmp$gene_name
  samples<-colnames(count_mat)
  total_count<-colSums(count_mat)
  nm_fpkm <- t(do.call( rbind,
                        lapply(1:length(total_count),
                               function(i){
                                 10^9*count_mat[,i]/lengths/total_count[i]
                                 #lengths向量自动遍历
                               }) ))
  colnames(nm_fpkm)<-paste(samples,"FPKM",sep = "_")
  return(nm_fpkm)
}

kasumi_fpkm<-fpkm_cal(Kasumi$counts,anno)
NB4_fpkm<-fpkm_cal(NB4$counts,anno)

#sampleHeatmap(HepG2_MHCC_fpkm[,1:9])
#p1<-sampleHeatmap(HepG2_MHCC_fpkm[,10:18])
pdf("Kasumi_samp_heatmap.pdf",width = 10,height = 8)
sampleHeatmap(kasumi_fpkm)
dev.off()

pdf("NB4_samp_heatmap.pdf",width = 10,height = 8)
sampleHeatmap(NB4_fpkm)
dev.off()

#PCA

sample_pca<-function(logCPM,group){
pca <- prcomp(t(logCPM),center = TRUE,scale. = TRUE)
df <- pca$x
df <- as.data.frame(df)
summ <- summary(pca)
xlab <- paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab <- paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
library(ggplot2)
ggplot(data = df,aes(x = PC1,y = PC2,color=group))+
  geom_point(size=4)+
  labs(x = xlab,y = ylab,color = "Samples")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
}
#sample_pca(HepG2_logCPM,group1[1:9])
#p2<-sample_pca(MHCC_logCPM,group1[10:18])
pdf("Kasumi_samp_PCA.pdf",width = 8,height = 6)
sample_pca(Kasumi_logCPM,group1)
dev.off()

pdf("NB4_samp_PCA.pdf",width = 8,height = 6)
sample_pca(NB4_logCPM,group2)
dev.off()

#boxplot overview
library(tidyr)
library(reshape2)
library(ggplot2)
sample_boxplot<-function(logCPM){
logcpm_df<-as.data.frame(logCPM)
logcpm_long<-gather(logcpm_df,S,value)
ggplot(data=logcpm_long,aes(x=S,y=log2(value+1),fill=factor(S)))+
  geom_boxplot()+
  labs(x="Samples", y = "Log2FPKM")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
}

#sample_boxplot(HepG2_MHCC_fpkm[,1:9])
#p3<-sample_boxplot(HepG2_MHCC_fpkm[,10:18])
pdf("Kasumi_samp_boxplot.pdf",width = 10,height = 8)
sample_boxplot(kasumi_fpkm)
dev.off()

pdf("NB4_samp_boxplot.pdf",width = 10,height = 8)
sample_boxplot(NB4_fpkm)
dev.off()

#--------------------edgeR analysis---------------------------
design1 <- model.matrix(~0+group1)
my_contrasts1 <- makeContrasts(Kasumi_1sg1VSLACZ=group1KASUMI_1sg1-group1KASUMI_1LACZ_1, 
                               KKasumi_1sg4VSLACZ=group1KKASUMI_1sg4-group1KASUMI_1LACZ_1,
                               levels=design1)


my_contrasts1 <- makeContrasts(Kasumi_1_PWsh1=group1Kasumi_1_PWsh1_dox-group1Kasumi_1_PWsh1, 
                               Kasumi_1_PWsh4=group1Kasumi_1_PWsh4_dox-group1Kasumi_1_PWsh4,
                               Kasumi_1_shEV=group1Kasumi_1_shEV_dox-group1Kasumi_1_shEV,
                               KasumiPWsh1VSEV_dox=group1Kasumi_1_PWsh1_dox-group1Kasumi_1_shEV_dox,
                               KasumiPWsh4VSEV_dox=group1Kasumi_1_PWsh4_dox-group1Kasumi_1_shEV_dox,
                              levels=design1)
design2 <- model.matrix(~0+group2)
my.contrasts2 <- makeContrasts(NB4_PWsh1=group2NB4_PWsh1_dox-group2NB4_PWsh1,
                               NB4_PWsh4=group2NB4_PWsh4_dox-group2NB4_PWsh4,
                               NB4_shEV=group2NB4_shEV_dox-group2NB4_shEV,
                               NB4PWsh1VSEV_dox=group2NB4_PWsh1_dox-group2NB4_shEV_dox,
                               NB4PWsh4VSEV_dox=group2NB4_PWsh4_dox-group2NB4_shEV_dox,
                               levels=design2)

library(biomaRt)
library(curl)

my_mart<-useMart("ensembl")
db<-listDatasets(my_mart)
my_dataset<-useDataset("hsapiens_gene_ensembl",mart = my_mart)

DGEcal<-function(dge,design,my.contrasts,anno){
dge <- estimateDisp(dge,design)
fit <- glmQLFit(dge, design)
for (i in 1:dim(my.contrasts)[2]) {
  #tabname<-pasete("qlf.",colnames(my.contrasts)[i])
  #str1<-paste(tabname," <- glmQLFTest(fit,contrast=my.contrasts[,i])",sep = "")
  #eval(parse(text=str1))
  test_results<-glmQLFTest(fit,contrast=my.contrasts[,i])
  ordered_tags <- topTags(test_results, n=100000)
  allDiff=ordered_tags$table
  allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
  allDiff[,"gene_name"]<-rownames(allDiff)
  #allDiff<-merge(allDiff,anno2[,c(10,11)],by.x="ID",by.y="gene_id",all.x=T)
  #allDiff[,"transcript_id_new"]<-gsub("\\.\\d*","",allDiff$transcript_id)
  #allDiff<-allDiff[!duplicated(allDiff$ID),]
  #rawid<-rownames(allDiff)
  #rawid<-stringr::str_split(rawid,"\\.",simplify = T)[,1]
  my_result<-getBM(attributes = c("external_gene_name",
                                  "description"),
                   filters ="external_gene_name",
                   values=allDiff$gene_name,
                   mart = my_dataset)
  my_result<-my_result[!duplicated(my_result$external_gene_name),]
  allDiff<-merge(allDiff,my_result,by.x="gene_name",by.y="external_gene_name",all.x=T)
  allDiff<-allDiff[order(allDiff$FDR),]
  #Genediff<-allDiff[allDiff$transcript_biotype=="protein_coding",]
  #Genediff<-na.omit(Genediff)
  #otherdiif<-allDiff[allDiff$transcript_biotype!="protein_coding",]
  #otherdiif<-na.omit(otherdiif)
  #savename1<-paste(colnames(my.contrasts)[i],"_gene.csv")
  #savename2<-paste(colnames(my.contrasts)[i],"_others.csv")
  nm_fpkm<-fpkm_cal(dge$counts,anno)%>%as.data.frame()
  nm_fpkm[,"gene_name"]<-rownames(nm_fpkm)
  allDiff<-merge(allDiff,nm_fpkm,by="gene_name")
  savename<-paste(colnames(my.contrasts)[i],"_edgeR_DEG.csv",sep = "")
  write.csv(allDiff, savename,row.names = F)
  #write.csv(Genediff, savename1,row.names = F)
  #write.csv(otherdiif, savename2,row.names = F)
  }
}
DGEcal(kasumi_norm,design1,my_contrasts1,anno)
DGEcal(NB4_norm,design2,my.contrasts2,anno)

write.csv(Kasumi$counts,"Kasumi_RawCounts.csv")
write.csv(NB4$counts,"NB4_RawCounts.csv")
write.csv(kasumi_fpkm,"Kasumi_FPKM.csv")
write.csv(NB4_fpkm,"NB4_FPKM.csv")
#--------------------DGE results visualization----------------------------
#火山图
library(ggrepel)

valcano_data<-function(annodata,alldata){
  alldata<-merge(alldata,annodata[,c(2,10)],by="ID",all.x=T)
  alldata$labels <- ifelse(abs(alldata$logFC) > 1 & alldata$FDR < 0.05, "Both", 
                         ifelse(alldata$FDR < 0.05, "Significant",
                                ifelse(abs(alldata$logFC) > 1, "log2FC > 1", "None")))
  return(alldata)
}

valcano_plot<-function(valcano_data){
  # draw the basic structure of volcano plot
  p <- ggplot(valcano_data, aes(logFC, -1*log10(FDR)))+
    geom_point(aes(color=labels),
                      #alpha=0.4,
                      size=1.75)+
    scale_color_manual(values =c("red", "blue", "grey","orange"))+
    labs(x=expression(log[2]("Fold Change")),
         y=expression(-log[10]("p adjusted value")))+
    geom_vline(xintercept=c(-1,1), lty=4, col="black", lwd=0.6)+
    geom_hline(yintercept = -log10(0.05), lty=4, col="black", lwd=0.6)+
    geom_text_repel(data =na.omit(subset(valcano_data, labels == "Both")),
              aes(label =external_gene_name),
              size=3,force = 1,
              segment.color = "#cccccc", segment.size = 0.5)+
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black")) +
    theme(
    # Change the color, size and the face of main title, x and y axis labels
    plot.title = element_text(size=15, face="bold", hjust = 0.5), #centered the title
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    # set legend 
    legend.position="right",
    legend.title = element_blank())
  return(p)
}

annodata<-read.csv("Kasumi_1_PWsh1_edgeR_DEG.csv")
#alldata<-read.csv("HepGsh1vsEV _raw.csv")
volcano_data<-valcano_data(annodata)
p<-valcano_plot(volcano_data)
p

library(Glimma)
annodata<-read.csv("Kasumi_1sg4VSLACZ_edgeR_DEG.csv")
interact_volcano_plot<-function(dge,design,my.contrast,annodata,plotname){
  dge <- estimateDisp(dge,design)
  fit <- glmQLFit(dge,design)
  re<- glmQLFTest(fit,contrast=my.contrast)
  rownames(annodata)<-annodata$gene_name
  annodata<-annodata[,c(2:7)]
  colnames(annodata)<-c(colnames(annodata)[1:5],"Description")
  #volcano_data<-volcano_data[order(volcano_data$FDR),]
  re$table<-annodata
  glimmaVolcano(re,anno = annodata[,5:6],html = plotname,main = strsplit(plotname,split = ".",fixed = T)[[1]][1])
}

interact_volcano_plot(kasumi_norm,design1,my_contrasts1[,2],annodata,
                      paste(colnames(my_contrasts1)[2],".html",sep = ""))


rownames(volcano_data)<-volcano_data$ID
test_results$table<-volcano_data
glimmaVolcano(test_results,anno = volcano_data[,6:8],html = "hepG2sh1.html")


#热图
reg_pattern_hm<-function(norm_mat){
  gene_matrix<-as.data.frame(norm_mat)
  gene_matrix[,"ID"]<-rownames(norm_mat)
  gene_matrix<-merge(annodata[,c(2,3,7:8,10)],gene_matrix,by="ID",all.x=T)
  gene_matrix<-gene_matrix[gene_matrix$FDR<=0.05,]
  gene_matrix<-gene_matrix[abs(gene_matrix$logFC)>=1,]
  gene_matrix[gene_matrix$logFC<=0,"labels"]<-"down"
  gene_matrix[gene_matrix$logFC>0,"labels"]<-"up"
  gene_matrix<-gene_matrix[order(gene_matrix$FDR),]
  rownames(gene_matrix)<-gene_matrix$ID
  annotation_row<-as.data.frame(factor(gene_matrix$labels))
  rownames(annotation_row)<-rownames(gene_matrix)
  colnames(annotation_row)<-"Regulation"
  #upreg<-gene_matrix[gene_matrix$labels=="up",c(6:14)]
  p<-pheatmap(as.matrix(gene_matrix[,6:14]),
              scale="row",
              border_color=NA,
              cutree_rows = 2,
              treeheight_col = 2,
              treeheight_row = 2,
              show_rownames=F,
              annotation_row = annotation_row,
              annotation_names_row=F,
              angle_col=45)
}
p<-reg_pattern_hm(HepG2_logCPM)



for (i in 1:20){
  i=20
  name<-paste(strsplit(GSEA_BP@result$ID[i],split = ":",fixed = T)[[1]][2],"pdf",sep = ".")
  pdf(name,width = 15,height = 10)
  gseaplot2(GSEA_BP, geneSetID = i, title = GSEA_CC$Description[i])
  dev.off()
}

library(clusterProfiler)
library(org.Hs.eg.db)

go_ORA_dotplot<-function(genelist,savename){
  ORA_BP<- enrichGO(gene  = names(geneList),
                    #universe      = names(geneList),
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'SYMBOL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.2,
                    readable      = TRUE)
  ORA_MF <- enrichGO(gene  = names(geneList),
                     #universe      = names(geneList),
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'SYMBOL',
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.2,
                     readable      = TRUE)
  ORA_CC <- enrichGO(gene  = names(geneList),
                     #universe      = names(geneList),
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'SYMBOL',
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.2,
                     readable      = TRUE)
  dotplot_BP<-dotplot(ORA_BP,showCategory=20)
  dotplot_MF<-dotplot(ORA_MF,showCategory=20)
  dotplot_CC<-dotplot(ORA_CC,showCategory=20)
  
  ORA_BP<-as.data.frame(ORA_BP)
  ORA_MF<-as.data.frame(ORA_MF)
  ORA_CC<-as.data.frame(ORA_CC)
  
  ORA_BP[,"ontology"]<-rep("BP",dim(ORA_BP)[1])
  ORA_MF[,"ontology"]<-rep("MF",dim(ORA_MF)[1])
  ORA_CC[,"ontology"]<-rep("CC",dim(ORA_CC)[1])
  ORA_ALL<-rbind(ORA_BP,ORA_MF,ORA_CC)
  write.csv(ORA_ALL,paste(rootpath,savename,"_go_ORA.csv",sep = ""),row.names = F)
  
  ggsave(paste(newlist[i],"_go_ORA.pdf",sep=""),plot_grid(dotplot_BP,
            dotplot_MF,
            dotplot_CC, 
            labels = c("BP", "MF","CC"),nrow = 1),height = 10,width = 20)
  
}

go_GSEA_dotplot<-function(genelist,savename){
  GSEA_BP<- gseGO(geneList     = genelist,
                  OrgDb        = org.Hs.eg.db,
                  keyType       = 'SYMBOL',
                  ont          = "BP",
                  minGSSize    = 10,
                  maxGSSize    = 500,
                  pvalueCutoff = 1,
                  verbose      = FALSE)
  
  GSEA_MF <- gseGO(geneList     = genelist,
                   OrgDb        = org.Hs.eg.db,
                   keyType       = 'SYMBOL',
                   ont          = "MF",
                   minGSSize    = 10,
                   maxGSSize    = 500,
                   pvalueCutoff = 1,
                   verbose      = FALSE)
  
  GSEA_CC <- gseGO(geneList     = genelist,
                   OrgDb        = org.Hs.eg.db,
                   keyType       = 'SYMBOL',
                   ont          = "CC",
                   minGSSize    = 10,
                   maxGSSize    = 500,
                   pvalueCutoff = 1,
                   verbose      = FALSE)
  
  dotplot_BP<-dotplot(GSEA_BP,showCategory=20)
  dotplot_MF<-dotplot(GSEA_MF,showCategory=20)
  dotplot_CC<-dotplot(GSEA_CC,showCategory=20)
  
  GSEA_BP<-as.data.frame(GSEA_BP)
  GSEA_MF<-as.data.frame(GSEA_MF)
  GSEA_CC<-as.data.frame(GSEA_CC)
  
  GSEA_BP[,"ontology"]<-rep("BP",dim(GSEA_BP)[1])
  GSEA_MF[,"ontology"]<-rep("MF",dim(GSEA_MF)[1])
  GSEA_CC[,"ontology"]<-rep("CC",dim(GSEA_CC)[1])
  GSEA_ALL<-rbind(GSEA_BP,GSEA_MF,GSEA_CC)
  write.csv(GSEA_ALL,paste(rootpath,savename,"_go_GSEA.csv",sep = ""),row.names = F)
  
  ggsave(paste(newlist[i],"_go_GSEA.pdf",sep=""),
         plot_grid(dotplot_BP,
                   dotplot_MF,
                   dotplot_CC, 
                   labels = c("BP", "MF","CC"),nrow = 1),
         height = 10,width = 20)
}


#------------------------------GO --------------------------------------------

newlist<-c()
for (i in 1:2){
  file_list<-strsplit(dir(),split = "_",fixed = T)
  newlist<-c(newlist,paste(file_list[[i]][1],
                           file_list[[i]][2],
                           #file_list[[i]][3],
                           sep = "_"))
}
newlist<-unique(newlist)

for (i in 1:length(newlist)){
  annodata<-read.csv(paste(newlist[i],"_edgeR_DEG.csv",sep = ""))
  geneList<-as.vector(annodata[,2])
  names(geneList)<-as.vector(annodata[,1])
  geneList<-sort(geneList,decreasing = T)
  
  #pdf(paste(newlist[i],"_go_ORA.pdf",sep=""),height = 10,width = 20)
  go_ORA_dotplot(geneList,newlist[i])
  go_GSEA_dotplot(geneList,newlist[i])
  #dev.off()
}

#---------------------KEGG-----------------------------------------------------

KEGG_GSEA_dotplot<-function(geneList,savename){
  kk_gsea<- gseKEGG(geneList= geneList,
                  organism     = 'hsa',
                  minGSSize    = 10,
                  pvalueCutoff = 1,
                  verbose      = FALSE)
  kk_gsea<-setReadable(kk_gsea,OrgDb = org.Hs.eg.db, keyType = "ENTREZID" )
  ggsave(paste(rootpath,savename,"_KEGG_GSEA.pdf",sep = ""),dotplot(kk_gsea,showCategory=20),height = 10,width = 10)
  write.csv(as.data.frame(kk_gsea),paste(rootpath,savename,"_KEGG_GSEA.csv",sep = ""),row.names = F)
}

KEGG_ORA_dotplot<-function(geneList,savename){
  kk_ORA<- enrichKEGG(gene = gene.df[,2],
                    organism     = 'hsa',
                    pvalueCutoff = 0.05)
  kk_ORA<-setReadable(kk_ORA,OrgDb = org.Hs.eg.db, keyType = "ENTREZID" )
  ggsave(paste(rootpath,savename,"_KEGG_ORA.pdf",sep = ""),dotplot(kk_ORA,showCategory=20),height = 10,width = 10)
  write.csv(as.data.frame(kk_ORA),paste(rootpath,savename,"_KEGG_ORA.csv",sep = ""),row.names = F)
}

newlist<-c()
for (i in 1:10){
  file_list<-strsplit(dir(),split = "_",fixed = T)
  newlist<-c(newlist,paste(file_list[[i]][1],file_list[[i]][2],file_list[[i]][3],sep = "_"))
}
newlist<-unique(newlist)

for (i in 1:length(newlist)){
  annodata<-read.csv(paste(newlist[i],"_edgeR_DEG.csv",sep = ""))
  geneList<-as.vector(annodata[,2])
  names(geneList)<-as.vector(annodata[,1])
  geneList<-sort(geneList,decreasing = T)
  gene.df <- bitr(names(geneList), fromType ="SYMBOL" ,
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)
  names(geneList)<-gene.df[,2]
  
  #pdf(paste(newlist[i],"_go_ORA.pdf",sep=""),height = 10,width = 20)
  KEGG_ORA_dotplot(geneList,newlist[i])
  KEGG_GSEA_dotplot(geneList,newlist[i])
  #dev.off()
}










#for (i in 5:20){
#pathview(gene.data  = geneList,
                   #pathway.id = kk_gsea@result$ID[i],
                   #species    = "hsa",
#gene.idtype = "SYMBOL",
#kegg.native=F,
         #kegg.dir="./Kasami_NB4_KEGG_GSEA",
         #imit      = list(gene=max(abs(geneList)), cpd=1))
#}



for (i in 1:20){
pathview(gene.data  = geneList,
                   pathway.id = kk_ORA@result$ID[i],
                   species    = "hsa",
         kegg.dir="./Kasami_NB4_KEGG_ORA",
         imit      = list(gene=max(abs(geneList)), cpd=1))
}

