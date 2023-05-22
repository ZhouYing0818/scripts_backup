#DEseq2
library(DESeq2)
HepG2sh1_counts<-counts_mat[,c(1:7)]
rownames(HepG2sh1_counts)<-HepG2sh1_counts$Gene
HepG2sh1_counts<-HepG2sh1_counts[,-1]
condition <- factor(c(rep("HepG2EV",3), rep("HepG2sh1",3)))
coldata <- data.frame(row.names = colnames(HepG2sh1_counts), condition)
dds<-DESeqDataSetFromMatrix(countData=HepG2sh1_counts, colData=coldata, design=~condition)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
dds1 <- DESeq(dds)
res <- results(dds1)
res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds1, normalized=TRUE)),by="row.names",sort=FALSE)

resdata=resdata[is.na(resdata$padj)==FALSE,]
resdata[,"ID"]<-resdata$Row.names
resdata<-resdata[,-1]
resdata<-merge(resdata,anno2[,c(10,11)],by.x="ID",by.y="gene_id",all.x=T)
resdata[,"transcript_id_new"]<-gsub("\\.\\d*","",resdata$transcript_id)
resdata<-resdata[!duplicated(resdata$ID),]
#rawid<-rownames(resdata)
#rawid<-stringr::str_split(rawid,"\\.",simplify = T)[,1]
my_result<-getBM(attributes = c("ensembl_transcript_id",
                                "transcript_biotype",
                                "external_gene_name",
                                "description"),
                 filters ="ensembl_transcript_id",
                 values=resdata$transcript_id_new,
                 mart = my_dataset)
my_result<-my_result[!duplicated(my_result$ensembl_transcript_id),]
resdata<-merge(resdata,my_result,by.x="transcript_id_new",by.y="ensembl_transcript_id",all.x=T)

Genediff<-resdata[resdata$transcript_biotype=="protein_coding",]
Genediff<-na.omit(Genediff)
write.csv(Genediff, "HepG2sh1_DESeq2.csv",row.names = F)

install.packages("VennDiagram")
library(VennDiagram)
DESeq2_DEG<-Genediff[Genediff$padj<0.2,]$external_gene_name
edgeR_DEG<-annodata[annodata$FDR<0.2,]$external_gene_name
p<-venn.diagram(list(DESeq2_0.2=DESeq2_DEG,edgeR_0.2=edgeR_DEG),
                col=c('red','green'),fill=c('red','green'),cat.col=c('red','green'),
             filename=NULL)
grid.draw(p)

HepG2_MHCC_TMM<-calcNormFactors(HepG2_MHCC, method = 'TMM')
f<-HepG2_MHCC_TMM$samples$norm.factors
lib_size<-HepG2_MHCC_TMM$samples$lib.size
HepG2_MHCC_TMM_mat<-t(t(HepG2_MHCC_TMM$counts) / (lib_size * f)) * 1e6
