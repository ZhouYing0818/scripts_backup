#读取gtf文件
library(rtracklayer)
gtf<-import("E:/data/anno/gencode.v36.annotation.gtf")%>%as.data.frame()
anno<-gtf[gtf$type=="gene",]
anno<-anno[anno$gene_type=="protein_coding",c(10,12)]
counts_mat<-merge(anno,counts_mat,by.x="gene_id",by.y="Gene")
counts_mat<-counts_mat[which(rowSums(counts_mat[,-c(1:2)])>3),]
dup_name<-counts_mat[duplicated(counts_mat$gene_name),]$gene_name
for (i in dup_name){
  tmp<-counts_mat[counts_mat$gene_name==i,]
  counts_mat<-counts_mat[counts_mat$gene_name!=i,]
  tmp1<-as.data.frame(c(tmp[1,1:2],colSums(tmp[,-c(1:2)])))
  counts_mat<-rbind(counts_mat,tmp1)
}





my_mart<-useMart("ensembl")
db<-listDatasets(my_mart)
my_dataset<-useDataset("hsapiens_gene_ensembl",mart = my_mart)
test<-getBM(attributes = c('ensembl_transcript_id',
                           "transcript_length",
                           "transcript_biotype",
                           "external_gene_name",
                           "description"),mart = my_dataset)
anno_mat1<-merge(counts_mat2,gtf[,c(10:11)],by="gene_id",all.x=T)
anno_mat1<-anno_mat1[!duplicated(anno_mat1$gene_id),]
anno_mat1[,"transcript_id_new"]<-gsub("\\.\\d*","",anno_mat1$transcript_id)
anno_mat1<-merge(anno_mat1,test,by.x="transcript_id_new",by.y="ensembl_transcript_id",all.x=T)
write.csv(na.omit(anno_mat1[anno_mat1$transcript_biotype=="protein_coding",]),"Kasumi_NB4_RawCounts.csv",row.names = F)

lengths<-anno_mat1$transcript_length
names(lengths)<-anno_mat1$transcript_id_new
nm_count<-anno_mat1[,c(3:38)]
rownames(nm_count)<-anno_mat1$transcript_id_new
samples<-colnames(nm_count)
colnames(nm_count)<-paste(samples,"RawCounts",sep = "_")
total_count<-colSums(anno_mat1[,c(3:38)])
nm_fpkm <- t(do.call( rbind,
                      lapply(1:length(total_count),
                             function(i){
                               10^9*nm_count[,i]/lengths/total_count[i]
                               #lengths向量自动遍历
                             }) ))
colnames(nm_fpkm)<-paste(samples,"FPKM",sep = "_")
nm_tpm <- log2((t(t(nm_fpkm)/colSums(nm_fpkm,na.rm = T))*10^6)+1)
colnames(nm_tpm)<-paste(samples,"log2(TPM+1)",sep = "_")


nm_fpkm<-as.data.frame(nm_fpkm)
nm_fpkm[,"transcript_id_new"]<-rownames(nm_fpkm)
tmp<-merge(anno_mat1[,c(1:2,39:43)],nm_fpkm,by="transcript_id_new",all.x=T)
tmp<-na.omit(tmp[tmp$transcript_biotype=="protein_coding",])
write.csv(tmp,"Kasumi_NB4_FPKM.csv",row.names = F)

nm_tpm<-as.data.frame(nm_tpm)
nm_tpm[,"transcript_id_new"]<-rownames(nm_tpm)
tmp<-merge(anno_mat1[,c(1:2,39:43)],nm_tpm,by="transcript_id_new",all.x=T)
tmp<-na.omit(tmp[tmp$transcript_biotype=="protein_coding",])
write.csv(tmp,"Kasumi_NB4_log2(TPM+1).csv",row.names = F)

HepG2_MHCC_TMM_mat<-as.data.frame(HepG2_MHCC_TMM_mat)
HepG2_MHCC_TMM_mat[,"gene_id"]<-rownames(HepG2_MHCC_TMM_mat)
HepG2_MHCC_TMM_mat<-merge(HepG2_MHCC_TMM_mat,gtf[,c(10:11)],by="gene_id",all.x=T)
HepG2_MHCC_TMM_mat<-HepG2_MHCC_TMM_mat[!duplicated(HepG2_MHCC_TMM_mat$gene_id),]
HepG2_MHCC_TMM_mat[,"transcript_id_new"]<-gsub("\\.\\d*","",HepG2_MHCC_TMM_mat$transcript_id)
tmp<-merge(anno_mat1[,c(1:2,39:43)],HepG2_MHCC_TMM_mat,by="transcript_id_new",all.x=T)
tmp<-na.omit(tmp[tmp$transcript_biotype=="protein_coding",])
write.csv(tmp,"HepG2_MHCC_TMM.csv",row.names = F)
