input_dir<-"/alldata/zhouying/RNAseq/Fan_20230201/X101SC22084125-Z01-J006/04_align_star"
out_dir<-"/alldata/zhouying/RNAseq/Dapars_test"
setwd(input_dir)
summary_data<-read.csv("sumamry_mapped_stat.csv")
PE_Mapped_reads<-summary_data[,c(1,6)]
PE_Mapped_reads_new<-data.frame()
for (i in 1:dim(PE_Mapped_reads)[1]){
  tmp_row<-PE_Mapped_reads[i,]
  tmp<-strsplit(PE_Mapped_reads[i,2],split = "(",fixed = T)[[1]][1]
  tmp_row$PE_Mapped_reads<-tmp
  tmp_row$Sample<-paste(tmp_row$Sample,"_dapars.wig",sep = "")
  PE_Mapped_reads_new<-rbind(PE_Mapped_reads_new,tmp_row)
}
write.table(PE_Mapped_reads_new,paste(out_dir,"all_samples_mapped.txt",sep = "/"),
            row.names = F,col.names = F,sep = "\t")