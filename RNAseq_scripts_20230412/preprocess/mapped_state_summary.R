# extracting alignstat information
align_dir<-"/alldata/zhouying/RNAseq/Fan_20230412/X101SC22084125-Z01-J009/04_align_star/"
bamstat_dir<-"/alldata/zhouying/RNAseq/Fan_20230412/X101SC22084125-Z01-J009/04_align_star/bam_stat/"
stattable<-data.frame()
maptable<-data.frame()
file_list<-dir(bamstat_dir)
for (bamstat_file in file_list){
  id<-strsplit(bamstat_file,split = ".",fixed = T)[[1]][1]
  file_path<-paste(bamstat_dir,bamstat_file,sep = "")
  tmp<-readLines(file_path)
  tmp<-strsplit(tmp,split = ":",fixed = T)
  tmp1<-c()
  tmp1_names<-c()
  for (i in c(6,8,9,11:12,14:22)){
    tmp1<-c(tmp1,as.numeric(tmp[[i]][2]))
    tmp1_names<-c(tmp1_names,tmp[[i]][1])
  }
  names(tmp1)<-tmp1_names
  tmp1<-t(as.data.frame(tmp1))
  rownames(tmp1)<-id
  stattable<-rbind(stattable,tmp1)
  file_path<-paste(align_dir,id,"_sortLog.final.out",sep = "")
  tmp<-readLines(file_path)
  tmp<-strsplit(tmp,split = "|",fixed = T)
  tmp2<-c()
  tmp2_names<-c()
  for (i in 1:length(tmp)){
    tmp2<-c(tmp2,tmp[[i]][2])
    tmp2_names<-c(tmp2_names,tmp[[i]][1])
  }
  names(tmp2)<-tmp2_names
  tmp2<-t(as.data.frame(tmp2))
  rownames(tmp2)<-id
  maptable<-rbind(maptable,tmp2)
}

mapstat_results<-data.frame(Sample=rownames(stattable),
           Valid_reads=as.numeric(maptable[,6])*2,
           Mapped_reads=paste(stattable[,1],"(",round(as.numeric(stattable[,1])/(as.numeric(maptable[,6])*2)*100,2),"%)",sep = ""),
           Unique_Mapped_reads=paste(as.numeric(stattable[,6]),"(",round(as.numeric(stattable[,6])/(as.numeric(maptable[,6])*2)*100,2),"%)",sep = ""),
           Multi_Mapped_reads=paste(as.numeric(stattable[,5]),"(",round(as.numeric(stattable[,5])/(as.numeric(maptable[,6])*2)*100,2),"%)",sep = ""),
           PE_Mapped_reads=paste(as.numeric(stattable[,13]),"(",round(as.numeric(stattable[,13])/(as.numeric(maptable[,6])*2)*100,2),"%)",sep = ""),
           Read_map_to_sense_strand=paste(as.numeric(stattable[,9]),"(",round(as.numeric(stattable[,9])/(as.numeric(maptable[,6])*2)*100,2),"%)",sep = ""),
           Read_map_to_antisense_strand=paste(as.numeric(stattable[,10]),"(",round(as.numeric(stattable[,10])/(as.numeric(maptable[,6])*2)*100,2),"%)",sep = ""),
           Non_splice_reads=paste(as.numeric(stattable[,11]),"(",round(as.numeric(stattable[,11])/(as.numeric(maptable[,6])*2)*100,2),"%)",sep = ""),
           Splice_reads=paste(as.numeric(stattable[,12]),"(",round(as.numeric(stattable[,12])/(as.numeric(maptable[,6])*2)*100,2),"%)",sep = "")
           )
write.csv(mapstat_results,paste(align_dir,"sumamry_mapped_stat.csv",sep=""),row.names = F)
