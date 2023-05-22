#! /bin/bash/
rootdir=/alldata/zhouying/RNAseq/Fan_20230412/X101SC22084125-Z01-J009/04_align_star
mkdir -p ${rootdir}/bam_stat
ls ${rootdir}|grep _sortAligned.sortedByCoord.out.bam|grep -v "bai"|while read id;
do
	echo ${id}' running...'
	#echo ${id%_sortA*}
	bam_stat.py -i ${rootdir}/${id}>${rootdir}/bam_stat/${id%_sortA*}.bamstat.txt
done
