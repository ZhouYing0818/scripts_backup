#! /bin/bash/
rootdir=/alldata/zhouying/RNAseq/Fan_20230412/X101SC22084125-Z01-J009/04_align_star
ls ${rootdir}|grep _sortAligned.sortedByCoord.out.bam|while read id;
do
	echo ${id}' running...'
	samtools index ${rootdir}/${id}
done
