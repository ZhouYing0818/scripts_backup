#! /bin/bash/
rootdir=/alldata/zhouying/RNAseq/Fan_20230412/X101SC22084125-Z01-J009/04_align_star
mkdir -p ${rootdir}/readsDist_stat/
ls ${rootdir}|grep _sortAligned.sortedByCoord.out.bam|grep -v "bai"|while read id;
do
	echo ${id}' running...'
	#echo ${id%_sortA*}
	read_distribution.py -i ${rootdir}/${id} -r /alldata/zhouying/ref_data/anno/GRCh38.p13/hg38_GENCODE.v38.bed>${rootdir}/readsDist_stat/${id%_sortA*}.readsDistribution.txt
done
