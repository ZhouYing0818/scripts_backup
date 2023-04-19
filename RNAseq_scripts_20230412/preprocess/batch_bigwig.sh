#! /bin/bash/
rootdir=/alldata/zhouying/RNAseq/Fan_20230110/X101SC22084125-Z01-J003
outdir=${rootdir}/06_bigwig
indir=${rootdir}/04_align_star
ls ${indir}|grep _sortAligned.sortedByCoord.out.bam|grep -v "bai"|while read id;
do
	echo ${id}' running...'
	#echo ${id%_sortA*}
	#bam_stat.py -i ${rootdir}/${id}>${rootdir}/bam_stat/${id%_sortA*}.bamstat.txt
	bamCoverage -p 20  --bam ${indir}/${id} --binSize 10 --normalizeUsing RPKM -o ${outdir}/${id%_sortA*}.bw 
done
