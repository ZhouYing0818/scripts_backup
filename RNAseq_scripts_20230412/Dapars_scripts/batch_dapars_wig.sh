#! /bin/bash/
rootdir=/alldata/zhouying/RNAseq/Fan_20230201/X101SC22084125-Z01-J006
outdir=/alldata/zhouying/RNAseq/Dapars_test
indir=${rootdir}/04_align_star
ls ${indir}|grep _sortAligned.sortedByCoord.out.bam|grep -v "bai"|while read id;
do
	echo ${id}' running...'
	#echo ${id%_sortA*}
	#bam_stat.py -i ${rootdir}/${id}>${rootdir}/bam_stat/${id%_sortA*}.bamstat.txt
	#bamCoverage -p 20  --bam ${indir}/${id} --binSize 10 --normalizeUsing RPKM -o ${outdir}/${id%_sortA*}.bw 
	bedtools genomecov -ibam ${indir}/${id} -bga -split -trackline>${outdir}/${id%_sortA*}_dapars.wig
done
