#! /bin/bash/
rootdir=/alldata/zhouying/Ribo-seq/ribo-TISH_test_data/05_align_star
ls ${rootdir}|grep _RiboAligned.sortedByCoord.out.bam|while read id;
do
	echo ${id}' running...'
	samtools index ${rootdir}/${id}
done
