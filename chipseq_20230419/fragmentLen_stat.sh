#! /bin/bash/
rootdir=/alldata/zhouying/CUT_and_TAG/test_demo
index=/alldata/zhouying/ref_data/genome/Ecoli/index/bowtie2/E.coli_K12_MG1655
outdir=${rootdir}/04_RemoveDup
indir=${rootdir}/04_RemoveDup
cat ${rootdir}/01_RawData/download_list.txt|while read id;
do
	echo ${id}' running...'
	samtools view -F 0x04 ${indir}/${id}_hg38_sort.rmdup.bam|awk -F '\t' 'function abs(x){return ((x < 0.0) ? -x : x)}{print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}'>${outdir}/${id}.fragmentLen.txt
done
