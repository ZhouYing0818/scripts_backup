#! /bin/bash/
rootdir=/alldata/zhouying/CUT_and_TAG/Jing_20230224
#index=/alldata/zhouying/ref_data/genome/Ecoli/index/bowtie2/E.coli_K12_MG1655
outdir=${rootdir}/06_calling-peaks
indir=${rootdir}/05_RemoveDup
ls ${rootdir}/01_RawData|while read id;
do
	echo ${id}' running...'
	samtools view ${indir}/${id}_hg38_sort.rmdup.bam|awk -F '\t' '$9<150{print $0}'|samtools view -@8>${outdir}/${id}_filted_len_150.bam
done
