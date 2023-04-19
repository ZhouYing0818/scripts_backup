#! /bin/bash/
rootdir=/alldata/zhouying/CUT\&RUN/GSE136332_SALL4
index=/alldata/zhouying/ref_data/genome/hg38/index/bwa/hg38.fa
outdir=${rootdir}/05_align
indir=${rootdir}/03_clean_fastq
ls ${rootdir}/01.RawData|grep SRR|while read id;
do
	echo ${id}' running...'
	read1=${indir}/${id}_1_clean.fq.gz
	read2=${indir}/${id}_2_clean.fq.gz
	#echo $read1
	#echo $read2

	bwa aln -t 8 ${index} ${read1} > ${outdir}/${id}_1.sai
	bwa aln -t 8 ${index} ${read2} > ${outdir}/${id}_2.sai
	bwa sampe ${index} ${outdir}/${id}_1.sai ${outdir}/${id}_2.sai ${read1} ${read2}|samtools view -Sb - > ${outdir}/${id}_hg38.bam

done
