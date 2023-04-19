#! /bin/bash/
rootdir=/alldata/zhouying/chipseq/ES2_20221209
index=/alldata/zhouying/ref_data/genome/hg38/index/bwa/hg38.fa
outdir=${rootdir}/04_align
indir=${rootdir}/00_CleanData
ls ${indir}|while read id;
do
	echo ${id}' running...'
	read1=${indir}/${id}/${id}_1.clean.fq.gz
	read2=${indir}/${id}/${id}_2.clean.fq.gz
	#echo $read1
	#echo $read2

	bwa aln -t 8 ${index} ${read1} > ${outdir}/${id}_1.sai
	bwa aln -t 8 ${index} ${read2} > ${outdir}/${id}_2.sai
	bwa sampe ${index} ${outdir}/${id}_1.sai ${outdir}/${id}_2.sai ${read1} ${read2}|samtools sort -@ 20 -o ${outdir}/${id}_sort.bam
	rm ${outdir}/*.sai
done
