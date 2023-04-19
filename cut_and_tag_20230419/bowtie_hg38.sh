#! /bin/bash/
rootdir=/alldata/zhouying/CUT_and_TAG/Jing_20230224
index=/alldata/zhouying/ref_data/genome/hg38/index/bowtie2/hg38
outdir=${rootdir}/04_align_bowtie2
indir=${rootdir}/02_CleanData
ls ${rootdir}/01_RawData|while read id;
do
	echo ${id}' running...'
	#read1=${indir}/${id}_1.fastq.gz
	#read2=${indir}/${id}_2.fastq.gz
	read1=${indir}/${id}_1_clean.fq.gz
	read2=${indir}/${id}_2_clean.fq.gz
	#echo $read1
	#echo $read2

	#bwa aln -t 8 ${index} ${read1} > ${outdir}/${id}_1.sai
	#bwa aln -t 8 ${index} ${read2} > ${outdir}/${id}_2.sai
	#bwa sampe ${index} ${outdir}/${id}_1.sai ${outdir}/${id}_2.sai ${read1} ${read2}|samtools view -Sb - > ${outdir}/${id}_hg38.bam

	bowtie2 --end-to-end \
		--very-sensitive \
		--no-mixed \
		--no-discordant \
		--phred33 -I 10 -X 700 -p 20 \
		-x ${index} \
		-1 ${read1} -2 ${read2}|samtools sort -O bam -@ 20 -o ${outdir}/${id}_sort.bam

done
