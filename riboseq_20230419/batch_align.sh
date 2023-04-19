#! /bin/bash/
#filtering rRNA and tRNA as well as mapping and counting using STAR
rootdir=/alldata/zhouying/Ribo-seq/ribo-TISH_test_data
gtf=/alldata/zhouying/ref_data/anno/GRCh38.p13/gencode.v36.annotation.gtf
index_star=/alldata/zhouying/ref_data/genome/GRCh38.p13/index/star
index_rRNA=/alldata/zhouying/ref_data/genome/GRCh38.p13/index/bowtie2/GRCh38.p13_rRNA
index_tRNA=/alldata/zhouying/ref_data/genome/GRCh38.p13/index/bowtie2/GRCh38.p13_tRNA
outdir=${rootdir}/04_rmrRNA_rmtRNA_bowtie2
indir=${rootdir}/02_CleanData
STAR_align_output=${rootdir}/05_align_star
ls ${rootdir}/01_RawData|grep SRR|while read id;
do
	echo ${id}' running...'
	read1=${indir}/${id}_clean.fq.gz
	#read2=${indir}/${id}_2_clean.fq.gz
	#echo $read1
	#echo $read2

	#bwa aln -t 8 ${index} ${read1} > ${outdir}/${id}_1.sai
	#bwa aln -t 8 ${index} ${read2} > ${outdir}/${id}_2.sai
	#bwa sampe ${index} ${outdir}/${id}_1.sai ${outdir}/${id}_2.sai ${read1} ${read2}|samtools view -Sb - > ${outdir}/${id}_hg38.bam
	
	#echo ${id}' removing rRNA...'
	#bowtie2 -p 20 -x ${index_rRNA} \
	#	--un-gz ${outdir}/${id}.rmrRNA.fq.gz \
	#	-U ${read1}
	#echo ${id}' removing tRNA...'
	#bowtie2 -p 20 -x ${index_tRNA} \
	#	--un-gz ${outdir}/${id}.rmrRNA.rmtRNA.fq.gz \
	#	-U ${outdir}/${id}.rmrRNA.fq.gz
	#rm ${outdir}/${id}.rmrRNA.fq.gz
	echo ${id}' STAR mapping...'
	STAR --twopassMode Basic \
		--quantMode TranscriptomeSAM GeneCounts \
		--runThreadN 20 \
		--genomeDir ${index_star} \
		--alignIntronMin 20 \
		--alignIntronMax 20000 \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMattrRGline ID:${id}_riboseq SM:${id}_riboseq PL:ILLUMINA \
		--outFilterMismatchNmax 2 \
		--outSJfilterReads Unique \
		--outSAMmultNmax 1 \
		--outFileNamePrefix ${STAR_align_output}/${id}_Ribo \
		--outSAMmapqUnique 60 \
		--readFilesCommand gunzip -c \
		--outSAMattributes All \
		--readFilesIn ${outdir}/${id}.rmrRNA.rmtRNA.fq.gz
done
