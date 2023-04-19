#! /bin/bash/
rootdir=/alldata/zhouying/RNAseq/Fan_20230412/X101SC22084125-Z01-J009
index=/alldata/zhouying/ref_data/genome/GRCh38.p13/index/star
gtf=/alldata/zhouying/ref_data/anno/GRCh38.p13/gencode.v36.annotation.gtf
outdir=${rootdir}/04_align_star
indir=${rootdir}/02_CleanData
ls ${rootdir}/02_CleanData|while read id;
do
	echo ${id}' running...'
	read1=${indir}/${id}/${id}_1.clean.fq.gz
	read2=${indir}/${id}/${id}_2.clean.fq.gz
	#echo $read1
	#echo $read2
	STAR --twopassMode Basic \
		--quantMode TranscriptomeSAM GeneCounts \
		--runThreadN 20 \
		--genomeDir ${index} \
		--alignIntronMin 20 \
		--alignIntronMax 50000 \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMattrRGline ID:${id} SM:${id} PL:ILLUMINA \
		--outFilterMismatchNmax 2 \
		--outSJfilterReads Unique \
		--outSAMmultNmax 1 \
		--outFileNamePrefix ${outdir}/${id}_sort \
		--outSAMmapqUnique 60 \
		--readFilesCommand gunzip -c \
		--readFilesIn ${read1} ${read2}
done
