#! /bin/bash/
rootdir=/alldata/zhouying/CUT_and_TAG/Jing_20230224
#index=/alldata/zhouying/ref_data/genome/Ecoli/index/bowtie2/E.coli_K12_MG1655
outdir=${rootdir}/05_RemoveDup
indir=${rootdir}/04_align_bowtie2
ls ${rootdir}/01_RawData|while read id;
do
	echo ${id}' running...'
	java -jar /alldata/zhouying/evnconfig/picard.jar MarkDuplicates \
		REMOVE_DUPLICATES=true \
		I=${indir}/${id}_sort.bam \
		O=${outdir}/${id}_hg38_sort.rmdup.bam \
		VALIDATION_STRINGENCY=SILENT \
		M=${outdir}/${id}_hg38_picard.txt

	#java -jar /alldata/zhouying/evnconfig/picard.jar MarkDuplicates \
	#	REMOVE_DUPLICATES=true \
	#	I=${indir}/${id}_Ecoli_sort.bam \
	#	O=${outdir}/${id}_Ecoli_sort.rmdup.bam \
	#	VALIDATION_STRINGENCY=SILENT \
	#	M=${outdir}/${id}_Ecoli_picard.txt
done
