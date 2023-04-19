#! /bin/bash/
rootdir=/alldata/zhouying/Ribo-seq/ribo-TISH_test_data
outdir=${rootdir}/02_CleanData
cd ../01_RawData
ls ./|grep SRR|while read id;
do
	echo ${id}' running...'
	read1=`find $id -name '*.fastq.gz'`
	#$read2=`find $id -name '*_2.fastq.gz'`
	#echo $read1
	#echo $read2
	#echo ${outdir}
	fastp -i $read1 -o ${outdir}/${id}_clean.fq.gz\
	#	-I $read2 -O ${outdir}/${id}_2_clean.fq.gz\
		-h ${outdir}/${id}.html\
		-j ${outdir}/${id}.json
		#-z 5
	#cutadapt -j 10 \
	#	 -a CTGTAGGCACCATCAATTCGTATGCCGTCTT \
	#	 -O 6 -m 20 -M 45 --discard-untrimmed \
	#	 -o ${outdir}/${id}_clean.fq.gz ${read1}
done
