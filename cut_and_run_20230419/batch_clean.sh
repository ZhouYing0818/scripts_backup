#! /bin/bash/
outdir=/alldata/zhouying/CUT\&RUN/GSE136332_SALL4/03_clean_fastq
cd ../01.RawData
ls ./|grep SRR|while read id;
do
	echo ${id}' running...'
	read1=`find $id -name '*_1.fastq.gz'`
	read2=`find $id -name '*_2.fastq.gz'`
	#echo $read1
	#echo $read2
	#echo ${outdir}
	fastp -i $read1 -o ${outdir}/${id}_1_clean.fq.gz\
		-I $read2 -O ${outdir}/${id}_2_clean.fq.gz\
		-h ${outdir}/${id}.html\
		-j ${outdir}/${id}.json\
		-z 5
done
