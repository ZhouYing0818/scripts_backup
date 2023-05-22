#! /bin/bash/
outdir=/alldata/zhouying/RNAseq/X101SC22084125-Z01-J001-B1-22/02.clean_data
cd ../01.RawData
ls ./|while read id;
do
	echo ${id}' running...'
	read1=`find $id -name '*_1.fq.gz'`
	read2=`find $id -name '*_2.fq.gz'`
	#echo $read1
	#echo $read2
	fastp -i $read1 -o ${outdir}/${id}_1_clean.fq.gz\
		-I $read2 -O ${outdir}/${id}_2_clean.fq.gz\
		-h ${outdir}/${id}.html\
		-j ${outdir}/${id}.json\
		-z 5 -f 15 -F 15
done
