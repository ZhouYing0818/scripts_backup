#! /bin/bash/
rootdir=/alldata/zhouying/Ribo-seq/ribo-TISH_test_data
outdir=${rootdir}/03_qc_results/clean_qc
cd ../02_CleanData
ls ./|grep SRR|while read id;
do
	echo ${id}' running...'
	#read1=`find $id -name '*.fastq.gz'`
	#read2=`find $id -name '*_2.fastq.gz'`
	#echo $read1
	#echo $read2
	#echo ${outdir}
	fastqc -t 12 -q -o  ${outdir} ${id}
done

#cd ${outdir}
#multiqc -n GSE93830_ssb_raw -o ./ ./
