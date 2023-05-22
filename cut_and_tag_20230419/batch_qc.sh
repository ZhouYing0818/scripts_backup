#! /bin/bash/
rootdir=/alldata/zhouying/CUT_and_TAG/Jing_20230224
index=/alldata/zhouying/ref_data/genome/GRCh38.p13/index/star
gtf=/alldata/zhouying/ref_data/anno/GRCh38.p13/gencode.v36.annotation.gtf
outdir=${rootdir}/03_qc_results/RawData
indir=${rootdir}/02_RawData
ls ${indir}|while read id;
do
	echo ${id}' running...'
	id1=`echo ${id}|awk -F '_' '{print $3}'`_`echo ${id}|awk -F '_' '{print $4}'`
	id2=`echo ${id}|awk -F '_' '{print $5}'`
	read1=${indir}/${id}/${id1}_${id2##*-}_1.fq.gz
	read2=${indir}/${id}/${id1}_${id2##*-}_2.fq.gz
	#echo $read1
	#echo $read2
	fastqc -t 20 -o ${outdir} ${read1} ${read2}
done
