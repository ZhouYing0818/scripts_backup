#! /bin/bash/
rootdir=/alldata/zhouying/CUT_and_TAG/Jing_20230224
#index=/alldata/zhouying/ref_data/genome/Ecoli/index/bowtie2/E.coli_K12_MG1655
outdir=${rootdir}/06_calling-peaks/histone_broad
indir=${rootdir}/06_calling-peaks/histone_broad
ls ${rootdir}/01_RawData|while read id;
do
	echo ${id}' running...'
	#macs2 callpeak -t ${indir}/${id}_hg38_sort.rmdup.bam -f BAMPE -B -g hs -n ${id} -q 0.05 --outdir ${outdir}
	bedGraphToBigWig ${indir}/${id}_treat_pileup.bdg /alldata/zhouying/ref_data/genome/hg38/hg38.fa.fai ${id}.bw
done
