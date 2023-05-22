#! /bin/bash/
#ribocode find P-sites
rootdir=/alldata/zhouying/Ribo-seq/ribo-TISH_test_data
RiboCode_annot=/alldata/zhouying/ref_data/RiboCode_anno
OUT=${rootdir}/ribocode
indir=${rootdir}/05_align_star

mkdir -p ${OUT}
mkdir -p ${OUT}/metaplots
mkdir -p ${OUT}/ORF
mkdir -p ${OUT}/ORF_count

ls ${rootdir}/01_RawData|grep SRR|while read id;
do	
	echo ${id}' running metaplots...'
	metaplots -a ${RiboCode_annot} \
		-r ${indir}/${id}_RiboAligned.toTranscriptome.out.bam \
		-o ${OUT}/metaplots/${id}_metaplots \
		-m 26 -M 50 -s yes -pv1 1 -pv2 1
	
	echo ${id}' ORF calling...'
	config=${OUT}/metaplots/${id}_metaplots_pre_config.txt
	
	RiboCode -a ${RiboCode_annot} -l no -g -c ${config} -o ${OUT}/ORF/${id}_ORF
	
	echo ${id}' ORF counting...'
	ORFcount -g ${OUT}/ORF/${id}_ORF.gtf \
		-r ${indir}/${id}_RiboAligned.sortedByCoord.out.bam \
		-f 15 -l 5 -e 100 -m 24 -M 35 -s yes \
		-o ${OUT}/ORF_count/${id}_ORF_count_all.txt
done
#plot_orf_density -a ${RiboCode_annot} -c ${config} -t ENST00000414273.1 -s 274 -e 309 --start-codon STARTCODON -o ${OUT}/ORF/ORF
