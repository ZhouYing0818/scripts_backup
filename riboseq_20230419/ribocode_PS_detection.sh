#! /bin/bash/
#ribocode find P-sites
rootdir=/alldata/zhouying/Ribo-seq/ribo-TISH_test_data
RiboCode_annot=/alldata/zhouying/ref_data/RiboCode_anno
OUT=${rootdir}/ribocoda
indir=${rootdir}/05_align_star

mkdir -p ${OUT}
mkdir -p ${OUT}/metaplots

metaplots -a ${RiboCode_annot} \
-r ${indir}/SRR6327777_RiboAligned.toTranscriptome.out.bam \
-o ${OUT}/metaplots/SRR6327777_metaplots_ \
-m 26 -M 50 -s yes -pv1 1 -pv2 1
