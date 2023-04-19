#! /bin/bash/
ref_dir=/alldata/zhouying/evnconfig/DaPars2-master/Dapars2_Test_Dataset
outdir=/alldata/zhouying/RNAseq/Dapars_test
indir=/alldata/zhouying/RNAseq/Dapars_test

cd ${outdir}

echo "Annotated_3UTR=${ref_dir}/RefSeq_hg38_3UTR_annotation.bed" >> dapars_config.txt
echo "Aligned_Wig_files=` cut -f1 all_samples_mapped.txt|sed '#g'|paste -s -d ","`" >> dapars_config.txt
echo "Output_directory=HepG2" >> dapars_config.txt
echo "Output_result_file=HepG2" >> dapars_config.txt
echo "Coverage_threshold=10" >> dapars_config.txt
echo "Num_Threads=30" >> dapars22_config.txt
echo "sequencing_depth_file=all_samples_mapped.txt" >> dapars_config.txt
