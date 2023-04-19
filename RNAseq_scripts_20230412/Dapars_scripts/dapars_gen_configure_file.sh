#! /bin/bash/
ref_dir=/alldata/zhouying/evnconfig/DaPars2-master/Dapars2_Test_Dataset
outdir=/alldata/zhouying/RNAseq/Dapars_test
indir=/alldata/zhouying/RNAseq/Dapars_test

cd ${outdir}

echo "Annotated_3UTR=${ref_dir}/RefSeq_hg38_3UTR_annotation.bed" > dapars_config.txt
echo "Group1_aligned_Wig=` cut -f1 all_samples_mapped.txt|sed -n '1,3p'|paste -s -d ","`" >> dapars_config.txt
echo "Group2_aligned_Wig=` cut -f1 all_samples_mapped.txt|sed -n '10,12p'|paste -s -d ",    "`" >> dapars_config.txt
echo "Num_least_in_group1=1" >> dapars_config.txt
echo "Num_least_in_group2=2" >> dapars_config.txt
echo "Coverag_cutoff=30" >> dapars_config.txt
echo "FDR_cutoff=0.05" >> dapars_config.txt
echo "PDUI_cutoff=0.4" >> dapars_config.txt
echo "Fold_change_cutoff=0.59" >> dapars_config.txt
