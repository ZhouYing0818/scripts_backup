input_file='GSM4859382_HEK293_siBAHD1_FLAG_WT_1_cell_PE_CR.bw GSM4859383_HEK293_siBAHD1_FLAG_WT_2_cell_PE_CR2.bw GSM4859387_HEK293_siBAHD1_WT_cell_H3K27ac_ChIP.bw GSM4859388_HEK293_siBAHD1_W667G_mutant_cell_H3K27ac_ChIP.bw GSM4859391_HEK293_siBAHD1_WT_cell_H3K27me3_ChIP.bw GSM4859392_HEK293_siBAHD1_W667G_mutant_cell_H3K27me3_ChIP.bw'
Reference=/alldata/zhouying/ref_data/anno/hg19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz

computeMatrix reference-point \
	--referencePoint TSS \
	-b 5000 -a 5000 \
	-R ${Reference} \
	-S ${input_file} \
	--skipZeros \
	-o BAHD1_H3K27ac_H3K27me3.gz

plotHeatmap -m BAHD1_H3K27ac_H3K27me3.gz \
	-out BAHD1_H3k27ac_H3K27me3_TSS.pdf \
	--colorMap Reds \--whatToShow 'heatmap and colorbar' \
	--zMin 0 --zMax 3 \--kmeans 3 \
	--samplesLabel BAHD1-1 BAHD1-2 H3K27ac-WT H3K27ac-W667G H3K27me3-WT H3K27me3-W667G
