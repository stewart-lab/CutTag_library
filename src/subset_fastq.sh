# get seqtk to modeify fasta files
git clone https://github.com/lh3/seqtk.git;
cd seqtk; make
# add to .bash_profile
export PATH=~/seqtk/:$PATH
# subsample (using same random seed to keep pairing)
# grep '^@' 240226_EO_Ax1_K27me3_1_S26_L005_R1_001_trimmed.fastq | wc -l = 41166436
# 25% of reads
seqtk sample -s100 240226_EO_Ax1_K27me3_1_S26_L005_R1_001_trimmed.fastq.gz 10291609 > 240226_EO_Ax1_K27me3_1_S26_L005_R1_001_trimmed_sub25.fastq
seqtk sample -s100 240226_EO_Ax1_K27me3_1_S26_L005_R2_001_trimmed.fastq.gz 10291609 > 240226_EO_Ax1_K27me3_1_S26_L005_R2_001_trimmed_sub25.fastq
# for 50% of reads
seqtk sample -s100 240226_EO_Ax1_K27me3_1_S26_L005_R1_001_trimmed.fastq.gz 20583218 > 240226_EO_Ax1_K27me3_1_S26_L005_R1_001_trimmed_sub50.fastq
seqtk sample -s100 240226_EO_Ax1_K27me3_1_S26_L005_R2_001_trimmed.fastq.gz 20583218 > 240226_EO_Ax1_K27me3_1_S26_L005_R2_001_trimmed_sub50.fastq
# for 10% of reads
seqtk sample -s100 240226_EO_Ax1_K27me3_1_S26_L005_R1_001_trimmed.fastq.gz 4116644 > 240226_EO_Ax1_K27me3_1_S26_L005_R1_001_trimmed_sub10.fastq
seqtk sample -s100 240226_EO_Ax1_K27me3_1_S26_L005_R2_001_trimmed.fastq.gz 4116644 > 240226_EO_Ax1_K27me3_1_S26_L005_R2_001_trimmed_sub10.fastq

# for next mark
# zcat 240226_EO_Ax1_K36me3_1_S27_L005_R1_001_trimmed.fastq.gz | grep '^@' | wc -l 
# 25% of reads
seqtk sample -s100 240226_EO_Ax1_K36me3_1_S27_L005_R1_001_trimmed.fastq 9614105 > 240226_EO_Ax1_K36me3_1_S27_L005_R1_001_trimmed_sub25.fastq
seqtk sample -s100 240226_EO_Ax1_K36me3_1_S27_L005_R2_001_trimmed.fastq.gz 9614105 > 240226_EO_Ax1_K36me3_1_S27_L005_R2_001_trimmed_sub25.fastq
# 50% of reads
seqtk sample -s100 240226_EO_Ax1_K36me3_1_S27_L005_R1_001_trimmed.fastq 19228210 > 240226_EO_Ax1_K36me3_1_S27_L005_R1_001_trimmed_sub50.fastq
seqtk sample -s100 240226_EO_Ax1_K36me3_1_S27_L005_R2_001_trimmed.fastq.gz 19228210 > 240226_EO_Ax1_K36me3_1_S27_L005_R2_001_trimmed_sub50.fastq
# 10% of reads
seqtk sample -s100 240226_EO_Ax1_K36me3_1_S27_L005_R1_001_trimmed.fastq 3845642 > 240226_EO_Ax1_K36me3_1_S27_L005_R1_001_trimmed_sub10.fastq
seqtk sample -s100 240226_EO_Ax1_K36me3_1_S27_L005_R2_001_trimmed.fastq.gz 3845642 > 240226_EO_Ax1_K36me3_1_S27_L005_R2_001_trimmed_sub10.fastq