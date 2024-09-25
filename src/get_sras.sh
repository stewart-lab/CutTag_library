# Droplet_Paired-Tag_H3K27ac_rep1_DNA; Mus musculus;
~/programs/sratoolkit.3.0.5-centos_linux64/bin/fastq-dump --split-files --gzip SRR23343778
# Droplet_Paired-Tag_H3K27ac_rep2_DNA; Mus musculus; 
~/programs/sratoolkit.3.0.5-centos_linux64/bin/fastq-dump  --split-files --gzip SRR23343774

# note on ATAC fastqs:
# ATAC FASTQs
# [Sample Name]S1_L00[Lane Number][Read Type]_001.fastq.gz

# Where Read Type is one of:

# I1: Dual index i7 read (optional)
# R1: Read 1
# R2: Dual index i5 read
# R3: Read 2
# Cell Ranger ARC will also accept ATAC FASTQs in this format:

# I1: Dual index i7 read (optional)
# R1: Read 1
# I2: Dual index i5 read
# R2: Read 2
