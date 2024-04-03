# get genome from ncbi

# install datasets via conda
conda create -n ncbi_datasets
conda activate ncbi_datasets
conda install -c conda-forge ncbi-datasets-cli
# download axolotle genome
datasets download genome accession GCA_002915635.3 --include gff3,rna,cds,protein,genome,seq-report
# download ecoli genome
datasets download genome accession GCF_000005845.2 --include gff3,rna,cds,protein,genome,seq-report