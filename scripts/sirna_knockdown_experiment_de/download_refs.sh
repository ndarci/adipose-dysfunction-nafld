#!/bin/bash

cd ../../data/

# download GRCh37 reference and annotation from GENCODE

# need genome + annotation to do mapping

# ref = 1-22, X, Y, M
# main (CHR) = ref only
# primary assembly (PRI) = ref + scaffolds
# everything (ALL) = ref + scaffolds + "assembly patches" + "alternate loci (haplotypes)"

# back-mapped CHR genome for 37
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
gunzip -f GRCh37.primary_assembly.genome.fa.gz

# original CHR annotation for 37 
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gunzip -f gencode.v19.annotation.gtf.gz

# refFlat
wget https://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz
