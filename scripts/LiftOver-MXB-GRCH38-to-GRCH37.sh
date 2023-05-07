#!/bin/bash

# Valeria Añorve-Garibay | 2023
# LiftOver was made using CrossMap (https://crossmap.sourceforge.net/)
# Reference Fasta (https://data.broadinstitute.org/snowman/hg19/Homo_sapiens_assembly19.fasta) was downloaded from The Broad Data Resources
# Chain File (https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh38_to_GRCh37.chain.gz) was downloaded from Ensembl FTP

chrom=${SLURM_ARRAY_TASK_ID}

# LiftOver
CrossMap.py vcf GRCh38_to_GRCh37.chain.gz MXB_6k.chr$chrom.b151.numeric.vcf.gz Homo_sapiens_assembly19.fasta \
  hg19/hg19.MXB_6k.chr$chrom.b151.numeric.vcf --compress

# Sort positions
bcftools sort --max-mem 5G --temp-dir /users/vgaribay/scratch/temp/ hg19/hg19.MXB_6k.chr$chrom.b151.numeric.vcf -Oz > \
  sort/hg19.MXB_6k.chr$chrom.b151.numeric.sort.vcf.gz

# Index
tabix -p vcf sort/hg19.MXB_6k.chr$chrom.b151.numeric.sort.vcf.gz
