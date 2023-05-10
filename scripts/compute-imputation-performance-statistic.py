__author__ = 'vgaribay'

import allel, argparse, numpy as np, pandas as pd
from scipy.stats import pearsonr

USAGE = """

imputation-performance-statistic.py -v1 -v2 --chr --out

"""

parser = argparse.ArgumentParser()

parser.add_argument('-v1', '--vcf1', action="store", dest="vcf_1", help="", default=None, type=str)
parser.add_argument('-v2', '--vcf2', action="store", dest="vcf_2", help="", default=None, type=str)
parser.add_argument('-c', '--chr', action="store", dest="chrom", help="", default=None, type=int)
parser.add_argument('-o', '--out', action="store", dest="out_file", help="", default=None, type=str)
args = parser.parse_args()

args = parser.parse_args()

# read VCFs
imputed = allel.read_vcf(args.vcf_1)
whole_genome = allel.read_vcf(args.vcf_2)

# extract genotype calls
imputed_gt = allel.GenotypeArray(imputed['calldata/GT'])
whole_genome_gt = allel.GenotypeArray(whole_genome['calldata/GT'])

# SNPs positions
positions = whole_genome['variants/POS']

# Compute statistics:
# 1) genotype correlation
# 2) genotype accuracy
# 3) heterozygotic precision
# 4)homozygotic alternative precision

genotype_correlation = []
genotype_accuracy = []
heterozygotic_precision = []
homozygotic_precision = []

for pos in range(0, len(positions)):
    
    is_missing = whole_genome_gt[pos].is_missing()
    
    whole_genome_genotypes = whole_genome_gt[pos][~is_missing].to_n_alt()
    imputed_genotypes = imputed_gt[pos][~is_missing].to_n_alt()
    
    # genotype correlation: pearson correlation coefficient between genotype vectors
    # correlation coefficient is undefined if an input array is constant (if a snp is fixed in either vector)
    
    genotype_correlation.append(pearsonr(whole_genome_genotypes, imputed_genotypes)[0])
    
    # genotype accuracy: the ratio of the number of correctly called genotypes to the number of called genotypes
    # a genotype is correctly called if it is the value in the imputed data is equal to the real data
    genotype_accuracy.append(sum(whole_genome_genotypes == imputed_genotypes)/len(whole_genome_genotypes))
    
    is_het = (whole_genome_genotypes == 1)
    
    if sum(is_het) > 0:    
        heterozygotic_precision.append(sum(imputed_genotypes[is_het] == 1)/sum(is_het))
    
    else:
        heterozygotic_precision.append(np.nan)
    
    is_hom = (whole_genome_genotypes == 2)
    
    if sum(is_hom) > 0:
        homozygotic_precision.append(sum(imputed_genotypes[is_hom] == 2)/sum(is_hom))
    
    else:
        homozygotic_precision.append(np.nan)
        
out = pd.DataFrame({'chrom': [args.chrom] * len(positions), 'position': positions,
                        'genotype_correlation': genotype_correlation,
                        'genotype_accuracy': genotype_accuracy,
                        'heterozygotic_precision': heterozygotic_precision,
                        'homozygotic_precision': homozygotic_precision})

out.to_csv(args.out_file, sep='\t', index=False, header=True)
