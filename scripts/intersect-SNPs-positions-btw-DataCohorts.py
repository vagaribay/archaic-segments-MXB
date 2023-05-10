__author__ = 'vgaribay'

import argparse, pandas as pd

USAGE = """
find_shared_sites_btw_data_cohorts.py -v1
                -v2
                -o
"""

parser = argparse.ArgumentParser()

parser.add_argument('-v1', '--vcf1', action="store", dest="data_1", help="", default=None, type=str)
parser.add_argument('-v2', '--vcf2', action="store", dest="data_2", help="", default=None, type=str)
parser.add_argument('-o', '--out', action="store", dest="out_file", help="", default=None, type=str)

args = parser.parse_args()

data_1 = pd.read_csv(args.data_1, sep = "\t", header = None, usecols = [0,3,4,5],
                     names = ["CHROM", "POS", "REF", "ALT"])

data_2 = pd.read_csv(args.data_2, sep = "\t", header = None, usecols = [0,3,4,5],
                  names = ["CHROM", "POS", "REF", "ALT"])

intersection = pd.merge(data_2, data_1, how='inner', on=['CHROM','POS', 'REF', 'ALT'])

intersection.to_csv(args.out_file, sep='\t', index=False, header=False)
