#!/usr/bin/env python3

from sys import argv, stdout
import numpy as np
import scipy.stats as stats
from pandas import read_table, merge, Index
from random import sample

np.random.seed(1)  # So that I get the same result every time.
                   # Is this good practice?

data = read_table(argv[1])
genes = data.columns[2:]
shuffled = data.copy()
shuffled.index = Index(sample(list(shuffled.index), len(shuffled)))
paired = merge(data, shuffled,
               left_index=True, right_index=True, suffixes=("_A", "_B"))
paired['delta_log_16S'] = np.log2(paired['copies16S_A']) - \
                          np.log2(paired['copies16S_B'])

stdout.write("gene\tpearson\tspearman\tkendall\n")
for gene in genes:
    delta_log16S = paired['delta_log_16S']
    delta_gene = paired[gene + "_A"] - paired[gene + "_B"]
    pearson = stats.pearsonr(delta_log16S, delta_gene)[0]
    spearman = stats.spearmanr(delta_log16S, delta_gene)[0]
    kendall = stats.kendalltau(delta_log16S, delta_gene)[0]

    # if not np.isnan(correlation):
    stdout.write("%s\t%f\t%f\t%f\n" % (gene, pearson, spearman, kendall))
