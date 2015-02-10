#!/usr/bin/env python3

from sys import argv, stdout
import numpy as np
from pandas import read_table, merge
from random import shuffle

np.random.seed(1)  # So that I get the same result every time.
                   # Is this good practice?

data = read_table(argv[1])
genes = data.columns[2:]
paired = merge(data, data.apply(np.random.permutation),
               left_index=True, right_index=True, suffixes=("_A", "_B"))
paired['delta_log_16S'] = np.log2(paired['copies16S_A']) - \
                      np.log2(paired['copies16S_B'])

correlation = {}
for gene in genes:
    correlation = np.corrcoef(paired['delta_log_16S'],
                              paired[gene + "_A"] - paired[gene + "_B"]
                             )[0,1]
    if not np.isnan(correlation):
        stdout.write("%s\t%f\n" % (gene, correlation))
