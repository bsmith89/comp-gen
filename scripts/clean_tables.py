#!/usr/bin/env python3

from sys import argv, stdout
from pandas import read_csv

data = read_csv(argv[1])
data[data.columns[1:]].to_csv(stdout, '\t', index=False)
