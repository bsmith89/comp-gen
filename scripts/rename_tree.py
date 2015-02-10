#!/usr/bin/env python3

from sys import argv, stdout
from Bio.Phylo import read, write

tree_path = argv[1]
names_path = argv[2]

tree = read(tree_path, 'newick')
with open(names_path) as names_handle:
    names_map = dict(tuple(line.split()) for line in names_handle)
for leaf in tree.get_terminals():
    leaf.name = names_map[leaf.name]

write([tree], stdout, 'newick')
