#!/usr/bin/env python3

from sys import argv, stdout
from Bio.Phylo import read as read_tree
from pandas import read_table, DataFrame
from itertools import chain, permutations
from scipy.stats import pearsonr, spearmanr, kendalltau
from numpy import log2

def differ(a, b, by_all=False):
    if not by_all:
        return not a == b
    else:
        return not any([x == y for x, y in zip(a, b)])

def mixed_clades(clade, key, by_all=False):
    """Search for phylogenetically independent mixed clades.

    Find the maximum number of evolutionarily independent taxon pairs.
    See Maddison 2000, J. theor. Bio.

    The function is applied recursively to a node, returning and mixed
    clades above that node, and the key value of that node.

    The key value of a node is the key value of all taxa above that
    node which are not already paired.

    """
    if clade.is_terminal():
        return ([], key(clade))

    paired_above = []
    subclade_values = []
    for subclade in clade.clades:
        paired, value = mixed_clades(subclade, key, by_all=by_all)
        paired_above.extend(paired)
        if value:
            subclade_values.append(value)

    if not subclade_values:
        return (paired_above, None)
    elif len(subclade_values) == 1:
        return (paired_above, subclade_values[0])
    elif differ(subclade_values[0], subclade_values[1], by_all=by_all):
        terminals = clade.get_terminals()
        # Remove any taxa which have already been paired above.
        pair = list(set(terminals) - set(chain(*paired_above)))
        return (paired_above + [pair], None)
    else:
        return (paired_above, subclade_values[0])
    # TODO: When you are looking for pairs differing in both traits,
    # how do you account for the loss of some combinations of values in this
    # step.

def choose_pairs(taxa, key, by_all=False):
    """Choose pairs of taxa which differ by the key.

    If *by_all* then each element of the key value must be different.

    Will return the first pair which satisfy this test.

    """
    for pair in permutations(taxa, 2):
        if differ(*[key(taxon) for taxon in pair], by_all=by_all):
            return pair

def indep_pairs(clade, key, by_all=False):
    return [choose_pairs(indep_sets, key, by_all=by_all)
            for indep_sets in mixed_clades(clade, key, by_all=by_all)[0]]

def contrast_pairs(pairs, characters, **keys):
    """Contrast each taxon for each key value in keys."""
    return DataFrame(dict(("delta_%s" % trait,
                           [keys[trait](tax1) - keys[trait](tax2)
                            for tax1, tax2 in pairs])
                           for trait in keys))

if __name__ == "__main__":
    tree = read_tree(argv[1], 'newick', rooted=True)
    characters = read_table(argv[2])
    characters.index = characters.treelabel
    characters["log_copies16S"] = log2(characters.copies16S)
    char1, char2 = "log_copies16S", argv[3]

    key1 = lambda x: characters[char1][x.name]
    key2 = lambda x: characters[char2][x.name]
    key = lambda x: (key1(x), key2(x))
    pairs = indep_pairs(tree.clade, key)
    data = contrast_pairs(pairs, characters,
                          **{char1: key1, char2: key2})

    covar1 = data['delta_%s' % char1]
    covar2 = data['delta_%s' % char2]
    pr = pearsonr(covar1, covar2)
    sr = spearmanr(covar1, covar2)
    kt = kendalltau(covar1, covar2)

    stdout.write("(%s) x (%s)\n"
                 "---------------------------\n"
                 "no. pairs: %d\n"
                 "pearson: %f (%f)\n"
                 "spearman: %f (%f)\n"
                 "kendall: %f (%f)\n"
                 %
                 (char1, char2, len(pairs),
                 pr[0], pr[1], sr[0], sr[1], kt[0], kt[1]))
