#!/usr/bin/env python3

from sys import argv, stdout
from Bio.Phylo import read as read_tree
from pandas import read_table, DataFrame
from itertools import permutations
from scipy.stats import pearsonr, spearmanr, kendalltau
import numpy as np
from numpy import log2
import random


def _pair_that_differ(mapping):
    """Select two keys in a dictionary which differ for all values.

    *mapping* is a dictionary of labels mapped to a sequence of values.
    
    Find two keys in dictionary *mapping* for which each index in the
    associated sequence are not equal.
    
    """
    if len(set(mapping.values())) == 1:
        # If only one value tuple in all of the leaves, there's no mixed pair.
        return None
    for key_series in zip(*mapping.values()):
        if len(set(key_series)) == 1:
            # If any value series has only one values, then
            # there are no mixed pairs.
            return None
    # If these things aren't true, then let's search for a mixed pair.
    # The keys have to be sorted here, or else the result is non-deterministic.
    for key1, key2 in permutations(sorted(mapping.keys()), 2):
        seq1, seq2 = mapping[key1], mapping[key2]
        if np.all([val1 != val2 for val1, val2 in zip(seq1, seq2)]):
            return key1, key2
    # TODO: Refactor so that the pairs are returned with their values,
    # therefore values will only have to be looked up once.

def _recur_mixed_pairs(clade, *key_funcs):
    """Worker function for *mixed_pairs*.

    The function is applied recursively to a node, returning pairs and
    independent branches from above that node.

    """
    if clade.is_terminal():
        value = tuple(func(clade.name) for func in key_funcs)  # TODO: Speed this up.
        # It's a leaf.  There are no pairs above.  Return the leaf's
        # name mapping to the leaf's value for each key as a tuple.
        return [], {clade.name: value}

    # Collect the pairs already found above this node and any free leaves
    # above this node.
    pairs_above = []
    free_leaves_above = {}
    for subclade in clade.clades:
        pairs, free_leaves = _recur_mixed_pairs(subclade, *key_funcs)
        pairs_above.extend(pairs)
        free_leaves_above.update(free_leaves)

    if not free_leaves_above:
        # The node above must have formed a pair.
        # There's no way you're going to get another pair at this node.
        # Just return all of the pairs formed above this point and
        # the empty dictionary.
        return pairs_above, free_leaves_above
    else:
        found_pair = _pair_that_differ(free_leaves_above)
        if found_pair:
            # You found a new pair at this node.  Return it
            # at the end of the list of pairs above.  There are now no
            # free leaves above (since they wouldn't be evolutionarily
            # independent from the leaves in this pair), so return an
            # empty dictionary.
            return pairs_above + [found_pair], {}
        else:
            # Well, there's no mixed pair of leaves above.  Return the
            # pairs list and all of the free leaves above this node.
            return pairs_above, free_leaves_above

def mixed_pairs(clade, *key_funcs):
    """Search for phylogenetically independent mixed clades.

    Very shallow wrapper around *_recur_mixed_pairs*.

    Find the maximum number of evolutionarily independent taxon pairs.
    See Maddison 2000, J. theor. Bio.

    If multiple key functions are given, all of them must be unequal for nodes
    to be considered mixed.
    If you would instead like any difference to indicate a mixed pair,
    feed a single key which returns a tuple of values.

    `mixed_pairs(clade, lambda x: data[trait1][x],
                        lambda x: data[trait2][x])`
    will give a different result than
    `mixed_pairs(clade, lambda x: (data[trait1][x],
                                   data[trait2][x]))`
    since the second key returns a tuple of values which are compared
    as a whole, rather than being unpacked by the comparison routine.

    """
    return _recur_mixed_pairs(clade, *key_funcs)[0]


def contrast_pairs(trait_table, pairs):
    tax1s, tax2s = zip(*pairs)
    data = {}
    for trait in trait_table.columns:
        tax1_values = trait_table[trait].ix[list(tax1s)].values
        tax2_values = trait_table[trait].ix[list(tax2s)].values
        data[trait] = tax1_values - tax2_values
    return DataFrame(data, index=pairs)

def trait_corr(trait_table, clade, mixed='both', corrfunc = kendalltau):
    """Calculate the correlation between two traits based on independent pairs.
    
    Takes a DataFrame of two traits, indexed by the names of leaves on
    the tree and a Biopython tree object.
    
    """
    trait1, trait2 = trait_table.columns
    make_key_func = lambda trait: (lambda x: trait_table[trait][x])
    key_func1 = make_key_func(trait1)
    key_func2 = make_key_func(trait2)

    if mixed == 'both':
        pairs = mixed_pairs(clade, key_func1, key_func2)
    elif mixed == 'first':
        pairs = mixed_pairs(clade, key_func1)
    elif mixed == 'second':
        pairs = mixed_pairs(clade, key_func2)
    else:
        key_func = make_key_func(mixed)
        try:
            pairs = mixed_pairs(clade, key_func)
        except Exception as err:
            raise err

    if len(pairs) == 0:
        return pairs, float('nan'), float('nan')

    contrast = contrast_pairs(trait_table, pairs)
    if np.any(contrast.abs().sum() == 0):
        return pairs, float('nan'), float('nan')
    else:
        corr, p_value = corrfunc(contrast[trait1], contrast[trait2])
        return pairs, corr, p_value

def permute_data(table, columns):
    """Permute particular columns of a DataFrame."""
    out = table.copy()
    out[columns] = table.apply(np.random.permutation)[columns]
    return out

def null_corr(trait_table, clade, cols=[2]):
    """Pull from the null distribution of correlations.
    
    Performs a permutation on traits to get at the null distribution of
    correlation scores.

    *cols* is the index of the column to shuffle:
    0: the labels
    1: the first column
    2: the second column

    # TODO: Fix up this function.
    
    """
    original_columns = trait_table.columns
    trait_table = trait_table.reset_index()
    trait_table = permute_data(trait_table, cols)
    trait_table.index = trait_table.treelabel
    trait_table = trait_table[original_columns]
    return trait_corr(trait_table, clade)

def main():
    tree = read_tree(argv[1], 'newick', rooted=True)
    characters = read_table(argv[2])
    characters.index = characters.treelabel
    characters['log_copies16S'] = np.log2(characters.copies16S)
    trt1 = 'log_copies16S'
    traits = characters.columns[2:-1]
    pairs = mixed_pairs(tree.clade, lambda taxon: characters[trt1][taxon])

    print("trait1", "trait2", "kendall_t", "kendall_p", sep="\t", file=stdout)
    for trt2 in traits:
        contrast = contrast_pairs(characters[[trt1, trt2]], pairs)
        kendall_t, kendall_p = kendalltau(contrast[trt1], contrast[trt2])
        print(trt1, trt2, kendall_t, kendall_p, sep="\t", file=stdout)

if __name__ == '__main__':
    main()
