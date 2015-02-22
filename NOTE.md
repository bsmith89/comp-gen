---
title: Comparative Genomics in Bacteria
author: Byron Smith & Ben Roller
...

# Introduction #
Here we will analyze the genomic signature of a growth-rate/yield trade-off.

16S rRNA copy number is an indicator of position along a life history axis.
Here we will explore the genomic content which goes along with that indicator.

# Requirements and Environment #
## Required Applications ##
-  [`python`](http://www.python.org/) (3.4+)
-  GNU Make
-  BASH

# Notebook #
## 2015-02-06 ##
### Prototyping a Sister-pairs Analysis ###
I'm going to try and test each gene for a relationship with copy number.
Usually this would require accounting for phylogenetic relationship,
so instead I'm going to use sister-paired taxa, to control for this.

## 2015-02-09 ##
### Prototyping a Sister-pairs Analysis ###
I removed lines 15 and 16 and renamed the RAxML tree to `tre/bacteria.nex`.
This allows the file to be opened with figtree.

## 2015-02-10 ##
### Sister-pairs Analysis ###
Unlikely the all-by-all potential pairings that I have worked with up to this
point, paired analyses actually require that there not be overlap between
the evolutionary histories of the pairs being compared (because this
represents non-independence/pseudoreplication).
Looking at Maddison 2000, I believe that I can implement their algorithm for
selecting a maximal number of taxon pairs which differ in both characters.
I intend to use this to devise paired sets for each character contrast.
I will also have to calculate a null distribution with a permutation test
BEFORE matching pairs.

## 2015-02-11 ##
### Independent pairs ###
I managed to get the algorithm working.
Now I just need to package it into a python script which will do the full
analysis.

## 2015-02-12 ##
### Phylogenetically Independent Correlations ###
I'd like to use evolutionarily independent pairing to test the correlations
between various genes and 16S copy number (actually _log_ copy number.
Ultimately I might want to compare it pairwise by genes).

## 2015-02-15 ##
### Cont. ###
I can run all of the basic correlations for each K0 in about 10 minutes.
M0 take far less time.
All results have been saved to `res/*0.paired_corr.tsv`.
I get a few hundred results for K0s with more than 50 pairs considered and
p-values less than 0.005.
The top two hits are 16S and 23S rRNA genes, as you might expect.
Hits further down the list are more interesting.  What I should now do is
compare these correlations with a null distribution based on a permutation
test.


## 2015-02-19 ##


### Fixed the non-deterministic behavior of `scripts/pairwise_log16S_corr.py` ###
I found that each run of the pairwise correlation script gave a different output.
This didn't make sense, and scared me.

Turns out that I was assuming my dictionary keys would be iterated back in the
same order each time.
This was incorrect.


### Phylogenetically Independent Correlations ###
I re-did some of the analysis that had been done over the last few days.
This has allowed me to clean up a few figures, and think critically about
these results.
