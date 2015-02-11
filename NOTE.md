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
