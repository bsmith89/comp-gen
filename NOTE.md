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
