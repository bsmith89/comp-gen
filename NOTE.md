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

## Glossary ##
[Kendall's $\tau$](https://en.wikipedia.org/wiki/Kendall_tau_rank_correlation_coefficient)
:   A non-parametric measure of the correlation between two vectors.
    It is very similar to Spearman's $\rho$, but has some slightly superior
    statistical properties.

Focal trait
:   A genomic trait which is compared to 16S copy number in the search
    for a statistical correlation.
    In this project, these traits are either K0 orthologous groups,
    or M0 modules.
    The trait value can either be binary (presence/absence) or positive
    integers (number of copies).
    For example, the 16S ribosomal RNA gene (K01977) itself can have any value
    between 0 and 15.

Contrast table
:   A table of two vectors encoding the arithmetic differences between
    paired species in two traits.
    For example if two pairs (4 taxa) have the following trait values
    for two traits:

    Pair    Trait 1 (A)    Trait 1 (B)    Trait 2 (A)    Trait 2 (B)
    ----    -----------    -----------    -----------    -----------
    A/B     1              2              5              3
    C/D     2              1              4              5

    then the contrast table will be:

    Pair    Trait 1    Trait 2
    ----    -------    -------
    A/B     -1         2
    C/D     1          -1
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


#### Q1 ####
What are the relationships between number of pairs, $\tau$, and p-value?

![](static/2015-02-19_fig1.png)

A few things stick out.

1. The number of pairs (I'll also refer to it as $n$) ranges from 0 to
    about 350.
2. Clearly $\tau$ is closely related to the number of pairs found.
    -  At low numbers of pairs, a great range of $\tau$'s are possible,
        including many with high values.
        These are probably worthless.
    -  The distribution of $\tau$ doesn't change much at higher values of
        $n$ (pairs).
    -  There is a subtle but consistent skew towards values of $\tau$ above 0.
        This especially stands out at intermediate number of pairs, although
        this could be the result of thinning numbers of points making
        it more visible.
3. A relatively small minority of traits pass both a $n$ and p-value cutoff.
4. Modules generally have a lower $n$ than genes.


#### Q2 ####
What's the relationship between the absolute number of species with a gene,
and the number of pairs that were found?

![](static/2015-02-19_fig2.png)

1. There is no obvious relationship between the number of genomes with a trait
    and p-value.
2. The number of pairs increases as the frequency of the trait increases,
    although this trend is not perfect.
3. $\tau$ _may_ increase with frequency, but it's not clear.
4. Many high $\tau$ traits are at high frequency
    -  These may represent traits which mostly vary in number,
        not presence/absence. (e.g. tRNAs)


#### Q3 ####
How does the distribution of $\tau$ values change as you adjust the p-value
cutoff?

I'm really proud of this figure:

![](static/2015-02-19_fig3.png)

I call it a fire plot.
The point is that you can trace a histogram at an arbitrary cutoff value,
in this case p-value, by tracing the 'iso-hue' line.

This figure shows that the shape of the distribution changes a _lot_
with p-value.
As you might expect, it become multimodal, since small absolute values of
$\tau$ don't lend themselves to statistical significance.


## 2015-02-20 ##
#### Q4 ####
Is there bias in the estimated correlations introduced by the picking of
pairs which vary in both traits?


![](static/2015-02-20_fig1.png)

Here I've used various permutation methods to get at the null distribution of
p-values.
When I permute the vector of values for K00001 (which does not show any
statistical significance), and then re-pick mixed pairs, contrast the pairs,
and calculate a p-value for Kendall's test, I still get a higher than
expected frequency of _low_ p-values.
(The expectation for p-values is uniform.)
This means that my picking method is introducing bias, since I _know_ that
there is no evolutionary correlation between copy number and the permuted
trait values.
I also know that this bias is not the result of the Kendall test on
vectors of contrasts (which, by definition have no 0's), since
I compare the p-value distribution to that obtained by permuting the
contrast table;
the result of that is a fairly uniform distribution.

It therefore seems possible that I am getting some correlations with
p-values under 0.05 purely due to this bias.  In fact, in this permutation run,
8.7% of my p-values came back under 0.05 and 1 in 200 came back $< 0.001$.

This seems problematic, since I want to be able to identify genes
with a p-value threshold.
I'd also like to be able to account for multiple testing
by doing a Bonferonni correction, but bias will invalidate this approach.

...

What fraction of the time might my hits under some p-value threshold
be accounted for by this bias?
I can be a liiiitle comforted by the next result:

![](static/2015-02-20_fig2.png)

This figure shows that for traits which passed a 0.001 ("highly significant")
p-value threshold, the majority of them are no longer under that
threshold when we permute. (2.7% have p-values under this cutoff.
We would have expected 0.1%).
