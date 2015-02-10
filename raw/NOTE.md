### RAxML_3473_1_flux_trim4_10.19.14.nex ###
From an email Ben Roller sent me on 2015-02-06, subject: KEGG orthology,
module, and tree data.

This file is a phylogenetic tree in the NEXUS format.  Figtree cannot read
the file correctly unless I've removed two of the directives.
According to that email:

```
The nexus formatted phylogenetic tree is the file
"RAxML_3473_1_flux_trim4_10.19.14.nex". There are 1147 bacteria, labeled by
their arb accession number.
```

Ben tells me that he made the file with his own decisions about what parts of
the 16S alignment should be masked.
He did not use the SILVA tree, although I think that would be valid.

Rooting is informative.  Ben apparently used an Archaea out-group, rooted,
and then pruned off that clade.

### 2.6.15_KEGG_K0.csv ###
From the same email as 'RAxML...'.

According to that email:

```
The K0 data file includes the abundance of each K0 (columns) within each of
1147 genomes (rows) of distinct species. The M0 data file includes the
presence/absence of each KEGG Module (columns labeled M0) within each of the
genomes. Each file also contains two accession number columns used to connect
it to the phylogenetic tree (rowname and column named 'treelabel') and an rrn
copy number column (named 'copies16S'). FYI I left in all K0 and M0 columns
that correspond to the rrn operon (i.e. 16S, 23S, 5S and tRNA genes).
Finally, the M0 file also contains two columns I personally made to help
identify autotrophs (named 'autmodscurated') and oxygenic photoautotrophs
(named 'photoautmodscurated').
```

### 2.6.15_KEGG_M0.csv ###
Same email:

```
The K0 data file includes the abundance of each K0 (columns) within each of
1147 genomes (rows) of distinct species. The M0 data file includes the
presence/absence of each KEGG Module (columns labeled M0) within each of the
genomes. Each file also contains two accession number columns used to connect
it to the phylogenetic tree (rowname and column named 'treelabel') and an rrn
copy number column (named 'copies16S'). FYI I left in all K0 and M0 columns
that correspond to the rrn operon (i.e. 16S, 23S, 5S and tRNA genes).
Finally, the M0 file also contains two columns I personally made to help
identify autotrophs (named 'autmodscurated') and oxygenic photoautotrophs
(named 'photoautmodscurated').
```
