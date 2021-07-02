# TargetSupercontigs
Phylogenomics: extract supercontigs from (draft) genome assemblies for more effective target capture sequencing in **CLOSELY related species**

## Rationale
Phylogenomic projects in non-model organisms typically use exons as their targets, often even if their study group consists of very recently diverged, and therefore phylogenetically difficult, taxa. This results in significantly reduced coverage and poor assembly in intronic or flanking regions, which are consistently shown to improve phylogenetic inference due to their faster mutation rates. 

Nowadays whole-genome shotgun sequencing is relatively cheap, making draft genome assembly possible for us poor phylogeneticists working on non-model organisms. 

TargetSupercontigs aims to do what it says on the box by locating target exons in the draft genome, identifying putative paralogs and missing genes, and extracting the full exon + intron for genes passing these filters.

It can also be used to identify genes with huge intronic regions (which are often present e.g in many mammals [[[1]](#1)] and should definitely **NOT** be assumed to be in linkage equilibrium).

# TODO:
Make script to group exons from the same gene into separate target loci based on inferred intronic regions.

## References
<a id="1">[1]</a> 
Gatesy, J. & Springer, M. S. (2014). 
Phylogenetic analysis at deep timescales: Unreliable gene trees, bypassed hidden support, and the coalescence/concatalescence conundrum,
Molecular Phylogenetics and Evolution, 80, 231-266.
