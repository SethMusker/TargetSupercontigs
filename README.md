# TargetSupercontigs
Phylogenomics: extract supercontigs from (draft) genome assemblies for more effective target capture sequencing in **CLOSELY related species**.

## Rationale
Phylogenomic projects in non-model organisms typically use exons as their targets, often even if their study group consists of very recently diverged, and therefore phylogenetically difficult, taxa. This results in significantly reduced coverage and poor assembly of intronic or flanking regions, which are consistently shown to improve phylogenetic inference due to their faster mutation rates. 

At the same time, whole-genome shotgun sequencing is now relatively cheap, making draft genome assembly possible for us poor phylogeneticists working on non-model organisms. 

TargetSupercontigs aims to locate target exons in the draft genome, identify putative paralogs and missing genes, and extract the full exon + intron for genes passing these filters. These can then be used as targets instead of (or in combination with) the original exons, provided all the species you're working with are relatively closely related (at deep time scales introns will become too divergent to be useful as targets; exactly how deep the time scale limit is is very hard to say).

It can also be used to identify genes with huge intronic regions (which are often present e.g. in many mammals [[[1]](#1)] and should definitely **NOT** be assumed to be in linkage equilibrium).

## TODO:
Make script to group exons from the same gene into separate target loci based on inferred intronic regions.
Make sure to flag as paralogs regions on the same chromosome matching the same target. Plotting option might be to label each segment individually, arrange vertically based on position on chromosome after stratifying chromosomes

## References
<a id="1">[1]</a> 
Gatesy, J. & Springer, M. S. (2014). 
Phylogenetic analysis at deep timescales: Unreliable gene trees, bypassed hidden support, and the coalescence/concatalescence conundrum,
Molecular Phylogenetics and Evolution, 80, 231-266.
