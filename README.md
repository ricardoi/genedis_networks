# Genetic Distance Network

Phythophtora genomes were sequenced, and used to generate VCFs files containing SNPs.
The SNPs were used to calculate the genetic distances of _P. ramorum_ for lineages EU1 and NA1. The genetic distances were used to construct phylogenetic trees, however, due to its homogeneity  .

We can visualize this relationships with epidemic networks using genetic data. The potential cutoff for your analysis can be determined by identifying a threshold that best differentiates a bimodal distribution of the genetic distances that are typically present in your sequences [MicroTrace Documentation](https://raw.githubusercontent.com/CDCgov/MicrobeTrace/master/docs/MicrobeTrace_Manual.pdf).  Example: The frequency distribution of genetic distances representing genetic distances from known pathogen transmission cases, and the red boxes representing genetic distances from cases without evidence of pathogen transmission.

Variant calling sequences (VCFs)
|Genome reference | Lineage|
|--|--|
|PR-102_v3 | NA1 |
|PR-15-019 | EU1 |

## Genetic Distances
Genetic distances of _P. ramorum_ calulated by N. Carleson
![Genetic Distances](https://github.com/ricardoi/genedis_networks/blob/main/Figures/Genetic_distances_subset10_.png)

## Phylogenetic tree: NJ
NJ tree using the original genetic Distances
![NJ tree](https://github.com/ricardoi/genedis_networks/blob/main/Figures/EU1-NA1_trees_NJ.png)

## Grafen rescaling trees
Trees rescaled using the Grafen method
![Grafen tree](https://github.com/ricardoi/genedis_networks/blob/main/Figures/EU1-NA1_Grafen_trees.png)

## Branch lengths from Grafen tree
Visualizing the branch length distribution of _P. ramorum_ trees
![Branch length](https://github.com/ricardoi/genedis_networks/blob/main/Figures/EU1-NA1_branchlenghts_Grafen_tree.png)

## Phylogenetic network
The phylogenetic tree was visualized as a network, tips and nodes are indicated in circles.
![Phylonet](https://github.com/ricardoi/genedis_networks/blob/main/Figures/EU1-NA1_phylogenetic_network_full.png)

## Phylogenetic network prunned
Removing branches with distances less than 0.08
![prunned network]

> Notes:
`$PATH` to folders `/dfs/Grunwald_Lab/lab_collabs/diagnostic_assays/genome_data/`
