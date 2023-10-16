# genetic-diversity-in-legionella
Compairing genomes of different strains of legionella to count the diversity of this bacteria

## Introduction

This documentation provides an overview of the genetic diversity analysis project for Legionella strains. In this project, I performed various analyses, including multiple sequence alignment, phylogenetic tree construction, and visualization of genetic variations.

## Multiple Sequence Alignment (MSA)

I used multiple sequence alignment to align genome sequences of Legionella strains. The alignment was visualized in this heatmap:

![Nucleotide MSA Heatmap](./data_images/nucleotide_msa_heatmap.png)


For this, I decided to create my own algorythm instead of `MUSCLEv5` or `Clustalw2` which were just too slow for me.

## Phylogenetic Tree

A phylogenetic tree was constructed to visualize the evolutionary relationships among the Legionella strains:

![Phylogenetic Tree](./data_images/phylogenetic_tree.png)

## Visualization of Genetic Variations

To visualize genetic variations, I created another heatmap:

![Legionella Strains Heatmap](./data_images/vis_legionella_strains_heatmap.png)

We also generated a nucleotide diversity plot to quantify the level of genetic diversity:

![Nucleotide Diversity Plot](./data_images/vis_nucleotide_diversity_plot.png)

## Conclusion

The genetic diversity analysis of Legionella strains has provided valuable insights into their evolutionary history and genetic variations.

There will be some more tools added and taken away, to make this project even better.