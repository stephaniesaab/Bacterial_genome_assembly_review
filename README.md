# Bacterial_genome_assembly_review
Assignment 1 for Genomic Methods of Bioinformatics, assembling a bacterial genome for Salmonella enterica

# Assignment 1

## Introduction

## Proposed Methods

## References


# DRAFT and NOTES
Challenges with sequence alignment from: Simpson & Pop 2015
- Genome assembly is the computational process of assembly together fragments of a genome based on sequenced DNA reads, in order to make one contiguous sequence.
- Challenge is that it's very complex: shortest common superstring problem --> computationally not viable and may have exponential number of solutions
- Challenge with complex and vertebrate genomes: many repeats that are almost identical, strays from parsimony (principle that the simplest explanation with the fewest evolutionary changes is most likely correct to explain observed data)
- 

Proposed method 
QC with FastQC --> Look at N50, CG content, 
Long read assembly with FLYE
Alignment with Minimap2
