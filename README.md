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

More Notes from readings: - What is genome assembly?
    
    - “[Genome assembly](https://www.google.com/search?q=Genome+assembly&oq=genome+assembly&gs_lcrp=EgZjaHJvbWUyCQgAEEUYORiABDIHCAEQABiABDIHCAIQABiABDIHCAMQABiABDIHCAQQABiABDIGCAUQRRg8MgYIBhBFGDwyBggHEEUYPNIBCDM5MjJqMGo3qAIAsAIA&sourceid=chrome&ie=UTF-8&ved=2ahUKEwjjmc28546SAxVVm4kEHV0dJvIQgK4QegYIAQgAEAM) is **the computational process of piecing together millions of short, fragmented DNA sequences (reads) from a genome into a complete, contiguous sequence”**
        

- `Why is genome assembly hard?` 
    
    - Reconstructing an entire genome means gluing together small fragments, in a particular order
        
- `Why are repeats problematic?` 
    
    - Repeats make it difficult to align sequences
        
    - Need the reads to be longer than the repeats to solve the assembly, if the reads are shorter than the reads than there are exponential number of genomes
        
    - Greedy strategy is based on locally optimal joining, so it can’t handle repeated genomic regions (doesn’t consider the whole genome order at once). May collapse repeats
        

- `Why long reads help (and hurt)?` 
    

- `What metrics matter (N50, accuracy, contiguity)?` 
    

Challenge of assembly: Ukkonen and others studied the computational complexity of the assembly problem, formalized as an instance of the shortest common superstring problem (61), the problem of finding the shortest string that encompasses all the reads as substrings. This Occam’s razor formulation assumes that the genome being reconstructed is the most parsimonious explanation of the set of reads. Their work showed that genome sequence assembly is computationally intractablefinding the correct solution may require exploring an exponential number of possible solutions” ([Simpson and Pop, 2015, p. 155](zotero://select/library/items/Y3IVKHPE))

“complex (e.g., vertebrate) genomes: repeats. Most genomes contain DNA segments that are repeated in nearly identical form, thus straying from parsimony” ([Simpson and Pop, 2015, p. 155](zotero://select/library/items/Y3IVKHPE)) ([pdf](zotero://open-pdf/library/items/QKT9JUUS?page=3&annotation=NUDJPGPT))

Parsimony in bio: “In biology, parsimony is **the principle that the simplest explanation requiring the fewest evolutionary changes (like mutations or trait developments) to explain observed data is the most likely to be correct**,”

“An explanation for this all-too-familiar gap between theory and practice is that theoretical intractability results are based on worst-case scenarios, which rarely occur in practice.” ([Simpson and Pop, 2015, p. 155](zotero://select/library/items/Y3IVKHPE)) ([pdf](zotero://open-pdf/library/items/QKT9JUUS?page=3&annotation=T6QTNCGH)) \--> Depends on ratio between size of sequence reads and size of repeats (Need reads to be longer than repeats to make a genome)

Outline below: 
# Bacterial_genome_assembly_review
Assignment 1 for Genomic Methods of Bioinformatics, assembling a bacterial genome for Salmonella enterica

# Assignment 1
# Table of contents

## 1. Introduction
### 1.1 Biological Background
Salmonella enterica is a Gram-negative bacterium that can cause food-borne illness in humans (Haendiges et al. 2019). They cause a high number of infections worldwide, with an estimated 1.2 million illnesses (Haendiges et al. 2019). Domestic animals can act as reservoirs for the food-borne spread of these pathogens. These put major strains on the population and economy as it was estimated that costs of food-borne diseases in the US range from 4.8 to 23 billion dollars, with salmonella being a major contributor. It has a circular genome of about 4.8 Mbp (McClelland, 2001). 
Genome assembly is the computational process of piecing together millions of short, fragmented DNA sequences (reads) from a genome into a complete, contiguous sequence (Simpson & Pop, 2015). It is a foundational tool for biological research and understanding the evolution of species, their physiological processes, and gives insights into the genetic components of diseases. An (almost) complete genome assembly can give a map of an organism's genetic build, which can serve as a foundation for detecting variants in nucleotides, genome structure, which are key for progression of biotechnology, evolutionary biology, precision medicine, and more (CD Genomics, n.d.). Long-read sequencing (LRS) by Oxford Nanopore Technologies (ONT) provides longer contigs of sequence reads of tens of kilobases or longer, but is more error prone (Boostrom et al., 2022). Here, we used LRS by ONT using R10 chemistry (expected Q20+, N50: 5-15 kb) and propose a methodology to computationally assembly the S. enterica genome.


- Salmonella enterica is gram negative bacteria
  - Genome size about 4.8 Mbp
  - Pathogenic, causes human gastroenteritis
     Important in public health, found in foods and domestic animals
- Why genome assembly is useful?
  -  Can detect variants, SNPs, can detect structural diff, genomic distances, assess quality of sequencing
- LRS used to sequence genome (long read sequencing) Oxford Nanopore, R10 chemistry, expected Q20+ (N50: 5-15kb). Oxford nanopore enables near-complete assemblies
### 1.2 Challenges of genome assembly
Genome assembly involves reconstructing an entire genome means gluing together small fragments, in a particular order. There are many challenges with genome assembly. One challenge is repeat sequences as these make it difficult to align sequences based on overlapping segments. For algorithms to align repeats, the reads need to be longer than the repeats to solve the assembly, if the reads are shorter than the reads than there are exponential number of genomes possible and the problem is not computationally tractable (Simpson & Pop 2015). Most of the developments in assembly algorithms addressed challenges with short reads, but long reads can address these problems, especially for small genomes like bacterial, and make complex algorithms unnecessary. Analysis of genomic variants is better with longer reads, long sequence reads make assembling a genome sequence easier. Challenge High-error reads associated with third-generation sequencing technologies (longer reads) --> longer reads but more errors. Another challenge with assembly is metagenomics and that there may be mixing of genomes in a sample such as microbiome data, or a mix of tumor cells. Also, as more sequence data become available, especially for larger complex genomes become available, such as vertebrates, there is a greater need for assemblers to be able to efficiently handle Big Data and scale with the amount of sequence reads produced. 

### 1.3 Comparison (justify approach)
When choosing a genome assembly tool, we must consider several factors including the nature of the sequencing data, the different algorithms used in various assemblers, and their performance on previous bacterial genomes. 
Firstly, the data we are working with is LRS from ONT. Therefore, while hybrid assemblers that use both short-read sequencing and LRS have been shown to perform better than LRS alone, we must choose from a LRS assembler (Boostrum et al., 2022).
Different genome assemblers use various algorithms which may be overlapping-based. 
Users use multiple assemblers and select the best according to multiple metrics, not just one (e.g. N50) 
 - Tradeoffs between read-length vs accuracy, need read length > repeat lengths, contiguity vs correctness
 - Assembly strategies: De Bruijn graphs rely on exact overlaps, repeat graphs allow for some error so they tolerate more --> used in Flye (cite 2019 paper). Overlap based strategies
   - Cite paper comparing Flye, Canu, Raven
Boostrum et al. (2022) compared several LRS assemblers to assemble LRS produced by ONT for several bacteria. They determined that Flye produced most similar metrics to the hybrid assembly methods, which performed the best. Flye produced the smallest genomic distance, indicating it had the most accurate assembly to the reference genome (0.000769-0.000930). Additionally, most of the LR assemblies they produced with < 100 SNVs were Flye based, and flye produced the fewest number of assembles with > 2000 SNVs. Overall, based on the 3Cs of assembly quality (contiguity, completeness, correctness), Flye performed the best when compared to other LRS assemblers (Canu, Raven, Miniasm). 

## 2. Proposed methods
### 2.1 Data & QC


 - ONT, R10, FastQ reads
 - QC using tolls like Nanoplot
 - Metrics assessed: Read length, N50, select reads with length > X kb, quality distribution
### 2.2 Assembly & Polishing
 - LRS assembly with Flye version 2.9.6 (latest)
 - Parameters --> default
 - R10 ONT parameters (3% error) set as default for --nano-hq --> see Flye github
### 2.3 Reference alignment & Visualization
 - Reference genome downloaded from ___ (NCBI)
 - Alignment with minimap2
 - Visualization with ____ (IGV)


## 3. References (check Zotero)

