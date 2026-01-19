# Bacterial_genome_assembly_review
Assignment 1 for Genomic Methods of Bioinformatics, assembling a bacterial genome for *Salmonella enterica*, aligning it to a reference genome, and visualizing variants.
## Table of Contents
* [1.Introduction](#1-introduction)
* [1.1 Biological Background](#11-biological-background)
* [1.2 Challenges of genome assembly](#12-challenges-of-genome-assembly)
* [1.3 Comparison](#13-comparison)
* [2. Proposed methods](#2-proposed-methods)
* [2.1 QC](#21-qc)
* [2.2 Read Processing](#22-read-processing)
* [2.3 *De novo* genome assembly](#23-de-novo-genome-assembly)
* [2.4 Assembly polishing](#24-assembly-polishing)
* [2.5 Reference alignment & Visualization](#25-reference-alignment--visualization)
* [3. References](#3-references)


## 1. Introduction 

### 1.1 Biological Background 

*Salmonella enterica* is a gram-negative bacterium and a major cause of food-borne illness in humans (Haendiges et al., 2019). *S. enterica* has been extensively studied at the genomic level and has a circular genome of about 4.8 Mbp (McClelland et al., 2001). Genome assembly is the computational process of piecing together fragmented DNA sequences from a genome into a complete, contiguous sequence (Simpson & Pop, 2015). It maps an organism's genetic blueprint, which enables identification of nucleotide and structural variants relevant to microbial genomics and public health (Simpson & Pop, 2015). Long-read sequencing (LRS) by Oxford Nanopore Technologies (ONT) provides longer contigs of sequence reads of tens of kilobases or longer but it is more error prone at about 85-95% accuracy (Boostrom et al., 2022, DeCoster et al., 2018). Here, we used long reads generated with R10 chemistry (expected Q20+, N50: 5-15 kb) and propose a methodology to computationally assemble the *S. enterica* genome. 

### 1.2 Challenges of genome assembly 

Genome assembly of bacteria is complicated by several challenges. One difficulty is repeat sequences as these create ambiguities when organizing reads based on overlapping segments. If the reads are shorter than the repeats, then there are an exponential number of possible genomes, and the problem is not computationally tractable (Simpson & Pop 2015). LRS can address this problem to make genome assembly feasible with overlapping-based algorithms (Simpson & Pop, 2015). Also, as more sequence data become available, there is a greater need for assemblers to be able to efficiently handle big data and scale with the volume of data produced. With LRS data, there is a growing need for assemblers to tolerate sequencing errors as well (Simpson & Pop, 2015). 

### 1.3 Comparison 

When choosing a genome assembly tool, we must consider several factors including the nature of the sequencing data, the assembly algorithm, and its performance on previous bacterial genomes. While hybrid assemblers that use both short-read sequencing and LRS have been shown to perform better than LRS alone, the available data necessitates the use of a long-read assembler (Boostrom et al., 2022). Flye uses repeat graphs, which are built on approximate sequence matches so they can tolerate the higher noise of LRS (Kolmogorov et al., 2019). Boostrom et al. (2022) compared several LRS assemblers (Canu, Raven, Miniasm, Flye) for *de novo* bacterial genome assembly. They determined that Flye produced most similar metrics to the hybrid assembly methods, which performed the best. Flye produced the smallest genomic distance (0.000769-0.000930), indicating that it had the most accurate assembly to the reference genome. Additionally, Flye produced the assemblies with the fewest single-nucleotide variants (SNVs). Overall, Flye is a robust assembler well-suited for long-read sequences produced by ONT (Boostrom et al., 2022). 

While LRS enables highly contiguous assemblies and improves the resolution of repetitive regions, it is associated with a higher per-base error rate compared to short-read sequencing (DeCoster et al., 2018). This necessitates post-assembly polishing to address sequencing errors (Boostrom et al., 2022). Additionally, aligning a *de novo* assembly to a reference genome allows for validation of the assembly accuracy and identification of SNVs and structural variants. However, the alignment may introduce reference bias as *S. enterica* comprises multiple serotypes, and the available reference genome may not be an exact match for the sequenced sample (McClelland et al., 2001). A mismatch may complicate interpretation of observed variants as differences may reflect different serotypes, rather than SNVs or assembly/sequencing errors. The use of high-quality ONT R10 chemistry and a robust long-read assembler such as Flye provides an effective balance between contiguity and accuracy for bacterial genome assembly of the available *S. enterica* long-read sequences (Boostrom et al., 2022, Simpson & Pop, 2015).

## 2. Proposed methods 

The proposed bioinformatics pipeline for this project is: quality control (QC) of ONT sequences using NanoPlot (version 1.46.2, DeCoster et al., 2018), optional filtering and trimming using NanoFilt (version 2.8.0, DeCoster et al., 2018), *de novo* genome assembly using Flye (version 2.9.6, Kolmogorov et al., 2019), assembly polishing using Medaka (version 2.2.0, Oxford Nanopore Technologies Ltd, 2018), reference alignment using Minimap2 (version 2.30, Li, 2018), variant calling with Bcftools (version 1.23, Danecek et al., 2021) and visualization using IGV (desktop application, Robinson et al., 2011).  

### 2.1 QC 

QC of ONT LRS data will be done with NanoPlot (version 1.46.2) to assess the data quality and whether further processing is necessary. Summary metrics will include read length distribution, N50, total number of reads, and mean read quality score (Q score). NanoPlot is appropriate for this data as high quality reads (+Q20) are expected, so a more complex tool (e.g. LongQC, Fukasawa et al. 2020) is unnecessary. A read length histogram will be used to confirm the presence of long reads and identify excess short fragments, while a bivariate plot of read length verses mean read quality will be used to detect low-quality read clusters.  

### 2.2 Read Processing 

Data processing will be done in NanoFilt (version 2.8.0) if indicated by the QC results. Filtering will be applied if mean read quality falls below the expected Q20, read length N50 is insufficient for Flye assembly (< 1kb), there are extremely low-quality clusters, or the total sequencing yield is too low for Flye (< 30 – 50x coverage) (Kolmogorov et al., 2019). Filtering parameters will include a minimum read length threshold (~1kb) and removal of lowest-quality reads (10%). These parameters will be used unless deviations are required based on QC results from NanoPlot. Aggressive trimming will be avoided to prevent unnecessary loss of coverage. 

### 2.3 *De novo* genome assembly 

*De novo* assembly of the long-read sequencing data will be performed using Flye (version 2.9.6), a long-read assembler optimized for Oxford Nanopore long-read data. Assembly will be conducted using default parameters appropriate for bacterial genomes, producing a draft consensus assembly. Based on the Flye manual, R10 ONT parameters (3% error) are default for ```--nano-hq``` (Kolmogorov et al., 2019). Additional parameters that will be used apart from the defaults will specify the output directory and threads used (half of the available CPU resources). 

### 2.4 Assembly Polishing 

The draft assembly will be polished using Medaka (version 2.2.0) to improve base-level accuracy by correcting sequencing errors associated with LRS. Specific parameters will be set for the input base calls, the assembly FASTA, the output directory, the bacteria flag, and the threads. A bacterial-specific Medaka consensus model compatible with R10 ONT chemistry will be used to optimize accuracy for the *S. enterica* genome (Oxford Nanopore Technologies Ltd., 2018). 

### 2.5 Reference alignment & Visualization 

A reference genome for *S. enterica* will be downloaded from NCBI (RefSeq assembly GCF_000006945.2). The polished assembly will be aligned to the reference genome using Minimap2 (version 2.30). Variant calling will be performed using Bcftools (version 1.23) to identify differences between the assembled genome and the reference. For Minimap2, default parameters will be used for alignment of entire genomes. The parameters will also include the asm5 flag as the data are high-quality alignments and this flag aligns assemblies with less than 5% divergence (Li, 2018). Variant calling will be performed with Bcftools to identify genetic variants to the reference genome. The parameters will set the ploidy as haploid (1) as it is a bacterial genome. Variant files will be compressed, sorted, and indexed for visualization. The assembly structure, alignments, and identified variants will be visualized in IGV (desktop application) to assess the assembly quality and validate variant calls. The shade-by-base quality parameter will be set to distinguish true variants from potential sequencing errors (Robinson et al., 2011). 

 

## 3. References 

Boostrom, I., Portal, E. A. R., Spiller, O. B., Walsh, T. R., & Sands, K. (2022). Comparing Long-Read Assemblers to Explore the Potential of a Sustainable Low-Cost, Low-Infrastructure Approach to Sequence Antimicrobial Resistant Bacteria With Oxford Nanopore Sequencing. Frontiers in Microbiology, 13. https://doi.org/10.3389/fmicb.2022.796465 

Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008 

De Coster, W., D’Hert, S., Schultz, D. T., Cruts, M., & Van Broeckhoven, C. (2018). NanoPack: Visualizing and processing long-read sequencing data. Bioinformatics, 34(15), 2666–2669. https://doi.org/10.1093/bioinformatics/bty149 

Fukasawa, Y., Ermini, L., Wang, H., Carty, K., & Cheung, M.-S. (2020). LongQC: A Quality Control Tool for Third Generation Sequencing Long Read Data. G3: Genes|Genomes|Genetics, 10(4), 1193–1196. https://doi.org/10.1534/g3.119.400864 

Haendiges, J., Gonzalez-Escalona, N., Miller, J. D., & Hoffmann, M. (2019). Complete Genome Sequences of Four *Salmonella enterica* Strains Associated with Pistachios Assembled Using a Combination of Short- and Long-Read Sequencing. Microbiology Resource Announcements, 8(38), 10.1128/mra.00975-19. https://doi.org/10.1128/mra.00975-19 

Kolmogorov, M., Yuan, J., Lin, Y., & Pevzner, P. A. (2019). Assembly of long, error-prone reads using repeat graphs. Nature Biotechnology, 37(5), 540–546. https://doi.org/10.1038/s41587-019-0072-8 

Li, H. (2018). Minimap2: Pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191 

McClelland, M., Sanderson, K. E., Spieth, J., Clifton, S. W., Latreille, P., Courtney, L., Porwollik, S., Ali, J., Dante, M., Du, F., Hou, S., Layman, D., Leonard, S., Nguyen, C., Scott, K., Holmes, A., Grewal, N., Mulvaney, E., Ryan, E., … Wilson, R. K. (2001). Complete genome sequence of *Salmonella enterica* serovar Typhimurium LT2. Nature, 413(6858), 852–856. https://doi.org/10.1038/35101614 

Oxford Nanopore Technologies Ltd. 2018. Medaka: Sequence correction provided by ONT Research. https://github.com/nanoporetech/medaka.

Robinson, J. T., Thorvaldsdóttir, H., Winckler, W., Guttman, M., Lander, E. S., Getz, G., & Mesirov, J. P. (2011). Integrative genomics viewer. Nature Biotechnology, 29(1), 24–26. https://doi.org/10.1038/nbt.1754 

Simpson, J. T., & Pop, M. (2015). The Theory and Practice of Genome Sequence Assembly. Annual Review of Genomics and Human Genetics, 16(1), 153–172. https://doi.org/10.1146/annurev-genom-090314-050032 
