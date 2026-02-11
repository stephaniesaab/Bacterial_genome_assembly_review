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
* [3. Results & Discussion](#3-results--discussion)
* [3.1 Assembly Summary](#31-assembly-summary)
* [3.2 Variant Calling](#32-variant-calling)
* [3.3 Conclusion](#33-conclusion)
* [4. References](#4-references)


## 1. Introduction 

### 1.1 Biological Background 

*Salmonella enterica* is a gram-negative bacterium and a major cause of food-borne illness in humans (Haendiges et al., 2019). *S. enterica* has been extensively studied at the genomic level and has a circular genome of about 4.8 Mbp (McClelland et al., 2001). Genome assembly is the computational process of piecing together fragmented DNA sequences from a genome into a complete, contiguous sequence (Simpson & Pop, 2015). It maps an organism's genetic blueprint, which enables identification of nucleotide and structural variants relevant to microbial genomics and public health (Simpson & Pop, 2015). Long-read sequencing (LRS) by Oxford Nanopore Technologies (ONT) provides longer contigs of sequence reads of tens of kilobases or longer but it is more error prone at about 85-95% accuracy (Boostrom et al., 2022, DeCoster et al., 2018). Here, we used long reads generated with R10 chemistry (expected Q20+, N50: 5-15 kb) and propose a methodology to computationally assemble the *S. enterica* genome. 

### 1.2 Challenges of genome assembly 

Genome assembly of bacteria is complicated by several challenges. One difficulty is repeat sequences as these create ambiguities when organizing reads based on overlapping segments. If the reads are shorter than the repeats, then there are an exponential number of possible genomes, and the problem is not computationally tractable (Simpson & Pop, 2015). LRS can address this problem to make genome assembly feasible with overlapping-based algorithms (Simpson & Pop, 2015). Also, as more sequence data become available, there is a greater need for assemblers to be able to efficiently handle big data and scale with the volume of data produced. With LRS data, there is a growing need for assemblers to tolerate sequencing errors as well (Simpson & Pop, 2015). 

### 1.3 Comparison 

When choosing a genome assembly tool, we must consider several factors including the nature of the sequencing data, the assembly algorithm, and its performance on previous bacterial genomes. While hybrid assemblers that use both short-read sequencing and LRS have been shown to perform better than LRS alone, the available data necessitates the use of a long-read assembler (Boostrom et al., 2022). Flye uses repeat graphs, which are built on approximate sequence matches so they can tolerate the higher noise of LRS (Kolmogorov et al., 2019). Boostrom et al. (2022) compared several LRS assemblers (Canu, Raven, Miniasm, Flye) for *de novo* bacterial genome assembly. They determined that Flye produced most similar metrics to the hybrid assembly methods, which performed the best. Flye produced the smallest genomic distance (0.000769-0.000930), indicating that it had the most accurate assembly to the reference genome. Additionally, Flye produced the assemblies with the fewest single-nucleotide variants (SNVs). Overall, Flye is a robust assembler well-suited for long-read sequences produced by ONT (Boostrom et al., 2022). 

Since LRS has a higher per-base error rate compared to short-read sequencing, it necessitates post-assembly polishing (DeCoster et al., 2018, Boostrom et al., 2022). Additionally, aligning a *de novo* assembly to a reference genome allows for validation of the assembly accuracy and identification of SNVs and structural variants. However, the alignment may introduce reference bias as *S. enterica* comprises multiple serotypes, and the available reference genome may not be an exact match for the sequenced sample (McClelland et al., 2001). A mismatch may complicate the interpretation of variants called as they may reflect different serotypes, rather than SNVs or assembly/sequencing errors. Despite these challenges, the use of high-quality ONT R10 chemistry and a robust long-read assembler such as Flye provides an effective balance between contiguity and accuracy for bacterial genome assembly of the available *S. enterica* long-read sequences (Boostrom et al., 2022, Simpson & Pop, 2015).

## 2. Proposed methods 

The bioinformatics pipeline completed for this project has the steps: quality control (QC) of ONT sequences using NanoPlot (version 1.46.2, DeCoster et al., 2018), filtering and trimming using NanoFilt (version 2.8.0, DeCoster et al., 2018), *de novo* genome assembly using Flye (version 2.9.6, Kolmogorov et al., 2019), two rounds of assembly polishing using Racon (version 1.5.0, Vaser et al., 2017) and Medaka (version 2.2.0, Oxford Nanopore Technologies Ltd, 2018), reference alignment using Minimap2 (version 2.30, Li, 2018), variant calling with Bcftools (version 1.23, Danecek et al., 2021), indexing with Samtools (Li et al., 2009) and visualization using IGV (desktop application, Robinson et al., 2011). The assembly process is based on those outlined by Kumar et al.(2025). 

### 2.1 QC 

QC of ONT long-read sequencing data was completed with NanoPlot (version 1.46.2) to assess the data quality and whether further processing is necessary. The metrics investigated include read length distribution, N50, total number of reads, and mean read quality score (Q score). NanoPlot is appropriate for this data as high quality reads are expected (Q20) so a more complex tool (e.g. LongQC, Fukasawa et al. 2020) is unnecessary. Read lengths histograms were used to confirm the presence of long reads, identify excess short fragments, while a bivariate plot of read length verses mean read quality will be used to detect low-quality read clusters. The --loglength command was used to enable logarithmic scaling of read lengths in the generated plots. From the report given by NanoPlot, it was found that there were several long low-quality reads (Q <= 8), and several short reads (< 2kb). 

### 2.2 Read Processing 

Data processing was done in NanoFilt if indicated by the QC results. Since there were found to be sufficient reads of long length and high quality (600-700Mb) it was determined that filtering could have stricter options for bacterial genomes, so that less polishing would be necessary after assembly. Filtering used the flags q –15 to select for reads of high quality (Q > = 15) and –l 3000 to select for reads of sufficient length (>= 3kb), as this maintains sufficient reads and coverage for Flye assembly (< 30 – 50x coverage) (Kolmogorov et al., 2019). Aggressive trimming was avoided to prevent unnecessary loss of coverage. Overall, 519140 reads were kept from the original 784124 reads. A quality control was done on the filtered reads with NanoPlot and showed that 628 Mb were kept from the original 809 Mb, and the mean quality score increased from 18.4 to 21.8.  

### 2.3 *De novo* genome assembly 

*De novo* assembly of the long-read sequencing data was performed using Flye (version 2.9.6), a long-read assembler optimized for Oxford Nanopore long-read data. Assembly was conducted using default parameters appropriate for bacterial genomes, producing a draft consensus assembly. The flags included ```--nano-hq``` as, based on the Flye manual, R10 ONT parameters (3% error) are default for --nano-hq (Kolmogorove et al., 2019). Additional parameters included the genome size at about 5Mb (```--genome-size 5m```), the output directory, and the threads used (8 threads). 

### 2.4 Assembly Polishing 

The draft assembly was polished using two rounds of Racon (version 1.5.0) and one round of Medaka (version 2.2.0) to improve base-level accuracy by correcting sequencing errors associated with LRS. Firstly, the filtered reads mapped to the Flye-generated assembly using Minimap2 (version 2.30). Specific parameters were set to make the output a SAM file, to use a preset optimized for ONT reads, and specify the threads to use (half of what is available on the CPU). This full command was: ```minimap2 -t 6 -ax map-ont flye_out/assembly.fasta filtered.fastq.gz > flye.sam```

For the assembly polishing with racon, the threads were increased to the entire CPU availability (12 threads). Default parameters were used for racon polishing, and the input files included the filtered reads and the SAM file produced by minimap2. This was done twice to produce a final polished genome racon2.fasta.  

For polishing with Medaka (one run) a bacterial-specific Medaka consensus, r1041_e82_400bps_sup_v5.2.0, model compatible with R10 ONT chemistry was be used to optimize accuracy for the S. enterica genome (Oxford Nanopore Technologies Ltd., 2018). This was chosen based on the default tools available for consensus genomes, which were listed by running ```medaka tools list_models```. All threads were used (12) and the input reads were set to the filtered reads, the draft was racon2.fasta. Once medaka finished, it was checked that it produced a genome assembly of appropriate size (approximately 5Mb), it was found that it produced an assembly of 4.9 Mb. 
### 2.5 Reference alignment & Visualization 

A reference genome and genome annotation file for *S. Enterica* was downloaded from NCBI (Accession GCF_000006945.2_ ASM694v2). The polished assembly was aligned to the reference genome using minimap2 (version 2.30) as previously described. Variant calling was performed using Bcftools (version 1.23) to identify differences between the assembled genome and the reference. For minimap2, default parameters were used for alignment of entire genomes, and the ```asm5``` flag as they are high-quality alignments, this flag aligns assemblies with less than 5% divergence (Li, 2018). The output was sorted by genomic location and indexed using samtools. Variant calling was performed with Bcftools to identify genetic variants to the reference genome. The parameters set the ploidy as haploid (1) as it is a bacterial genome, multiallelic variant calling, and compressed the file into a VCF for visualization. The assembly structure, alignments, and identified variants were visualized in IGV to assess assembly the assembly quality and validate variant calls. To identify high-quality variants only, the variants were filtered using bcftools to exclude low-quality calls (low Q and low depth) using ```view -i 'QUAL>30 && DP>10'``` to filter for variants with high quality and high depth, to distinguish true variants from potential sequencing errors in one or a few reads.

## 3. Results & Discussion
### 3.1 Assembly Summary

In this study, Oxford Nanopore (R10 chemistry) long-read sequencing raw reads were assembled to create a de novo genome assembly of *Salmonella enterica*. The assembly process made use of several tools and followed steps to quality control the raw reads, filter the raw reads, assemble the genome, polish the assembly (two rounds), align the genome to a reference, conduct variant calling of raw reads against the reference, and graph the assembly against the reference genome to assess its coverage and depth. The assembly produced three contigs that make up a total size of 4.9Mb, which is similar in size to reported *S. enterica* genomes (5.1 Mb) (Stevens et al., 2017). The results of this genome analysis include a de novo genome assembly, a variant calling that produced a number of SNPs, indels, transitions and transversions. One gene of interest identified to have many SNPs was the *STM2913* gene, which encodes a putative permease protein similar to the *E. coli* gluconate permease (McClelland et al., 2001). The SNPs highlighted include missense mutations that may alter the bacteria's ability to uptake gluconate from the environment and hinder energy metabolism under certain environmental conditions. 

A sample representation of the aligned reads in contig 2 with the reference genome are shown in IGV visualization in *Figure 1*, highlighting the adequate coverage of the assembly against the reference genome. The BAM files were validated in IGV to visualize the assembly. IGV highlights the efficacy and accuracy of the assembly. The upper reads_vs_ref track highlights the alignment of the raw reads against a section of the reference genome, the asm_vs_ref track shows the assembled genome against the same sample. The *de novo* assembly shows only one contig (contig 2), which aligns well with the reference genome and the coverage track shows that there is sufficiently uniform coverage across the genome. There is a slight dip in coverage that aligns with some gaps in alignment amongst the reads. Overall, the assembled contigs show high collinearity with the reference, indicating a contiguous and accurate genome assembly. 
<img width="1536" height="806" alt="igv_snapshot" src="https://github.com/user-attachments/assets/afa9fed2-6e2c-4ef9-83db-8133463722d2" />
**Figure 1. IGV Visualization of sample region of raw Oxford Nanopore reads (R10) aligned to the polished *Salmonella enterica* *de novo* genome assembly.** The sample contig (18kb) includes a coverage track showing multiple reads aligning across the large contig of the reference genome (NC_003197.2 contig, coordinates 2,775,187-2,793,401), indicating a high depth. There is a middle region with more white space indicating gaps in the reads alignments, which also aligns with a decline in read coverage as seen in the coverage track. Overall, the read alignment depth is sufficient across the majority of the contig and shows a high-quality assembly. This contig represents the large circular DNA of the bacteria. 

### 3.2 Variant Calling
The variant calling produced 10149 SNPs, of which 5081 transitions and 5068 transversions, and 76 indels. Variants plotted on IGV for selection of important variants were filtered for quality variants that had high quality (quality >30) and read depth (depth > 10) to be more confident that the variants identified were true mutations and not due to errors in LRS in one or a few reads.

One limitation with variant calling was that it could not filter for important variants that may have biological significance or variants in important genes, as bcftools requires a very specific, strict format for parsing gene, transcript, and exon hierarchies that is different from standard NCBI GFF outputs. Therefore, the selection of important variants was done by hand. In the future, a post-processing tool for variant call format (VCF) files such as Ensemble Variant Effector Predictor (VEP) can be useful for annotating the functional consequence of variants on genes (McLaren et al., 2016). 

A cluster of five SNPs were identified in the coding region of the *STM2913* gene. The genome annotation indicated that this gene encodes a putative permease of 498 amino acids, similar to the *E. coli* putative membrane transport protein (AAC75782.1). The mutations were found to be at amino acid positions: 235F>I, 241K>E, 367Q, 429T>A, 456H>R. All but one were an adenine to guanine transition, while that at position 235 was a threonine to adenine transversion. Four of these mutations are nonsynonymous and thus likely have a biological impact on the gene. The phenylalanine to isoleucine mutation is conservative as they are both bulky non-polar amino acids, and both are frequently found in transmembrane helices (De Marothy et al., 2015). The threonine to alanine mutation likely affects the protein structure as threonine is known to contribute hydrogen bonding to transmembrane helices, while alanine cannot contribute to hydrogen bonding (De Marothy et al., 2015). Lysine and arginine are both charged amino acids that are less likely to be found in transmembrane helices, they may alter dipole interactions and hinder the protein from embedding into phospholipids. The distribution and predicted effects of SNPs within STM2913 are visualized in *Figure 2* using a lollipop plot. 

Studies in *Pseudomonas stutzeri* produced and purified STM2913 from *S. enterica* and produced a Gluconate:H+ symporter, which is a member of the GntP family of transporter. These transporters utilize the electrochemical gradiant of protons to facilitate the uptake of gluconate into the cell (Sommer et al., 2017). Studies in *E. coli* explain how this transporter can uptake gluconate, a salt of gluconic acid derived from glucose, which can be used in the Entner-Doudoroff pathway as an alternative pathway to glycolisis, in which glucose is converted into pyruvate by bacteria (Eisenberg & Dobrogosz, 1967). The use of gluconate is dependent on the environmental conditions, growth rate, oxygen availability, and intracellular metabolite levels. Variants affecting the transmembrane transporter of gluconate likely alter the carbon utilization efficiency of S. enterica under certain environmental conditions. 
<img width="1443" height="675" alt="SNP_lollipop" src="https://github.com/user-attachments/assets/350f0591-ed10-4d39-b04b-52a306aae268" />
**Figure 2. Lollipop plot showing the distribution of SNPs across the STM2913 coding region of *Salmonella enterica*.** Each lollipop represents a single SNP positioned along the protein sequence. Although multiple variants occur within the gene, colours represent whether mutations are missense (nonsynonymous) or synonymous (no amino acid change), four are nonsynonymous and likely to affect protein function, one is silent and unlikely to affect the protein structure. The labels above each dot indicate the amino acid position of the SNP, which were calculated by hand based on the genome visualization on IGV.

One other gene that had multiple variants and may have altered function in the sample is *STM1009* which encodes Gifsy-2 prophage exodeoxyribonuclease. 

### 3.3 Conclusion
Overall, several variants in bases were found in the raw reads of *S. enterica*. Some of these variants likely contribute to the efficiency of the bacteria to utilize alternative carbon sources like gluconate under specific environmental conditions. Future assessments in variant calling should utilize alternative tools that are able to integrate genome annotations to filter variants for likely biologically-relevant SNPs and indels. 

## 4. References 

Boostrom, I., Portal, E. A. R., Spiller, O. B., Walsh, T. R., & Sands, K. (2022). Comparing Long-Read Assemblers to Explore the Potential of a Sustainable Low-Cost, Low-Infrastructure Approach to Sequence Antimicrobial Resistant Bacteria With Oxford Nanopore Sequencing. *Frontiers in Microbiology*, 13. https://doi.org/10.3389/fmicb.2022.796465 

Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008 

De Coster, W., D’Hert, S., Schultz, D. T., Cruts, M., & Van Broeckhoven, C. (2018). NanoPack: Visualizing and processing long-read sequencing data. *Bioinformatics*, 34(15), 2666–2669. https://doi.org/10.1093/bioinformatics/bty149 

De Marothy, M. T., & Elofsson, A. (2015). Marginally hydrophobic transmembrane α-helices shaping membrane protein folding. *Protein Science : A Publication of the Protein Society*, 24(7), 1057–1074. https://doi.org/10.1002/pro.2698

Eisenberg, R. C., & Dobrogosz, W. J. (1967). Gluconate Metabolism in Escherichia coli1. *Journal of Bacteriology*, 93(3), 941–949. https://doi.org/10.1128/jb.93.3.941-949.1967

Fukasawa, Y., Ermini, L., Wang, H., Carty, K., & Cheung, M.-S. (2020). LongQC: A Quality Control Tool for Third Generation Sequencing Long Read Data. *G3: Genes|Genomes|Genetics*, 10(4), 1193–1196. https://doi.org/10.1534/g3.119.400864 

Haendiges, J., Gonzalez-Escalona, N., Miller, J. D., & Hoffmann, M. (2019). Complete Genome Sequences of Four *Salmonella enterica* Strains Associated with Pistachios Assembled Using a Combination of Short- and Long-Read Sequencing. *Microbiology Resource Announcements*, 8(38), 10.1128/mra.00975-19. https://doi.org/10.1128/mra.00975-19 

Kolmogorov, M., Yuan, J., Lin, Y., & Pevzner, P. A. (2019). Assembly of long, error-prone reads using repeat graphs. *Nature Biotechnology*, 37(5), 540–546. https://doi.org/10.1038/s41587-019-0072-8 

Kumar, A., Häggblom, M. M., & Kerkhof, L. J. (2025). A Step-by-Step Guide to Sequencing and Assembly of Complete Bacterial Genomes Using the Oxford Nanopore MinION. In V. J. Carabetta & O. Akintunde (Eds.), High Throughput Gene Screening: Methods and Protocols (pp. 31–43). *Springer US*. https://doi.org/10.1007/978-1-0716-4192-7_2 

Li, H. (2018). Minimap2: Pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191 

Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., & Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352 

McClelland, M., Sanderson, K. E., Spieth, J., Clifton, S. W., Latreille, P., Courtney, L., Porwollik, S., Ali, J., Dante, M., Du, F., Hou, S., Layman, D., Leonard, S., Nguyen, C., Scott, K., Holmes, A., Grewal, N., Mulvaney, E., Ryan, E., … Wilson, R. K. (2001). Complete genome sequence of *Salmonella enterica* serovar Typhimurium LT2. *Nature*, 413(6858), 852–856. https://doi.org/10.1038/35101614 

McLaren, W., Gil, L., Hunt, S. E., Riat, H. S., Ritchie, G. R. S., Thormann, A., Flicek, P., & Cunningham, F. (2016). The Ensembl Variant Effect Predictor. *Genome Biology*, 17(1), 122. https://doi.org/10.1186/s13059-016-0974-4

Oxford Nanopore Technologies Ltd. 2018. Medaka: Sequence correction provided by ONT Research. https://github.com/nanoporetech/medaka.

Robinson, J. T., Thorvaldsdóttir, H., Winckler, W., Guttman, M., Lander, E. S., Getz, G., & Mesirov, J. P. (2011). Integrative genomics viewer. *Nature Biotechnology*, 29(1), 24–26. https://doi.org/10.1038/nbt.1754 

Simpson, J. T., & Pop, M. (2015). The Theory and Practice of Genome Sequence Assembly. *Annual Review of Genomics and Human Genetics*, 16(1), 153–172. https://doi.org/10.1146/annurev-genom-090314-050032 

Sommer, M., Xie, H., & Michel, H. (2017). Pseudomonas stutzeri as an alternative host for membrane proteins. *Microbial Cell Factories*, 16, 157. https://doi.org/10.1186/s12934-017-0771-0 

Vaser, R., Sović, I., Nagarajan, N., & Šikić, M. (2017). Fast and accurate de novo genome assembly from long uncorrected reads. *Genome Research*, 27(5), 737–746. https://doi.org/10.1101/gr.214270.116 
