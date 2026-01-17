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
