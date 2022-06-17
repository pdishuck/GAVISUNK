# GAVISUNK
Genome Assembly Validation via Inter-SUNK distances in nanopore reads

Setup source files in config.yaml and ont.tsv

config.yaml requires:
- a manifest of ONT data for validation "ONT.tsv",
- and the kmer-length to use for validation "SUNK_len" (default: 20)
- bed: .bed file of regions to visualize

The ONT manifest is a .tsv file "ont.tsv" with the following columns:
- sample: must be unique for each file
- hap[1,2]\_ONT : location of ONT reads (gzipped or raw .fa or .fq)
- hap[1,2]\_asm : location of haplotype assembly (FASTA)
- hap[1,2]\_paf : PAF format alignment of haplotype to reference (optional)

I recommend a single input file for each haplotype, and to haplotype phase with canu and parental illumina (not currently incorporated into this pipeline)

Run snakefile locally 

Troubleshooting tips: 
To get snakemake conda envs to work correctly, you may need to deactivate your local conda env ($PATH issues as in https://github.com/snakemake/snakemake/issues/883)
