# GAVISUNK
Genome Assembly Validation via Inter-SUNK distances in nanopore reads

(WIP)

Setup source files in config.yaml and ont.tsv

config.yaml requires:
- a haplotype-phase HiFi assembly "asm",
- a manifest of ONT data for validation "ONT_manifest",
- and the kmer-length to use for validation "SUNK_len"

The ONT manifest is a .tsv file "ont.tsv" with the following columns:
- sampname: must be unique for each file,
- path: location of ONT reads (gzipped or raw .fa or .fq)
- hap: haplotype identifier (not currently used)
- bed: .bed file of regions to visualize with these reads

I recommend a single input file for each haplotype, and to haplotype phase with canu and parental illumina (not currently incorporated into this pipeline)

For now, run Snakefile locally on high memory machine

Troubleshooting: 
may need to deactivate local conda env to get snakemake conda envs to work correctly ($PATH issues as in https://github.com/snakemake/snakemake/issues/883)
