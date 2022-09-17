# GAVISUNK
Genome Assembly Validation via Inter-SUNK distances in nanopore reads

Setup source files in config.yaml and ont.tsv

config.yaml requires:
- a manifest of ONT data for validation "ont.tsv",
- and the kmer-length to use for validation "SUNK_len" (default: 20)

The ONT manifest is a .tsv file "ont.tsv" with the following columns:
- sample: must be unique for each file
- hap[1,2]\_ONT : location of ONT reads (gzipped or raw .fa or .fq)
- hap[1,2]\_asm : location of haplotype assembly (FASTA, must also have .fai index)
- hap[1,2]\_bed : BED format regions of assembly to visualize (optional) 
- hap[1,2]\_colotrack : BED format color track to include in visualizations (optional, no headers, columns Chrom Start End Color) 

I recommend a single input file for each haplotype, and to haplotype phase with canu and parental illumina (not currently incorporated into this pipeline)

Prerequisites: Snakemake (tested with versions 6.12.1, 7.8.2)

To run snakefile locally on the provided test cases (AMY locus of CHM13/1 pseudodiploid and HG02723), generating optional image outputs, clone this repo and execute:
```
snakemake -R --use-conda --cores 8 --configfile .test/config.yaml --resources load=1000
```
For Centromeres 1st run above step and then execute
```
snakemake -R --use-conda --cores 8 --configfile .test/config.yaml --resources load=1000 checkGaps
```

.BED results are found in the `results/[sample]/final_outs` directory

Automated visualizations of validation gaps are found in the `results/pngs/gaps/[sample] directory

For the included test cases, the following output should be generated: `results/gaps/AMY_HG02723/AMY_HG02723_hap1_AMY_h1_84861_524275.png`, corresponding to the region displayed in Figure 1 of the manuscript:
![plot](./.test/data/HG02723/AMY_HG02723_hap1_AMY_h1_84861_524275.png)




Troubleshooting tips: 
To get snakemake conda envs to work correctly, you may need to deactivate your local conda env ($PATH issues as in https://github.com/snakemake/snakemake/issues/883)
