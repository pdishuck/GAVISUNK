include: "gatherSplits.smk"
rule combine_asm_haps:
  input:
    hap1_asm = lambda wildcards: manifest_df.at[wildcards.sample, "hap1_asm"],
    hap2_asm = lambda wildcards: manifest_df.at[wildcards.sample, "hap2_asm"]
  output:
    combined = temp("results/{sample}/mrsfast/ref.fa")
  resources:
    mem = 8,
    load = 100
  threads: 1
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/combine_asm_haps.log"
  shell:
    """
    zcat -f {input.hap1_asm} {input.hap2_asm} > {output.combined}
    samtools faidx {output.combined}
    """

#Get the kmers from the assembly
rule jellyfish_count:
  input:
    asm = rules.combine_asm_haps.output.combined,
  output:
    counts = "results/{sample}/db/jellyfish.counts"
  resources:
    mem = 1,
    load = 900
  threads:32
  conda:
    "../envs/viz.yaml"
  params:
    sunk_len = config["SUNK_len"]
  log:
    "logs/{sample}/jellyfish_count.log"
  shell:
    """
    jellyfish count -m {params.sunk_len} -s 10000000 -t {threads} -C -c 1 -U 1 {input.asm} -o {output.counts}
    """

#Filter to kmers seen once in the assembly i.e. SUNKs
rule define_SUNKs:
  input:
    counts = rules.jellyfish_count.output.counts
  output:
    db = "results/{sample}/db/jellyfish.db",
    fa = "results/{sample}/db/jellyfish.fa"
  resources:
    mem=1,
    load = 900
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/define_SUNKs.log"
  threads: 32
  shell:
    """
    jellyfish dump -c -t {input.counts}  | awk '{{print $1}}' > {output.db}
    awk '{{ print ">"$0"\\n"$0 }}' {output.db} > {output.fa}
    """

#Create index of assembly for mrsfast mapping
rule mrsfast_index:
  input:
    asm = rules.combine_asm_haps.output.combined
  output:
    index = temp("results/{sample}/mrsfast/ref.fa.index"),
  resources:
    mem=8, #1.1G used
    load = 100
  threads: 1
  conda:
    "../envs/viz.yaml" 
  log:
    "logs/{sample}/mrsfast_index.log"    
  shell:
    """
    mrsfast --ws 14 --index {input.asm}
    """

#Map SUNKs back to assembly
rule mrsfast_search:
  input:
    ref = rules.combine_asm_haps.output.combined,
    index = rules.mrsfast_index.output.index,
    db = rules.define_SUNKs.output.fa
  output:
    sam = temp("results/{sample}/mrsfast/kmer_map.sam"),
    bam = temp("results/{sample}/mrsfast/kmer_map.bam"),
  resources:
      mem=2, #73 G VMS
      load=900
  threads: 16
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/mrsfast_search.log"  
  shell:
    """
    mrsfast --search {input.ref} --threads {threads} --mem 16 --seq {input.db} -o {output.sam} --disable-nohits -e 0 
    samtools sort -@ {threads} {output.sam} -o {output.bam} 
    """

#Convert mrsfast bam file to bed file
rule bed_convert:
  input:
    sam = rules.mrsfast_search.output.sam,
    bam = rules.mrsfast_search.output.bam,
  output:
    bed = temp("results/{sample}/mrsfast/kmer_map.bed"),
    bedmerge = temp("results/{sample}/mrsfast/kmer_map.merge.bed"),
    locs = "results/{sample}/mrsfast/kmer.loc",
  resources:
      mem=16,
      load=100,
  threads: 2
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/bed_convert.log"  
  shell:
    """
    bedtools bamtobed -i {input.bam} > {output.bed} && \
    bedtools merge -i {output.bed} > {output.bedmerge} && \
    bedtools intersect -a {output.bed} -b {output.bedmerge} -wo | cut -f 1,2,4,8 > {output.locs}
    """