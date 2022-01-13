# tested with Snakemake v6.7.0
import pandas as pd

configfile: "config.yaml"

ENV = os.path.dirname(workflow.snakefile) + "/env.cfg"
shell.prefix("source {ENV}; set -eo pipefail; ")

onttsv = config['ONT_manifest']
file_df = pd.read_csv(onttsv,sep="\t",)
file_df.set_index("sampname",drop=False,inplace=True)

rule all:
  input: 
    sunkannot = expand("sunkpos/{sample}.sunkpos", sample=file_df.index),  
    rlens = expand("sunkpos/{sample}.rlen", sample=file_df.index),  
    plots = expand('plots/{sample}_done',sample=file_df.index),  
      
    # sampname needs to be unique identifier
    # sunkannot = "sunkpos/ONT.sunkpos"


rule define_SUNKs:
  input:
    asm = config['asm'],
  output:
    counts = "db/jellyfish.counts",
    db = "db/jellyfish.db",
    fa = "db/jellyfish.fa"
  resources:
    mem=1 #28gb used
  threads: 32
  benchmark:
    "benchmarks/define_SUNKs.txt"
  shell:
    """
    jellyfish count -m {config[SUNK_len]} -s 10000000 -t 32 -C -c 1 -U 1 {input.asm} -o {output.counts}
    jellyfish dump -c -t {output.counts}  | awk '{{print $1}}' > {output.db}
    awk '{{ print ">"$0"\\n"$0 }}' {output.db} > {output.fa}
    """

rule mrsfast_index:
  input:
    asm = config['asm']
  output:
      index = "mrsfast/ref.fa.index"
  resources:
    mem=8 #1.1G used
  threads: 1
  benchmark:
    "benchmarks/mrsfast_index.txt"     
  # create softlinked ref in mrsfast folder
  shell:
    """
    cd mrsfast && \
    ln -s -f {input.asm} ./ref.fa && \
    mrsfast --ws 14 --index  ./ref.fa
    """

rule mrsfast_search:
  input:
    index = rules.mrsfast_index.output.index,
    db = rules.define_SUNKs.output.fa
  output:
    sam = "mrsfast/kmer_map.sam",
    bam = "mrsfast/kmer_map.bam",
  resources:
      mem=2 #73 G VMS
  threads: 64
  benchmark:
    "benchmarks/mrsfast_search.txt"  
  shell:
    """
    mrsfast --search mrsfast/ref.fa --threads 64 --mem 96 --seq {input.db} -o {output.sam} --disable-nohits -e 0 
    samtools sort -@ 64 -m 1.5G {output.sam} -T $TMPDIR -o {output.bam} 
    """

rule bed_convert:
  input:
    sam = rules.mrsfast_search.output.sam,
    bam = rules.mrsfast_search.output.bam,
  output:
    bed = "mrsfast/kmer_map.bed",
    bedmerge = "mrsfast/kmer_map.merge.bed",
    locs = "mrsfast/kmer.loc",
  resources:
      mem=16
  threads: 2
  benchmark:
    "benchmarks/bed_convert.txt"  
  shell:
    """
    bedtools bamtobed -i {input.bam} > {output.bed} && \
    bedtools merge -i {output.bed} > {output.bedmerge} && \
    bedtools intersect -a {output.bed} -b {output.bedmerge} -wo | cut -f 1,2,4,8 > {output.locs}
    """



def get_reads(wildcard):
  return file_df.loc[wildcard]['path']

rule SUNK_annot:
  input:
    locs = rules.bed_convert.output.locs,
    db = rules.define_SUNKs.output.fa,
    ONT = get_reads,
  output:
    'sunkpos/{sample}.sunkpos',
  resources:
      mem=8
  threads: 2
  benchmark:
    "benchmarks/SUNK_annot_{sample}.txt",      
  shell:
    """
    scripts/kmerpos_annot3 {input.ONT} {input.db} {input.locs} {output}
    """



rule read_lengths:
  input:
    ONT = get_reads,
  output:
    'sunkpos/{sample}.rlen',
  resources:
      mem=8
  threads: 2
  benchmark:
    "benchmarks/read_lengths_{sample}.txt",      
  shell:
    """
    scripts/rlen {input.ONT} {output}
    """

def get_bed(wildcard):
  return file_df.loc[wildcard]['bed']

rule viz:
  input:
    bed = get_bed,
    locs = rules.bed_convert.output.locs,
    rlen = 'sunkpos/{sample}.rlen',
    SUNKPOS = 'sunkpos/{sample}.sunkpos',
  output:
    'plots/{sample}_done',    
  resources:
      mem=160
  threads: 2    
  benchmark:
    "benchmarks/viz_{sample}.txt",        
  conda:
    "envs/viz.yaml"
  shell:
    """
    python scripts/viz.py {input.bed} {input.locs} {input.SUNKPOS} {input.rlen} &&
    touch {output}
    """