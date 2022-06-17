# tested with Snakemake v6.7.0
import pandas as pd
import os

configfile: "config.yaml"

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)


MANIFEST = config["ONT_manifest"]
BED = config['bed']

manifest_df = pd.read_csv(MANIFEST, sep='\t', header=0, index_col='sample')
print(manifest_df)
bed_df = pd.read_csv(BED, sep='\t', header=None, names=['contig','start','stop','gene'])
bed_df.set_index(['gene'], inplace=True, drop=True)

scattergather:
    split=config.get("nchunks", 10),

wildcard_constraints:
    hap="hap1|hap2"

def splitlocs(hap_type,group_df, dirname, exten):
  for name, contig_group in group_df:
    contig_name = name.replace("#", "_")
    #contig_name_final = contig_name.replace(".", "_")
    output_path = f'{dirname}/{contig_name}_{hap_type}.{exten}'
    print(output_path)
    contig_group.to_csv(output_path, header=None, index=None, sep="\t")
  return(True)

def getONT(wildcards):
  return manifest_df.at[wildcards.sample,wildcards.hap+"_ONT"]

def getasms(wildcards):
    return(manifest_df.at[wildcards.sample, wildcards.hap+"_asm"])

def getParAsm(wildcards):
  return {'hap1_asm': manifest_df.at[wildcards.sample, "hap1_asm"], 'hap2_asm': manifest_df.at[wildcards.sample, "hap2_asm"]}

def gethapFai(wildcards):
  return manifest_df.at[wildcards.sample,wildcards.hap+"_asm"]+".fai"

def getParFai(wildcards):
  return {'hap1_fai': manifest_df.at[wildcards.sample, "hap1_asm"]+".fai", 'hap2_fai': manifest_df.at[wildcards.sample, "hap2_asm"]+".fai"}

def gatherAsmBeds(wildcards):
  CONTIGS = glob_wildcards("results/{sample}/breaks/{contigs}_{hap}.sunkpos".format(sample=wildcards.sample, hap=wildcards.hap, contigs='{contigs}')).contigs
  return expand(rules.process_by_contig.output.bed, sample=wildcards.sample, hap=wildcards.hap, contigs=CONTIGS)

def getSunkLocs(wildcards):
  CONTIGS = glob_wildcards("results/{sample}/breaks/{contigs}_{hap}.loc".format(sample=wildcards.sample, hap=wildcards.hap, contigs='{contigs}')).contigs
  return expand("results/{sample}/breaks/{contigs}_{hap}.loc", sample=wildcards.sample, hap=wildcards.hap, contigs=CONTIGS)

def getInterOut(wildcards):
  CONTIGS = glob_wildcards("results/{sample}/inter_outs/{contigs}_{hap}.tsv".format(sample=wildcards.sample, hap=wildcards.hap, contigs='{contigs}')).contigs
  return expand("results/{sample}/inter_outs/{contigs}_{hap}.tsv", sample=wildcards.sample, hap=wildcards.hap, contigs=CONTIGS)

def getSunkPos(wildcards):
  CONTIGS = glob_wildcards("results/{sample}/breaks/{contigs}_{hap}.sunkpos".format(sample=wildcards.sample, hap=wildcards.hap, contigs='{contigs}')).contigs
  return expand("results/{sample}/breaks/{contigs}_{hap}.sunkpos", sample=wildcards.sample, hap=wildcards.hap, contigs=CONTIGS)

def getViz(wildcards):
  return expand("results/pngs/{gene}/{sample}/{sample}_{hap}.done", gene=bed_df.index, sample=wildcards.sample, hap=wildcards.hap)

def getPaf(wildcards):
  if pd.isnull(manifest_df.at[wildcards.sample, wildcards.hap+"_paf"]):
    return rules.align_to_ref.output.paf
  else:
    return manifest_df.at[wildcards.sample, wildcards.hap+"_paf"]

rule all:
  input: 
    expand('results/{sample}/sunkpos/{hap}.sunkpos', sample=manifest_df.index, hap=['hap1','hap2']),
    expand('results/{sample}/sunkpos/{hap}.rlen', sample=manifest_df.index, hap=['hap1','hap2']),
    expand('results/{sample}/sunkpos/bad_sunks.txt', sample=manifest_df.index),
    expand("results/{sample}/breaks/{hap}_splits_pos.done",sample=manifest_df.index, hap=['hap1','hap2']),
    expand("results/{sample}/final_out/{hap}.bed", sample=manifest_df.index, hap=['hap1','hap2']), 

rule split_ONT: # accept FOFN
  input:
    reads=getONT
  output:
    reads=temp(scatter.split("temp/{{sample}}/reads/{{hap}}_{scatteritem}.fq.gz")),
  resources:
    mem=2,
    load = 100
  threads: 8
  conda:
    "envs/viz.yaml"
  shell:
    """
    zcat -f $(cat {input.reads} ) | seqtk seq -F '#' | rustybam fastq-split {output.reads}
    """

rule combine_asm_haps:
  input:
    unpack(getParAsm)
  output:
    combined = temp("results/{sample}/mrsfast/ref.fa")
  resources:
    mem = 8,
    load = 100
  threads: 1
  conda:
    "envs/viz.yaml"
  shell:
    """
    zcat -f {input.hap1_asm} {input.hap2_asm} > {output.combined}
    samtools faidx {output.combined}
    """

rule jellyfish_count:
  input:
    asm = rules.combine_asm_haps.output.combined,
  output:
    counts = "results/{sample}/db/jellyfish.counts"
  resources:
    mem = 1,
    load = 900
  threads:32
  params:
    sunk_len = config["SUNK_len"]
  shell:
    """
    jellyfish count -m {params.sunk_len} -s 10000000 -t {threads} -C -c 1 -U 1 {input.asm} -o {output.counts}
    """

rule define_SUNKs:
  input:
    counts = rules.jellyfish_count.output.counts
  output:
    db = "results/{sample}/db/jellyfish.db",
    fa = "results/{sample}/db/jellyfish.fa"
  resources:
    mem=1,
    load = 900
  threads: 32
  shell:
    """
    jellyfish dump -c -t {input.counts}  | awk '{{print $1}}' > {output.db}
    awk '{{ print ">"$0"\\n"$0 }}' {output.db} > {output.fa}
    """

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
    "envs/viz.yaml"     
  shell:
    """
    mrsfast --ws 14 --index {input.asm}
    """

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
  threads: 32

  conda:
    "envs/viz.yaml"  
  shell:
    """
    mrsfast --search {input.ref} --threads {threads} --mem 96 --seq {input.db} -o {output.sam} --disable-nohits -e 0 
    samtools sort -@ {threads} -m 2G {output.sam} -T $TMPDIR -o {output.bam} 
    """

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
    "envs/viz.yaml"  
  shell:
    """
    bedtools bamtobed -i {input.bam} > {output.bed} && \
    bedtools merge -i {output.bed} > {output.bedmerge} && \
    bedtools intersect -a {output.bed} -b {output.bedmerge} -wo | cut -f 1,2,4,8 > {output.locs}
    """

rule SUNK_annot:
  input:
    locs = rules.bed_convert.output.locs,
    db = rules.define_SUNKs.output.fa,
    ONT = "temp/{sample}/reads/{hap}_{scatteritem}.fq.gz",
  output:
    temp('results/{sample}/sunkpos/{hap}_{scatteritem}.sunkpos')
  resources:
      mem=8,
      load=100
  threads: 2    
  shell:
    """
    scripts/kmerpos_annot3 {input.ONT} {input.db} {input.locs} {output}
    """

rule read_lengths:
  input:
    ONT = "temp/{sample}/reads/{hap}_{scatteritem}.fq.gz",
  output:
    temp('results/{sample}/sunkpos/{hap}_{scatteritem}.rlen'),
  resources:
      mem=8,
      load=100
  threads: 2   
  shell:
    """
    scripts/rlen {input.ONT} {output}
    """

rule diag_filter_step:
  input:
    ONT_pos = rules.SUNK_annot.output,
    fai =  gethapFai 
  output:
    ONT_pos_diag = temp('results/{sample}/sunkpos/{hap}_{scatteritem}_diag.sunkpos')
  resources:
    mem=8,
    load=100
  threads: 1
  shell:
    """
    scripts/diag_filter_v3 {input.ONT_pos} {input.fai} > {output.ONT_pos_diag}
    """

rule diag_filter_final:
  input:
    ONT_pos = rules.SUNK_annot.output,
    ONT_pos_diag = rules.diag_filter_step.output.ONT_pos_diag
  output:
    ONT_pos_diag_final = temp('results/{sample}/sunkpos/{hap}_{scatteritem}_diag2.sunkpos')
  resources:
    mem=8,
    load=100
  threads: 1
  shell:
    """
    scripts/diag_filter_step2 {input.ONT_pos} {input.ONT_pos_diag} > {output.ONT_pos_diag_final}
    """

rule combine_ont:
  input:
    gather_ONT_pos = gather.split("results/{{sample}}/sunkpos/{{hap}}_{scatteritem}_diag2.sunkpos"),
    gather_ONT_len = gather.split("results/{{sample}}/sunkpos/{{hap}}_{scatteritem}.rlen")
  output:
    ONT_pos = 'results/{sample}/sunkpos/{hap}.sunkpos',
    ONT_len = 'results/{sample}/sunkpos/{hap}.rlen'
  resources:
    mem=8,
    load=100
  threads: 1
  shell:
    """
    cat {input.gather_ONT_pos} > {output.ONT_pos}
    cat {input.gather_ONT_len} > {output.ONT_len}
    """
rule bad_sunks:
  input:
    unpack(getParFai),
    hap1_sunkpos = "results/{sample}/sunkpos/hap1.sunkpos",
    hap2_sunkpos =  "results/{sample}/sunkpos/hap2.sunkpos"
  output:
    badsunks = "results/{sample}/sunkpos/bad_sunks.txt"
  resources:
    mem=25,
    load=100
  threads: 1
  conda:
    "envs/viz.yaml"
  shell:
    """
    python scripts/badsunks_AR.py {input.hap1_fai} {input.hap2_fai} {input.hap1_sunkpos} {input.hap2_sunkpos} {output.badsunks}
    """

checkpoint split_sunkpos:
  input:
    ONT_pos = rules.combine_ont.output.ONT_pos,
    kmer_loc = rules.bed_convert.output.locs
  output:
    flag = touch("results/{sample}/breaks/{hap}_splits_pos.done")
  resources:
    mem=50,
    load=100
  threads:1
  benchmark:
    "benchmarks/split_sunkpos_{sample}_{hap}.done"
  run:
    ONT_pos_df = pd.read_csv(input.ONT_pos, header=None, sep="\t", names=['read','read_pos','contig','contig_start','contig_stop'])
    kmer_loc = pd.read_csv(input.kmer_loc, header=None, sep="\t", names=['contig','pos1','SUNK','pos2'])
    contig_list = ONT_pos_df['contig'].unique().tolist()
    kmer_loc_subset = kmer_loc[kmer_loc['contig'].isin(contig_list)]
    dir_name = os.path.dirname(output.flag)  
    group_ONT = ONT_pos_df.groupby(['contig'])
    group_loc = kmer_loc_subset.groupby(['contig'])
    output_complete = splitlocs(wildcards.hap,group_ONT, dir_name,"sunkpos")
    output_complete_2 = splitlocs(wildcards.hap,group_loc, dir_name,"loc")
 
rule process_by_contig:
  input:
    sunk_pos_contig = "results/{sample}/breaks/{contigs}_{hap}.sunkpos",
    locs_contig = "results/{sample}/breaks/{contigs}_{hap}.loc",
    rlen = rules.combine_ont.output.ONT_len,
    bad_sunks = rules.bad_sunks.output.badsunks,
  output:
    outputdf = "results/{sample}/inter_outs/{contigs}_{hap}.tsv",
    bed = "results/{sample}/bed_files/{contigs}_{hap}.bed"
  resources:
    mem=64,
    load=250,
  threads: 1
  conda:
    "envs/viz.yaml"
  shell:
    """
    python scripts/process-by-contig_lowmem_AR.py {input.locs_contig} {input.sunk_pos_contig} {input.rlen} {input.bad_sunks} {output.outputdf} {output.bed}
    touch {output.bed}
    """

rule gather_process_by_contig:
  input:
    bed = gatherAsmBeds,
    flag = rules.split_sunkpos.output.flag
  output:
    allbeds = touch("results/{sample}/final_out/{hap}.bed")
  resources:
    mem=8,
    load=100

rule align_to_ref:
  input:
    fasta = getasms,
    ref = config.get('ref','/net/eichler/vol26/eee_shared/assemblies/CHM13/T2T/v1.1/chm13.v1.1.fasta')
  output:
    paf = 'results/{sample}/align_to_ref/{sample}_{hap}_to_ref.paf'
  resources:
    mem=8,
    load=100
  threads: 12
  conda:
    "envs/viz.yaml"
  shell:
    '''
    module load minimap2/2.24
    minimap2 -t {threads} --eqx -c -x asm20 --secondary=no -s 25000 {input.ref} {input.fasta} > {output.paf}
    '''

rule subset_paf:
  input:
    bed = config['bed'],
    paf = getPaf
  output:
    subset_paf = 'results/{sample}/align_to_ref/{sample}_{hap}_to_ref_subset.paf',
    inter_bed = temp('results/{sample}/align_to_ref/inter_{sample}_{hap}_to_ref_subset.bed'),
    final_bed = 'results/{sample}/align_to_ref/{sample}_{hap}_to_ref_subset.bed'
  resources:
    mem=8,
    load=100,
  threads: 1
  conda:
    "envs/viz.yaml"
  shell:
    '''
    module load rustybam/0.1.27
    rustybam liftover -l --bed {input.bed} {input.paf} > {output.subset_paf}
    awk '{{OFS="\\t"; print $6, $8, $9, $1, $3, $4}}' {output.subset_paf}| sort -k1,1 -k2,2n > {output.inter_bed}
    bedtools intersect -a {output.inter_bed} -b {input.bed} -wo | cut -f4,5,6,10 | sort -k1,1 -k2,2n | bedtools merge -i - -d 175000 -c 4 -o distinct > {output.final_bed}
    '''
rule subset_bed:
  input:
    bed = rules.subset_paf.output.final_bed
  output:
    flag = touch('results/{sample}/align_to_ref/gene_beds/{sample}_{hap}_to_ref_{gene}.bed')
  resources:
    mem=8,
    load=100,
  threads: 1
  run:
    bed_df = pd.read_csv(input.bed, header=None, sep="\t", names=['contig','start','stop','region_name'])
    bed_df_group = bed_df.groupby(['region_name'])
    for (splitno, split) in bed_df_group:
      split.to_csv(f'results/{wildcards.sample}/align_to_ref/gene_beds/{wildcards.sample}_{wildcards.hap}_to_ref_{splitno}.bed', header=None, sep="\t", index=None )
      #print(f'results/{wildcards.sample}/align_to_ref/gene_beds/{wildcards.sample}_{wildcards.hap}_to_ref_{splitno}.bed')
rule viz_contigs:
  input:
    flag = rules.subset_bed.output.flag,
    bed = 'results/{sample}/align_to_ref/gene_beds/{sample}_{hap}_to_ref_{gene}.bed',
    rlen = rules.combine_ont.output.ONT_len,
    interout = getInterOut,
    pos_locs = rules.split_sunkpos.output.flag
  output:
    flag_gene = touch('results/pngs/{gene}/{sample}/{sample}_{hap}.done')
  resources:
    mem = 160,
    load=200
  threads: 2
  conda:
    "envs/viz.yaml"
  shell:
    '''
    python scripts/viz_AR.py {input.bed} {input.rlen} $(dirname {output.flag_gene}) $(dirname {input.interout[0]}) $(dirname {input.pos_locs}) $(dirname {input.pos_locs}) {wildcards.sample} {wildcards.hap} 
    '''
