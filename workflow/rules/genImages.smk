rule align_to_ref:
  input:
    fasta = lambda wildcards: manifest_df.at[wildcards.sample, f"{wildcards.hap}_asm"],
    ref = config.get('ref')
  output:
    paf = 'results/{sample}/align_to_ref/{sample}_{hap}_to_ref.paf'
  resources:
    mem=8,
    load=100
  threads: 12
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/{hap}_align_to_ref.log"
  shell:
    '''
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
    "../envs/viz.yaml"
  log:
    "logs/{sample}/{hap}_subset_paf.log"
  shell:
    '''
    rustybam liftover -l --bed {input.bed} {input.paf} > {output.subset_paf}
    awk '{{OFS="\\t"; print $6, $8, $9, $1, $3, $4}}' {output.subset_paf}| sort -k1,1 -k2,2n > {output.inter_bed}
    bedtools intersect -a {output.inter_bed} -b {input.bed} -wo | cut -f4,5,6,10 | sort -k1,1 -k2,2n | bedtools merge -i - -d 175000 -c 4 -o distinct > {output.final_bed}
    '''
rule subset_bed:
  input:
    bed = rules.subset_paf.output.final_bed
  output:
    flag = touch('results/{sample}/align_to_ref/gene_beds/{sample}_{hap}_to_ref.done')
  resources:
    mem=8,
    load=100,
  threads: 1
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/{hap}_subset_bed."
  script:
    "../scripts/split_beds.py"

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
    "../envs/viz.yaml"
  log:
    "logs/{sample}/{gene}_{hap}_viz_contigs.log"
  shell:
    '''
    python ../scripts/viz_AR.py {input.bed} {input.rlen} $(dirname {output.flag_gene}) $(dirname {input.interout[0]}) $(dirname {input.pos_locs}) $(dirname {input.pos_locs}) {wildcards.sample} {wildcards.hap} 
    '''

rule generate_images:
  input:
    expand('results/pngs/{gene}/{sample}/{sample}_{hap}.done', gene=bed_df.index, sample=manifest_df.index, hap=['hap1', 'hap2'] )
