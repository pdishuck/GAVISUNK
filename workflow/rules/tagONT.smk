include: "gatherSplits.smk"
rule split_ONT: # accept FOFN
  input:
    reads=lambda wildcards: manifest_df.at[wildcards.sample,f"{wildcards.hap}_ONT"]
  output:
    reads=scatter.split("temp/{{sample}}/reads/{{hap}}_{scatteritem}.fq.gz"),
  resources:
    mem=2,
    load = 100
  threads: 8
  conda:
    "../envs/viz.yaml"
  log:
    scatter.split("logs/{{sample}}/{{hap}}_{scatteritem}_split_ONT.log")
  shell:
    """
    cat {input.reads} | seqtk seq -F '#' | rustybam fastq-split {output.reads}
    """
rule SUNK_annot:
  input:
    locs = rules.bed_convert.output.locs,
    db = rules.define_SUNKs.output.db,
    ONT = "temp/{sample}/reads/{hap}_{scatteritem}.fq.gz",
  output:
    'results/{sample}/sunkpos/{hap}_{scatteritem}.sunkpos'
  resources:
    mem=8,
    load=100
  threads: 2 
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/{hap}_{scatteritem}_SUNK_annot.log"  
  shell:
    """
    workflow/scripts/kmerpos_annot3 {input.ONT} {input.db} {input.locs} {output}
    """

rule read_lengths:
  input:
    ONT = "temp/{sample}/reads/{hap}_{scatteritem}.fq.gz",
  output:
    'results/{sample}/sunkpos/{hap}_{scatteritem}.rlen',
  resources:
    mem=8,
    load=100
  threads: 2
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/{hap}_{scatteritem}_read_lengths.log"  
  shell:
    """
    workflow/scripts/rlen {input.ONT} {output}
    """

rule diag_filter_step:
  input:
    ONT_pos = rules.SUNK_annot.output,
    fai =  lambda wildcards: manifest_df.at[wildcards.sample,f"{wildcards.hap}_asm"]+".fai"
  output:
    ONT_pos_diag = 'results/{sample}/sunkpos/{hap}_{scatteritem}_diag.sunkpos'
  resources:
    mem=8,
    load=100
  threads: 1
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/{hap}_{scatteritem}_diag_filter_step.log"
  shell:
    """
    workflow/scripts/diag_filter_v3 {input.ONT_pos} {input.fai} > {output.ONT_pos_diag}
    """
rule diag_filter_final:
  input:
    ONT_pos = rules.SUNK_annot.output,
    ONT_pos_diag = rules.diag_filter_step.output.ONT_pos_diag
  output:
    ONT_pos_diag_final = 'results/{sample}/sunkpos/{hap}_{scatteritem}_diag2.sunkpos'
  resources:
    mem=8,
    load=100
  threads: 1
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/{hap}_{scatteritem}_diag_filter_final.log"
  shell:
    """
    workflow/scripts/diag_filter_step2 {input.ONT_pos} {input.ONT_pos_diag} > {output.ONT_pos_diag_final}
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
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/{hap}_combine_ont.log"
  shell:
    """
    cat {input.gather_ONT_pos} > {output.ONT_pos}
    cat {input.gather_ONT_len} > {output.ONT_len}
    """
rule bad_sunks:
  input:
    hap1_fai = lambda wildcards: manifest_df.at[wildcards.sample, "hap1_asm"]+".fai",
    hap2_fai = lambda wildcards: manifest_df.at[wildcards.sample, "hap2_asm"]+".fai",
    hap1_sunkpos = "results/{sample}/sunkpos/hap1.sunkpos",
    hap2_sunkpos =  "results/{sample}/sunkpos/hap2.sunkpos"
  output:
    badsunks = "results/{sample}/sunkpos/bad_sunks.txt"
  resources:
    mem=10,
    load=100
  threads: 1
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/bad_sunks.log"
  shell:
    """
    python workflow/scripts/badsunks_AR.py {input.hap1_fai} {input.hap2_fai} {input.hap1_sunkpos} {input.hap2_sunkpos} {output.badsunks}
    """

checkpoint split_sunkpos:
  input:
    ONT_pos = rules.combine_ont.output.ONT_pos,
    kmer_loc = rules.bed_convert.output.locs
  output:
    flag = touch("results/{sample}/breaks/{hap}_splits_pos.done")
  resources:
    mem=10,
    load=100
  threads:1
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/{hap}_split_sunkpos.log"
  script:
    "../scripts/split_locs.py"
 
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
    mem=16,
    load=250,
  threads: 1
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/{contigs}_{hap}_process_by_contig.log"
  shell:
    """
    python workflow/scripts/process-by-contig_lowmem_AR.py {input.locs_contig} {input.sunk_pos_contig} {input.rlen} {input.bad_sunks} {output.outputdf} {output.bed}
    touch {output.bed}
    """

rule gather_process_by_contig:
  input:
    bed = gatherAsmBeds,
    flag = rules.split_sunkpos.output.flag
  output:
    allbeds = "results/{sample}/final_out/{hap}.valid.bed"
  resources:
    mem=8,
    load=100
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/{hap}_gather_proccess_by_contig.log"
  shell:
    '''
    cat {input.bed} > {output.allbeds}
    '''

rule get_gaps:
  input:
    hap1_fai = lambda wildcards: manifest_df.at[wildcards.sample, "hap1_asm"]+".fai",
    hap2_fai = lambda wildcards: manifest_df.at[wildcards.sample, "hap2_asm"]+".fai",
    allbeds_hap1 = "results/{sample}/final_out/hap1.valid.bed",
    allbeds_hap2 = "results/{sample}/final_out/hap2.valid.bed",
  output:
    gaps1 = "results/{sample}/final_out/hap1.gaps.bed",
    gaps2 = "results/{sample}/final_out/hap2.gaps.bed",
  resources:
    mem=10,
    load=100
  threads: 1
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/get_gaps.log"
  shell:
    """
    python workflow/scripts/get_gaps.py {input.hap1_fai} {input.hap2_fai}  {wildcards.sample} results/{wildcards.sample}/bed_files/  results/{wildcards.sample}/final_out/
    """

rule slop_gaps:
  input:
    gaps = "results/{sample}/final_out/{hap}.gaps.bed",
    fai =  lambda wildcards: manifest_df.at[wildcards.sample,f"{wildcards.hap}_asm"]+".fai"
  output:
    gaps_slop = "results/{sample}/final_out/{hap}.gaps.slop.bed",
  resources:
    mem=10,
    load=100
  threads: 1
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/get_gaps_{hap}.log"
  shell:
    """
    bedtools slop -i {input.gaps} -g {input.fai}  -b 200000 > {output.gaps_slop}
    """

