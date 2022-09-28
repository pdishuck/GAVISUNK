include: "gatherSplits.smk"
rule confirm_out:
  input:
    final_bed = expand("results/{sample}/final_out/{hap}.valid.bed", sample=manifest_df.index, hap=['hap1','hap2']),
    interout = getInterOut,
  output:
    flag = touch("results/{sample}/inter_outs/{hap}.done")

rule viz_contigs_detailed:
  input:
    unpack(vizInputsDetailed)
  output:
    flag = touch('results/pngs/{runmode}/{sample}/{sample}_{hap}_detailed.done')
  resources:
    mem = 160,
    load = 200,
  threads: 2
  conda:
    "../envs/viz.yaml"
  log:
     "logs/{runmode}_{sample}_{hap}_viz_contigs_detailed.log"
  script:
    "../scripts/viz_detailed.py"
    
rule viz_contigs:
  input:
    unpack(vizInputs)
  output:
    flag = touch('results/pngs/{runmode}/{sample}/{sample}_{hap}.done')
  resources:
    mem = 160,
    load = 200,
  threads: 2
  conda:
    "../envs/viz.yaml"
  log:
     "logs/{runmode}_{sample}_{hap}_viz_contigs.log"
  script:
    "../scripts/viz.py"


rule covprob:
  input:
    unpack(covprobInputs),
    fai =  lambda wildcards: manifest_df.at[wildcards.sample,f"{wildcards.hap}_asm"]+".fai",
    # SUNK_len = str(config['SUNK_len']),
  output:
    tsv = 'results/{sample}/final_out/{hap}.gaps.covprob.tsv',
  resources:
    mem = 160,
    load = 200,
  threads: 2
  benchmark:
    "benchmarks/{sample}_{hap}_covprob.bench"
  conda:
    "../envs/viz.yaml"
  log:
     "logs/{sample}_{hap}_covprob.log"
  script:
    "../scripts/covprob.py"


