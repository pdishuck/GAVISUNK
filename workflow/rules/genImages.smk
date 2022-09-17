include: "gatherSplits.smk"
rule confirm_out:
  input:
    final_bed = expand("results/{sample}/final_out/{hap}.valid.bed", sample=manifest_df.index, hap=['hap1','hap2']),
    interout = getInterOut,
  output:
    flag = touch("results/{sample}/inter_outs/{hap}.done")

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


rule viz_check_gaps:
  input:
    unpack(vizInputsCheck)
  output:
    flag = touch('results/pngs/check_gaps/{sample}/{sample}_{hap}.done')
  resources:
    mem = 160,
    load = 200,
  threads: 2
  conda:
    "../envs/viz.yaml"
  log:
     "logs/{sample}_{hap}_viz_contigs.log"
  script:
    "../scripts/viz.py"

    
rule checkGaps:
  input:
     expand('results/pngs/check_gaps/{sample}/{sample}_{hap}.done', sample=manifest_df.index, hap=['hap1','hap2'])
