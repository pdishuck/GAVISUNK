include: "gatherSplits.smk"

checkpoint split_gaps:
  input:
    gaps_bed = "results/{sample}/final_out/{hap}.gaps.bed"
  output:
    flag = touch("results/{sample}/check_gaps/{hap}.gaps.done"),
    define_path = touch("results/{sample}/check_gaps/{hap}.determine_path.done")
  resources:
    mem=10,
    load=100
  threads: 1
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/get_gaps_{hap}.log"
  script:
    "../scripts/split_beds.py"

rule subset_sunkpos:
  input:
    split_bed = "results/{sample}/check_gaps/{contgap}_{hap}.gaps.bed",
    sunkpos = 'results/{sample}/sunkpos/{hap}_{range}-of-{lens}.sunkpos'
  output:
    subset_sunkpos="results/{sample}/check_gaps/subset/{contgap}_{hap}_{range}-of-{lens}.subset_region.sunkpos"
  resources:
    mem=10,
    load=100,
  threads: 1
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/sunk_sunkpos_{contgap}_{hap}_{range}-of-{lens}.log"  
  shell:
    '''
    cat {input.split_bed} | while read line; do
        contigstart=$(echo ${{line}} | awk '{{print $2}}')
        contigstop=$(echo ${{line}} | awk '{{print $3}}')
        contig=$(echo ${{line}} | awk '{{print $1}}')
        grep -P "${{contig}}\\t" {input.sunkpos}  | awk -v var="$contigstart" '($4>=var)' |  awk -v var="$contigstop" '($4<=var)' >> {output.subset_sunkpos}
    done
    '''

rule subset_diag2:
  input:
    split_bed = "results/{sample}/check_gaps/{contgap}_{hap}.gaps.bed",
    diag2 = 'results/{sample}/sunkpos/{hap}_{range}-of-{lens}_diag2.sunkpos'
  output:
    subset_diag="results/{sample}/check_gaps/subset/{contgap}_{hap}_{range}-of-{lens}.subset_diag2.sunkpos"
  resources:
    mem=10,
    load=100,
  threads: 1
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/sunk_diag2_{contgap}_{hap}_{range}-of-{lens}.log"  
  shell:
    '''
    cat {input.split_bed} | while read line; do
        contigstart=$(echo ${{line}} | awk '{{print $2}}')
        contigstop=$(echo ${{line}} | awk '{{print $3}}')
        contig=$(echo ${{line}} | awk '{{print $1}}')
        grep -P "${{contig}}\\t" {input.diag2}  | awk -v var="$contigstart" '($4>=var)' |  awk -v var="$contigstop" '($4<=var)' >> {output.subset_diag}
    done
    '''
rule cat_sunkpos:
  input: 
    all_sunkpos = lambda wildcards: expand("results/{sample}/check_gaps/subset/{contgap}_{hap}_{range}-of-{lens}.subset_region.sunkpos",sample=wildcards.sample, hap=wildcards.hap, contgap=wildcards.contgap, range=range(1,splitsize+1), lens=splitsize),
    all_diag2 = lambda wildcards: expand("results/{sample}/check_gaps/subset/{contgap}_{hap}_{range}-of-{lens}.subset_diag2.sunkpos",sample=wildcards.sample, hap=wildcards.hap, contgap=wildcards.contgap, range=range(1,splitsize+1), lens=splitsize),
  output:
    all_sunkpos = 'results/{sample}/check_gaps/subset/{contgap}_{hap}_subset_region.sunkpos',
    all_diag2 = 'results/{sample}/check_gaps/subset/{contgap}_{hap}_subset_diag2.sunkpos'
  resources:
    mem = 12,
    hrs = 12
  threads: 1
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/{contgap}_{hap}.cat_sunkpos.log"
  shell:
    '''
    cat {input.all_sunkpos} > {output.all_sunkpos}
    cat {input.all_diag2} > {output.all_diag2}
    '''
rule get_lists:
  input:
    sunkpos = rules.cat_sunkpos.output.all_sunkpos,
    diag2 = rules.cat_sunkpos.output.all_diag2,
  output:
    sunkpos = 'results/{sample}/check_gaps/subset/{contgap}_{hap}_subset_region.sunkpos.list',
    diag2 = 'results/{sample}/check_gaps/subset/{contgap}_{hap}_subset_region_diag2.sunkpos.list'
  resources:
    mem = 12,
    hrs = 12
  threads: 1
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/{contgap}_{hap}.get_lists.log"
  shell:
    '''
    cut -f1 {input.sunkpos} | sort -u > {output.sunkpos}
    cut -f1 {input.diag2} | sort -u > {output.diag2}
    '''
rule full_extract:
  input: 
    all_sunkpos = 'results/{sample}/sunkpos/{hap}_{range}-of-{lens}.sunkpos',
    all_diag2 = 'results/{sample}/sunkpos/{hap}_{range}-of-{lens}_diag2.sunkpos',
    region_reads = rules.get_lists.output.sunkpos,
    diag2_region_reads = rules.get_lists.output.diag2
  output:
    subset_sunkpos = 'results/{sample}/check_gaps/subset/{contgap}_{hap}_{range}-of-{lens}_subset_reads.sunkpos',
    subset_diag2 = 'results/{sample}/check_gaps/subset/{contgap}_{hap}_{range}-of-{lens}_subset_reads_diag2.sunkpos'
  params:
    contig = lambda wildcards: wildcards.contgap.replace("_","#")
  resources:
    mem = 12,
    hrs = 12
  threads: 1
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/{contgap}_{hap}_{range}_{lens}.full_extact.log"
  shell:
    '''
    grep -P "{params.contig}\\t" {input.all_sunkpos} | grep -f {input.region_reads} > {output.subset_sunkpos} || true
    grep -P "{params.contig}\\t" {input.all_diag2} | grep -f {input.diag2_region_reads} > {output.subset_diag2} || true
    '''
rule combine_subset:
  input:
    all_sunkpos = lambda wildcards: expand("results/{sample}/check_gaps/subset/{contgap}_{hap}_{range}-of-{lens}_subset_reads.sunkpos", sample=wildcards.sample, hap=wildcards.hap, contgap=wildcards.contgap, range=range(1,splitsize+1), lens=splitsize),
    all_diag2 = lambda wildcards: expand("results/{sample}/check_gaps/subset/{contgap}_{hap}_{range}-of-{lens}_subset_reads_diag2.sunkpos", sample=wildcards.sample, hap=wildcards.hap, contgap=wildcards.contgap, range=range(1,splitsize+1), lens=splitsize)
  output:
    all_sunkpos = 'results/{sample}/check_gaps/subset/{contgap}_{hap}_subset_reads.sunkpos',
    all_diag2 = 'results/{sample}/check_gaps/subset/{contgap}_{hap}_subset_reads_diag2.sunkpos'
  resources:
    mem =12,
    hrs = 12,
  threads: 1
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/{contgap}_{hap}.combine_subset.log"
  shell:
    '''
    cat {input.all_sunkpos} > {output.all_sunkpos}
    cat {input.all_diag2} > {output.all_diag2}
    '''
rule only_not_placed:
  input:
    all_sunkpos_list = rules.get_lists.output.sunkpos,
    all_sunkpos = rules.combine_subset.output.all_sunkpos,
    all_diag = lambda wildcards: expand('results/{sample}/sunkpos/{hap}_{range}-of-{lens}_diag.sunkpos',sample=wildcards.sample, hap=wildcards.hap, contgap=wildcards.contgap, range=range(1,splitsize+1), lens=splitsize)
  output:
    placed_list = 'results/{sample}/check_gaps/subset/{contgap}_{hap}_subset.sunkpos_placed.list',
    not_placed_list = 'results/{sample}/check_gaps/subset/{contgap}_{hap}_subset.sunkpos_unplaced.list',
    not_placed_sunkpos = 'results/{sample}/check_gaps/subset/{contgap}_{hap}_subset.sunkpos_unplaced.sunkpos'
  resources:
    mem =12,
    hrs = 12,
  threads: 1
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/{contgap}_{hap}.only_not_placed.log"
  shell:
    '''
    grep -h -f {input.all_sunkpos_list} {input.all_diag} | cut -f1 | sort -u > {output.placed_list} || true
    grep -v -f {output.placed_list} {input.all_sunkpos} > {output.not_placed_sunkpos} || true
    cut -f1 {output.not_placed_sunkpos} | sort -u > {output.not_placed_list}
    '''

rule combine_diag:
  input:
    not_placed = rules.only_not_placed.output.not_placed_sunkpos,
    diag2 = rules.combine_subset.output.all_diag2
  output:
    combined = 'results/{sample}/check_gaps/subset/{contgap}_{hap}.sunkpos'
  resources:
    mem =12,
    hrs = 12,
  threads: 1
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/{contgap}_{hap}.combine_diag.log"
  shell:
    '''
    cat {input.not_placed} {input.diag2} > {output.combined}
    '''

rule process_by_contig_join_reads:
  input:
    sunk_pos_contig = rules.combine_diag.output.combined,
    locs_contig = rules.bed_convert.output.locs,
    rlen = rules.combine_ont.output.ONT_len,
    bad_sunks = rules.bad_sunks.output.badsunks ,
  output:
    outputdf = "results/{sample}/check_gaps/inter_outs/{contgap}_{hap}_test.tsv",
    bed = "results/{sample}/check_gaps/bed_files/{contgap}_{hap}.bed"
  resources:
    mem=16,
    load=250,
  threads: 1
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/{contgap}_{hap}.pbcjr.log"
  shell:
    """
    python workflow/scripts/process-by-contig_lowmem_AR.py {input.locs_contig} {input.sunk_pos_contig} {input.rlen} {input.bad_sunks} {output.outputdf} {output.bed}
    touch {output.bed}
    """
rule add_reads:
  input:
    sunkpos_orig="results/{sample}/breaks/{contgap}_{hap}.sunkpos",
    sunkpos_new=rules.only_not_placed.output.not_placed_sunkpos,
    tsv_orig="results/{sample}/inter_outs/{contgap}_{hap}.tsv",
    tsv_new=rules.process_by_contig_join_reads.output.outputdf,
  output:
    sunkpos_combine = "results/{sample}/check_gaps/breaks/{contgap}_{hap}.sunkpos",
    tsv_combine = "results/{sample}/check_gaps/inter_outs/{contgap}_{hap}_gaps.tsv"
  resources:
    mem=16,
    load=250,
  threads: 1
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/{contgap}_{hap}_add_reads.log"
  shell:
    '''
    cat {input.sunkpos_orig} {input.sunkpos_new} > {output.sunkpos_combine}
    cat {input.tsv_orig} {input.tsv_new} > {output.tsv_combine}
    '''
rule filter_no_reads:
  input:
    tsv = rules.add_reads.output.tsv_combine
  output:
    correct_df = 'results/{sample}/check_gaps/inter_outs/{contgap}_{hap}.tsv'
  shell:
    '''
    grep -v "no_reads" {input.tsv} > {output.correct_df}
    '''
rule gather_checks:
  input:
    sunkpos = getSunkPos2,
    tsv = getInterOutsGaps,
    flag = rules.split_gaps.output.flag
  output:
    sunkpos_flag = touch('results/{sample}/check_gaps/breaks/{hap}.done'),
    tsv_flag = touch('results/{sample}/check_gaps/inter_outs/{hap}_inter.done')
  log:
    "logs/{sample}/{hap}_gather_checks.log"

rule merge_beds:
  input:
    bed_orig='results/{sample}/bed_files/{contgap}_{hap}.bed',
    bed_new=rules.process_by_contig_join_reads.output.bed
  output:
    bed_combine = "results/{sample}/check_gaps/bed_files/{contgap}_{hap}.merged.valid.bed"
  resources:
    mem=16,
    load=250,
  threads: 1
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/{contgap}_{hap}.merge_beds.log"
  shell:
    '''
    cat {input.bed_orig} {input.bed_new} | bedtools sort -i stdin | bedtools merge -i stdin > {output.bed_combine}
    '''

rule gather_process_by_contig_check:
  input:
    bed = getGaps,
    flag = rules.split_gaps.output.flag
  output:
    allbeds = "results/{sample}/check_gaps/final_out/{hap}.merged.valid.bed"
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
rule get_gaps_check:
  input:
    hap1_fai = lambda wildcards: manifest_df.at[wildcards.sample, f'{wildcards.hap}_asm']+".fai",
    allbeds_hap1 = "results/{sample}/check_gaps/final_out/{hap}.merged.valid.bed",
  output:
    gaps1 = "results/{sample}/check_gaps/final_out/{hap}.merged.gaps.bed",
  resources:
    mem=10,
    load=100
  threads: 1
  conda:
    "../envs/viz.yaml"
  log:
    "logs/{sample}/{hap}_get_gaps.log"
  shell:
    """
    python workflow/scripts/get_gaps_2.py {input.hap1_fai} {wildcards.sample} {wildcards.hap} results/{wildcards.sample}/check_gaps/bed_files/  results/{wildcards.sample}/check_gaps/final_out/
    """


rule slop_gaps_check:
  input:
    gaps = "results/{sample}/check_gaps/final_out/{hap}.merged.gaps.bed",
    fai =  lambda wildcards: manifest_df.at[wildcards.sample,f"{wildcards.hap}_asm"]+".fai",
  output:
    gaps_slop = "results/{sample}/check_gaps/final_out/{hap}.merged.gaps.slop.bed",
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

rule no_gaps:
  input:
    flag = rules.split_gaps.output.flag
  output:
    no_gaps = touch('results/check_gaps/{sample}_{hap}_no_gaps.txt')
