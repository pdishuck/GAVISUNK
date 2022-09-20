def gatherAsmBeds(wildcards):
  CONTIGS = glob_wildcards("results/{sample}/breaks/{contigs}_{hap}.sunkpos".format(sample=wildcards.sample, hap=wildcards.hap, contigs='{contigs}')).contigs
  return expand(rules.process_by_contig.output.bed, sample=wildcards.sample, hap=wildcards.hap, contigs=CONTIGS)

def getSunkLocs(wildcards):
  CONTIGS = glob_wildcards("results/{sample}/breaks/{contigs}_{hap}.loc".format(sample=wildcards.sample, hap=wildcards.hap, contigs='{contigs}')).contigs
  return expand("results/{sample}/breaks/{contigs}_{hap}.loc", sample=wildcards.sample, hap=wildcards.hap, contigs=CONTIGS)

def getInterOut(wildcards):
  CONTIGS = glob_wildcards("results/{sample}/inter_outs/{contigs}_{hap}.tsv".format(sample=wildcards.sample, hap=wildcards.hap, contigs='{contigs}')).contigs
  # print("results/{sample}/inter_outs/{contigs}_{hap}.tsv")
  return expand("results/{sample}/inter_outs/{contigs}_{hap}.tsv", sample=wildcards.sample, hap=wildcards.hap, contigs=CONTIGS)

def getSunkPos(wildcards):
  CONTIGS = glob_wildcards("results/{sample}/breaks/{contigs}_{hap}.sunkpos".format(sample=wildcards.sample, hap=wildcards.hap, contigs='{contigs}')).contigs
  return expand("results/{sample}/breaks/{contigs}_{hap}.sunkpos", sample=wildcards.sample, hap=wildcards.hap, contigs=CONTIGS)

def vizInputs(wildcards):
  # print(wildcards.runmode)
  if wildcards.runmode == 'user_bed':
    bed = manifest_df.at[wildcards.sample, f'{wildcards.hap}_bed']
  elif wildcards.runmode == 'gaps':
    bed = rules.slop_gaps.output.gaps_slop
  input_dict = {'bed': bed, 'rlen': rules.combine_ont.output.ONT_len, 'interout': rules.confirm_out.output.flag, 'pos_locs': rules.split_sunkpos.output.flag}
  if not pd.isnull(manifest_df.at[wildcards.sample, f'{wildcards.hap}_colortrack']):
    input_dict['colorbed'] = manifest_df.at[wildcards.sample, f'{wildcards.hap}_colortrack']
  #print(input_dict)
  return input_dict
