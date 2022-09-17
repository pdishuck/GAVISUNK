include: "gatherSplits.smk"
def gatherAsmBeds(wildcards):
  CONTIGS = glob_wildcards("results/{sample}/breaks/{contigs}_{hap}.sunkpos".format(sample=wildcards.sample, hap=wildcards.hap, contigs='{contigs}')).contigs
  return expand(rules.process_by_contig.output.bed, sample=wildcards.sample, hap=wildcards.hap, contigs=CONTIGS)

def getSunkLocs(wildcards):
  CONTIGS = glob_wildcards("results/{sample}/breaks/{contigs}_{hap}.loc".format(sample=wildcards.sample, hap=wildcards.hap, contigs='{contigs}')).contigs
  return expand("results/{sample}/breaks/{contigs}_{hap}.loc", sample=wildcards.sample, hap=wildcards.hap, contigs=CONTIGS)

def getInterOut(wildcards):
  CONTIGS = glob_wildcards("results/{sample}/inter_outs/{contigs}_{hap}.tsv".format(sample=wildcards.sample, hap=wildcards.hap, contigs='{contigs}')).contigs
  print("results/{sample}/inter_outs/{contigs}_{hap}.tsv")
  return expand("results/{sample}/inter_outs/{contigs}_{hap}.tsv", sample=wildcards.sample, hap=wildcards.hap, contigs=CONTIGS)

def getSunkPos(wildcards):
  CONTIGS = glob_wildcards("results/{sample}/breaks/{contigs}_{hap}.sunkpos".format(sample=wildcards.sample, hap=wildcards.hap, contigs='{contigs}')).contigs
  return expand("results/{sample}/breaks/{contigs}_{hap}.sunkpos", sample=wildcards.sample, hap=wildcards.hap, contigs=CONTIGS)

def getSunkPos2(wildcards):
  CONTIGS = glob_wildcards("results/{sample}/check_gaps/{contgap}_{hap}.gaps.bed".format(sample=wildcards.sample, hap=wildcards.hap, contgap='{contgap}')).contgap
  print(expand("results/{sample}/check_gaps/subset/{contgap}_{hap}_subset_region.sunkpos", sample=wildcards.sample, hap=wildcards.hap, contgap=CONTIGS))
  return expand(rules.add_reads.output.sunkpos_combine, sample=wildcards.sample, hap=wildcards.hap, contgap=CONTIGS)

def getInterOutsGaps(wildcards):
  CONTIGS = glob_wildcards("results/{sample}/check_gaps/{contgap}_{hap}.gaps.bed".format(sample=wildcards.sample, hap=wildcards.hap, contgap='{contgap}')).contgap
  return expand(rules.filter_no_reads.output.correct_df, sample=wildcards.sample, hap=wildcards.hap, contgap=CONTIGS)

def getGaps(wildcards):
  CONTIGS=glob_wildcards("results/{sample}/check_gaps/{contgap}_{hap}.gaps.bed".format(sample=wildcards.sample, hap=wildcards.hap, contgap='{contgap}')).contgap
  return expand("results/{sample}/check_gaps/bed_files/{contgap}_{hap}.merged.valid.bed", sample=wildcards.sample, hap=wildcards.hap, contgap=CONTIGS)

def vizInputsCheck(wildcards):
  file_check=checkpoints.split_gaps.get(sample=wildcards.sample, hap=wildcards.hap).output[1]
  if os.stat(file_check).st_size != 0:
    bed = rules.no_gaps.output.no_gaps
    interout = rules.no_gaps.output.no_gaps
    pos = rules.no_gaps.output.no_gaps
  else:
    bed = rules.slop_gaps_check.output.gaps_slop
    interout = rules.gather_checks.output.tsv_flag
    pos = rules.gather_checks.output.sunkpos_flag
  input_dict = {'bed': bed, 'rlen': rules.combine_ont.output.ONT_len, 'interout': interout, 'pos': pos, 'locs': rules.split_sunkpos.output.flag}
  print(input_dict)
  return input_dict

def vizInputs(wildcards):
  print(wildcards.runmode)
  if wildcards.runmode == 'user_bed':
    bed = manifest_df.at[wildcards.sample, f'{wildcards.hap}_bed']
  elif wildcards.runmode == 'gaps':
    bed = rules.slop_gaps.output.gaps_slop
  input_dict = {'bed': bed, 'rlen': rules.combine_ont.output.ONT_len, 'interout': rules.confirm_out.output.flag, 'pos': rules.split_sunkpos.output.flag, 'locs': rules.split_sunkpos.output.flag}
  if not pd.isnull(manifest_df.at[wildcards.sample, f'{wildcards.hap}_colortrack']):
    input_dict['colorbed'] = manifest_df.at[wildcards.sample, f'{wildcards.hap}_colortrack']
  return input_dict