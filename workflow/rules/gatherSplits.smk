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

