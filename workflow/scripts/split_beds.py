import pandas as pd

bed_df = pd.read_csv(snakefile.input.bed, header=None, sep="\t", names=['contig','start','stop','region_name'])
bed_df_group = bed_df.groupby(['region_name'])
for (splitno, split) in bed_df_group:
  split.to_csv(f'results/{wildcards.sample}/align_to_ref/gene_beds/{wildcards.sample}_{wildcards.hap}_to_ref_{splitno}.bed', header=None, sep="\t", index=None )
