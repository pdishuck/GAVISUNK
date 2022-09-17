import pandas as pd

bed_df = pd.read_csv(snakemake.input.gaps_bed, header=None, sep="\t", names=['contig','start','stop'])
if len(bed_df) == 0:
  pd.DataFrame(list(['no_gaps'])).to_csv(snakemake.output.define_path,header=False,sep="\t",index=False)
  exit()

bed_df_group = bed_df.groupby(['contig'])
flag = pd.DataFrame()
for (splitno, split) in bed_df_group:
  split.to_csv(f'results/{snakemake.wildcards.sample}/check_gaps/{splitno}_{snakemake.wildcards.hap}.gaps.bed', header=None, sep="\t", index=None )
print(snakemake.output.flag)
#flag.to_csv(snakemake.output.flag, header=None, sep="\t", index=None)

