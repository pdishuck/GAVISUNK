import pandas as pd
import os


def splitlocs(hap_type,group_df, dirname, exten):
  for name, contig_group in group_df:
    contig_name = name.replace("#", "_")
    output_path = f'{dirname}/{contig_name}_{hap_type}.{exten}'
    print(output_path)
    contig_group.to_csv(output_path, header=None, index=None, sep="\t")
  return(True)


ONT_pos_df = pd.read_csv(snakemake.input.ONT_pos, header=None, sep="\t", names=['read','read_pos','contig','contig_start','contig_stop'])
kmer_loc = pd.read_csv(snakemake.input.kmer_loc, header=None, sep="\t", names=['contig','pos1','SUNK','pos2'])
contig_list = ONT_pos_df['contig'].unique().tolist()
kmer_loc_subset = kmer_loc[kmer_loc['contig'].isin(contig_list)]
dir_name = os.path.dirname(snakemake.output.flag)  
group_ONT = ONT_pos_df.groupby(['contig'])
group_loc = kmer_loc_subset.groupby(['contig'])
output_complete = splitlocs(snakemake.wildcards.hap,group_ONT, dir_name,"sunkpos")
output_complete_2 = splitlocs(snakemake.wildcards.hap,group_loc, dir_name,"loc")