import pandas as pd, numpy as np
import os
from sympy.solvers import solve, solveset
from sympy import Symbol, Reals
import pyranges as pr

import argparse

print(snakemake.input)



# genome_size=3.055e6
fai = pd.read_csv(snakemake.input.fai,sep='\t',header=None,names=['contig','length','cumlen','3','4'])
genome_size = fai['length'].sum()/1000
# print("assembly size: ", genome_size)
df = pd.read_csv(snakemake.input.bed, header=None, names="Chromosome Start End".split(), sep="\t")
df['type'] = "gap"
# df2 = pd.read_csv(breakdir + "hap" + hap + ".nodata.bed", header=None, names="Chromosome Start End".split(), sep="\t")
# df2['type'] = 'nodata'
# df3 = pd.concat([df,df2])
df3 = df
df3.reset_index(inplace=True)
breaks = pr.PyRanges(df3)




# breaks

# df3.eval("End - Start").sort_values()

# kmermergefile = os.path.dirname(snakemake.input.pos_locs) +"/"+ cont2 + "_" + snakemake.wildcards.hap + ".loc"

kmermerge = pd.read_csv(snakemake.input.locs,sep="\t",header=None,names=['chrom','loc','kmer','ID'])# this could be split too
kmermerge = kmermerge.drop_duplicates(subset='kmer')
kmermerge['ID2'] = kmermerge['chrom'].astype(str) +":"+ kmermerge['ID'].astype(str)
kmer_groups = kmermerge.drop_duplicates(subset='ID2') #  
kmer_groups['dist'] = kmer_groups.ID.diff().fillna(kmer_groups.ID).astype(int) #rows)
kmer_groups['dist'] = kmer_groups['dist'].map( lambda x: max(x,0) ) 


lendf = pd.read_csv(snakemake.input.rlen,sep="\t",header=None,names=['rname','len'],dtype={'rname':'string','len':'uint32'})
lendf.drop_duplicates(inplace=True)
lendf.set_index("rname",drop=True,inplace=True)
lendf.sort_values(by="len",ascending=False,inplace=True)
lendf['bp']=lendf['len'].cumsum()
# lendf['cov'] = lendf['bp']/3.055e9

lendf['kbp'] = lendf['len']/1000
lendf['kbp'] = lendf['kbp'].astype(int)

kbp_cnt = lendf['kbp'].value_counts(sort=False)


def cov_prob(intsize,readlen,readcnt,basequal=1): # everything in kbp ?
    if intsize > readlen:
#         print("no help here")
        return 1
    return  ( (1-   (  ((readlen-intsize) /genome_size) * (basequal**2) )     )**(readcnt) )


x = Symbol('x')
p=0.94 # read quality
q=1-p
r=int(snakemake.config['SUNK_len']) # kmer size

roots2 = solveset( 1 - x + q*((p)**(r)) * ((x)**(r+1)), x, domain=Reals)
pinv = 1/p
roots2 = [x for x in roots2  if x > 0]
roots2.sort(key=lambda e: abs(pinv - e))
assert(len(roots2) == 2)
roots3  = roots2[1]

print("Should be unequal: ",roots3, pinv)



n=30 # group size
qn = ( ( 1- p*roots3)/(q*(r + 1 - r*roots3))) * (1 / (roots3**(n+1)))
pn = float(1-qn)
# print(qn, pn) # odds perfect 20mer not seen in kmergroup length 22



covprobs = []
for intsize in range(1,3500):
    cum_prob = 1
#     counter = 0
    for i,x in kbp_cnt.iteritems():
        prob = cov_prob(intsize, i, x,pn) 
        cum_prob = cum_prob * prob
    covprobs.append((intsize,1-cum_prob))
#         counter += 1    
covs,probs = zip(*covprobs)



covprobsdict = {k:v for k,v in covprobs}
covprobsdict[0] = 1.0

for cov, prob in covprobs:
    if prob <=0.99: 
        print(cov,prob)
        break



max_gaps = []
max_gaps2 = []
for x in breaks.df.itertuples():
    ix,ix2,contig,regstart,regend,types = x
    xmin = regstart
    xmax = regend    
    kmers = kmer_groups.query("(chrom == @contig) & (ID > (@xmin-2)) & (ID < (@xmax+2))")
#     print(len(kmers))
    kmerlist = list(pd.unique(kmers.ID))       
    max_gaps.append(kmers['dist'].max())    
    max_gaps2.append((contig,regstart,regend,kmers['dist'].max(),kmers.iloc[kmers['dist'].argmax()].kmer)) 


# breaks

# kmer_groups.dist.max()



df1 = breaks.df
df1['max_gap'] = max_gaps
df1['covprob'] = df1['max_gap'].apply(lambda x: covprobsdict[int(x/1000)] )
# df1.sort_values(by='covprob').head(25)

df1.drop(columns='index').to_csv(snakemake.output.tsv,index=False,sep="\t")
