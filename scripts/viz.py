import networkx as nx
# import graph_tool as gt
import timeit

from scipy.spatial.distance import pdist
import itertools as it

import numpy as np
import pandas as pd, numpy as np
from matplotlib import (pyplot as plt,
                        lines)
import seaborn as sns
import timeit

from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt

import seaborn as sns
import numpy as np
 
from matplotlib import (pyplot as plt,
                        lines)
import matplotlib.ticker as mtick
import argparse

print(pd.__version__)


parser = argparse.ArgumentParser()

parser.add_argument("bed", help=".bed file of regions to visualize")
parser.add_argument("SUNKs", help=".loc file of SUNKs aligned to assembly")
parser.add_argument("sunkpos", help=".sunkpos file of SUNK locations on ONT reads")
parser.add_argument("rlen", help=".rlen file of ONT read lengths")
parser.add_argument("plotdir", help="output directory for plots")
parser.add_argument("--minlen", help="minimum length filter for reads to visualize",default=10000,type=int)
parser.add_argument("--opt_filt", help="apply optional filter",action='store_true')

args = parser.parse_args()


# variables
#asm=
#ont=
#hap=

kmermergefile = args.SUNKs
kmermerge = pd.read_csv(kmermergefile,sep="\t",header=None,names=['chrom','loc','kmer','ID'])# 
kmermerge = kmermerge.drop_duplicates(subset='kmer')
plotdir = args.plotdir
sunkposfile = args.sunkpos
sunkposcat = pd.read_csv(sunkposfile,sep="\t",header=None,names=['rname','pos','chrom','start','ID'],dtype={'rname':'string','pos':'uint32','chrom':'category','start':'uint32','ID':'uint32'})# 

lendf = pd.read_csv(args.rlen,sep="\t",header=None,names=['rname','len'],dtype={'rname':'string','len':'uint32'})
lendf.set_index("rname",drop=True,inplace=True)
minlen = args.minlen
opt_filter = args.opt_filt
if opt_filter:
    pct80 = 0.3333333333333333

# Dupmask file
# winnowmap depth input


# import region bed
regbedfile = args.bed
bedreg = pd.read_csv(regbedfile, delimiter="\t",encoding='utf-8',header=None)
header = ['chrom','start','end','orient','label']
bedreg.columns = header + [''] * (len(bedreg.columns) - len(header))

# LOOP OVER BED REGIONS
for row in bedreg.itertuples():
    contig = row.chrom
    regstart = row.start
    regend = row.end
    xmin = regstart
    xmax = regend

    kmer_groups = kmermerge.query("chrom == @contig").drop_duplicates(subset='ID')
    kmer_groups['dist'] = kmer_groups.ID.diff().fillna(kmer_groups.ID).astype(int) #rows)
    
    print("processing " + str(contig)+"_"+str(regstart)+"_"+str(regend))
    # contig = bedreg1['chrom'].values[0]
    # xmin = bedreg1['start'].values[0]
    # xmax = bedreg1['end'].values[0]

    multisunk = sunkposcat.query("chrom == @contig & start >= @xmin & start <= @xmax").groupby('rname',as_index=False).agg({'ID':'nunique'}).sort_values(by='ID').query('ID > 1')
    print("multisunk len: ", len(multisunk))
    # if len(multisunk) == 0: continue
    multisunkset = set(multisunk.rname.tolist())


    ### RESTRICT TO REGION
    sunkpos3 = sunkposcat.query('rname in @multisunkset ').sort_values(by=['rname','start']) # reads with multiple SUNK group matches # & start >= @xmin & start <= @xmax
    #         sunkpos3 = sunkpos2.drop_duplicates(subset=['rname','ID'],keep='first') # one representation per read * SUNK group # done in kmerpos_annot script
    print(len(sunkpos3))
    print(len(pd.unique(sunkpos3.rname)))

    # keep only reads where plurality of sunkpos on contig of interest
    counts = sunkpos3.groupby(['rname','chrom'],as_index=False,dropna=False).size().sort_values(by=['rname','size'],ascending=False)
    maxset = set(counts.drop_duplicates(subset=['rname'],keep='first').query('chrom == @contig').rname.tolist())
    print(len(maxset))
    sunkpos3b = sunkposcat.query('rname in @maxset and chrom==@contig').sort_values(by=['rname','start']) # reads with multiple SUNK group matches
#     del(sunkpos3)
    # del(sunkposcat)
    print(len(sunkpos3b))
    print(len(pd.unique(sunkpos3b.rname)))


    # part2 bookmark

    minlenset = set(lendf.query("len >= @minlen").index.tolist())


    ### keeep pos
    start = timeit.default_timer()
    sub = sunkpos3b.query("rname in @minlenset") #sunkpos3
    print(len(sub))
    sub['start'] = sub['start'].astype('int64')
    sub['pos'] = sub['pos'].astype('int64')
    startarray = sub.groupby("rname").apply(lambda g: pdist(g[['start']]))
    posarray = sub.groupby("rname").apply(lambda g: pdist(g[['pos']]))
    signarray = sub.groupby("rname").apply(lambda g: pdist(g[['pos']], lambda u,v: u[0]>v[0] ))
    idarray = sub.groupby("rname").apply(lambda g: list(it.combinations(g['ID'],2)))
    locarray = sub.groupby("rname").apply(lambda g: list(it.combinations(g['pos'],2))) # keep track of ID-Pos corresponence
    startpos = pd.concat([startarray,posarray,idarray,locarray,signarray],axis=1,names=['exp','dist','IDs','Poss','distsign'])
    startpos['diff'] = startpos[1]/startpos[0]
    startexp = startpos.explode([0,1,2,3,4,'diff'])
    startexp.rename(columns={0:'exp',1:'dist',2:'IDs',3:'Poss',4:'distsign'},inplace=True)
    startexp['rname'] = startexp.index
    end = timeit.default_timer()
    print(end-start)


    distdf2 = startexp.query('not ((diff >= 1.1) | (diff <= 0.9)) ')
    distdf2[['ID', 'ID2']] = pd.DataFrame(distdf2['IDs'].tolist(), index=distdf2.index)
    distdf2[['pos1', 'pos2']] = pd.DataFrame(distdf2['Poss'].tolist(), index=distdf2.index)
    print(len(distdf2))


    # mergedf = pd.concat( [
    #     distdf2[['rname','ID','pos1']].rename(columns={'pos1':'pos'}),
    #     distdf2[['rname','ID2','pos2']].rename(columns={'ID2':'ID','pos2':'pos'})
    #     ]).drop_duplicates().query("ID > @xmin & ID < @xmax")


    # OPTIONAL FILTER
    if opt_filter:
        print("APPLYING OPTIONAL FILTER")
        cands = pd.concat([distdf2[['rname','ID']],distdf2[['rname','ID2']].rename(columns={'ID2':'ID'})]).drop_duplicates()
        cands2 = pd.concat([distdf.merge(cands,on=['rname','ID'],how='inner'),distdf.merge(cands.rename(columns={'ID':'ID2'}),on=['rname','ID2'],how='inner')]).drop_duplicates()
        def gooddiff(x):
           return sum([v > 1.02 or v < 0.98 for v in x])/len(x)
        diffpct = cands2.groupby('rname',as_index=False).agg({'diff':lambda x: gooddiff(x)})
        goodreads_alt = set(diffpct.query('diff <= @pct80')['rname']) # pct of 
        distdf2 = distdf2.query("rname in @goodreads_alt")

    # get longest fully connected subgraph for each read

    start = timeit.default_timer()
    outputs = []
    skip=0
    # df=distdf2
    rnamelist = list(pd.unique(distdf2.rname))
    # rnamelist = ['00060b17-6288-4b65-a67e-93e7498e0633']
    # rnamelist = ['fffdf958-07a2-4bc1-95aa-7c6b05cc70c6']
    print(len(rnamelist))
    for rname in rnamelist:
        sub = distdf2.query("rname == @rname")

        trueOrient = sub.distsign.value_counts().index[0]
        sub2 = sub.query("distsign == @trueOrient")    
        
        G3 = nx.convert_matrix.from_pandas_edgelist(sub2,'ID','ID2',['pos1','pos2'])
    #     S3 = G3.subgraph(max(nx.connected_components(G3), key=len))
        S3 = G3.subgraph(max(nx.clique.find_cliques(G3), key=len)) # fully connected subgraphs only
    #     [g for g in nx.clique.find_cliques(G3) if len(g) > 2]
        if len(S3) < 2: 
            skip+=1
            continue
        S3pd = nx.convert_matrix.to_pandas_edgelist(S3,'ID','ID2')
        S3pd2 = pd.concat([S3pd[['ID']],S3pd[['ID2']].rename(columns={'ID2':"ID"})]).drop_duplicates()
        S3pd2['rname'] = rname
        outputs.append(S3pd2)
    outputsdf = pd.concat(outputs)
    print(len(outputsdf))
    print("skipped ", skip)
    print(len(pd.unique(outputsdf['rname'])))
    end = timeit.default_timer()
    print(end-start) # 



    # PLOT

    df = outputsdf.query("ID > @xmin & ID < @xmax")

    kmerlist = list(pd.unique(kmer_groups.query("ID > @xmin & ID < @xmax").ID))       
    kmerset = set(kmerlist)
    never_seen2 = kmerset - set(df.ID)
    print(len(never_seen2))

    constart = df.ID.min()
    conend = df.ID.max()
    # plt.figure(figsize=(25,15))
    # plt.grid(False)
    maxes={}
    first=True

    fig,ax = plt.subplots(figsize=(20,9))
    fig.subplots_adjust(right=0.75)
    ax.set_xlim(regstart,regend)
    # twin1 = ax.twinx()
    # fig.yticks([])
    # ax.set_yticks([])
    # twin2 = ax.twinx()

    # twin1.spines.right.set_position(('axes', 1.2))
    counter =0

    for rname in list(pd.unique(df.sort_values(by='ID')['rname'])):
        counter +=1

        try: rlen = lendf.loc[rname].len
        except:
            print("no rlen for: ", rname)
            continue
    #     if counter > 5: continue
        sub = df.query('rname == @rname')
        idlist = set(sub.ID)    
        possub = sunkpos3.query('rname == @rname & ID in @idlist').sort_values(by='ID')
    #     possuball = sunkposcat.query('rname == @rname ').drop_duplicates(subset=['chrom','ID'])
    #     badlocs_same = possuball.query(" chrom==@contig & (not ID in @idlist) ")    
    #     badlocs_diff = possuball.query(" chrom != @contig ")    
        posfirst = possub.iloc[0]
        poslast = possub.iloc[-1]
        forward = True
        badlocs_trans = []

        if posfirst.pos > poslast.pos: forward=False


        if forward:
            posstart = posfirst.start - posfirst.pos
            posend =  poslast.start + (rlen - poslast.pos)  
    #         badlocs_same_trans = [ (posfirst.start + (x - posfirst.pos)) for x in badlocs_same.pos.tolist() ]
    #         badlocs_diff_trans = [ (posfirst.start + (x - posfirst.pos)) for x in badlocs_diff.pos.tolist() ]
        else:
            posend = poslast.pos + poslast.start
            posstart =  posfirst.start - (rlen - posfirst.pos)
    #         badlocs_same_trans = [ (posfirst.start - (x - posfirst.pos)) for x in badlocs_same.pos.tolist() ]
    #         badlocs_diff_trans = [ (posfirst.start - (x - posfirst.pos)) for x in badlocs_diff.pos.tolist() ]
    #     print(posstart,posend)

        start = posstart
        end = posend

        if first:
            rownum=0
            maxes[rownum] = end
            first=False
        else:
            for x in range(len(maxes)+1):
                if x in maxes.keys():
                    prevmax = maxes[x]
    #                 print(contig,x,start,prevmax)
                    if start>prevmax:
                        rownum=x
                        maxes[x] = end
                        break
                else:
                    maxes[x] = end
                    rownum=x
    #     print("X ",contig,rownum,start,maxes[rownum])




        ax.plot([posstart,posend],[rownum+1,rownum+1],color='lightgray',linewidth=9,solid_capstyle='butt',zorder=1)
        ax.scatter(data=sub,x='ID',y=[rownum+1]*len(sub),color='k',s=55,marker="|",facecolors=None,linewidths=0.4,zorder=4)
    #     ax.scatter(x=badlocs_same_trans,y=[rownum]*len(badlocs_same_trans),color='r',s=55,marker="|",facecolors=None,linewidths=0.4,zorder=3)
    #     ax.scatter(x=badlocs_diff_trans,y=[rownum]*len(badlocs_diff_trans),color='g',s=55,marker="|",facecolors=None,linewidths=0.4,zorder=2)
        rownum = 0
    never_seen3 = [x for x in never_seen2 if not ((x<constart) | (x>conend))]
    ax.scatter(x=never_seen3,y=[-2.5]*len(never_seen3),color='k',s=55,marker="|",facecolors=None,linewidths=0.4,zorder=2) 
    ax.scatter(x=kmerlist,y=[-1.5]*len(kmerlist),color='k',s=55,marker="|",facecolors=None,linewidths=0.4,zorder=2)
    # # plt.scatter(data=graph3,x='ID',y='counter',color='red',s=2)
    fmt = '{x:,.0f}'
    tick = mtick.StrMethodFormatter(fmt)
    ax.xaxis.set_major_formatter(tick) 


    ax.set_ylabel("ONT Read Depth")
    ax.set_xlabel("Contig coordinate")

    plt.savefig(plotdir+contig+"_"+str(regstart)+"_"+str(regend)+"_v3.png",dpi=300,pad_inches=0,bbox_inches='tight')
    plt.savefig(plotdir+contig+"_"+str(regstart)+"_"+str(regend)+"_v3.pdf",dpi=300,pad_inches=0,bbox_inches='tight')
    plt.savefig(plotdir+contig+"_"+str(regstart)+"_"+str(regend)+"_v3.svg",dpi=300,pad_inches=0,bbox_inches='tight')  

    print("Plotting complete for " +contig+"_"+str(regstart)+"_"+str(regend))
