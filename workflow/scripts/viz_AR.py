#!/usr/bin/python
import pandas as pd
from matplotlib import (pyplot as plt,
                        lines)
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle, Wedge, Polygon, Rectangle
import seaborn as sns
import numpy as nptwi
import matplotlib.ticker as mtick
import argparse
print(pd.__version__)
import os
from ncls import NCLS

parser = argparse.ArgumentParser()

parser.add_argument("bed", help=".bed file of regions to visualize")
parser.add_argument("rlen", help=".rlen file of ONT read lengths")
parser.add_argument("plotdir", help="output directory for plots")
parser.add_argument("interoutdir", help="output directory for plots")
parser.add_argument("posdir", help="output directory for plots")
parser.add_argument("locsdir", help="output directory for plots")
parser.add_argument("sampname")
parser.add_argument("hap")
parser.add_argument("--minlen", help="minimum length filter for reads to visualize",default=10000,type=int)
parser.add_argument("--opt_filt", help="apply optional filter",action='store_true')
parser.add_argument("--colorbed", help="color track for output")



args = parser.parse_args()

regbedfile = args.bed
bedreg1 = pd.read_csv(regbedfile, delimiter="\t",encoding='utf-8',header=None)
header = ['chrom','start','end','name']
bedreg1.columns = header + [''] * (len(bedreg1.columns) - len(header))  

lendf = pd.read_csv(args.rlen,sep="\t",header=None,names=['rname','len'],dtype={'rname':'string','len':'uint32'})
lendf.drop_duplicates(inplace=True)
lendf.set_index("rname",drop=True,inplace=True)
minlen = args.minlen
opt_filter = args.opt_filt

if opt_filter:
    pct80 = 0.3333333333333333
bedreg1['chrom'].value_counts()
contigs = list(bedreg1['chrom'].value_counts().index)
print(len(contigs))
if args.colorbed == "NULL": args.colorbed = None
if args.colorbed is not None:
    dm = pd.read_csv(args.colorbed, delim_whitespace=True, names=['chr','chrStart','chrEnd','color'], header=None,dtype = {'chr':'string','chrStart':'int','chrEnd':'int','color':'string'})
    dm = dm[dm.chrStart >= 0]
    dm["y"] = 1
    dm["func"] = ""
    
for contig in contigs:
    cont2 = contig.replace("#","_")
    print(contig)
    
    bedreg= bedreg1.query('chrom == @contig')
    try: outputsdf= pd.read_csv(args.interoutdir + "/" + cont2 + "_" + args.hap + ".tsv" ,sep="\t",names=['ID','rname'],dtype = {'ID':'int','rname':'string'})
    except:
        print("outputsdf not found for", contig)
    #     continue
    outputsdf.set_index("rname",inplace=True,drop=True)
    if len(outputsdf) <2:
        print("NO READS FOR", contig)
    #     continue

    kmermergefile = args.locsdir +"/"+ cont2 + "_" + args.hap + ".loc"
    kmermerge = pd.read_csv(kmermergefile,sep="\t",header=None,names=['chrom','loc','kmer','ID'])# this could be split too
    kmermerge = kmermerge.drop_duplicates(subset='kmer')
    plotdir = args.plotdir + "/"
    sunkposfile = args.posdir+"/"+ cont2 + "_" + args.hap + ".sunkpos"
    sunkposcat = pd.read_csv(sunkposfile,sep="\t",header=None,names=['rname','pos','chrom','start','ID'],dtype={'rname':'string','pos':'int64','chrom':'category','start':'int64','ID':'int64'})


    kmermerge['ID2'] = kmermerge['chrom'].astype(str) +":"+ kmermerge['ID'].astype(str)
    
    sunkposcat['ID2'] = sunkposcat['chrom'].astype(str) +":"+ sunkposcat['ID'].astype(str)
    
    print("sunkposcat len: ", len(sunkposcat))
    #     sunkposcat = sunkposcat.query("ID2 not in @badsunkin")
    #     print("sunkposcat len: ", len(sunkposcat))    
    
    
    kmer_groups = kmermerge.drop_duplicates(subset='ID') #  
    kmer_groups['dist'] = kmer_groups.ID.diff().fillna(kmer_groups.ID).astype(int) #rows)
    
    multisunk = sunkposcat.groupby('rname',as_index=False).agg({'ID':'nunique'}).sort_values(by='ID').query('ID > 1')
    print("Filter 1 multisunk len: ", len(multisunk))
    if len(multisunk) == 0:
      print("No viable reads for this region: "+str(contig)+"_"+str(regstart)+"_"+str(regend))
    #   continue
    
    multisunkset = set(multisunk.rname.tolist())
    
    ### RESTRICT TO REGION
    sunkpos3 = sunkposcat.query('rname in @multisunkset ').sort_values(by=['rname','start']) # reads with multiple SUNK group matches # & start >= @xmin & start <= @xmax
    #         sunkpos3 = sunkpos2.drop_duplicates(subset=['rname','ID'],keep='first') # one representation per read * SUNK group # done in kmerpos_annot script
    print(len(sunkpos3))
    print(len(pd.unique(sunkpos3.rname)))
    print(len(pd.unique(sunkposcat.rname)))
    sunkpos3.set_index("rname",drop=True,inplace=True)

    for row in bedreg.itertuples():
        contig = row.chrom
        regstart = row.start
        regend = row.end
        xmin = regstart
        xmax = regend
        df = outputsdf.query("ID > @xmin & ID < @xmax")
        if len(df) == 0:
            print("NO HITS FOR", row)
            continue
        kmerlist = list(pd.unique(kmer_groups.query("ID > @xmin & ID < @xmax").ID))
        kmerset = set(kmerlist)
        never_seen2 = kmerset - set(df.ID)
        print(len(never_seen2))
        constart = df.ID.min()
        conend = df.ID.max()
        maxes={}
        first=True
        fig,ax = plt.subplots(figsize=(20,9))
        fig.subplots_adjust(right=0.75)
        ax.set_xlim(regstart,regend)
        intervals = []
        poss = []
        for rname in list(pd.unique(df.index)):
            try: rlen = lendf.loc[rname].len
            except:
                print("no rlen for: ", rname)
                continue
            sub = df.loc[rname]
            try: idlist = set(sub.ID)
            except: idlist = set([sub.ID])
            possub = sunkpos3.loc[rname].query("ID in @idlist")
            posfirst = possub.iloc[0]
            poslast = possub.iloc[-1]
            forward = True
            badlocs_trans = []
            if posfirst.pos > poslast.pos: forward=False
            if forward:
                posstart = posfirst.start - posfirst.pos
                posend =  poslast.start + (rlen - poslast.pos)  
            else:
                posend = poslast.pos + poslast.start
                posstart =  posfirst.start - (rlen - posfirst.pos)
            start = posstart
            end = posend
            intervals.append((start,end))
            poss.append(list(idlist))

        interval_df = pd.DataFrame(intervals,columns=['Start','End'],)
        ncls = NCLS(interval_df.Start,interval_df.End,interval_df.index)
        row_asn = dict()
        interval_df.sort_values(by='Start',inplace=True)

        rows_used = set([0])
        for i,s,e in interval_df.itertuples():
            overlap = ncls.find_overlap(s,e)
            rows_overlap = set()
            for i2,s2,e2 in overlap:
                rows_overlap.add( row_asn.get(e2,0))
            rows_avail = rows_used - rows_overlap
            if len(rows_avail) == 0:
                row_pick = sorted(rows_used)[-1]+1
                rows_used.add(row_pick)
                row_asn[i] = row_pick
            else:
                row_asn[i] = sorted(rows_avail)[0]
        interval_df['row_asn'] = interval_df.index.map(lambda x: row_asn[x])
        patches = []
        ticks = []
        for i,s,e,r in interval_df.itertuples():
            patches.append(Rectangle((s, r), e-s, 0.7))
            tick = [Rectangle((x, r), 1, 0.7)  for x in poss[i]]
            ticks = ticks + tick
        never_seen3 = [x for x in never_seen2 if not ((x<constart) | (x>conend))]

        if args.colorbed is not None:
            patches2=[]
            for r in dm[['chrStart','chrEnd','color','chr']].query('chr==@contig & chrEnd > @xmin & chrStart < @xmax ').itertuples():
                polygon = Polygon([[r[1],-4],[r[1],-5],[r[2],-5],[r[2],-4]], True, color=r[3],alpha=1)
                patches2.append(polygon)
            p = PatchCollection(patches2, alpha=1, match_original=True)
            ax.add_collection(p)
        ax.add_collection(PatchCollection(patches,color='lightgray',zorder=1,aa=True,edgecolor='w',linewidth=0.01))
        ax.add_collection(PatchCollection(ticks,color='k',zorder=7,linewidth=0.5,edgecolor='k',aa=True))
        ax.scatter(x=never_seen3,y=[-2.5]*len(never_seen3),color='k',s=55,marker="|",facecolors=None,linewidths=0.4,zorder=2)
        ax.scatter(x=kmerlist,y=[-1.5]*len(kmerlist),color='k',s=55,marker="|",facecolors=None,linewidths=0.4,zorder=2)
        fmt = '{x:,.0f}'
        tick = mtick.StrMethodFormatter(fmt)
        ax.xaxis.set_major_formatter(tick) 
        ax.set_ylabel("ONT Read Depth")
        ax.set_xlabel("Contig coordinate")
        plt.rcParams['svg.fonttype'] = 'none'
        print(args.plotdir + "/" + args.sampname +"_" + args.hap + "_" +cont2+"_"+str(regstart)+"_"+str(regend)+".svg")
        plt.savefig(args.plotdir + "/" + args.sampname +"_" + args.hap + "_" +cont2+"_"+str(regstart)+"_"+str(regend)+".svg",format="svg",pad_inches=0,bbox_inches='tight')
        plt.savefig(args.plotdir + "/" + args.sampname +"_" + args.hap + "_" +cont2+"_"+str(regstart)+"_"+str(regend)+".png",pad_inches=0,bbox_inches='tight',dpi=300)
        plt.savefig(args.plotdir + "/" + args.sampname +"_" + args.hap + "_" +cont2+"_"+str(regstart)+"_"+str(regend)+".pdf",pad_inches=0,bbox_inches='tight')
        plt.close()
        print("Plotting complete for " +contig+"_"+str(regstart)+"_"+str(regend))
