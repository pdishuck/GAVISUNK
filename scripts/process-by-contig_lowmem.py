#! /usr/bin/python
import networkx as nx
import graph_tool as gt
from graph_tool.topology import extract_largest_component

import timeit

from scipy.spatial.distance import pdist
import itertools as it

import numpy as npt
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

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("contig", help="contig to process")
    parser.add_argument("SUNKs", help=".loc file of SUNKs aligned to assembly")
    parser.add_argument("sunkpos", help=".sunkpos file of SUNK locations on ONT reads")
    parser.add_argument("rlen", help=".rlen file of ONT read lengths")
    # parser.add_argument("badsunks", help="list of SUNKs to exclude")
    parser.add_argument("plotdir", help="output directory for plots")
    parser.add_argument("outputfile", help="output file for intermediate")

    parser.add_argument("--minlen", help="minimum length filter for reads to visualize",default=10000,type=int)
    parser.add_argument("--opt_filt", help="apply optional filter",action='store_true')

    args = parser.parse_args()
    contig = args.contig
    ofile = args.outputfile


    lendf = pd.read_csv(args.rlen,sep="\t",header=None,names=['rname','len'],dtype={'rname':'string','len':'uint32'})
    lendf.set_index("rname",drop=True,inplace=True)
    minlen = args.minlen


    kmermergefile = args.SUNKs
    kmermerge = pd.read_csv(kmermergefile,sep="\t",header=None,names=['chrom','loc','kmer','ID'])# this could be split too
    kmermerge = kmermerge.drop_duplicates(subset='kmer')
    kmermerge['ID2'] = kmermerge['chrom'].astype(str) +":"+ kmermerge['ID'].astype(str)
    plotdir = args.plotdir + "/"
    sunkposfile = args.sunkpos
    sunkposcat = pd.read_csv(sunkposfile,sep="\t",header=None,names=['rname','pos','chrom','start','ID'],dtype={'rname':'string','pos':'uint32','chrom':'category','start':'uint32','ID':'uint32'})# 



    ## get list of bad sunks
    with open("badsunks.txt", "r") as f:
        badsunkin = set(f.read().splitlines())
    sunkposcat['ID2'] = sunkposcat['chrom'].astype(str) +":"+ sunkposcat['ID'].astype(str)

    print("sunkposcat len: ", len(sunkposcat))
    sunkposcat = sunkposcat.query("ID2 not in @badsunkin")
    print("sunkposcat len: ", len(sunkposcat))


    # contig = row.chrom
    # regstart = row.start
    # regend = row.end
    # xmin = regstart
    # xmax = regend

    # contig = 'chr20'
    kmer_groups = kmermerge.drop_duplicates(subset='ID') #  
    kmer_groups['dist'] = kmer_groups.ID.diff().fillna(kmer_groups.ID).astype(int) #rows)

    print("processing " + str(contig))
    # contig = bedreg1['chrom'].values[0]
    # xmin = bedreg1['start'].values[0]
    # xmax = bedreg1['end'].values[0]

    #Filter 1: only keep reads with at least two SUNKs within target region
    multisunk = sunkposcat.groupby('rname',as_index=False).agg({'ID':'nunique'}).sort_values(by='ID').query('ID > 1')
    print("Filter 1 multisunk len: ", len(multisunk))
    if len(multisunk) == 0:
      print("No viable reads for this region: "+str(contig))
      pd.DataFrame(list([contig])).to_csv(ofile,header=False,sep="\t",index=False)
      exit()
    #   continue

    multisunkset = set(multisunk.rname.tolist())

    ### RESTRICT TO REGION
    sunkpos3 = sunkposcat.query('rname in @multisunkset ').sort_values(by=['rname','start']) # reads with multiple SUNK group matches # & start >= @xmin & start <= @xmax
    #         sunkpos3 = sunkpos2.drop_duplicates(subset=['rname','ID'],keep='first') # one representation per read * SUNK group # done in kmerpos_annot script
    print(len(sunkpos3))
    print(len(pd.unique(sunkpos3.rname)))
    print(len(pd.unique(sunkposcat.rname)))


    sunkpos3b = sunkpos3.drop_duplicates()

    minlen = 10000
    minlenset = set(lendf.query("len >= @minlen").index.tolist())
    sub = sunkpos3b.query("rname in @minlenset") #sunkpos3





    # part2 bookmark

    # minlenset = set(lendf.query("len >= @minlen").index.tolist())


    ### keeep pos
    start = timeit.default_timer()

    print(len(sub))
    sub['start'] = sub['start'].astype('int64')
    sub['pos'] = sub['pos'].astype('int64')

    # rnamelist = list(pd.unique(distdf2.rname))
    # rnamelist = ['00060b17-6288-4b65-a67e-93e7498e0633']
    # rnamelist = ['fffdf958-07a2-4bc1-95aa-7c6b05cc70c6']
    print("rnamelist len: ", len(pd.unique(sub.rname)))
    outputs= []
    counter=0
    subgraphs=[]
    start2=timeit.default_timer()
    times=[]
    sub['start'] = sub['start'].astype('int64')
    sub['pos'] = sub['pos'].astype('int64')
    grouped = sub.groupby("rname") # .drop(columns='rname')
    for rname, g in grouped:
    #     if counter > 10:break
        counter +=1 
        start=timeit.default_timer()
        startarray = pdist(g[['start']])
        posarray = pdist(g[['pos']])
        signarray = pdist(g[['pos']], lambda u,v: u[0]>v[0] )
        idarray = list(it.combinations(g['ID'],2))
        locarray = list(it.combinations(g['pos'],2)) # keep track of ID-Pos corresponence
    
        diffarray = posarray/startarray
        mask = np.logical_not((diffarray >= 1.1) + (diffarray <= 0.9))
        if sum(mask) < 1: continue
        # only calculate on masked version (distdf2 equiv)
        signarray = signarray.astype(int)
        vals,counts = np.unique(signarray[mask], return_counts = True)
        trueOrient = vals[np.argmax(counts)]
        mask2 = np.logical_and(mask,(signarray == trueOrient))
    
        stack = np.column_stack([idarray,locarray])[mask2,:]
        sub2 = pd.DataFrame(np.column_stack([idarray,locarray])[mask2,:])
        sub2.columns = ['ID','ID2','pos1','pos2']
        mergedf = pd.DataFrame(np.row_stack([stack[:,[0,2]],stack[:,[1,3]]]))
        mergedf.columns = ['ID','pos']
    
        mergedf = pd.concat( [
            sub2[['ID','pos1']].rename(columns={'pos1':'pos'}),
            sub2[['ID2','pos2']].rename(columns={'ID2':'ID','pos2':'pos'})
            ]) # .drop_duplicates().query("ID > @xmin & ID < @xmax")
    
    
        multipos = mergedf.groupby("ID")['pos'].agg(nunique='nunique').query("nunique>1")
    
        if len(multipos)>=1:
            badrowlist = []
    #     multipos
    
            for i in multipos.index:
        #         i=38570942
                goodpos = mergedf.query("ID == @i")['pos'].value_counts().index[0] # any checks that this is good?
                badrowlist += list(sub2.query("(ID == @i) & not ((pos1 == @goodpos)|(pos2 == @goodpos))").index)
    
    
    
            sub2b = sub2.drop(badrowlist)
        else: sub2b = sub2
        
    
        mergedf2 = pd.concat( [
            sub2b[['ID','pos1']].rename(columns={'pos1':'pos'}),
            sub2b[['ID2','pos2']].rename(columns={'ID2':'ID','pos2':'pos'})
            ]) # .drop_duplicates().query("ID > @xmin & ID < @xmax")
    
        g = gt.Graph(directed=False)
        # g.vp.name = g.new_vertex_property("int64_t")
        name = g.add_edge_list(sub2b[['ID','ID2']].values,hashed=True, hash_type='int64_t') # eprops = [a_pos1,a_pos2]
        g.vp.name = name
        glarge = gt.topology.extract_largest_component(g)
    
        glarge.vp.name.get_array()
        df=pd.DataFrame(g.vp.name.get_array()[glarge.get_vertices()],columns=['ID'])
    #     print("num of SUNKs: ", len(df))
        df['rname'] = rname
        outputs.append(df)
        end=timeit.default_timer()
        times.append(end-start)
        if counter % 100000 == 0: print(counter,rname,end-start)

    if len(outputs) == 0:
      print("No good reads for this contig: "+str(contig))
      pd.DataFrame(list([contig])).to_csv(ofile,header=False,sep="\t",index=False)
      exit()
    else: 
      outputsdf = pd.concat(outputs)
      outputsdf.to_csv(ofile,header=False,sep="\t",index=False)
    end2=timeit.default_timer()
    print('overall time: ',end2-start2)


































#     startarray = sub.groupby("rname").apply(lambda g: pdist(g[['start']]))
#     posarray = sub.groupby("rname").apply(lambda g: pdist(g[['pos']]))
#     signarray = sub.groupby("rname").apply(lambda g: pdist(g[['pos']], lambda u,v: u[0]>v[0] ))
#     idarray = sub.groupby("rname").apply(lambda g: list(it.combinations(g['ID'],2)))
#     locarray = sub.groupby("rname").apply(lambda g: list(it.combinations(g['pos'],2))) # keep track of ID-Pos corresponence
#     startpos = pd.concat([startarray,posarray,idarray,locarray,signarray],axis=1)
# 
#     startpos['diff'] = startpos[1]/startpos[0]
#     startexp = startpos.explode([0,1,2,3,4,'diff'])
#     startexp.rename(columns={0:'exp',1:'dist',2:'IDs',3:'Poss',4:'distsign'},inplace=True)
#     startexp['rname'] = startexp.index
#     end = timeit.default_timer()
#     print(end-start)
# 
# 
#     distdf2 = startexp.query('not ((diff >= 1.1) | (diff <= 0.9)) ')
#     distdf2[['ID', 'ID2']] = pd.DataFrame(distdf2['IDs'].tolist(), index=distdf2.index)
#     distdf2[['pos1', 'pos2']] = pd.DataFrame(distdf2['Poss'].tolist(), index=distdf2.index)
#     print(len(distdf2))
# 
# 
# 
# 
# 
# 
#     ###
# 
# 
#     rnamelist = list(pd.unique(distdf2.rname))
#     # rnamelist = ['00060b17-6288-4b65-a67e-93e7498e0633']
#     # rnamelist = ['fffdf958-07a2-4bc1-95aa-7c6b05cc70c6']
#     print("rnamelist len: ", len(rnamelist))
#     outputs= []
#     counter=0
#     subgraphs=[]
#     start2=timeit.default_timer()
#     times=[]
#     grouped = distdf2.drop(columns='rname').groupby("rname")
#     for rname, sub in grouped:
#     #     if counter > 10:break
#         counter+=1
#         start=timeit.default_timer()
#     #     sub = distdf2.query("rname == @rname")
#         trueOrient = sub.distsign.value_counts().index[0]
#         sub2 = sub.query("distsign == @trueOrient")    
#         sub2.reset_index(inplace=True,drop=True)
# 
#         mergedf = pd.concat( [
#             sub2[['ID','pos1']].rename(columns={'pos1':'pos'}),
#             sub2[['ID2','pos2']].rename(columns={'ID2':'ID','pos2':'pos'})
#             ]) # .drop_duplicates().query("ID > @xmin & ID < @xmax")
# 
# 
#         multipos = mergedf.groupby("ID")['pos'].agg(nunique='nunique').query("nunique>1")
# 
#         if len(multipos)>=1:
#             badrowlist = []
#     #     multipos
# 
#             for i in multipos.index:
#         #         i=38570942
#                 goodpos = mergedf.query("ID == @i")['pos'].value_counts().index[0] # any checks that this is good?
#                 badrowlist += list(sub2.query("(ID == @i) & not ((pos1 == @goodpos)|(pos2 == @goodpos))").index)
# 
# 
#             #print("rows to remove: ",len(badrowlist))
# 
#             sub2b = sub2.drop(badrowlist)
#         else: sub2b = sub2
#         
# 
#         mergedf2 = pd.concat( [
#             sub2b[['ID','pos1']].rename(columns={'pos1':'pos'}),
#             sub2b[['ID2','pos2']].rename(columns={'ID2':'ID','pos2':'pos'})
#             ]) # .drop_duplicates().query("ID > @xmin & ID < @xmax")
# 
#         g = gt.Graph(directed=False)
#         # g.vp.name = g.new_vertex_property("int64_t")
#         name = g.add_edge_list(sub2b[['ID','ID2']].values,hashed=True, hash_type='int64_t') # eprops = [a_pos1,a_pos2]
#         g.vp.name = name
#         glarge = gt.topology.extract_largest_component(g)
# 
#         glarge.vp.name.get_array()
#         df=pd.DataFrame(g.vp.name.get_array()[glarge.get_vertices()],columns=['ID'])
#     #     print("num of SUNKs: ", len(df))
#         if len(df) != len(pd.unique(mergedf2['ID'])): 
#             #print("size diff:! ",len(df),len(pd.unique(mergedf2['ID'])))
#             subgraphs.append(rname)
#         df['rname'] = rname
#         outputs.append(df)
#         end=timeit.default_timer()
#         times.append(end-start)
#         if counter % 10000 == 0: print(counter,rname,end-start)
#     outputsdf = pd.concat(outputs)
#     end2=timeit.default_timer()
#     print('overall time: ',end2-start2)
# 




#    outputsdf.to_csv("splits/" + contig + ".outputsdf.gz",header=False,sep="\t",index=False)



    # made outputdf

    ### make contig-wide graph
    graphtimes=[]
    start = timeit.default_timer()
    graphtimes.append(start)
    idarray = pd.DataFrame(outputsdf.groupby("rname").apply(lambda g: list(it.combinations(g['ID'],2)))).explode([0])
    end = timeit.default_timer()
    graphtimes.append(end)
    print(start-end)
    idarray.rename(columns={0:'IDs',},inplace=True)
    idarray[['ID', 'ID2']] = pd.DataFrame(idarray['IDs'].tolist(), index=idarray.index)
    graphtimes=[]
    end = timeit.default_timer()
    print(start-end)
    graphtimes.append(end)

    g = gt.Graph(directed=False)
    # g.vp.name = g.new_vertex_property("int64_t")
    name = g.add_edge_list(idarray[['ID','ID2']].values,hashed=True, hash_type='int64_t') # eprops = [a_pos1,a_pos2]
    g.vp.name = name

    end = timeit.default_timer()
    print(start-end)
    graphtimes.append(end)

    # glarge = gt.topology.extract_largest_component(g)
    c = gt.topology.label_components(g)[0]

    end = timeit.default_timer()
    print(start-end)
    graphtimes.append(end)


    print("confirmed regions")
    regions = []
    for s in list(pd.unique(c.a)):
        color='lightgray'
        offset=counter
        sunks = g.vp.name.get_array()[gt.GraphView(g, vfilt=c.a == s).get_vertices()]
        if len(sunks) <= 2: continue
        start = int(sunks.min())
        end = int(sunks.max())
        span = end-start
        print(start,end,span,len(sunks))
        regions.append((start,end,span,len(sunks)))
    #     ax.plot([start,end],[rownum+offset,rownum+offset],color=color,linewidth=9,solid_capstyle='butt',zorder=1)    
        end = timeit.default_timer()
        counter +=1 

    regdf = pd.DataFrame(regions,columns=['start','end','span','sunks']).sort_values(by='start')    

    regdf['contig'] = contig
    regdf[['contig','start','end']].to_csv("breaks/" + contig + ".bed",header=False,sep="\t",index=False)

main()
