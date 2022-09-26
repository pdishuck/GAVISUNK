import pandas as pd
from matplotlib import (pyplot as plt,
                        lines)
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle, Wedge, Polygon, Rectangle
import seaborn as sns
import numpy as nptwi
import matplotlib.ticker as mtick
import argparse
#print(pd.__version__)
import os
from ncls import NCLS

print(snakemake.input)

runmode=snakemake.wildcards.runmode
opt_keys = list(snakemake.input.keys())

regbedfile = snakemake.input.bed
if os.stat(regbedfile).st_size == 0:
    print("empty gap file")
    exit()

bedreg1 = pd.read_csv(regbedfile, delimiter="\t",encoding='utf-8',header=None)
header = ['chrom','start','end']
bedreg1.columns = header + [''] * (len(bedreg1.columns) - len(header))  

lendf = pd.read_csv(snakemake.input.rlen,sep="\t",header=None,names=['rname','len'],dtype={'rname':'string','len':'uint32'})
lendf.drop_duplicates(inplace=True)
lendf.set_index("rname",drop=True,inplace=True)


bedreg1['chrom'].value_counts()
contigs = list(bedreg1['chrom'].value_counts().index)
# print(len(contigs))

if 'colorbed' in opt_keys:
    dm = pd.read_csv(snakemake.input.colorbed, delim_whitespace=True, names=['chr','chrStart','chrEnd','color'], header=None,dtype = {'chr':'string','chrStart':'int','chrEnd':'int','color':'string'})
    dm = dm[dm.chrStart >= 0]
    dm["y"] = 1
    dm["func"] = ""





# sunkposfile = os.path.dirname(snakemake.input.pos_locs)+"/"+ cont2 + "_" + snakemake.wildcards.hap + ".sunkpos"
sunkposfile = snakemake.input.pos_locs
sunkposcat = pd.read_csv(sunkposfile,sep="\t",header=None,names=['rname','pos','chrom','start','ID'],dtype={'rname':'string','pos':'int64','chrom':'category','start':'int64','ID':'int64'})

sunkposcat.set_index("rname",drop=True,inplace=True)    
sunkposcat['ID2'] = sunkposcat['chrom'].astype(str) +":"+ sunkposcat['ID'].astype(str)

for contig in contigs:

    cont2 = contig.replace("#","_")
    bedreg= bedreg1.query('chrom == @contig')

    #print(type(snakemake.input.interout))
    #my_str=os.path.dirname(snakemake.input.interout)
    #print(my_str)

    try: outputsdf= pd.read_csv(os.path.dirname(snakemake.input.interout) + "/" + cont2 + "_" + snakemake.wildcards.hap + ".tsv" ,sep="\t",names=['ID','rname'],dtype = {'ID':'int','rname':'string'})
    except:
        print("outputsdf not found for", contig)

    outputsdf.set_index("rname",inplace=True,drop=True)

    if len(outputsdf) <2:
        print("NO READS FOR", contig)


    kmermergefile = os.path.dirname(snakemake.input.splits) +"/"+ cont2 + "_" + snakemake.wildcards.hap + ".loc"
    kmermerge = pd.read_csv(kmermergefile,sep="\t",header=None,names=['chrom','loc','kmer','ID'])# this could be split too
    kmermerge = kmermerge.drop_duplicates(subset='kmer')
    plotdir = os.path.dirname(snakemake.output.flag) + "/"
    #print(plotdir)

    kmermerge['ID2'] = kmermerge['chrom'].astype(str) +":"+ kmermerge['ID'].astype(str)

  
    kmer_groups = kmermerge.drop_duplicates(subset='ID') #  
    kmer_groups['dist'] = kmer_groups.ID.diff().fillna(kmer_groups.ID).astype(int) #rows)
    
    sunkpos_sub = sunkposcat.query("chrom == @contig")
    multisunk = sunkpos_sub.groupby(sunkpos_sub.index,as_index=True).agg({'ID':'nunique'}).sort_values(by='ID').query('ID > 1')
    
    if len(multisunk) == 0:
      print("No viable reads for this region: "+str(contig)+"_"+str(regstart)+"_"+str(regend))

    multisunklist = list(multisunk.index)
    multisunkset = set(multisunklist)
    
    ### RESTRICT TO REGION
    # sunkpos3 = sunkposcat.loc[multisunklist] # reads with multiple SUNK group matches # & start >= @xmin & start <= @xmax
# orig sunkposfile stuff

    sunkposfile_split = os.path.dirname(snakemake.input.splits)+"/"+ cont2 + "_" + snakemake.wildcards.hap + ".sunkpos"
    sunkposcat_split = pd.read_csv(sunkposfile_split,sep="\t",header=None,names=['rname','pos','chrom','start','ID'],dtype={'rname':'string','pos':'int64','chrom':'category','start':'int64','ID':'int64'})
    sunkposcat_split['ID2'] = sunkposcat_split['chrom'].astype(str) +":"+ sunkposcat_split['ID'].astype(str)
    sunkpos3 = sunkposcat_split.query('rname in @multisunkset ').sort_values(by=['rname','start']) # reads with multiple SUNK group matches # & start >= @xmin & start <= @xmax
    sunkpos3.set_index("rname",drop=True,inplace=True)



    # sunkpos3.set_index("rname",drop=True,inplace=True)

    for row in bedreg.itertuples():
        contig = row.chrom
        regstart = row.start
        regend = row.end
        xmin = regstart
        xmax = regend
        breakloc = int((xmin + xmax)/2)
        rightset = set(outputsdf.query("ID > @breakloc").index)
        # print("skipped ", skip)
        #         print(len(pd.unique(outputsdf['rname'])))
        #     end = timeit.default_timer()
        #     print(end-start) # 


        sunklines=True

        # PLOT

        df = outputsdf.query("ID > @xmin & ID < @xmax")


        if len(df) == 0:
            print("NO HITS FOR", row)
        #     continue
        #         df.set_index("rname",inplace=True,drop=True)


        kmerlist = list(pd.unique(kmer_groups.query("ID > @xmin & ID < @xmax").ID))
        kmerset = set(kmerlist)
        never_seen2 = kmerset - set(df.ID)


        constart = df.ID.min()
        conend = df.ID.max()
        # plt.figure(figsize=(25,15))
        # plt.grid(False)
        maxes={}
        first=True

        fig,ax = plt.subplots(figsize=(20,9))
        fig.subplots_adjust(right=0.75)
        ax.set_xlim(regstart,regend)

        intervals = []
        valid_intervals = []
        poss = []
        poss_bad = []

        for rname in list(pd.unique(df.index)):
            try: rlen = lendf.loc[rname].len
            except:
                print("no rlen for: ", rname)
                if runmode == 'user_bed':
                    continue
            sub = df.loc[rname]
            
            try: idlist = set(sub.ID)
            except: idlist = set([sub.ID])

            possub = sunkpos3.loc[rname].query("ID in @idlist")
            posfirst = possub.iloc[0]
            poslast = possub.iloc[-1]
            forward = True
            badlocs_trans = []

            possuball = sunkposcat.loc[rname].drop_duplicates(subset=['chrom','ID'])
            badlocs_same = possuball.query(" chrom==@contig & (not ID in @idlist) ")    
        #     badlocs_diff = possuball.query(" chrom != @contig ")    
            badlocs_diff = possuball.query(" chrom != @contig  ")  
            badlocs_diff_most = badlocs_diff.chrom.value_counts().index[0]
            badlocs_diff = badlocs_diff.query("chrom == @badlocs_diff_most") 
        #     if len(badlocs_diff) > 0 : break


            if posfirst.pos > poslast.pos: forward=False


            if forward:
                posstart = posfirst.start - posfirst.pos
                posend =  poslast.start + (rlen - poslast.pos)
                badlocs_same_trans = [ (posfirst.start + (x - posfirst.pos)) for x in badlocs_same.pos.tolist() ]
                badlocs_diff_trans = [ (posfirst.start + (x - posfirst.pos)) for x in badlocs_diff.pos.tolist() ]
            else:
                posend = poslast.pos + poslast.start
                posstart =  posfirst.start - (rlen - posfirst.pos)
                badlocs_same_trans = [ (posfirst.start - (x - posfirst.pos)) for x in badlocs_same.pos.tolist() ]
                badlocs_diff_trans = [ (posfirst.start - (x - posfirst.pos)) for x in badlocs_diff.pos.tolist() ]
        #     print(posstart,posend)

            start = posstart
            end = posend
            
            if runmode == 'gaps':
                plotgroup=1
                columns = ['Start','End','Plotgroup']
                if rname in rightset: plotgroup=0
                intervals.append((start,end,min(idlist),max(idlist),plotgroup))
            elif runmode == 'user_bed':
                intervals.append((start,end,min(idlist),max(idlist),plotgroup))
                columns = ['Start','End']
            poss_bad.append(badlocs_diff_trans)
        #     valid_intervals.append((min(idlist),max(idlist),plotgroup))
            
            poss.append(list(idlist))

        interval_df = pd.DataFrame(intervals,columns=['Start','End','valid_Start','valid_End','Plotgroup'],)
        # valid_interval_df = pd.DataFrame(valid_intervals,columns=['Start','End','Plotgroup'],)


        ncls = NCLS(interval_df.Start,interval_df.End,interval_df.index)

        row_asn = dict()


        # interval_df.sort_values(by=['Plotgroup','Start'],inplace=True)
        interval_df.sort_values(by=['Start'],inplace=True)

        rows_used = set([0])
        if runmode == 'user_bed':
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
        elif runmode == 'gaps':
            for pgi in [0,1]:
                intv_sub = interval_df.query("Plotgroup == @pgi")
                intv_sub['End'] = intv_sub['End'] + 500 # Minimum separation between reads
                if pgi==1:
                    intv_sub = intv_sub.sort_values(by=['End'],ascending=False)
                for i,s,e,sv,ev,pg in intv_sub.itertuples():

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
                    maxrow_prev = max(row_asn.values())

        interval_df['row_asn'] = interval_df.index.map(lambda x: row_asn[x])

        patches = []
        valid_patches = []                     
        ticks = []
        bad_ticks = []
        for i,s,e,sv,ev,pg,r in interval_df.itertuples():
            patches.append(Rectangle((s, r+0.2), e-s, 0.4))
            valid_patches.append(Rectangle((sv, r), ev-sv, 0.7))                         
            tick = [Rectangle((x, r), 1, 0.7)  for x in poss[i]]
            bad_tick = [Rectangle((x, r+0.2), 1, 0.4)  for x in poss_bad[i]]
            bad_ticks = bad_ticks + bad_tick
            ticks = ticks + tick

        never_seen3 = [x for x in never_seen2 if not ((x<constart) | (x>conend))]

        if 'colorbed' in opt_keys:
            # print("COLORBED!!!")
            patches2=[]
            for r in dm[['chrStart','chrEnd','color','chr']].query('chr==@contig & chrEnd > @xmin & chrStart < @xmax ').itertuples():
                polygon = Polygon([[r[1],-4],[r[1],-5],[r[2],-5],[r[2],-4]], True, color=r[3],alpha=1)
                patches2.append(polygon)
            p = PatchCollection(patches2, alpha=1, match_original=True)
            # print("patchlen: ",len(patches2))
            ax.add_collection(p)
        
        ax.add_collection(PatchCollection(valid_patches,color='lightgray',zorder=2,aa=True,edgecolor='w',linewidth=0.01))
        ax.add_collection(PatchCollection(patches,color='paleturquoise',zorder=1,aa=True,edgecolor='w',linewidth=0.01))                     
        ax.add_collection(PatchCollection(ticks,color='k',zorder=7,linewidth=0.5,edgecolor='k',aa=True)) # aa=True
        ax.add_collection(PatchCollection(bad_ticks,color='r',zorder=8,linewidth=0.5,edgecolor='r',aa=True)) # aa=True

        ax.scatter(x=never_seen3,y=[-2.5]*len(never_seen3),color='k',s=55,marker="|",facecolors=None,linewidths=0.4,zorder=2)


        if sunklines:
            sl = [Rectangle((x, -1),0.5, (max(rows_used) + 1.5) )  for x in kmerlist]
            ax.add_collection(PatchCollection(sl,color='darkgray',zorder=0.9,linewidth=0.5,edgecolor='lightgray',aa=True,linestyle="--")) # aa=True

        # ax.scatter(x=kmerlist,y=[-1.5]*len(kmerlist),color='k',s=55,marker="|",facecolors=None,linewidths=0.4,zorder=2)
        # # plt.scatter(data=graph3,x='ID',y='counter',color='red',s=2)
        fmt = '{x:,.0f}'
        tick = mtick.StrMethodFormatter(fmt)
        ax.xaxis.set_major_formatter(tick)


        ax.set_ylabel("ONT Read Depth")
        ax.set_xlabel("Contig coordinate")

        plt.rcParams['svg.fonttype'] = 'none'
        plt.savefig(plotdir + "/" + snakemake.wildcards.sample +"_" + snakemake.wildcards.hap + "_" +cont2+"_"+str(regstart)+"_"+str(regend)+"_detailed.svg",format="svg",pad_inches=0,bbox_inches='tight')
        plt.savefig(plotdir + "/" + snakemake.wildcards.sample +"_" + snakemake.wildcards.hap + "_" +cont2+"_"+str(regstart)+"_"+str(regend)+"_detailed.png",pad_inches=0,bbox_inches='tight',dpi=300)
        plt.savefig(plotdir + "/" + snakemake.wildcards.sample +"_" + snakemake.wildcards.hap + "_" +cont2+"_"+str(regstart)+"_"+str(regend)+"_detailed.pdf",pad_inches=0,bbox_inches='tight')


        #plt.savefig(args.plotdir +row.name +"__" +  args.sampname +"_" + args.hap + "_" +contig+"_"+str(regstart)+"_"+str(regend)+".svg",format="svg",pad_inches=0,bbox_inches='tight')

        #         plt.savefig(plotdir+contig+"_"+str(regstart)+"_"+str(regend)+"_v4.png",dpi=300,pad_inches=0,bbox_inches='tight')
        # plt.savefig(plotdir+contig+"_"+str(regstart)+"_"+str(regend)+"_v4.pdf",dpi=300,pad_inches=0,bbox_inches='tight')
        # plt.savefig(plotdir+contig+"_"+str(regstart)+"_"+str(regend)+"_v4.svg",dpi=300,pad_inches=0,bbox_inches='tight')  
        # plt.savefig(args.plotdir +row.name +"__" +  args.sampname +"_" + args.hap + "_" +contig+"_"+str(regstart)+"_"+str(regend)+".png",dpi=300,pad_inches=0,bbox_inches='tight')

        #         plt.savefig(args.plotdir + args.sampname +"_" + args.hap + "_" + row.name +"_" +contig+"_"+str(regstart)+"_"+str(regend)+".png",dpi=300,pad_inches=0,bbox_inches='tight')
        # plt.savefig("breakplots/"+contig+"_"+str(regstart)+"_"+str(regend)+"_.pdf",dpi=300,pad_inches=0,bbox_inches='tight')
        # plt.savefig("breakplots/"+contig+"_"+str(regstart)+"_"+str(regend)+"_.svg",dpi=300,pad_inches=0,bbox_inches='tight')  
        plt.close()

        print("Plotting complete for " +contig+"_"+str(regstart)+"_"+str(regend)) 
