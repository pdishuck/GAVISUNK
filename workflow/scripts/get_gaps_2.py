import pandas as pd, numpy as np
import os
import argparse
import pyranges as pr

parser = argparse.ArgumentParser()
parser.add_argument("fai")

parser.add_argument("sample")
parser.add_argument("hap")
parser.add_argument('indir')
parser.add_argument('outdir')
args = parser.parse_args()
outdir = args.outdir
indir = args.indir
samp = args.sample
hap = args.hap

def main():
    print("hap1")
    fai = pd.read_csv(args.fai,sep='\t',header=None,names=['contig','length','cumlen','3','4'])
    contigs = fai.contig.tolist()
    allbreaks={}
    outbreaks = []
    bp1 = 0
    bp2 = 0
    nodata = []
    for x in contigs:
        # print(x)
        cont2 = x.replace("#","_")        
        contiglen = fai.query("contig == @x").length.values[0]
        try: bedin = pd.read_csv(indir + cont2 + "_"+hap+".merged.valid.bed",sep="\t",header=None,names=['contig','start','end'])
    #     try: bedin = pd.read_csv(x + ".bed",sep="\t",header=None,names=['contig','start','end'])        
        
        except:
            print("No bed for", x)
            nodata.append((x, 0, contiglen ))
            bp2 += contiglen
            continue
        if len(bedin) == 0:
            nodata.append((x, 0, contiglen ))
            bp2 += contiglen
            continue 
            
        bed = pr.PyRanges(bedin.rename(columns={'contig':'Chromosome','start':'Start','end':'End'}))
        bedin2 = bed.merge().df
        bedin2.columns = ['contig', 'start', 'end']
        bedin2['end_prev'] = bedin2['end'].shift(1,fill_value=np.nan)
        
    
        breaks = []
        first=True
        for row in bedin2.itertuples():
            if row.end < row.end_prev: 
                print("weird", x)
                break
            if first:
                first=False
                # print("confirmed region begins: ",row.start)
            else:
                # print(int(row.end_prev), row.start,)
                breaks.append((int(row.end_prev), row.start,))
                outbreaks.append((x, int(row.end_prev), row.start-1,))
                bp1 += (row.start - int(row.end_prev))

    #     print(row.end)
    #     print("unconfirmed bp at contig end: ", contiglen - row.end)
    # #     print(breaks)
        allbreaks[x] = breaks
    print("breaks: ", len(outbreaks), bp1)
    print("no data: ", len(nodata), bp2)
    pd.DataFrame(nodata).to_csv(outdir +hap+".nodata.bed",header=False,sep="\t",index=False)
    pd.DataFrame(outbreaks).to_csv(outdir +hap+".merged.gaps.bed",header=False,sep="\t",index=False)

'''
    fai = pd.read_csv(args.fai,sep='\t',header=None,names=['contig','length','cumlen','3','4'])
    contigs = fai.contig.tolist()
    print("hap2")
    allbreaks={}
    outbreaks = []
    bp1 = 0
    bp2 = 0
    nodata = []
    for x in contigs:
        cont2 = x.replace("#","_")        # print(x)
        contiglen = fai.query("contig == @x").length.values[0]
    #     try: bedin = pd.read_csv(x + "_hap2.bed",sep="\t",header=None,names=['contig','start','end'])
        try: bedin = pd.read_csv(indir + cont2 + "_hap2.bed",sep="\t",header=None,names=['contig','start','end'])
        except:
            print("No bed for", x)
            nodata.append((x, 0, contiglen ))
            bp2 += contiglen
            continue
        if len(bedin) == 0:
            nodata.append((x, 0, contiglen ))
            bp2 += contiglen
            continue 
            
        bed = pr.PyRanges(bedin.rename(columns={'contig':'Chromosome','start':'Start','end':'End'}))
        bedin2 = bed.merge().df
        bedin2.columns = ['contig', 'start', 'end']
        bedin2['end_prev'] = bedin2['end'].shift(1,fill_value=np.nan)
        
    
        breaks = []
        first=True
        for row in bedin2.itertuples():
            if first:
                first=False
                # print("confirmed region begins: ",row.start)
            else:
                # print(int(row.end_prev), row.start,)
                breaks.append((int(row.end_prev), row.start,))
                outbreaks.append((x, int(row.end_prev), row.start-1,))
                bp1 += (row.start - int(row.end_prev))

    #     print(row.end)
    #     print("unconfirmed bp at contig end: ", contiglen - row.end)
    # #     print(breaks)
        allbreaks[x] = breaks
    print("breaks: ", len(outbreaks), bp1)
    print("no data: ", len(nodata), bp2)
    pd.DataFrame(nodata).to_csv(outdir + "hap2.nodata.bed",header=False,sep="\t",index=False)
    pd.DataFrame(outbreaks).to_csv(outdir + "hap2.gaps.bed",header=False,sep="\t",index=False)
    '''
main()
