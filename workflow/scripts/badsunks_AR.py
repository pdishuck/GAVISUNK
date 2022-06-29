#! /usr/bin/python

import argparse
import pandas as pd
from matplotlib import (pyplot as plt,
                        lines)
import seaborn as sns
import os

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("fai1", help="list of contigs, hap1")
    parser.add_argument("fai2", help="list of contigs, hap2")
    parser.add_argument("sunkpos1", help=".sunkpos file of SUNK locations on hap1 ONT reads")
    parser.add_argument("sunkpos2", help=".sunkpos file of SUNK locations on hap2 ONT reads")
    parser.add_argument("outpath", help="output file")
    args = parser.parse_args()

    sunkposcat = pd.read_csv(args.sunkpos1,sep="\t",header=None,names=['rname','pos','chrom','start','ID'],dtype={'rname':'string','pos':'uint32','chrom':'category','start':'uint32','ID':'uint32'})# 


    sunkposcat['ID2'] = sunkposcat['chrom'].astype(str) +":"+ sunkposcat['ID'].astype(str)
    kmer_counts = sunkposcat.ID2.value_counts()
    kmer_counts2 = pd.DataFrame(kmer_counts,index=None)
    kmer_counts2.reset_index(inplace=True)
    kmer_counts2.rename(columns={'ID2':'count','index':'ID2'},inplace=True)

    kmer_counts2[['contig','ID']] =  kmer_counts2.ID2.str.split(":",1,expand=True)

    fai = pd.read_csv(args.fai1, sep='\t',header=None,names=['contig','length','cumlen','3','4'])
    contiglist = set(fai.contig.tolist())

    kmer_counts2['correct_hap'] = kmer_counts2['contig'].isin(contiglist)

    kmer_counts2.correct_hap.value_counts()

    fig,ax = plt.subplots(figsize=(22,22))
    g = sns.histplot(data = kmer_counts2.query('count < 150 & count > 1'), x='count',binwidth=1,hue='correct_hap')
#    plt.savefig("badsunks_hap1.png")

    # meancov = kmer_counts2[kmer_counts2['count'] > kmer_counts2['count'].median()/2]['count'].mode()
    meancov = kmer_counts2.query("correct_hap")['count'].mode().values[0]
    print("mean coverage", meancov)

    limit = meancov + (meancov)**0.5*4 #4 SDs above mean (though apparent mean is depressed by seq accuracy)

    badsunks1a = set(kmer_counts2.query("correct_hap").query("(count > @limit) or (count < 2)")['ID2'])
    print(len(badsunks1a))

    badsunks1b = set(kmer_counts2.query("not correct_hap").query("(count > @limit)")['ID2'])
    print(len(badsunks1b))


    # sunkposfile = args.sunkpos
    sunkposfile = args.sunkpos2

    sunkposcat = pd.read_csv(sunkposfile,sep="\t",header=None,names=['rname','pos','chrom','start','ID'],dtype={'rname':'string','pos':'uint32','chrom':'category','start':'uint32','ID':'uint32'})# 

    sunkposcat['ID2'] = sunkposcat['chrom'].astype(str) +":"+ sunkposcat['ID'].astype(str)
    kmer_counts = sunkposcat.ID2.value_counts()
    kmer_counts2 = pd.DataFrame(kmer_counts,index=None)
    kmer_counts2.reset_index(inplace=True)
    kmer_counts2.rename(columns={'ID2':'count','index':'ID2'},inplace=True)

    kmer_counts2[['contig','ID']] =  kmer_counts2.ID2.str.split(":",1,expand=True)



    fai = pd.read_csv(args.fai2 ,sep='\t',header=None,names=['contig','length','cumlen','3','4'])
    contiglist = set(fai.contig.tolist())

    kmer_counts2['correct_hap'] = kmer_counts2['contig'].isin(contiglist)

    kmer_counts2.correct_hap.value_counts()

    fig,ax = plt.subplots(figsize=(22,22))
    g = sns.histplot(data = kmer_counts2.query('count < 150 & count > 1'), x='count',binwidth=1,hue='correct_hap')
 #   plt.savefig("badsunks_hap2.png")

    # meancov = kmer_counts2[kmer_counts2['count'] > kmer_counts2['count'].median()/2]['count'].mode()
    meancov = kmer_counts2.query("correct_hap")['count'].mode().values[0]
    print("mean coverage", meancov)


    limit = meancov + (meancov)**0.5*4 #4 SDs above mean (though apparent mean is depressed by seq accuracy)


    badsunks2a = set(kmer_counts2.query("correct_hap").query("(count > @limit) or (count < 2)")['ID2'])
    print(len(badsunks2a))


    badsunks2b = set(kmer_counts2.query("not correct_hap").query("(count > @limit)")['ID2'])
    print(len(badsunks2b))


    badsunks_all = badsunks1a |badsunks1b | badsunks2a | badsunks2b
    print(len(badsunks_all))

    outfile = open(args.outpath, "w")
    for element in badsunks_all:
        outfile.write(element + "\n")
    outfile.close()

main()
