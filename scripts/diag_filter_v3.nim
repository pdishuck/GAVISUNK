## Philip Dishuck
# This version restricts considered "best" contigs to those in the target haplotype with new argument, fai
##


# import readfq #andreas-wilm, heng li
# import kmer #brentp
import os
import sets
import system
import tables
import parsecsv
import strutils
import streams

import arraymancer

proc main()=
  let args = commandLineParams()
  var
    sunkpos = args[0]
    fai = args[1] 
  # type
  #   kmerRec = tuple
  #     rname: string
  #     pos: int64
  #     chrom: string
  #     start: int64
  #     group: int64
  var # import contig list for this haplotype
    f: CsvParser
    s1 = newFileStream(fai, fmRead)
    contigs = initHashSet[string]()
  open(f, s1, fai, separator = '\t')
  while readrow(f):
    contigs.incl(f.row[0])
  f.close()

  var
    x: CsvParser
    # kmerRecs = initTable[uint64,kmerRec]()
    s = newFileStream(sunkpos, fmRead)
    bandwidth = 2500
  if s == nil:
    quit("cannot open the file" & args[0])

  open(x, s, sunkpos, separator = '\t')
  var rnamePrev = ""
  type
    posStart = tuple
      pos: seq[int64]
      start: seq[int64]
      group: seq[int64]
  var posStarts = initTable[string, posStart]()      
  # var poss = initTable[string, seq[int64]]()
  # var starts = initTable[string, seq[int64]]()
  # var i = 0
  var goodset = initHashSet[int64]()
  # var goodfrac: float
  var hitlen: int
  var maxgood: int
  var maxhitlen: int
  while readrow(x):
    if x.row[0] != rnamePrev:

    
      var median: int64

      maxgood = 0
      var bestcontig = ""
      var bestdir = ""
      maxhitlen = 0
      
      for k,v in posStarts.pairs():
        # echo k

        var p = v[0].toTensor()
        var s = v[1].toTensor()
        var g = v[2].toTensor()
        hitlen = g.shape[0]
        
        if  p.shape[0] == 1: continue
        # try reverse orientation of read to contig
        var n = p + s

        median = int(percentile(n,50))
        var good = 0
        clear(goodset)

        for coord, v in n:
          if (abs(v-median) < bandwidth): 
            # echo s.at(coord)[0], '\t', g.at(coord)[0],'\t', v-median
            goodset.incl(g.at(coord)[0]) 
            good += 1
        # if good > maxgood:
        # goodfrac = len(goodset)/hitlen
        # if len(goodset) > maxgood:
        if len(goodset) > maxgood:  
          maxgood = len(goodset)
          maxhitlen = hitlen
          bestcontig = k      
          bestdir = "-"
        if (len(goodset) == maxgood) and hitlen < maxhitlen:
          maxgood = len(goodset)
          maxhitlen = hitlen
          bestcontig = k      
          bestdir = "-"      

        # forward orientation 
        # echo "forward"
        n = s - p

        median = int(percentile(n,50))
        good = 0
        clear(goodset)
        for coord, v in n:
          # echo coord
          # echo v
          # echo p[coord]
          
          if (abs(v-median) < bandwidth): 
            # echo s.at(coord)[0], '\t', g.at(coord)[0], '\t', v-median
            goodset.incl(g.at(coord)[0])

            good += 1
        # if good > maxgood:
        # goodfrac = len(goodset)/hitlen
        # if len(goodset) > maxgood:
        if len(goodset) > maxgood:  
          maxgood = len(goodset)
          maxhitlen = hitlen
          bestcontig = k      
          bestdir = "+"
        if (len(goodset) == maxgood) and hitlen < maxhitlen:
          maxgood = len(goodset)
          maxhitlen = hitlen
          bestcontig = k      
          bestdir = "+"      


      if maxgood>1: echo rnamePrev, "\t", bestcontig, "\t", maxgood, "\t", bestdir   , '\t', maxgood

      posStarts.clear()

    if not posStarts.haskey(x.row[2]):
      posStarts[x.row[2]] = (@[],@[],@[]) # make new entry of empty seqs for contig
    if contigs.contains(x.row[2]):
      posStarts[x.row[2]].pos.add(parseInt(x.row[1])) # add pos and start to seqs for contigs
      posStarts[x.row[2]].start.add(parseInt(x.row[3]))
      posStarts[x.row[2]].group.add(parseInt(x.row[4]))

    rnamePrev = x.row[0]
  close(x)


  var median: int64
  maxgood = 0
  var bestcontig: string
  var bestdir: string
  maxhitlen = 0
        
  for k,v in posStarts.pairs():
    # echo k

    var p = v[0].toTensor()
    var s = v[1].toTensor()
    var g = v[2].toTensor()
    hitlen = g.shape[0]
    
    if  p.shape[0] == 1: continue
    # try reverse orientation of read to contig
    var n = p + s

    median = int(percentile(n,50))
    var good = 0
    clear(goodset)

    for coord, v in n:
      if (abs(v-median) < bandwidth): 
        # echo s.at(coord)[0], '\t', g.at(coord)[0],'\t', v-median
        goodset.incl(g.at(coord)[0]) 
        good += 1
    # if good > maxgood:
    # goodfrac = len(goodset)/hitlen
    # if len(goodset) > maxgood:
    if len(goodset) > maxgood:  
      maxgood = len(goodset)
      maxhitlen = hitlen
      bestcontig = k      
      bestdir = "-"
    if (len(goodset) == maxgood) and hitlen < maxhitlen:
      maxgood = len(goodset)
      maxhitlen = hitlen
      bestcontig = k      
      bestdir = "-"      

    # forward orientation 
    # echo "forward"
    n = s - p

    median = int(percentile(n,50))
    good = 0
    clear(goodset)
    for coord, v in n:
      # echo coord
      # echo v
      # echo p[coord]
      
      if (abs(v-median) < bandwidth): 
        # echo s.at(coord)[0], '\t', g.at(coord)[0], '\t', v-median
        goodset.incl(g.at(coord)[0])

        good += 1
    # if good > maxgood:
    # goodfrac = len(goodset)/hitlen
    # if len(goodset) > maxgood:
    if len(goodset) > maxgood:  
      maxgood = len(goodset)
      maxhitlen = hitlen
      bestcontig = k      
      bestdir = "+"
    if (len(goodset) == maxgood) and hitlen < maxhitlen:
      maxgood = len(goodset)
      maxhitlen = hitlen
      bestcontig = k      
      bestdir = "+"      


  if maxgood>1: echo rnamePrev, "\t", bestcontig, "\t", maxgood, "\t", bestdir   , '\t', maxgood



  # how many pos total?

    # maybe output multiple contigs

main()
