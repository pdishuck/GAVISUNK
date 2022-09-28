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
import zip/gzipfiles

proc main()=
  let args = commandLineParams()
  var
    sunkpos = args[0]
    hic_phase = args[1]
    bandwidth = parseInt(args[2])
    cntFile = args[3]
    contigTable = initTable[string, CountTable[char]()]()
    phase: char
  # type
  #   kmerRec = tuple
  #     rname: string
  #     pos: int64
  #     chrom: string
  #     start: int64
  #     group: int64

  var
    x: CsvParser
    y: CsvParser
    # kmerRecs = initTable[uint64,kmerRec]()
    sp = newFileStream(sunkpos, fmRead)
  if sp == nil:
    quit("cannot open the file" & args[0])
  var
    h = newGzFileStream(hic_phase, fmRead)
    phases = initTable[string,char]()

  open(y, h, hic_phase, separator = '\t')
  while readrow(y):
    phases[y.row[0]] = y.row[1][0]
  close(y)

  open(x, sp, sunkpos, separator = '\t')
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


      if maxgood>1:
        phase = phases[rnamePrev]
        echo rnamePrev, "\t", bestcontig, "\t", maxgood, "\t", bestdir   , '\t', maxgood, '\t', phase
        if not contigTable.haskey(bestcontig): contigTable[bestcontig] = initCountTable[char]()
        contigTable[bestcontig].inc(phase)

      posStarts.clear()

    if not posStarts.haskey(x.row[2]):
      posStarts[x.row[2]] = (@[],@[],@[]) # make new entry of empty seqs for contig
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


  if maxgood>1:
      phase = phases[rnamePrev]
      echo rnamePrev, "\t", bestcontig, "\t", maxgood, "\t", bestdir   , '\t', maxgood, '\t', phase
      if not contigTable.haskey(bestcontig): contigTable[bestcontig] = initCountTable[char]()
      contigTable[bestcontig].inc(phase)

  let f = open(cntFile, fmWrite)
  var maxphase: char
  var maxcnt: int
  for key,cnt in contigTable.pairs:
    maxphase = '0'
    maxcnt = 0
    for k,v in cnt:
      
#       if v > maxcnt and k != '0':
#         maxphase = k
#         maxcnt = v
        
      f.writeLine(key,"\t", k,"\t", v)
#     f.writeLine(key,"\t", maxphase,"\t", maxcnt)
  f.close()
  # how many pos total?

    # maybe output multiple contigs

main()
