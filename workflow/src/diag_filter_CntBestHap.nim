## Philip Dishuck
##


# import readfq #andreas-wilm, heng li
# import kmer #brentp
import os
import system
import tables
import parsecsv
import streams
import strUtils
import zip/gzipfiles

proc main()=
  let args = commandLineParams()
  var
    cntFile = args[0]
    cntFile_out = args[1]
    contigTable = initTable[string, CountTable[char]()]()
  # type
  #   kmerRec = tuple
  #     rname: string
  #     pos: int64
  #     chrom: string
  #     start: int64
  #     group: int64

  var
    x: CsvParser
    # kmerRecs = initTable[uint64,kmerRec]()
    cf = newFileStream(cntFile, fmRead)
  if cf == nil:
    quit("cannot open the file" & args[0])


  open(x, cf, cntFile, separator = '\t')
  while readrow(x):
    if not contigTable.haskey(x.row[0]): contigTable[x.row[0]] = initCountTable[char]()
    contigTable[x.row[0]].inc(x.row[1][0],parseInt(x.row[2]))
  close(x)
  

  var maxphase: char
  var maxcnt: int
  let f = open(cntFile_out, fmWrite)  
  for key,cnt in contigTable.pairs:
    maxphase = '0'
    maxcnt = 0
    for k,v in cnt:
      
      if v > maxcnt and k != '0':
        maxphase = k
        maxcnt = v
        
    #   f.writeLine(key,"\t", k,"\t", v)
#     f.writeLine(key,"\t", maxphase,"\t", maxcnt)
    f.writeLine(key,"\t", maxphase)
  f.close()
  # how many pos total?

    # maybe output multiple contigs

main()
