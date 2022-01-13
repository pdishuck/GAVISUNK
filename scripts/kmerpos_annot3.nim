import readfq #andreas-wilm, heng li
import kmer #brentp
import os
import sets
import system
import tables
import parsecsv
import strutils
import streams


proc main()=
  let args = commandLineParams()
  var
    reads = args[0]
    kmers = args[1]
    locs = args[2]
    ofile = args[3]
  # read in unique kmers
  var uniqueKmers = initHashSet[uint64]()
  let f = open(kmers)
  var line : string
  while f.read_line(line):
    uniqueKmers.incl(line.encode())
    # uniqueKmers.incl(line.reverse_complement)
  f.close()
  
  var First = true
  var k: int
  let g = open(kmers)
  while g.read_line(line):
    if First:
      k=line.len()
      First = false
    # echo k
    break
    # uniqueKmers.incl(line.reverse_complement)
  g.close()

  stderr.writeLine("FINISHED READING KMERS")    

  # loop through loc file, put in hash table with kmers as keys
  
  type
    coordloc = tuple
      chrom: string
      group: int64
  type
    coord = tuple
      chrom: string
      start: int64
      group: int64
  var
    x: CsvParser
    curcoord: coord
    coords = initTable[uint64,coord]()
    s = newFileStream(locs, fmRead)

  if s == nil:
    quit("cannot open the file" & args[0])

  open(x, s, locs, separator = '\t')
  
  while readrow(x):
    curcoord.chrom = x.row[0]
    curcoord.start = parseInt(x.row[1])
    curcoord.group = parseInt(x.row[3])
    coords[x.row[2].encode()] = curcoord
  close(x)
  

  stderr.writeLine("FINISHED READING LOCS")    
    
  
  # for x in uniqueKmers.items:
  #   echo x
  #   break
  
  # let k = 20
  # var r = 0
  let o = open(ofile, fmWrite)
  var
    curLoc: coordLoc
    prevLoc: coordLoc
  for rec in readfq(reads):
    var i = 0
    # var found = 0
    for x in rec.sequence.slide(k):
      if x[0] in uniqueKmers:
        curcoord = coords[x[0]] 
        curLoc = (curcoord.chrom,curcoord.group)
        if curLoc == prevloc: continue
        o.writeLine(rec.name, "\t", i, "\t", curcoord.chrom, "\t", curcoord.start, "\t", curcoord.group)
        prevLoc = curLoc
        # inc found
      inc i
  o.close()
    # stderr.writeLine("read#:",r, "\t", "kmercount:",found)    
    # inc r

main()



# #zstd output
#   import zstd/compress
#   import zstd/decompress

#   var source = readFile("tests/files/nixon.bmp")
#   var compressed = compress(source, level=3)
#   var decompressed = decompress(compressed)
#   check equalmem(decompressed[0].addr, source[0].addr, source.len)
