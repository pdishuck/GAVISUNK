import os
import system
import parsecsv
# import streams
# import sets
# from streams import newFileStream
from strutils import join
# import gzfile
import zip/gzipfiles
import tables


proc main() =
  let args = commandLineParams()
  
  var contigfile = args[1]
  var bestContigs = initTable[string,string]()
  var x: CsvParser
  x.open(newGZFileStream(contigfile), contigfile, separator = '\t')
  while x.readRow():
    bestContigs[x.row[0]] = x.row[1]
  x.close()

  # var contigs : HashSet[string]
  # let f = open(contigfile)
  # var line : string
  # while f.readline(line):
  #   contigs.incl(line)
  # f.close()
  
  # let ctgCol = 2
  type
    outLine = array[500_000,seq[string]]
  
  var outlines: outLine
  # var x: CsvParser
  x.open(newGZFileStream(args[0]), args[0], separator = '\t')

  # let workContig = "h1tg000116l"
  var bestContig = ""
  var i = 0
  # var match = false
  # columns: rname, seq, pos, contig, start
  var rname_prev = "string"
  while x.readRow():
    if x.row[0] != rname_prev: 
      # if match:
      for j in 0..<i:
        stdout.write(outLines[j].join("\t"),'\n')
      stdout.flushFile()
      bestContig = bestContigs.getOrDefault(x.row[0])
  #      output outLines from 0 to i-1
      i = 0
    # outLines[i] = x.row
    # inc i
    rname_prev = x.row[0]
    if x.row[2] == bestContig:
      outLines[i] = x.row
      inc i
  
  # if match:
  for j in 0..<i:
    stdout.write(outLines[j].join("\t"),'\n')
  stdout.flushFile()
  x.close()
main()
#  curcoord = coords[x.row[1]]
#  stdout.write(x.row.join("\t"),'\t',curcoord.chrom,'\t',curcoord.start,"\n")
#     stdout.write(coords[x.row[1]])



# import ./version

# let args = commandLineParams()
# 
#  
# var
# #  filtBam: Bam
#   rawBam: Bam
#   outBam: Bam
#   fasta: cstring
#   threads = 0
#   format = "b"
# 
# #open(filtBam, $args[0], threads=threads, fai=fasta)
# open(rawBam, $args[1], threads=threads, fai=fasta)
# open(outBam, $args[2], threads=threads, mode="wb", fai=fasta)
# 
# outBam.write_header(rawBam.hdr)
# # get set of all query names in filter bam
# 
# let filtTxt = open($args[0])
# for rec in filtTxt.lines:
#   names.incl(rec)
# filtTxt.close()
# 
# # for rec in filtBam:
# #   names.incl(rec.qname)
# #filtBam.close()
# 
# stderr.writeLine(names.len," qnames in filter")
# 
# #output raw queries to new file
# var writecnt = 0
# for rec in rawBam:
#   if names.contains(rec.qname):
#     outBam.write(rec)
#     inc(writecnt)
# stderr.writeLine(writecnt," records written")
# 
# rawBam.close()
# OutBam.close()
# 


# var First = true
# var k: int
# let g = open(kmers)
# while g.read_line(line):
#   if First:
#     k=line.len()
#     First = false
#   # echo k
#   break
#   # uniqueKmers.incl(line.reverse_complement)
# g.close()
# 
# 
# # for x in uniqueKmers.items:
# #   echo x
# #   break
# 
# # let k = 20
# var r = 0
# for rec in readfq(reads):
#   var i = 0
#   var found = 0
#   for x in rec.sequence.slide(k):
#     if x[0] in uniqueKmers: 
#       var n = newString(k)
#       x[0].decode(n)    
#       echo rec.name, "\t", n, "\t", i
#       inc found
#     inc i
#   stderr.writeLine("read#:",r, "\t", "kmercount:",found)    
#   inc r
