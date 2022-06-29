import readfq #andreas-wilm, heng li
import system
import os

proc main()=
  let args = commandLineParams()
  var
    reads = args[0]
    ofile = args[1] 
  # let k = 20
  # var r = 0
  let o = open(ofile, fmWrite)
  for rec in readfq(reads):
    o.writeLine(rec.name, "\t", len(rec.sequence))
  o.close()

main()



# #zstd output
#   import zstd/compress
#   import zstd/decompress

#   var source = readFile("tests/files/nixon.bmp")
#   var compressed = compress(source, level=3)
#   var decompressed = decompress(compressed)
#   check equalmem(decompressed[0].addr, source[0].addr, source.len)
