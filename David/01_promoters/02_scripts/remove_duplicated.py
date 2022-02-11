#!/usr/bin/env pythonm
import sys

lastChrom = None
lastStart = None
for line in open(sys.argv[1]):
    fields = line.strip().split()
    chrom = fields[0]
    start = int(fields[1])

    if lastStart is not None and start == lastStart: # this is a duplicated start point, a bug in my connect the dots script
        lastChrom = chrom
        lastStart = start
        continue

    print("\t".join(fields))
    lastChrom = chrom
    lastStart = start
