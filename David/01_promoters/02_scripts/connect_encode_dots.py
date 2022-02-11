#!/usr/bin/env python
import sys,numpy
#N=2
N=None
bedgraph = sys.argv[1]

lastChrom = None
lastStart = None
lastValue = None

for input_line_num, line in enumerate(open(bedgraph)):
    fields = line.strip().split()
    chrom = fields[0]
    start = int(fields[1])
    value = float(fields[3])

    if lastStart is None:
        lastChrom = chrom
        lastStart = start
        lastValue = value
        continue

    if chrom != lastChrom:
        # reset
        lastChrom = chrom
        lastStart = start
        lastValue = value
        continue

    if value == lastValue: # found a string of zeroes
        interp = [value] * 100
    else:
        try:
            interp = numpy.arange(lastValue, value, -(lastValue-value)/100)
        except:
            print(input_line_num, lastChrom, chrom, lastValue,value,file=sys.stderr)
            raise

    for i, val in enumerate(interp):
        print(lastChrom, lastStart + i,lastStart + i + 1, "%.6f" % val)

    lastChrom = chrom
    lastStart = start
    lastValue = value

    if N is not None:
        N = N - 1
        if N <= 0: 
            break
    

        
