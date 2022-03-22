#Usage : python printEdgeBedIntervals.py <genomesizefile> <L>
#prints two bed intervals per chromosome, the first L bases and the end L bases

import sys
import csv

endLengths=int(sys.argv[2])

with open(sys.argv[1], 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        contiglen = int(row[1])
        if contiglen < endLengths:
            print (row[0], "\t0\t", contiglen, sep="")
        else:
            print (row[0], "\t0\t", endLengths, sep="")
            print (row[0], "\t", contiglen - endLengths, "\t", contiglen, sep="")

