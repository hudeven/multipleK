#!/usr/bin/env python
import os
import sys
from Bio import SeqIO
from optparse import OptionParser
import time

parser = OptionParser()
parser.add_option('-o','--output', dest = "outputFile", help = "Name of the output file")
parser.add_option('-k','--klength', dest = "kLength", help = "Length of indexed kmer")
parser.add_option('-K','--Klength', dest = "KLength", help = "Length of query kmer")
parser.add_option('-b','--boxquery', dest = "boxquery", help = "Box query file name")

(options, args) = parser.parse_args(sys.argv[1:])
outputFilename = options.outputFile
boxqueryFilename = options.boxquery
k = int(options.kLength)
K = int(options.KLength)

analyzeFilename = outputFilename+".an"
subqueryFilename = boxqueryFilename + ".sub"
resultFilename = subqueryFilename + ".result"

if K > k:
    print 'Big K query'
    print '****** generating subqueries ******\n'
    subquery_file = open(subqueryFilename, 'w')
    for record in SeqIO.parse(boxqueryFilename, "fasta"):
	for i in range(0, K-k+1):
	    query = ">" + record.id + " dim="+ str(k) + "\n"
	    query += str(record.seq[i : i+k]) + "\n"
	    # print query
	    subquery_file.write(query)
    subquery_file.close();

    print '****** run ndtree ******\n'
    cmd = './ndTree --boxquery ' + subqueryFilename
    os.system(cmd)
    
    print '****** intersection result *******\n'
    cmd = 'python ../bin/intersect_result/fast_intersect_result.py --input '+resultFilename+' --output '+outputFilename
    os.system(cmd)
#
#    print '****** align and analyse\n'
    
#    cmd = 'python ../bin/align_kmer_read/fast_align_kmer_read.py --output '+analyzeFilename+' --readsfile ../data/cancer/read/0_reads.fa --resultsfile '+resultFilename+' --queryfile '+boxqueryFilename
 #   os.system(cmd)

elif K < k:
    print 'Small K query'

else:
    print 'Normal K query'


