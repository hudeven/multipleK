#!/usr/bin/env python
import os
import sys
from Bio import SeqIO
from optparse import OptionParser
import time
import timeit

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
UNI_ELEM = 'X'

analyzeFilename = outputFilename+".an"
subqueryFilename = boxqueryFilename + ".sub"
resultFilename = subqueryFilename + ".result"

start = timeit.default_timer();
print '\n****** Generating subqueries ******\n'
if K > k:
    print 'Big K query k = ' + str(k) + ' K = ' + str(K)
    subquery_file = open(subqueryFilename, 'w')
    for record in SeqIO.parse(boxqueryFilename, "fasta"):
	for i in range(0, K-k+1):
	    query = ">" + record.id + " dim="+ str(k) + "\n"
	    query += str(record.seq[i : i+k]) + "\n"
	    # print query
	    subquery_file.write(query)
    subquery_file.close();

    print '\n****** run ndtree ******\n'
    cmd = './ndTree --boxquery ' + subqueryFilename
    os.system(cmd)
    
    print '\n****** intersect result *******\n'
    cmd = 'python ../bin/intersect_result/fast_intersect_result.py --input '+resultFilename+' --output '+outputFilename + ' --action intersect'
    os.system(cmd)

elif K < k:
    print 'Small K query k = ' + str(k) + ' K = ' + str(K)
    subquery_file = open(subqueryFilename, 'w')
    subquery_file = open(subqueryFilename, 'w')
    max_num = k-K
    for record in SeqIO.parse(boxqueryFilename, "fasta"):
	for left in range(0, max_num+1):
	    query = ">" + record.id + " dim="+ str(k) + "\n"
	    for i in range(0, left):
		query += UNI_ELEM
	    query += str(record.seq)
	    right = max_num - left
	    for i in range(0, right):
		query += UNI_ELEM
	    query += '\n'
	    subquery_file.write(query)
    subquery_file.close();

    print '\n****** run ndtree ******\n'
    cmd = './ndTree --boxquery ' + subqueryFilename
    os.system(cmd)
 
    print '\n****** union result *******\n'
    cmd = 'python ../bin/intersect_result/fast_intersect_result.py --input '+resultFilename+' --output '+outputFilename + ' --action union'
    os.system(cmd)


else:
    print 'Normal K query k = ' + str(k) + ' K = ' + str(K)
    subquery_file = open(subqueryFilename, 'w')
    max_num = k-K
    for record in SeqIO.parse(boxqueryFilename, "fasta"):
	query = ">" + record.id + " dim="+ str(k) + "\n"
	query += str(record.seq) + '\n'
	subquery_file.write(query)
    subquery_file.close();

    print '\n****** run ndtree ******\n'
    cmd = './ndTree --boxquery ' + subqueryFilename
    os.system(cmd)
 
    print '\n****** union result *******\n'
    cmd = 'python ../bin/intersect_result/fast_intersect_result.py --input '+resultFilename+' --output '+outputFilename + ' --action union'
    os.system(cmd)

stop = timeit.default_timer();
print '###### query time = '+ str(stop - start)

   

