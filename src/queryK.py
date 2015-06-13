#!/usr/bin/env python
import os
import sys
from Bio import SeqIO
from optparse import OptionParser
import time
import timeit
import select_sub

parser = OptionParser()
parser.add_option('-o','--output', dest = "outputFile", help = "Name of the output file")
parser.add_option('-k','--klength', dest = "kLength", help = "Length of indexed kmer")
parser.add_option('-K','--Klength', dest = "KLength", help = "Length of query kmer")
parser.add_option('-b','--boxquery', dest = "boxquery", help = "Box query file name")
parser.add_option('-s','--subquery', dest = "subquery", help = "The amount of subqueries for each query")
parser.add_option('-t','--threshold', dest = "threshold", help = "read id set with size greater than threshold will be kick out!")


(options, args) = parser.parse_args(sys.argv[1:])
outputFilename = options.outputFile
boxqueryFilename = options.boxquery
k = int(options.kLength)
K = int(options.KLength)
threshold = '0'
if options.threshold != None:
    threshold = options.threshold

sub_num = 0
if options.subquery != None:
    sub_num = int(options.subquery)
print sub_num

UNI_ELEM = 'X'

subqueryFilename = boxqueryFilename + ".sub"
resultFilename = subqueryFilename + ".result"

start = timeit.default_timer();
opt = 0
query_num = 0

print '\n****** Generating subqueries ******\n'
if K > k:
    print 'Big K query k = ' + str(k) + ' K = ' + str(K) + " sub = " + str(sub_num)
    subquery_file = open(subqueryFilename, 'w')
    qset = []
    for i in range(0, K - k + 1):
	qset.append(i)

    for record in SeqIO.parse(boxqueryFilename, "fasta"):
	rset = select_sub.BTLT(qset, sub_num)
#	rset = select_sub.randomly(qset, sub_num)
#	rset = select_sub.sequentialy(qset, sub_num)
#	print "****** rset = " + str(rset)

        query_num += 1
        opt += len(rset)
	for i in rset:
	    query = ">" + record.id + " dim="+ str(k) + " start=" + str(i)+"\n"
	    query += str(record.seq[i : i+k]) + "\n"
	    # print query
	    subquery_file.write(query)
    subquery_file.close();

    print '\n****** run ndtree ******\n'
    cmd = './ndTree --boxquery ' + subqueryFilename
    os.system(cmd)
    
    print '\n****** intersect result *******\n'
    cmd = 'python ../bin/intersect_result/fast_intersect_result.py --input '+resultFilename+' --output '+outputFilename + ' --action intersect ' + ' --threshold ' + threshold
    os.system(cmd)

    print "optimization(reduces sub queries):"
    print opt/(query_num * len(qset))

elif K < k:
    print 'Small K query k = ' + str(k) + ' K = ' + str(K) + " sub = " + str(sub_num)
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
    cmd = 'python ../bin/intersect_result/fast_intersect_result.py --input '+resultFilename+' --output '+outputFilename + ' --action union ' + ' --threshold ' + threshold 
    os.system(cmd)


else:
    print 'Normal K query k = ' + str(k) + ' K = ' + str(K) + " sub = " + str(sub_num)
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
    cmd = 'python ../bin/intersect_result/fast_intersect_result.py --input '+resultFilename+' --output '+outputFilename + ' --action union ' + ' --threshold ' + threshold
    os.system(cmd)

stop = timeit.default_timer();
print '###### query time = '+ str(stop - start)

   

