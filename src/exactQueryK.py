#!/usr/bin/env python
import os
import sys
from Bio import SeqIO
from optparse import OptionParser
import time
import timeit
import select_sub
import subprocess
import itertools

def getProductions(distance, seq):
    alpha = ['A','T','G','C']
    prod = []
    for i in range(distance):
        prod.append(alpha)
    xSeqs = list(itertools.product(*prod))
    fullSeqs = []
    for i in range(distance+1):
        for xseq in xSeqs:
            fullSeqs.append(''.join(xseq[:i])+ seq + ''.join(xseq[i:]))
    return fullSeqs

def exactSmallQuery(rset, theta = -1):
    resultSet = set()
    lastSize = 0
    first = True
    qCount = 0
    
    for subquery in rset:
        qCount += 1
#        print subquery
        p = subprocess.Popen(["./ndTree", "-o", subquery], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, err = p.communicate() 
        p_status = p.wait()
        if ('No' in output):
   	    resultSet = []
	    break
        else:
            lines = output.split()
	    readset = lines[-1].split(',')
            if first:
                first = False
                resultSet = set(readset)
		lastSize = len(resultSet)
                #print "qCount=" + str(qCount) + "\t size=" + str(len(resultSet))
            else:
                resultSet = resultSet.union(set(readset))
		if len(resultSet) -lastSize <= theta:
                    break
            lastSize = len(resultSet)
    return (qCount, resultSet)


def optQuery(rset, fullMer, theta=0):
    resultSet = set()
    lastSize = 0
    first = True
    qCount = 0
    
    for i in rset:
        qCount += 1
        subquery = fullMer[i : i+k]
        #print "subquery= " + subquery
        #print '\n****** run ndtree ******\n'
        p = subprocess.Popen(["./ndTree", "-o", subquery], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, err = p.communicate() 
        p_status = p.wait()
        if ('No' in output):
   	    resultSet = []
	    break
        else:
            lines = output.split()
	    readset = lines[-1].split(',')
            if first:
                first = False
                resultSet = set(readset)
		lastSize = len(resultSet)
                #print "qCount=" + str(qCount) + "\t size=" + str(len(resultSet))
            else:
                if K - k >= 0:
                    resultSet = resultSet.intersection(set(readset))
		    if lastSize - len(resultSet) <= theta:
		        break
                else:
                    resultSet = resultSet.union(set(readset))
		    if len(resultSet) -lastSize <= theta:
                        break
            lastSize = len(resultSet)
    #print len(resultSet)
    return (qCount, resultSet)





parser = OptionParser()
parser.add_option('-o','--output', dest = "outputFile", help = "Name of the output file")
parser.add_option('-k','--klength', dest = "kLength", help = "Length of indexed kmer")
parser.add_option('-K','--Klength', dest = "KLength", help = "Length of query kmer")
parser.add_option('-b','--boxquery', dest = "boxquery", help = "Box query file name")
parser.add_option('-t','--theta', dest = "theta", help = "lastSize - len(resultSet) <= theta")
parser.add_option('-s','--strategy', dest = "strategy", help = "strategy to choose sub queries")
parser.add_option('-m','--mode', dest = "mode", help = "exact: do all equivalentexact query; others: do box query")


(options, args) = parser.parse_args(sys.argv[1:])
resultFilename = options.outputFile
boxqueryFilename = options.boxquery
k = int(options.kLength)
K = int(options.KLength)
theta = 0
if options.theta != None:
    theta = int(options.theta)
strategy = "bfs"
if options.strategy != None:
    strategy = options.strategy.lower()
print "Strategy is " + strategy

mode = 'box'
if options.mode != None:
    mode = options.mode


UNI_ELEM = 'X'

subqueryFilename = boxqueryFilename + ".sub"

start = timeit.default_timer();
opt = 0
query_num = 0

print '\n****** Generating subqueries ******\n'
print 'Multiple K query k = ' + str(k) + ' K = ' + str(K) + ' mode ' + mode 
totalQ = 0
usedQ = 0
validTotalQ = 0
validUsedQ = 0
result_file = open(resultFilename, 'w')
qset = []
distance = abs(K-k)
for i in range(0, distance + 1):
    qset.append(i)

for record in SeqIO.parse(boxqueryFilename, "fasta"):
    if strategy == "bfs":
	rset = select_sub.fullBFS(qset)
    elif strategy == "random":
        rset = select_sub.randomly(qset)
    elif strategy == "sequential":
        rset = select_sub.sequentialy(qset)
    else:
	rset = select_sub.fullBFS(qset)
      
	#print "****** rset = " + str(rset)
    if K - k >= 0:
        fullMer = str(record.seq)
        subq , resultSet = optQuery(rset, fullMer, theta)
    else:
	if mode.lower() == "exact":
            rset = getProductions(distance, str(record.seq))
            subq , resultSet = exactSmallQuery(rset, theta)
        else:
            fullMer = UNI_ELEM * distance + str(record.seq) + UNI_ELEM * distance 
            subq , resultSet = optQuery(rset, fullMer, theta)
            
    usedQ += subq
    totalQ += len(rset)

    if len(resultSet) != 0:
        #print str(subq) + " / " + str(len(rset))
        validUsedQ += subq
        validTotalQ += len(rset)
        outputStr = ">" + str(record.id) + "\n" + ",".join(resultSet) + "\n"
        result_file.write(outputStr)
    

print "valid optimization(reduces sub queries):"
if validTotalQ != 0:
    print validUsedQ * 1.0 / validTotalQ
print "overall optimization"
if totalQ != 0:
    print usedQ * 1.0 / totalQ
result_file.close()

stop = timeit.default_timer();
print '###### query time = '+ str(stop - start)

   


