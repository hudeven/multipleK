#!/usr/bin/env python
import sys
from Bio import SeqIO
from optparse import OptionParser
import random

parser = OptionParser()
parser.add_option('-d','--readsfile', dest = "readsFilename", help = "Name of the reads file")
parser.add_option('-r','--resultsfile', dest = "resultsFilename", help = "Name of the results file")
parser.add_option('-q','--queryfile', dest = "queryFilename", help = "Name of the original K query file")

(options, args) = parser.parse_args(sys.argv[1:])
readsFilename = options.readsFilename
if readsFilename[-1] == 'q':
    formatStr = "fastq"
else:
    formatStr = "fasta"
resultsFilename = options.resultsFilename
outputFilename = resultsFilename + ".ana"
queryFilename = options.queryFilename

output_file = open(outputFilename, 'w')
result_file = open(resultsFilename, 'r')
read_records = list(SeqIO.parse(readsFilename, formatStr))
query_records = list(SeqIO.parse(queryFilename, "fasta"))
outputStr = ""
sum_accurate = 0.0
BUFFER_SIZE = 128
count = 0 
for record in SeqIO.parse(resultsFilename, "fasta"):
    count += 1
    if count % BUFFER_SIZE == 0:
	output_file.write(outputStr);
	outputStr = ""
    kmerIndex = int(record.name)
    kmerStr = str(query_records[kmerIndex].seq)
    readid_list = str(record.seq).split(',')
    fcount = 0
    nfcount = 0
    for rid in readid_list: 
	read = read_records[int(rid)]
	readStr = str(read.seq)
	if(readStr.find(kmerStr) == -1):
	    outputStr += "Not query "+ str(kmerIndex)+" "  +kmerStr+" in read "+ rid +"\n"+readStr+"\n";
	    nfcount += 1;
        else:
	    outputStr += "Found query "+ str(kmerIndex)+ " " +kmerStr+" in read " + rid  +"\n"+readStr+"\n";
	    fcount += 1;
    sum_accurate += fcount * 1.0 / (fcount + nfcount)

if count == 0:
    print "Empty result file " + resultsFilename
else:
    rate = sum_accurate / count;
    print "###### Accuracy = "+str(rate)
    outputStr += "###### Accuracy = " + str(rate);
    output_file.write(outputStr);
output_file.close();
result_file.close();
