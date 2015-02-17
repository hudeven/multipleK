#!/usr/bin/env python
import sys
from Bio import SeqIO
from optparse import OptionParser
import random

parser = OptionParser()
parser.add_option('-d','--readsfile', dest = "readsFilename", help = "Name of the reads file")
parser.add_option('-r','--resultsfile', dest = "resultsFilename", help = "Name of the results file")
parser.add_option('-q','--queryfile', dest = "queryFilename", help = "Name of the original K query file")
parser.add_option('-f','--format', dest = "formatstr", help = "format of reads file fasta/fastq")


(options, args) = parser.parse_args(sys.argv[1:])
readsFilename = options.readsFilename
resultsFilename = options.resultsFilename
outputFilename = resultsFilename + ".ana"
queryFilename = options.queryFilename
formatStr = options.formatstr

output_file = open(outputFilename, 'w')
result_file = open(resultsFilename, 'r')
read_records = list(SeqIO.parse(readsFilename, formatStr))
query_records = list(SeqIO.parse(queryFilename, "fasta"))
outputStr = ""
fcount=0;
nfcount=0;

for record in SeqIO.parse(resultsFilename, "fasta"):
    kmerIndex = int(record.name)
    kmerStr = str(query_records[kmerIndex].seq)
    readid_list = str(record.seq).split(',')
    for rid in readid_list: 
	read = read_records[int(rid)]
	readStr = str(read.seq)
	if(readStr.find(kmerStr) == -1):
	    outputStr += "Not query "+ str(kmerIndex)+" "  +kmerStr+" in read "+ rid +"\n"+readStr+"\n";
	    nfcount += 1;
        else:
	    outputStr += "Found query "+ str(kmerIndex)+ " " +kmerStr+" in read " + rid  +"\n"+readStr+"\n";
	    fcount += 1;

total_count = fcount + nfcount;
if total_count == 0:
    print "total count is 0"
else:
    print 1.0*fcount / (fcount+nfcount);
    outputStr += str(1.0*fcount / (fcount+nfcount));
    output_file.write(outputStr);
    output_file.close();
    result_file.close();
