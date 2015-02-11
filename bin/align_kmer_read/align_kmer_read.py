#!/usr/bin/env python
import sys
from Bio import SeqIO
from optparse import OptionParser
import random

def letter2num(c):
     if (not cmp(c,'A')):
         elem = 0
     elif (not cmp(c,'T')):
         elem = 1
     elif (not cmp(c,'G')):
         elem = 2
     elif (not cmp(c,'C')):
         elem = 3
     else:
         elem = 4
 
     return elem

def num2letter(n):
     if (not cmp(n,'0')):
         elem = 'A'
     elif (not cmp(n,'1')):
         elem = 'T'
     elif (not cmp(n,'2')):
         elem = 'G'
     elif (not cmp(n,'3')):
         elem = 'C'
     else:
         elem = 'N'
 
     return elem



parser = OptionParser()
parser.add_option('-o','--output', dest = "outputFile", help = "Name of the output file")
parser.add_option('-f','--readsfile', dest = "readsFilename", help = "Name of the reads file")
parser.add_option('-r','--resultsfile', dest = "resultsFilename", help = "Name of the results file")


(options, args) = parser.parse_args(sys.argv[1:])
outputFilename = options.outputFile
readsFilename = options.readsFilename
resultsFilename = options.resultsFilename
print resultsFilename

output_file = open(outputFilename, 'w')
result_file = open(resultsFilename, 'r')
lines = result_file.readlines()
seq_records = list(SeqIO.parse(readsFilename,"fasta"))
outputStr = ""
i = 0
fcount=0;
nfcount=0;
   
end = len(lines)
while i<end:
    readids = lines[i].split()
    i+=1;
    kmer = lines[i].split()
    i+=1;
    kmerStr = ""
    for mer in kmer:
        kmerStr = kmerStr + num2letter(mer)
    for readid in readids:
	readidInt = int(readid)
        read = seq_records[readidInt-1000]
        readStr = ''.join(read)
        if(readStr.find(kmerStr) == -1):
	    outputStr += "Not found "+kmerStr+" in \n"+readStr+"\n";
	    nfcount += 1;
        else:
	    outputStr += "Found "+kmerStr+" in \n"+readStr+"\n";
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
