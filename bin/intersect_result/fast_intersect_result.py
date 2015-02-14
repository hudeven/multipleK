#!/usr/bin/env python
import sys
from optparse import OptionParser
from Bio import SeqIO

parser = OptionParser()
parser.add_option('-i', '--input', dest = 'inputFile', help = 'input file name')
parser.add_option('-o', '--output', dest ='outputFile', help = 'output file name')

(options, args) =parser.parse_args(sys.argv[1:])
inputFilename = options.inputFile
outputFilename = options.outputFile

result_file = open(outputFilename, "w")

old = -1
set_num =  0;
hist = {}
for record in SeqIO.parse(inputFilename, "fasta"):
    if record.id != old and old != -1:
	query = ">" + old + "\n";
	query += ','.join(key for key in hist if hist[key]==set_num)
	print query
	result_file.write(query+"\n");
	hist.clear()
	set_num = 0;
	
    set_num += 1

    seq_str = str(record.seq)
    li = seq_str.split(',')
    for key in li:
        if key in hist:
	    hist[key] += 1
    	else:
	    hist[key] = 1

    old = record.id

# write for the last record
query = ">" + old + "\n";
query += ','.join(key for key in hist if hist[key]==set_num)
print query
result_file.write(query+"\n");
result_file.close()
	    
