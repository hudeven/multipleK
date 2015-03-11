#!/usr/bin/env python
import sys
from optparse import OptionParser
from Bio import SeqIO

parser = OptionParser()
parser.add_option('-i', '--input', dest = 'inputFile', help = 'input file name')
parser.add_option('-o', '--output', dest ='outputFile', help = 'output file name')
parser.add_option('-a', '--action', dest ='action', help = 'intersect / union')

(options, args) =parser.parse_args(sys.argv[1:])
inputFilename = options.inputFile
outputFilename = options.outputFile
action = options.action

result_file = open(outputFilename, "w")

old_id = -1
set_num =  0;
hist = {}
for record in SeqIO.parse(inputFilename, "fasta"):
    if record.id != old_id and old_id != -1:
	query = ">" + str(old_id) + "\n";
	if action == 'union':
	    its = ','.join(key for key in hist)
	else:
	    # only save 1 read id
	    #its = ""
	    #for key in hist:
	    #	if hist[key] == set_num:
	    #	    its = key
	    #	    break
	    its = ','.join(key for key in hist if hist[key]==set_num)

	if its != "":
	    query += its
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

    old_id = record.id

# write for the last record
query = ">"+str(old_id) + "\n";
its = ','.join(key for key in hist if hist[key]==set_num)
if its != "":
    query += its
    result_file.write(query+"\n");
result_file.close()
	    
