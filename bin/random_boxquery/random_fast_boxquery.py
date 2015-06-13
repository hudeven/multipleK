#!/usr/bin/env python
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from optparse import OptionParser
import random
import patterns

def letter2num(c):
     if (not cmp(c,'A')):
         elem = '0'
     elif (not cmp(c,'T')):
         elem = '1'
     elif (not cmp(c,'G')):
         elem = '2'
     elif (not cmp(c,'C')):
         elem = '3'
     else:
         elem = '4'
 
     return elem

parser = OptionParser()
parser.add_option('-k', '--klength', dest = 'klength', help = "the length of kmer")
parser.add_option('-n', '--num', dest = 'num', help = "the amount of queries")
parser.add_option('-o', '--output', dest = 'output', help = "the filename of box query file")
parser.add_option('-r', '--read', dest = 'read', help = " the filename of read")
parser.add_option('-b', '--boxsize', dest = 'boxsize', help = "size of box query")
parser.add_option('-p', '--pattern', dest = 'pattern', help = "pattern id")

(options, args) = parser.parse_args(sys.argv[1:])
klength = int(options.klength)
num = int(options.num)
read_filename = options.read
if read_filename[-1] == 'q':
    formatStr = "fastq"
else:
    formatStr = "fasta"
query_filename = options.output
boxsize = int(options.boxsize)
patternId = options.pattern

alphabet = ['A', 'T', 'G', 'C']
query_file = open(query_filename, 'w')
patternDistributionFile = open("pattern_dist" + str(klength), 'w')
record_list = list(SeqIO.parse(read_filename,formatStr))
maxIndex = len(record_list) - 1
query_list = []
pHash = {}
i = 0
while i < num:
	# pick up a read randomly
	readIndex = random.randrange(0, maxIndex)
	record = record_list[readIndex]
	seqStr = str(record.seq[:])
	end =len(seqStr) - klength

	query_letter = '';
	query = '';
	r = random.randrange(0, end)
	kmer = seqStr[r : r+klength]   
	#print kmer
	if patternId != None:
	    pid = patterns.getPattern(kmer)
	    if pid in pHash:
		pHash[pid] += 1
	    else:
		pHash[pid] = 1
	    #print "############################ pid = " + str(pid)
	    if pid != int(patternId) and int(patternId) != 0:
                #i += 1  # result query number may not equal to num !!!
		continue
	box = random.randrange(r, r+klength)
	for p in range(r, r+klength):
	    if boxsize >= 2 and  p == box:
		query += "("
		for t in range(0, boxsize):
		    query += alphabet[t]
		query += ")"
	    else:
	        query += seqStr[p] 
	# id is very important for multiple K but name and desc is optional 
	# and only for user to check query info
	query_seq = SeqIO.SeqRecord(Seq(query,generic_dna), id = str(i), description="dim="+str(klength))
	query_list.append(query_seq)
	i += 1

SeqIO.write(query_list, query_file, "fasta");
query_file.close()

for key in pHash:
	patternDistributionFile.write(str(key) + " " + str(pHash[key]) + "\n")

patternDistributionFile.close()
#print aln_ref+'\n'+aln_sample

