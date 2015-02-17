#!/usr/bin/env python
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from optparse import OptionParser
import random

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
parser.add_option('-f', '--format', dest = 'formatstr', help = "format of input file fasta/fastq")

(options, args) = parser.parse_args(sys.argv[1:])
klength = int(options.klength)
num = int(options.num)
read_filename = options.read
query_filename = options.output
formatStr = options.formatstr
boxsize = int(options.boxsize)
alphabet = ['A', 'T', 'G', 'C']
query_file = open(query_filename, 'w')
record_list = list(SeqIO.parse(read_filename,formatStr))
maxIndex = len(record_list) - 1
query_list = []
for i in range(0, num):
	# pick up a read randomly
	readIndex = random.randrange(0, maxIndex)
	record = record_list[readIndex]
	seqStr = str(record.seq[:])
	end =len(seqStr) - klength

	query_letter = '';
	query = '';
	r = random.randrange(0, end)
	kmer = seqStr[r : r+klength]   
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

SeqIO.write(query_list, query_file, "fasta");
query_file.close()


#print aln_ref+'\n'+aln_sample

