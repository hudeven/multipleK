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
parser.add_option('-b', '--boxsize', dest = 'boxsize', help = "size of box query")
parser.add_option('-f', '--format', dest = 'formatstr', help = "format of input file fasta/fastq")

(options, args) = parser.parse_args(sys.argv[1:])
klength = int(options.klength)
num = int(options.num)
query_filename = options.output
formatStr = options.formatstr
boxsize = int(options.boxsize)
alphabet = ['A', 'T', 'G', 'C']
query_file = open(query_filename, 'w')
query_list = []
for i in range(0, num):
	query_letter = '';
	query = '';
	box = random.randrange(0, klength)
	for p in range(0, klength):
	    if boxsize >= 2 and  p == box:
		query += "("
		for t in range(0, boxsize):
		    query += alphabet[t]
		query += ")"
	    else:
	        query += alphabet[random.randrange(0,4)]

	query_seq = SeqIO.SeqRecord(Seq(query,generic_dna), id = str(i), description="dim="+str(klength))
	query_list.append(query_seq)

SeqIO.write(query_list, query_file, "fasta");
query_file.close()


#print aln_ref+'\n'+aln_sample

