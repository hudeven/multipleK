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
parser.add_option('-r', '--reference', dest = 'reference', help = " the filename of reference sequence")
parser.add_option('-b', '--boxsize', dest = 'boxsize', help = "size of box query")

(options, args) = parser.parse_args(sys.argv[1:])
klength = int(options.klength)
num = int(options.num)
ref_filename = options.reference
query_filename = options.output
boxsize = options.boxsize
query_letter_filename = options.output+".letter"

query_letter_file = open(query_letter_filename, 'w')
query_file = open(query_filename, 'w')
record = SeqIO.read(open(ref_filename),"fasta")
ref = str(record.seq[:])

end =len(ref) - klength
query_list = []
for i in range(0, num):
	query_letter = '';
	query = '';
	r = random.randrange(0, end)
	kmer = ref[r : r+klength]   
	for p in range(r, r+klength):
	    query_letter += '1' + ' ' + ref[p]+'\n' 
	    query += ref[p] 
		
        query_letter_file.write(str(i)+"\t"+kmer+'\n'+query_letter)
	query_seq = SeqIO.SeqRecord(Seq(query,generic_dna), id = record.id, name=record.name, description = record.description)
	query_list.append(query_seq)

SeqIO.write(query_list, query_file, "fasta");
query_file.close()
query_letter_file.close();


#print aln_ref+'\n'+aln_sample

