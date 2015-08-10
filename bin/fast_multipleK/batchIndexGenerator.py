#!/usr/bin/env python
import os
import sys
import time
import timeit
from optparse import OptionParser

startk = 8
endk = 8


WORKING_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)),"../../")

read_file = WORKING_DIR + 'data/sample/read.fasta'

#read_file = WORKING_DIR + 'data/multipleK/read_test.fasta'
#read_file = WORKING_DIR + 'data/multipleK/reads_ecoli.fasta'
#read_file = WORKING_DIR + 'data/multipleK/read_CP65671.fasta'
#read_file = WORKING_DIR + 'data/multipleK/read_CP65671_200.fasta'
#read_file = WORKING_DIR + 'data/multipleK/read_NC32798.fasta'
#read_file = WORKING_DIR + 'data/multipleK/random_read.fastq'
#read_file = WORKING_DIR + 'data/multipleK/read_NZ_CP007569.fasta'
#read_file = WORKING_DIR + 'data/multipleK/read_NC_007131.6.fasta'
#read_file = WORKING_DIR + 'data/multipleK/read_smallNC_000001.11.fasta'
#read_file = WORKING_DIR + 'data/multipleK/read_NC_000019.10.fasta'
#read_file = WORKING_DIR + 'data/multipleK/SRR606249_1.fastq'
#read_file = WORKING_DIR + 'data/multipleK/read_NC_005106.4.fasta'

kmer_file = WORKING_DIR + 'data/multipleK/kmer.fa'

prefix = "index_" + read_file[read_file.rindex('/')+1:-6] + "_"
print prefix
fast_format = 'fasta'
print "\n****** multiple K query ******\n"
for i in range(startk, endk+1):
        ikmer_file = kmer_file + '_' + str(i)
	print "\n****** generate kmers from reads ******\n"
        os.chdir(WORKING_DIR)
	# fasta_fastkmer.py is for kmer list of fasta format
	cmd = "python bin/reads2kmer/reads2kmer.py --output "+ ikmer_file +" --klength "+ str(i)  +" --readsfile " + read_file
	os.system(cmd)

	print "\n****** compile bond-tree ******\n"
	os.chdir(WORKING_DIR + 'src/')
	# modify kmer length in code file dim.h and then compile bondtree
	dim_file = open("dim.h", 'w')
	dim_file.write("const int DIM = " + str(i) + ";")
	dim_file.close()
	os.system("make");

	print "\n****** build index tree ******\n"
	cmd = './ndTree --data '+ ikmer_file + ' --mode rebuild ';
	os.system(cmd)

        indexPath = os.path.join(os.path.dirname(read_file),prefix + str(i))
        print "path is : " + indexPath
        os.makedirs(indexPath)
        os.system("mv ../data/index.bin ../data/record.typeid " + indexPath)
