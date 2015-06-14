#!/usr/bin/env python

klength = 12
startk = 7
endk = 11
querynum = 1000
theta = -1
patternId = '0'

import os
import sys
import time
import timeit
from optparse import OptionParser

WORKING_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)),"../../")
#read_file = WORKING_DIR + 'data/multipleK/homo_read.fastq'
#read_file = WORKING_DIR + 'data/multipleK/random_read.fastq'
read_file = WORKING_DIR + 'data/multipleK/read_CP65671.fasta'
#read_file = WORKING_DIR + 'data/multipleK/read_CP65671_200.fasta'
#read_file = WORKING_DIR + 'data/multipleK/read_NC32798.fasta '
#read_file = WORKING_DIR + 'data/multipleK/read_NZ_CP007569.fasta'
#read_file = WORKING_DIR + 'data/multipleK/reads_ecoli.fasta'

query_file = WORKING_DIR + 'data/multipleK/boxquery'
result_file = WORKING_DIR + 'data/multipleK/result'

print "\n****** compile bond-tree ******\n"
os.chdir(WORKING_DIR + 'src/')
# modify kmer length in code file dim.h and then compile bondtree
dim_file = open("dim.h", 'w')
dim_file.write("const int DIM = %s ;" % klength)
dim_file.close()
os.system("make");

for i in range(startk, endk+1):
    print "\n****** multiple K query k=%s K=%s ******\n" % (klength, i)
    print "\n******generate random box query******\n"
    os.chdir(WORKING_DIR)
    cmd = 'python bin/random_boxquery/random_fast_boxquery.py --num '+ str(querynum)  +'  --klength '+ str(i) +' --output '+ query_file + str(i) +' --read '+ read_file +' --boxsize 1 --pattern ' + patternId 
    os.system(cmd)

    print "\n****** do box query k=" + str(klength) + " K=" + str(i)
    os.chdir(WORKING_DIR + 'src/')
    #cmd = 'python queryK.py -k '+ klength + ' -K '+ str(i) + ' -b ' + query_file + str(i) + ' -o '+ result_file + str(i) + ' --threshold ' + threshold
    cmd = 'python optQueryK.py -k '+ str(klength) + ' -K '+ str(i) + ' -b ' + query_file + str(i) + ' -o '+ result_file + str(i) + " -t " + str(theta)
    os.system(cmd)

    print "\n****** alignment ******\n"
    cmd = 'python ../bin/align_kmer_read/fast_align_kmer_read.py --readsfile ' + read_file + ' --resultsfile ' + result_file + str(i) + ' --queryfile ' + query_file + str(i) 
    os.system(cmd)
    
    print "\n\n\n***** Done ******\n\n\n"
