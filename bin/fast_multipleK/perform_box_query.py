#!/usr/bin/env python
import os
import sys
import time
import timeit
from optparse import OptionParser

#readFilename = "read_CP65671.fasta"
readFilename = "read_NC32798.fasta"
#readFilename = "homo_read.fastq"
#readFilename = "random_read.fastq"
#readFilename = "reads_ecoli.fasta"
klength = 7
startk = 7
endk = 30
querynum = 1000

patternId = '0'
WORKING_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)),"../../")
read_file = WORKING_DIR + 'data/multipleK/' + readFilename
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
    os.chdir(WORKING_DIR)
    cmd = 'python bin/random_boxquery/random_fast_boxquery.py --num '+ str(querynum)  +'  --klength '+ str(i) +' --output '+ query_file + str(i) +' --read '+ read_file +' --boxsize 1 --pattern ' + patternId 
    os.system(cmd)

    os.chdir(WORKING_DIR + 'src/')
    
    # BFS 
    cmd = 'python optQueryK.py -k '+ str(klength) + ' -K '+ str(i) + ' -b ' + query_file + str(i) + ' -o '+ result_file + str(i) 
    os.system(cmd)
    
    cmd = 'python ../bin/align_kmer_read/fast_align_kmer_read.py --readsfile ' + read_file + ' --resultsfile ' + result_file + str(i) + ' --queryfile ' + query_file + str(i) 
    os.system(cmd)

    # Random 
    cmd = 'python optQueryK.py -k '+ str(klength) + ' -K '+ str(i) + ' -b ' + query_file + str(i) + ' -o '+ result_file + str(i) + " -s random"
    os.system(cmd)
    
    cmd = 'python ../bin/align_kmer_read/fast_align_kmer_read.py --readsfile ' + read_file + ' --resultsfile ' + result_file + str(i) + ' --queryfile ' + query_file + str(i) 
    os.system(cmd)

    # Sequential
    cmd = 'python optQueryK.py -k '+ str(klength) + ' -K '+ str(i) + ' -b ' + query_file + str(i) + ' -o '+ result_file + str(i) + " -s sequential"
    os.system(cmd)
    
    cmd = 'python ../bin/align_kmer_read/fast_align_kmer_read.py --readsfile ' + read_file + ' --resultsfile ' + result_file + str(i) + ' --queryfile ' + query_file + str(i)
    os.system(cmd)

    # Full 
    cmd = 'python optQueryK.py -k '+ str(klength) + ' -K '+ str(i) + ' -b ' + query_file + str(i) + ' -o '+ result_file + str(i) + ' -t -1'
    os.system(cmd)

    cmd = 'python ../bin/align_kmer_read/fast_align_kmer_read.py --readsfile ' + read_file + ' --resultsfile ' + result_file + str(i) + ' --queryfile ' + query_file + str(i) 
    os.system(cmd)
