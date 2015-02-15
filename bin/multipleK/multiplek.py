import os
import sys
import timeit
from optparse import OptionParser

#Modify it to your working dir
#WORKING_DIR = "/Users/stevenliu/workspace/new_ndtree/multipleK/"
WORKING_DIR = "/home/stevenliu/workspace/multipleK/"
#WORKING_DIR = "/media/psf/Home/MultipleK/bin/multipleK/"

start = timeit.default_timer()


parser = OptionParser()
parser.add_option('-r','--repeat', dest = "repeat", help = "repeat times")
parser.add_option('-k','--klength', dest = "klength", help = "length of kmer")
parser.add_option('-m','--multiplek', dest = "multiplek", help = "length of dynamic kmer")
(options, args) = parser.parse_args(sys.argv[1:])

klength = options.klength
multiplek = options.multiplek
repeat = int(options.repeat)
os.chdir(WORKING_DIR)

for i in range(0, repeat):
	print "******generate reads from references"
	cmd = "python bin/ref2reads/ref2reads.py --output data/multipleK/" + str(i) + "_reads.fa --readlength 100 --reference data/test.fa --coverage 1 --error 0"
	os.system(cmd)

	print "******generate kmers from reads******\n"
	# fasta_fastkmer.py is for kmer list of fasta format
	cmd = "python bin/reads2kmer/fasta_fastkmer.py --output data/multipleK/kmer.fa --klength "+ klength  +" --readsfile data/multipleK/" + str(i) + "_reads.fa"
	os.system(cmd)
 
	print "******generate random box query******\n"
	cmd = 'python bin/random_boxquery/random_fast_boxquery.py --num 10  --klength '+multiplek +' --output data/multipleK/boxquery --reference data/test.fa --boxsize 1'
	os.system(cmd)
	
	print "******run bond-tree ******\n"
	os.chdir(WORKING_DIR + "src/")

	# modify kmer length in code file dim.h and then compile bondtree
	dim_file = open("dim.h", 'w')
	dim_file.write("const int DIM = " + klength + ";")
	dim_file.close()
	os.system("make");
	cmd = './ndTree --dimension '+ klength  +' --data '+ '../data/multipleK/kmer.fa'+ ' --boxquery ../data/multipleK/boxquery --mode rebuild --record ../data/record --querydim '+multiplek;
	os.system(cmd)

	print "***** intersect result ******\n"
	os.chdir(WORKING_DIR)
	cmd = "python bin/intersect_result/fast_intersect_result.py --input data/multipleK/boxquery.result --output data/multipleK/result.fa"

	print "***** analyse accuracy ******\n"
	os.chdir(WORKING_DIR)
	cmd = 'python bin/align_kmer_read/align_kmer_read.py --output data/cancer/query/analyse.txt --readsfile data/cancer/read/0_reads.fa --resultsfile data/cancer/query/boxquery.result'
	os.system(cmd)

stop = timeit.default_timer()
print "***** running time *****\n"
print stop - start 

