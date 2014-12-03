import os
import sys
import timeit
from optparse import OptionParser

#Modify it to your working dir
#WORKING_DIR = "/Users/stevenliu/workspace/new_ndtree/multipleK/"
WORKING_DIR = "/home/stevenliu/workspace/multipleK/"

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
	cmd = "python bin/ref2reads/ref2reads.py --output data/cancer/read/" + str(i) + "_reads.fa --readlength 100 --reference data/test.fa --coverage 1 --error 0"
	os.system(cmd)

	print "******generate kmers from reads******\n"
	cmd = "python bin/reads2kmer/fasta_kmer.py --output data/cancer/kmer/all.kmer --klength "+ klength  +" --readsfile data/cancer/read/" + str(i) + "_reads.fa"
	os.system(cmd)
 
	print "******generate random box query******\n"
	cmd = 'python bin/random_boxquery/random_boxquery.py --num 1000  --klength '+multiplek +' --output data/cancer/query/boxquery --reference data/test.fa'
	os.system(cmd)
	
	print "******run bond-tree ******\n"
	os.chdir(WORKING_DIR + "src/")

	# modify kmer length in code file dim.h and then compile bondtree
	dim_file = open("dim.h", 'w')
	dim_file.write("const int DIM = " + klength + ";")
	dim_file.close()
	os.system("make");
	cmd = './ndTree  --index ../data/index --dimension '+ klength  +' --data '+ '../data/cancer/kmer/all.kmer'+ ' --boxquery ../data/cancer/query/boxquery --mode rebuild  --aux ../data/cancer/kmer/all.kmer.desc --record ../data/record --querydim '+multiplek;
	os.system(cmd)

	print "***** analyse accuracy ******\n"
	os.chdir(WORKING_DIR)
	cmd = 'python bin/align_kmer_read/align_kmer_read.py --output data/cancer/query/analyse.txt --readsfile data/cancer/read/0_reads.fa --resultsfile data/cancer/query/boxquery.result'
	os.system(cmd)

stop = timeit.default_timer()
print "***** running time *****\n"
print stop - start 

