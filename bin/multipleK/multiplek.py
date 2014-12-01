import os
import sys
from optparse import OptionParser

#Modify it to your working dir
WORKING_DIR = "/home/stevenliu/workspace/bondtree/"

parser = OptionParser()
parser.add_option('-n','--count', dest = "count", help = "number of random cancer dna")
parser.add_option('-r','--repeat', dest = "repeat", help = "repeat times")
parser.add_option('-k','--klength', dest = "klength", help = "length of kmer")
parser.add_option('-m','--multiplek', dest = "multiplek", help = "length of dynamic kmer")
(options, args) = parser.parse_args(sys.argv[1:])

count = int(options.count)
klength = options.klength
multiplek = options.multiplek
repeat = int(options.repeat)
os.chdir(WORKING_DIR)

for i in range(0, repeat):
	print "******generate cancers' sequences from reference******\n"
	for i in range(0, count):
	    cmd = "python bin/ref2cancer_ref/ref2cancer_ref.py --output data/cancer/ref/"+str(i)+"_ref.fa --reference data/test.fa --id "+ str(i) + " --desc "+ "cancer"+str(i)+" --rate 100"
	    os.system(cmd)

	print "******generate reads from references"
	#os.chdir(WORKING_DIR + "bin/ref2reads/")
	for i in range(0, count):
	    cmd = "python bin/ref2reads/ref2reads.py --output data/cancer/read/" + str(i) + "_reads.fa --readlength 100 --reference data/cancer/ref/" + str(i) + "_ref.fa --coverage 20 --error 0"
	    os.system(cmd)

	print "******generate kmers from reads******\n"
	for i in range(0, count):
	    cmd = "python bin/reads2kmer/fasta_kmer.py --output data/cancer/kmer/" + str(i)+".kmer --klength "+ klength  +" --readsfile data/cancer/read/" + str(i) + "_reads.fa"
	    os.system(cmd)

	print "****** combine all cancers' kmers to one file ******\n"
	os.chdir(WORKING_DIR + "data/cancer/kmer/")
	cmd1 = "cat "
	cmd2 = "cat "
	cmd3 = "cat "
	cmd4 = "cat "
	for i in range(0, count):
	    cmd1 += str(i)+".kmer "
	    cmd2 += str(i)+".kmer.id "
	    cmd3 += str(i)+".kmer.desc "
	    cmd4 += str(i)+".kmer.name "
	cmd1 += "> all.kmer"
	cmd2 += "> all.kmer.id"
	cmd3 += "> all.kmer.desc"
	cmd4 += "> all.kmer.name"
	os.system(cmd1)
	os.system(cmd2)
	os.system(cmd3)
	os.system(cmd4)
	 
	print "******generate random box query******\n"
	os.chdir(WORKING_DIR)

	cmd = 'python bin/random_boxquery/random_boxquery.py --klength '+multiplek +' --output data/cancer/query/boxquery --reference data/test.fa'
	os.system(cmd)
	
	print "******run bond-tree ******\n"
	os.chdir(WORKING_DIR + "src/")
	cmd = './ndTree  --index ../data/index_real --dimension '+ klength  +' --data '+ '../data/cancer/kmer/all.kmer'+ ' --boxquery ../data/cancer/query/boxquery --skip 0 --newtree 1 --aux ../data/cancer/kmer/all.kmer.desc --record ../data/record --querydim '+multiplek;
	os.system(cmd)

	print "***** analyse accuracy ******\n"
	os.chdir(WORKING_DIR)
	cmd = 'python bin/align_kmer_read/align_kmer_read.py --output data/cancer/query/analyse.txt --readsfile data/cancer/read/0_reads.fa --resultsfile data/cancer/query/boxquery.result'
	os.system(cmd)
