import os
import sys
import time
import timeit
from optparse import OptionParser

WORKING_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)),"../../")
#read_file = WORKING_DIR + 'data/multipleK/homo_read.fastq'
#read_file = WORKING_DIR + 'data/multipleK/random_read.fastq'
#read_file = WORKING_DIR + 'data/multipleK/reads_ecoli.fasta'
#read_file = WORKING_DIR + 'data/multipleK/read_CP65671.fasta'
#read_file = WORKING_DIR + 'data/multipleK/read_NC32798.fasta'
#read_file = WORKING_DIR + 'data/multipleK/read_NZ_CP007569.fasta'
#read_file = WORKING_DIR + 'data/multipleK/read_smallNC_000001.11.fasta'
read_file = WORKING_DIR + 'data/multipleK/read_NC_005106.4.fasta'

#pattern0 is pattern 4 !!!
query_file= WORKING_DIR + 'data/multipleK/patterns/pattern0/pattern0_'
query_file= WORKING_DIR + 'data/multipleK/patterns/pattern1/e1/pattern1_e1_'
#query_file = WORKING_DIR + 'data/multipleK/patterns/pattern1/e2/pattern1_e2_'
#query_file= WORKING_DIR + 'data/multipleK/patterns/pattern2/pattern2_'
#query_file= WORKING_DIR + 'data/multipleK/patterns/pattern3/pattern3_'
#query_file= WORKING_DIR + 'data/multipleK/patterns/pattern5/pattern5_'

result_file = WORKING_DIR + 'data/multipleK/result'

parser = OptionParser()
parser.add_option('-k','--klength', dest = "klength", help = "length of kmer")
parser.add_option('-s','--startk', dest = "startk", help = "minimum k")
parser.add_option('-e','--endk', dest = "endk", help = "maximum k")
parser.add_option('-q','--querynum', dest = "querynum", help = "number of query")
parser.add_option('-n','--subnum', dest = "subnum", help = "number of subqueries")
(options, args) = parser.parse_args(sys.argv[1:])

klength = options.klength
startk = int(options.startk)
endk = int(options.endk)
querynum = options.querynum
subnum = options.subnum

start = timeit.default_timer()

os.chdir(WORKING_DIR)

print "\n****** compile bond-tree ******\n"
os.chdir(WORKING_DIR + 'src/')
# modify kmer length in code file dim.h and then compile bondtree
dim_file = open("dim.h", 'w')
dim_file.write("const int DIM = " + klength + ";")
dim_file.close()
os.system("make");


print "\n****** multiple K query ******\n"
for i in range(startk, endk+1):
    print "\n****** do box query k=" + klength + " K=" + str(i)
    os.chdir(WORKING_DIR + 'src/')
    cmd = 'python queryK.py -k '+ klength + ' -K '+ str(i) + ' -b ' + query_file + str(i) + ' -o '+ result_file + str(i) + ' -s ' + subnum
    os.system(cmd)

    print "\n****** alignment ******\n"
    cmd = 'python ../bin/align_kmer_read/fast_align_kmer_read.py --readsfile ' + read_file + ' --resultsfile ' + result_file + str(i) + ' --queryfile ' + query_file + str(i) 
    os.system(cmd)
    
    print "\n\n\n***** Done ******\n\n\n"

stop = timeit.default_timer()
print "***** running time *****\n"
print stop - start 
