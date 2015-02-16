import os
import sys
import time
import timeit
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-k','--klength', dest = "klength", help = "length of kmer")
parser.add_option('-s','--startk', dest = "startk", help = "minimum k")
parser.add_option('-e','--endk', dest = "endk", help = "maximum k")
parser.add_option('-q','--querynum', dest = "querynum", help = "number of query")
parser.add_option('-w','--workingdir', dest = "workingdir", help = "working directory")

(options, args) = parser.parse_args(sys.argv[1:])

klength = options.klength
startk = int(options.startk)
endk = int(options.endk)
querynum = options.querynum
WORKING_DIR = options.workingdir

ref_file = WORKING_DIR + 'data/test.fa'
read_file = WORKING_DIR + 'data/multipleK/read.fa'
query_file = WORKING_DIR + 'data/multipleK/boxquery'
result_file = WORKING_DIR + 'data/multipleK/result'

start = timeit.default_timer()

print "\n****** batch box query ******\n"
for i in range(startk, endk+1):
    print "index k=" + klength + "\tquery K=" + str(i)
    print "\n******generate random box query******\n"
    os.chdir(WORKING_DIR)
    cmd = 'python bin/random_boxquery/random_fast_boxquery.py --num '+ querynum  +'  --klength '+ str(i) +' --output '+ query_file + str(i) +' --reference '+ ref_file +' --boxsize 1'
    os.system(cmd)
    print "output query file: \n" + query_file + str(i)

    print "\n****** do box query k=" + klength + " K=" + str(i)
    os.chdir(WORKING_DIR + 'src/')
    cmd = 'python queryK.py -k '+ klength + ' -K '+ str(i) + ' -b ' + query_file + str(i) + ' -o '+ result_file + str(i)
    os.system(cmd)

    print "\n****** alignment ******\n"
    cmd = 'python ../bin/align_kmer_read/fast_align_kmer_read.py --readsfile ' + read_file + ' --resultsfile ' + result_file + str(i) + ' --queryfile ' + query_file + str(i)
    os.system(cmd)
        
    print "\n\n\n***** Done ******\n\n\n"

stop = timeit.default_timer()
print "***** running time *****\n"
print stop - start 
