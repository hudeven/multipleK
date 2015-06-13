import sys
import select_sub

seq = "ATGCATGCATGCATGCATGCATGC"
qset = [0,1,2]
k = 4
num = 3
# rset = select_sub.sequential(qset, num)
# rset = select_sub.randomly(qset, num)
# rset = select_sub.binary(qset, num)
# rset = select_sub.skip(qset, num)
rset = select_sub.BTLT(qset, num)

print rset
if rset == None:
    sys.exit(0)

for start in rset:
    kmer = seq[start:start+k]
    print kmer
