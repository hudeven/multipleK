#!/usr/bin/env python
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from optparse import OptionParser
import random
 
parser = OptionParser()
parser.add_option('-o','--output', dest = "outputFile", help = "Name of the output file")
parser.add_option('-l','--readlength', dest = "readlength", help = "lengh of read")
parser.add_option('-r','--reference', dest = "reference", help = "Name of the reference file")
parser.add_option('-c','--coverage', dest = "coverage", help = "coverage for sequencer")
parser.add_option('-e','--error', dest = "error", help = "error rate for sequencer e/1000")
parser.add_option('-t','--type', dest = "error", help = "error rate for sequencer e/1000")


(options, args) = parser.parse_args(sys.argv[1:])
outFilename = options.outputFile
refFilename = options.reference
readlength = int(options.readlength)
coverage = int(options.coverage)
error = int(options.error)

#There should be one and only one record, the entire genome:
record = SeqIO.read(open(refFilename), "fasta")
output_handle = open(outFilename, "w")
seqNoN = str(record.seq).replace('N','')
frags=[]
limit=len(seqNoN)
readnum = coverage * limit / readlength
BUFFER_SIZE = 10000
count = 0
for i in range(0, readnum) :
    count += 1
    if count > BUFFER_SIZE:
        count = 0
        SeqIO.write(frags, output_handle, "fasta")
        frags = []

    start=random.randint(0,limit-readlength)
    end=start+readlength
    frag=seqNoN[start:end]
    readitem=SeqIO.SeqRecord(Seq(frag),id= str(i+1000),name=record.name, description=record.description)
    frags.append(readitem)
 
SeqIO.write(frags, output_handle, "fasta")
output_handle.close()
