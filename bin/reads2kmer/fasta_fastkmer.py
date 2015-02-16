#!/usr/bin/env python
import sys
from Bio import SeqIO
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-o','--output', dest = "outputFile", help = "Name of the output file")
parser.add_option('-k','--klength', dest = "kmerLength", help = "Length of kmer")
parser.add_option('-f','--readsfile', dest = "readsFilename", help = "Name of the reads file")
(options, args) = parser.parse_args(sys.argv[1:])
kmerFilename = options.outputFile
readsFilename = options.readsFilename
kmerLength = int(options.kmerLength)

kmer_file = open(kmerFilename, 'w')

kmer_list=""
buffer_size = 255
count = 0
for seq_record in SeqIO.parse(readsFilename,"fasta"):
    cur = 0
    cur_max = len(seq_record) - kmerLength
    for cur in range(0, cur_max):
	kmer_seq = str(seq_record.seq[cur:cur+kmerLength]);
	kmer = '>' + seq_record.id + '\n' + kmer_seq + '\n'
	# kmer = SeqIO.SeqRecord(kmer_seq, id=seq_record.id, description="")
	# kmer = SeqIO.SeqRecord(kmer_seq, id=seq_record.id, name=seq_record.name, description=seq_record.description)
	kmer_list += kmer
	count += 1;
        if count > buffer_size:
	    kmer_file.write(kmer_list);
	    count = 0
	    kmer_list = "";
	    # SeqIO.write(kmer_list, kmer_file, "fasta");
if count != 0:
    kmer_file.write(kmer_list)
kmer_file.close()

