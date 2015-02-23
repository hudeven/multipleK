#!/usr/bin/env python
import sys
from Bio import SeqIO
from optparse import OptionParser
import random

parser = OptionParser()
parser.add_option('-o','--output', dest = "outputFile", help = "Name of the output file")
parser.add_option('-r','--readlength', dest = "readLength", help = "Length of read")
parser.add_option('-n','--number', dest = "num", help = "number of read")
(options, args) = parser.parse_args(sys.argv[1:])
readFilename = options.outputFile
readLength = int(options.readLength)
read_id = -1;
num = int(options.num)

read_file = open(readFilename, 'w')
read_list=""
buffer_size = 255
count = 0
alphabet=['A','T','G','C']

read_quality = ""
for i in range(0,readLength):
    read_quality += '0';

for t in range(0,num):
    read_id += 1
    read_seq = ""
    for i in range(0,readLength):
	read_seq += alphabet[random.randrange(0,4)]
    read = '@' + str(read_id) + '\n' + read_seq + '\n+\n' + read_quality+'\n'
    read_list += read
    count += 1;
    if count > buffer_size:
	read_file.write(read_list);
	count = 0
	read_list = "";
if count != 0:
    read_file.write(read_list)
read_file.close()

