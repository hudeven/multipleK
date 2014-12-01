#!/bin/sh
sudo gdb --args ndTree --dimension 16 --index ../data/index_real --data ../data/cancer/kmer/all.kmer --boxquery ../data/cancer/query/boxquery --skip 0 --newtree 1 --aux ../data/cancer/kmer/all.kmer.desc --record ../data/record --querydim 6
