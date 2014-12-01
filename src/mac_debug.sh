#!/bin/sh
lldb -- ndTree --dimension 16 --index ../data/index_real --data ../data/data_random --boxquery ../data/query_random --skip 0 --newtree 1 --aux ../data/cancer/kmer/all.kmer.desc --record ../data/record --querydim 16
