CC=g++
CFLAGS=-g -std=c++11
ndTree: main.o logClass.o utility.o
	$(CC) $(CFLAGS) -o ndTree main.o -lz utility.o logClass.o 
main.o: main.cpp ndtreeHelper.h	logClass.h config.h Dir_entry.h  ND_tree.h kseq.h Dir_node.h  Leaf_entry.h Node.h  dim.h utility.h  Leaf_node.h 
	$(CC) $(CFLAGS) -c main.cpp
utility.o: utility.cpp utility.h
	$(CC) $(CFLAGS) -c utility.cpp
logClass.o: logClass.cpp logClass.h
	$(CC) $(CFLAGS) -c logClass.cpp
run:
	./ndTree --dimension 5 --index ../data/index_real --data ../data/cancer/kmer/all.kmer --boxquery ../data/cancer/query/boxquery --mode rebuild --aux ../data/cancer/kmer/all.kmer.desc --record ../data/record --querydim 6
mac:
	lldb -- ndTree --dimension 16 --index ../data/index --data ../data/cancer/kmer/all.kmer --boxquery ../data/cancer/query/boxquery --mode rebuild   --aux ../data/cancer/kmer/all.kmer.desc --record ../data/record  --querydim 18
gdb:
	gdb --args ndTree --dimension 6 --index ../data/index.bin --data ../data/data.fa --boxquery ../data/cancer/query/boxquery --mode rebuild   --aux ../data/cancer/kmer/all.kmer.desc --record ../data/record  --querydim 6
simplegdb:
	gdb --args ndTree --dimension 6 --mode rebuild  --querydim 7 --boxquery ../data/cancer/query/boxquery 
query:
	gdb --args ndTree --dimension 5 --index ../data/index --boxquery ../data/cancer/query/boxquery --record ../data/record  --querydim 6
clean:
	rm -f *.o

# Dir_node.h  Leaf_entry.h Node.h  dim.h utility.h  Leaf_node.h 
