CC=g++
CFLAGS=-g -c
ndTree: main.o logClass.o utility.o 
	$(CC) -o ndTree main.o utility.o logClass.o 
main.o: main.cpp logClass.h ndtreeHelper.h ndtreeHelper.h config.h Dir_entry.h Dir_node.h Leaf_entry.h Leaf_node.h logClass.h ND_tree.h Node.h utility.h
	$(CC) $(CFLAGS) main.cpp
utility.o: utility.cpp utility.h
	$(CC) $(CFLAGS) utility.cpp
logClass.o: logClass.cpp logClass.h
	$(CC) $(CFLAGS) logClass.cpp
run:
	./ndTree --dimension 5 --index ../data/index_real --data ../data/cancer/kmer/all.kmer --boxquery ../data/cancer/query/boxquery --skip 0 --newtree 1 --aux ../data/cancer/kmer/all.kmer.desc --record ../data/record --querydim 6
gdb:
	gdb --args ndTree --dimension 5 --index ../data/index --data ../data/cancer/kmer/all.kmer --boxquery ../data/cancer/query/boxquery --mode rebuild   --aux ../data/cancer/kmer/all.kmer.desc --record ../data/record  --querydim 6
query:
	gdb --args ndTree --dimension 5 --index ../data/index --boxquery ../data/cancer/query/boxquery --record ../data/record  --querydim 6
clean:
	rm -f *.o
