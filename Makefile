SHELL=bash
COMP=g++
CCOMP=gcc
PREFIX ?=/usr/local
CFLAGS = -fPIC
FLAGS=-std=c++11 --std=gnu++11 -fPIC -O3 -g
IFLAGS=-I$(PREFIX)/include
LFLAGS=-L$(PREFIX)/lib
ifeq ($(findstring cellbouncer, ${CONDA_PREFIX}), cellbouncer)
	IFLAGS += -I${CONDA_PREFIX}/include
	LFLAGS += -L${CONDA_PREFIX}/lib
endif
BC_LENX2 ?= 32
KX2 ?= 16

all: hash_bc unhash_bc lib/libhtswrapper.so lib/libhtswrapper.a

lib/libhtswrapper.so: build/bam.o build/bc.o build/bc_scanner.o build/umi.o build/edlib.o build/seq_fuzzy_match.o build/serialize.o build/gzreader.o build/khashtable.o build/kmsuftree.o
	$(COMP) -shared $(IFLAGS) $(LFLAGS) -o lib/libhtswrapper.so build/bam.o build/bc.o build/bc_scanner.o build/umi.o build/edlib.o build/seq_fuzzy_match.o build/serialize.o build/gzreader.o build/khashtable.o build/kmsuftree.o -lz -lhts

lib/libhtswrapper.a: build/bam.o build/bc.o build/bc_scanner.o build/umi.o build/edlib.o build/seq_fuzzy_match.o build/serialize.o build/gzreader.o build/khashtable.o build/kmsuftree.o
	ar rcs lib/libhtswrapper.a build/bam.o build/bc.o build/bc_scanner.o build/umi.o build/edlib.o build/seq_fuzzy_match.o build/serialize.o build/gzreader.o build/khashtable.o

hash_bc: src/bc.h build/bc.o build/gzreader.o
	$(COMP) $(IFLAGS) $(FLAGS) -DBC_LENX2=$(BC_LENX2) -DKX2=$(KX2) -o hash_bc src/hash_bc.cpp build/gzreader.o build/bc.o -lz

unhash_bc: src/bc.h build/bc.o build/gzreader.o
	$(COMP) $(IFLAGS) $(FLAGS) -DBC_LENX2=$(BC_LENX2) -DKX2=$(KX2) -o unhash_bc src/unhash_bc.cpp build/gzreader.o build/bc.o -lz 
wl_test: src/wl_test.cpp src/bc.h build/bc.o
	$(COMP) -DBC_LENX2=$(BC_LENX2) -DKX2=$(KX2) -g build/bc.o build/gzreader.o src/wl_test.cpp $(LFLAGS) -o wl_test -lz

khashtable_test: src/khashtable_test.cpp src/khashtable.h build/khashtable.o build/bc.o 
	$(COMP) $(CXXFLAGS) $(CXXIFLAGS) -DBC_LENX2=$(BC_LENX2) -DKX2=$(KX2) -g build/khashtable.o build/bc.o build/gzreader.o src/khashtable_test.cpp $(LFLAGS) -o khashtable_test -lz

kmsuftree_test: src/kmsuftree_test.cpp src/kmsuftree.h build/kmsuftree.o
	$(COMP) $(CXXFLAGS) $(CXXIFLAGS) build/kmsuftree.o src/kmsuftree_test.cpp $(LFLAGS) $(DEPS) -o kmsuftree_test -lz

build/bam.o: src/bam.cpp src/bam.h
	$(COMP) $(IFLAGS) $(FLAGS) -c -o build/bam.o src/bam.cpp

build/bc.o: src/bc.cpp src/bc.h
	$(COMP) $(IFLAGS) $(FLAGS) -DBC_LENX2=$(BC_LENX2) -DKX2=$(KX2) -c -o build/bc.o src/bc.cpp

build/bc_scanner.o: src/bc_scanner.cpp src/bc_scanner.h src/bc.h
	$(COMP) $(IFLAGS) $(FLAGS) -DBC_LENX2=$(BC_LENX2) -DKX2=$(KX2) -c -o build/bc_scanner.o src/bc_scanner.cpp

build/umi.o: src/umi.cpp src/umi.h src/bc.h
	$(COMP) $(IFLAGS) $(FLAGS) -DBC_LENX2=$(BC_LENX2) -DKX2=$(KX2) -c -o build/umi.o src/umi.cpp

build/seq_fuzzy_match.o: src/seq_fuzzy_match.cpp src/seq_fuzzy_match.h build/edlib.o
	$(COMP) $(IFLAGS) $(FLAGS) -c -o build/seq_fuzzy_match.o src/seq_fuzzy_match.cpp

build/edlib.o: src/edlib/edlib.cpp src/edlib/edlib.h
	$(COMP) $(IFLAGS) $(FLAGS) -c -o build/edlib.o src/edlib/edlib.cpp

build/serialize.o: src/serialize.cpp src/serialize.h
	$(COMP) $(IFLAGS) $(FLAGS) -c -o build/serialize.o src/serialize.cpp

build/gzreader.o: src/gzreader.cpp src/gzreader.h
	$(COMP) $(IFLAGS) $(FLAGS) -c -o build/gzreader.o src/gzreader.cpp

build/khashtable.o: src/khashtable.cpp src/khashtable.h $(DEPS)
	$(COMP) $(IFLAGS) $(FLAGS) -DBC_LENX2=$(BC_LENX2) -DKX2=$(KX2) -g src/khashtable.cpp -c -o build/khashtable.o

build/kmsuftree.o: src/kmsuftree.c src/kmsuftree.h
	$(CCOMP) $(CFLAGS) src/kmsuftree.c -c -o build/kmsuftree.o

clean:
	rm build/*.o
	rm lib/*.so
	rm lib/*.a

install: | $(PREFIX)/lib $(PREFIX)/include/htswrapper
	cp lib/*.so $(PREFIX)/lib
	cp lib/*.a $(PREFIX)/lib
	cp src/bam.h $(PREFIX)/include/htswrapper
	cp src/bc.h $(PREFIX)/include/htswrapper
	cp src/bc_scanner.h $(PREFIX)/include/htswrapper
	cp src/serialize.h $(PREFIX)/include/htswrapper
	cp src/gzreader.h $(PREFIX)/include/htswrapper
	cp src/khashtable.h $(PREFIX)/include/htswrapper
	cp src/kmsuftree.h $(PREFIX)/include/htswrapper
	cp src/umi.h $(PREFIX)/include/htswrapper
	cp src/seq_fuzzy_match.h $(PREFIX)/include/htswrapper
	cp src/robin_hood/robin_hood.h $(PREFIX)/include/htswrapper/robin_hood
	cp src/robin_hood/LICENSE $(PREFIX)/include/htswrapper/robin_hood
	cp src/edlib/LICENSE $(PREFIX)/include/htswrapper/edlib
	cp src/edlib/edlib.h $(PREFIX)/include/htswrapper/edlib

$(PREFIX)/lib:
	mkdir -p $(PREFIX)/lib

$(PREFIX)/include/htswrapper:
	mkdir -p $(PREFIX)/include/htswrapper
	mkdir -p $(PREFIX)/include/htswrapper/robin_hood
	mkdir -p $(PREFIX)/include/htswrapper/edlib
