SHELL=/bin/bash
COMP=g++
FLAGS=-std=c++11 --std=gnu++11 -fPIC
PREFIX ?=/usr/local

all: lib/libhtswrapper.so

lib/libhtswrapper.so: build/bam.o build/bc_hash.o
	$(COMP) -shared -o lib/libhtswrapper.so build/bam.o build/bc_hash.o -lhts

build/bam.o: src/bam.cpp src/bam.h
	$(COMP) $(FLAGS) -c -o build/bam.o src/bam.cpp

build/bc_hash.o: src/bc_hash.cpp src/bc_hash.h
	$(COMP) $(FLAGS) -c -o build/bc_hash.o src/bc_hash.cpp

clean:
	rm build/*.o
	rm lib/libmixturedist.so

install: | $(PREFIX)/lib $(PREFIX)/include/htswrapper
	cp lib/libhtswrapper.so $(PREFIX)/lib
	cp src/bam.h $(PREFIX)/include/htswrapper
	cp src/bc_hash.h $(PREFIX)/include/htswrapper

$(PREFIX)/lib:
	mkdir -p $(PREFIX)/lib

$(PREFIX)/include/htswrapper:
	mkdir -p $(PREFIX)/include/htswrapper
