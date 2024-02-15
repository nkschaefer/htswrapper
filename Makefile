SHELL=bash
COMP=g++
PREFIX ?=/usr/local
FLAGS=-std=c++11 --std=gnu++11 -fPIC
IFLAGS=-I$(PREFIX)/include
LFLAGS=-L$(PREFIX)/lib
BC_LEN ?= 16

all: lib/libhtswrapper.so lib/libhtswrapper.a

lib/libhtswrapper.so: build/bam.o build/bc.o build/serialize.o
	$(COMP) -shared $(IFLAGS) $(LFLAGS) -o lib/libhtswrapper.so build/bam.o build/bc.o build/serialize.o -lz -lhts

lib/libhtswrapper.a: build/bam.o build/bc.o build/serialize.o
	ar rcs lib/libhtswrapper.a build/bam.o build/bc.o build/serialize.o
 
build/bam.o: src/bam.cpp src/bam.h
	$(COMP) $(IFLAGS) $(FLAGS) -c -o build/bam.o src/bam.cpp

build/bc.o: src/bc.cpp src/bc.h
	$(COMP) $(IFLAGS) $(FLAGS) -DBC_LEN=$(BC_LEN) -c -o build/bc.o src/bc.cpp

build/serialize.o: src/serialize.cpp src/serialize.h
	$(COMP) $(IFLAGS) $(FLAGS) -c -o build/serialize.o src/serialize.cpp

clean:
	rm build/*.o
	rm lib/*.so

install: | $(PREFIX)/lib $(PREFIX)/include/htswrapper
	cp lib/*.so $(PREFIX)/lib
	cp lib/*.a $(PREFIX)/lib
	cp src/bam.h $(PREFIX)/include/htswrapper
	cp src/bc.h $(PREFIX)/include/htswrapper
	cp src/serialize.h $(PREFIX)/include/htswrapper

$(PREFIX)/lib:
	mkdir -p $(PREFIX)/lib

$(PREFIX)/include/htswrapper:
	mkdir -p $(PREFIX)/include/htswrapper
