SHELL=/bin/bash
COMP=g++
FLAGS=-std=c++11 --std=gnu++11 -fPIC
PREFIX ?=/usr/local

all: lib/libhtswrapper.so

lib/libhtswrapper.so: build/bam.o
	$(COMP) -shared -o lib/libhtswrapper.so build/bam.o -lhts

build/bam.o: src/bam.cpp src/bam.h
	$(COMP) $(FLAGS) -c -o build/bam.o src/bam.cpp

clean:
	rm build/*.o
	rm lib/libmixturedist.so

install: | $(PREFIX)/lib $(PREFIX)/include/htswrapper
	cp lib/libhtswrapper.so $(PREFIX)/lib
	cp src/bam.h $(PREFIX)/include/htswrapper

$(PREFIX)/lib:
	mkdir -p $(PREFIX)/lib

$(PREFIX)/include/htswrapper:
	mkdir -p $(PREFIX)/include/mixtureDist
