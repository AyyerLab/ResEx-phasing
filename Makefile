.PHONY: mkdir
CC=gcc
LDFLAGS=-lgsl -lgslcblas -lfftw3f_threads -lfftw3f -lm
CFLAGS=-fopenmp -O3 -Wall

src = $(wildcard src/*.c)
objects = $(patsubst src/%.c,bin/%.o,$(src)) 

utils_src = $(wildcard utils/src/*.c)
utils = $(patsubst utils/src/%.c,utils/%,$(utils_src))

#all: gen_data recon $(utils)
all: mkdir recon $(utils)

mkdir: bin

bin:
	mkdir -p bin

recon: $(filter-out bin/gen_data.o bin/hio.o, $(objects))
	$(LINK.c) $^ -o $@

$(objects): bin/%.o: src/%.c src/brcont.h
	$(CC) -c $< -o $@ $(CFLAGS)

gen_data: $(filter-out bin/recon.o, $(objects))
	$(LINK.c) $^ -o $@

$(utils): utils/%: utils/src/%.c
	$(CC) $< -o $@ $(LDFLAGS) $(CFLAGS)

clean:
	rm -f gen_data recon $(objects) $(utils)
