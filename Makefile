.PHONY: mkdir
CC=gcc
LDFLAGS=-lgsl -lgslcblas -lfftw3f_threads -lfftw3f -lm
CFLAGS=-std=gnu99 -fopenmp -O3 -Wall -g

recon_src = $(wildcard src/*.c)
recon_obj = $(patsubst src/%.c, bin/%.o, $(recon_src)) 
utils_src = $(wildcard utils/src/*.c)
utils = $(patsubst utils/src/%.c, utils/%, $(utils_src))
directories = data bin images results data/maps data/recon data/convert

all: mkdir recon $(utils) utils/calc_smoothness

mkdir: $(directories)
$(directories):
	mkdir -p $(directories)

recon: $(recon_obj)
	$(LINK.c) $^ -o $@
$(filter-out bin/recon.o, $(recon_obj)): bin/%.o: src/%.c src/%.h
	$(CC) -c $< -o $@ $(CFLAGS)
bin/recon.o: src/recon.c src/algorithm.h
	$(CC) -c $< -o $@ $(CFLAGS)

utils/gen_fdens: bin/fft.o
utils/gen_dens: bin/fft.o
$(utils): utils/%: utils/src/%.c
	$(CC) $^ -o $@ $(LDFLAGS) $(CFLAGS)

clean:
	rm -f recon $(recon_obj) $(utils)
