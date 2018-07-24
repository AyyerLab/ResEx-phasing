.PHONY: mkdir
CC=gcc
LDFLAGS=$(shell gsl-config --libs) -lfftw3f_threads -lfftw3f -lm
CFLAGS=$(shell gsl-config --cflags) -std=gnu99 -fopenmp -O3 -Wall -g

recon_src = $(wildcard src/*.c)
recon_obj = $(patsubst src/%.c, bin/%.o, $(recon_src)) 
utils_src = $(wildcard utils/src/*.c)
utils = $(patsubst utils/src/%.c, utils/%, $(utils_src))
directories = data bin images results data/maps data/recon data/convert

all: mkdir recon $(utils) utils/calc_smoothness

mkdir: $(directories)
$(directories):
	mkdir -p $(directories)

recon: $(filter-out bin/utils.o, $(recon_obj))
	$(LINK.c) $^ -o $@
$(filter-out bin/recon.o, $(recon_obj)): bin/%.o: src/%.c src/%.h
	$(CC) -c $< -o $@ $(CFLAGS)
bin/recon.o: src/recon.c src/algorithm.h
	$(CC) -c $< -o $@ $(CFLAGS)

utils/gen_fdens: bin/fft.o bin/volume.o
utils/band_limit utils/boost_cont utils/calc_fsc utils/combine utils/create_support utils/gen_dens utils/liquidize: bin/fft.o
utils/sharpen: bin/volume.o
$(utils): utils/%: utils/src/%.c bin/utils.o
	$(CC) $^ -o $@ $(LDFLAGS) $(CFLAGS)

clean:
	rm -f recon $(recon_obj) $(utils)
