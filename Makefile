.PHONY: mkdir
CC=gcc

LDFLAGS=-lfftw3f_threads -lfftw3f $(shell gsl-config --libs) -Wl,-rpath,$(shell gsl-config --prefix)/lib
CFLAGS=$(shell gsl-config --cflags) -std=gnu99 -fopenmp -O3 -Wall -Wno-unused-result -Wno-format-overflow

recon_src = $(wildcard src/*.c)
recon_obj = $(patsubst src/%.c, bin/%.o, $(recon_src)) 
utils_src = $(wildcard utils/src/*.c)
utils = $(patsubst utils/src/%.c, utils/%, $(utils_src))
directories = data bin images results data/maps data/recon data/convert

all: mkdir recon $(utils)

mkdir: $(directories)
$(directories):
	mkdir -p $(directories)

recon: $(filter-out bin/utils.o, $(recon_obj))
	$(CC) $^ -o $@ $(CFLAGS) $(LDFLAGS)
$(filter-out bin/recon.o, $(recon_obj)): bin/%.o: src/%.c src/%.h
	$(CC) -c $< -o $@ $(CFLAGS)
bin/recon.o: src/recon.c src/algorithm.h
	$(CC) -c $< -o $@ $(CFLAGS)

utils/gen_fdens: bin/fft.o bin/volume.o
utils/band_limit utils/boost_cont utils/calc_fsc utils/combine: bin/fft.o
utils/create_support utils/gen_dens utils/liquidize: bin/fft.o
utils/sharpen: bin/volume.o
$(utils): utils/%: utils/src/%.c bin/utils.o
	$(CC) $^ -o $@ $(LDFLAGS) $(CFLAGS)

clean:
	rm -f recon $(recon_obj) $(utils)
