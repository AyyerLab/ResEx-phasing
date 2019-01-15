.PHONY: mkdir
CC=gcc

FFTW_LIBS=$(shell pkg-config --libs fftw3f) -lfftw3f_threads -Wl,-rpath,$(shell pkg-config --libs-only-L fftw3f|cut -c 3-)
GSL_LIBS=$(shell gsl-config --libs) -Wl,-rpath,$(shell gsl-config --prefix)/lib
LDFLAGS=$(FFTW_LIBS) $(GSL_LIBS)
CFLAGS=$(shell gsl-config --cflags) $(shell pkg-config --cflags fftw3) -std=gnu99 -fopenmp -O3 -Wall

ifeq ($(CC), gcc)
	CFLAGS += -Wno-unused-result
endif

recon_src = $(wildcard src/*.c)
recon_obj = $(patsubst src/%.c, bin/%.o, $(recon_src)) 
utils_src = $(wildcard utils/src/*.c)
utils = $(patsubst utils/src/%.c, utils/%, $(utils_src))
directories = data bin images results data/recon data/convert data/merges

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
