.PHONY: mkdir
CC=gcc
LDFLAGS=-lgsl -lgslcblas -lfftw3f_threads -lfftw3f -lm
CFLAGS=-fopenmp -O3 -Wall

src = $(wildcard src/*.c)
objects = $(patsubst src/%.c,bin/%.o,$(src)) 

utils_src = $(wildcard utils/src/*.c)
utils = $(patsubst utils/src/%.c,utils/%,$(utils_src))

directories = bin images

all: mkdir gen_data recon $(utils)

mkdir: $(directories)

$(directories):
	mkdir -p $(directories)

recon: $(filter-out bin/gen_data.o bin/hio.o, $(objects))
	$(LINK.c) $^ -o $@

$(objects): bin/%.o: src/%.c src/brcont.h
	$(CC) -c $< -o $@ $(CFLAGS)

gen_data: bin/gen_data.o bin/utils.o bin/setup.o
	$(LINK.c) $^ -o $@

$(utils): utils/%: utils/src/%.c
	$(CC) $< -o $@ $(LDFLAGS) $(CFLAGS)

clean:
	rm -f gen_data recon $(objects) $(utils)
