CC=gcc
LDFLAGS=-lgsl -lgslcblas -lfftw3f_threads -lfftw3f -lm
CFLAGS=-fopenmp -O3 -Wall

objects = bin/diffmap.o bin/setup.o bin/utils.o
utils = utils/chop_bragg utils/assemble_hkl utils/gen_bmerge
autils = utils/forward utils/fstretch utils/bragg_gen

all: gen_data recon $(utils) $(autils)

bin/%.o: src/%.c src/brcont.h
	gcc -c $< -o $@ $(CFLAGS)

recon: bin/recon.o $(objects)
	$(LINK.c) $^ -o $@

gen_data: bin/gen_data.o $(objects)
	$(LINK.c) $^ -o $@

utils/%: utils/src/%.c
	gcc $< -o $@ $(LDFLAGS) $(CFLAGS)

clean:
	rm -f recon gen_data bin/gen_data.o bin/recon.o $(objects) $(utils) $(autils)
