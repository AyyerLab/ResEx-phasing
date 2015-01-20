CC=gcc
LDFLAGS=-lgsl -lgslcblas -lfftw3f_threads -lfftw3f -lm
CFLAGS=-fopenmp -O3 -Wall
UTILSCFLAGS=-lfftw3f -lm -O3 -Wall

objects = bin/recon.o bin/diffmap.o bin/setup.o bin/utils.o

all: recon $(objects)

bin/%.o: src/%.c src/brcont.h
	gcc -c $< -o $@ $(CFLAGS)

recon: $(objects)
	$(LINK.c) $^ -o $@

utils/%: src/utils/%.c
	gcc $< -o $@ $(UTILSCFLAGS)

clean:
#	rm -f recon $(objects) $(miscutils) $(datautils) $(supportutils) $(transformutils)
	rm -f recon $(objects)
