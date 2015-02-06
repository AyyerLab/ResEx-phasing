CC=gcc
LDFLAGS=-lgsl -lgslcblas -lfftw3f_threads -lfftw3f -lm
CFLAGS=-fopenmp -O3 -Wall

objects = bin/diffmap.o bin/setup.o bin/utils.o

all: gen_data recon bin/gen_data.o bin/recon.o $(objects)

bin/%.o: src/%.c src/brcont.h
	gcc -c $< -o $@ $(CFLAGS)

recon: bin/recon.o $(objects)
	$(LINK.c) $^ -o $@

gen_data: bin/gen_data.o $(objects)
	$(LINK.c) $^ -o $@

clean:
	rm -f recon gen_data bin/gen_data.o bin/recon.o $(objects)
