CFLAGS=-Wall -Wextra -Wpedantic -g -std=gnu99 -lfftw3_threads -lfftw3 -lm -lpthread -lgsl -lgslcblas -O2
.DEFAULT:= all
.PHONY: all

OBJECTS=main.o util.o estimator.o vector.o thread_pool.o transformations.o file_io.o

all: main


main: $(OBJECTS)
	 cc $(CFLAGS) $(OBJECTS) -o main 

clean:
	rm *.o main
