CFLAGS=-Wall -Wextra -Wpedantic -g -std=gnu99 -lfftw3_threads -lfftw3 -lm -lpthread -lgsl -lgslcblas -O2
.DEFAULT:= all
.PHONY: all

OBJECTS_MAIN=main.o file_io.o
OBJECTS_SHARED=util.o estimator.o vector.o thread_pool.o transformations.o
OBJECTS_TEST = lattice-fields-testing-ver/main.o lattice-fields-testing-ver/field_gen.o

all: main lattice-fields-test


main: $(OBJECTS_MAIN) $(OBJECTS_SHARED)
	 cc $(CFLAGS) $(OBJECTS_MAIN) $(OBJECTS_SHARED) -o main 


lattice-fields-test: $(OBJECTS_TEST) $(OBJECTS_SHARED)
	 cc $(CFLAGS) $(OBJECTS_TEST) $(OBJECTS_SHARED) -o lattice-fields-testing-ver/test 

clean:
	rm *.o main lattice-fields-testing-ver/*.o lattice-fields-testing-ver/test
