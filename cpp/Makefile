CC=g++
OBJECTS=ProblemParameters.o State.o Ring.o

%.o: %.c %.h
	$(CC) -c -o $@ $<

test: test.cpp $(OBJECTS)
	$(CC) -o $@ $^
	./tests

cuda:
	$(warning The CUDA-compliant version is still not ready.)

seq: main.cpp $(OBJECTS)
	$(CC) -o $@ $^

all: cuda seq

clear:
	rm -f *.o
	rm -f tests
	rm -f seq
	rm -f cuda
