CC=clang++
CFLAGS=-c -std=c++17
LFLAGS= -Xpreprocessor -fopenmp -lomp
PROGRAMFILE=mesh

all: mesh

mesh: main.o
	$(CC) $(LFLAGS) main.o -o $(PROGRAMFILE)

main.o:  mesh_main.cpp makefile #main.h
	$(CC) $(CFLAGS) mesh_main.cpp -o main.o

clean:
	rm -rf *.o $(PROGRAMFILE)
