CC=clang++
CFLAGS=-c -std=c++17
PROGRAMFILE=mesh

all: mesh

mesh: main.o
	$(CC) main.o -o $(PROGRAMFILE)

main.o:  mesh_main.cpp #main.h
	$(CC) $(CFLAGS) mesh_main.cpp -o main.o

clean:
	rm -rf *.o $(PROGRAMFILE)
