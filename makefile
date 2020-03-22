CC=clang++
CFLAGS=-c
PROGRAMFILE=mesh

all: mesh

mesh: main.o
	$(CC) main.o -o $(PROGRAMFILE)

main.o:  mesh_main.cpp #main.h
	$(CC) $(CFLAGS) mesh_main.cpp

clean:
	rm -rf *.o $(PROGRAMFILE)
