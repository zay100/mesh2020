CC?=mpicc

OPTFLAGS = -O2
INCLUDES = -I.
COPTIONS = -w

CFLAGS = $(COPTIONS) $(OPTFLAGS) $(INCLUDES)
LD = $(CC) -L. 


OBJS = rw_mesh.o  rw_mesh_header.o  rw_mesh_mesh.o  rw_mesh_off.o  rw_mesh_string.o  rw_mesh_vtk.o

.c.o:
	$(CC) $(CFLAGS) -c $*.c

rw_mesh.a: $(OBJS)
	ar rv $@ $(OBJS)
	ranlib $@

clean:
	rm -rf *.o rw_mesh.a
