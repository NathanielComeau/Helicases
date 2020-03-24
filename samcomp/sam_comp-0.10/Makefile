all: sam_comp

#CC=/nfs/users/nfs_r/rmd/opt/llvm/bin/clang++
#CFLAGS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -Wall -g

CC = g++
CXX = g++
#CFLAGS  = -g -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_FORTIFY_SOURCE=2 -fstack-protector-all
#CFLAGS  = -O3 -g -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
CFLAGS  = -O3 -g -fomit-frame-pointer -fstrict-aliasing -ffast-math -DNDEBUG -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -Wall

#CC=icpc
#CXX=icpc
#CFLAGS  = -O3 -g -fomit-frame-pointer -fstrict-aliasing -fast -DNDEBUG -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -Wall -m64 -mia32 -axAVX

LIBS=-lz

.c.o:
	$(CC) $(CFLAGS) -c $< 

.cpp.o:
	$(CXX) $(CFLAGS) -c $< 

SAM_COMP_OBJ = sam_comp.o bam.o hash_table.o pooled_alloc.o

sam_comp: $(SAM_COMP_OBJ)
	$(CXX) $(CFLAGS) -o $@ $(SAM_COMP_OBJ) $(LIBS)

clean:
	-rm *.o
	-rm sam_comp
