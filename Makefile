# Makefile
# GNU makefile for example Nucleation model code
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# includes
incdir = $(MMSP_PATH)/include

# compilers/flags
compiler = g++
pcompiler = mpic++
flags = -O3 -I $(incdir)
pflags = $(flags) -include mpi.h

# the program
nucleation: nucleation.cpp
	$(compiler) $(flags) $< -o $@ -lz

parallel: nucleation.cpp
	$(pcompiler) $(pflags) $< -o $@ -lz

clean:
	rm -f nucleation parallel
