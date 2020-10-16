# Makefile
# GNU makefile for Nucleation Benchmark code
# Questions/comments to trevor.keller@nist.gov (Trevor Keller)

# includes
incdir = $(MMSP_PATH)/include

# compilers/flags
compiler = g++
pcompiler = mpic++
flags = -O3 -I $(incdir)
pflags = $(flags) -include mpi.h

all: nuc8a nuc8b nuc8c nuc8d
.PHONY: all clean

# default program
debug: nucleation.cpp nucleation.hpp par8a.hpp
	$(compiler) $(flags) -O0 -g $< -include par8a.hpp -o $@ -lz

nucleation: nucleation.cpp nucleation.hpp
	$(pcompiler) $(pflags) $< -include par8a.hpp -o $@ -lz

# benchmark problems
nuc8a: nucleation.cpp nucleation.hpp
	$(pcompiler) $(pflags) $< -include par8a.hpp -o $@ -lz

nuc8b: nucleation.cpp nucleation.hpp
	$(pcompiler) $(pflags) $< -include par8b.hpp -o $@ -lz

nuc8c: nucleation.cpp nucleation.hpp
	$(pcompiler) $(pflags) $< -include par8c.hpp -o $@ -lz

nuc8d: nucleation.cpp nucleation.hpp
	$(pcompiler) $(pflags) $< -include par8d.hpp -o $@ -lz

# cleanup
clean:
	rm -vf debug nucleation nuc8a nuc8b nuc8c nuc8d
