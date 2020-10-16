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

all: nuc8.1a nuc8.1b nuc8.1c
.PHONY: all clean

# default program
debug: nucleation.cpp nucleation.hpp par8.1a.hpp
	$(compiler) $(flags) -O0 -g $< -include par8.1a.hpp -o $@ -lz

nucleation: nucleation.cpp nucleation.hpp
	$(pcompiler) $(pflags) $< -include par8.1a.hpp -o $@ -lz

# benchmark problems
nuc8.1a: nucleation.cpp nucleation.hpp
	$(pcompiler) $(pflags) $< -include par8.1a.hpp -o $@ -lz

nuc8.1b: nucleation.cpp nucleation.hpp
	$(pcompiler) $(pflags) $< -include par8.1b.hpp -o $@ -lz

nuc8.1c: nucleation.cpp nucleation.hpp
	$(pcompiler) $(pflags) $< -include par8.1c.hpp -o $@ -lz

nuc8.1d: nucleation.cpp nucleation.hpp
	$(pcompiler) $(pflags) $< -include par8.1d.hpp -o $@ -lz

# cleanup
clean:
	rm -vf debug nucleation nuc8.1a nuc8.1b nuc8.1c nuc8.1d
