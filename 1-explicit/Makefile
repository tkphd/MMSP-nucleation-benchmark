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

all: mmsp2norm nuc8.1a nuc8.1b nuc8.1c nuc8.1d1 nuc8.1d2 nuc8.1d3 nuc8.1d4 nuc8.1d5 nuc8.1d6 nuc8.1d7 nuc8.1d8
.PHONY: all clean

mmsp2norm: mmsp2norm.cpp nucleation.hpp
	$(compiler) $(flags) $< -o $@ -lz

debug: nucleation.cpp nucleation.hpp par8.1a.hpp
	$(compiler) $(flags) -O0 -g $< -include par8.1a.hpp -o $@ -lz

nucleation: nucleation.cpp nucleation.hpp
	$(pcompiler) $(pflags) $< -include par8.1a.hpp -o $@ -lz

# benchmark problems
nuc8.1a: nucleation.cpp nucleation.hpp par8.1a.hpp
	$(pcompiler) $(pflags) $< -include par8.1a.hpp -o $@ -lz

nuc8.1b: nucleation.cpp nucleation.hpp par8.1b.hpp
	$(pcompiler) $(pflags) $< -include par8.1b.hpp -o $@ -lz

nuc8.1c: nucleation.cpp nucleation.hpp par8.1c.hpp
	$(pcompiler) $(pflags) $< -include par8.1c.hpp -o $@ -lz

nuc8.1d1: nucleation.cpp nucleation.hpp par8.1d1.hpp
	$(pcompiler) $(pflags) $< -include par8.1d1.hpp -o $@ -lz

nuc8.1d2: nucleation.cpp nucleation.hpp par8.1d2.hpp
	$(pcompiler) $(pflags) $< -include par8.1d2.hpp -o $@ -lz

nuc8.1d3: nucleation.cpp nucleation.hpp par8.1d3.hpp
	$(pcompiler) $(pflags) $< -include par8.1d3.hpp -o $@ -lz

nuc8.1d4: nucleation.cpp nucleation.hpp par8.1d4.hpp
	$(pcompiler) $(pflags) $< -include par8.1d4.hpp -o $@ -lz

nuc8.1d5: nucleation.cpp nucleation.hpp par8.1d5.hpp
	$(pcompiler) $(pflags) $< -include par8.1d5.hpp -o $@ -lz

nuc8.1d6: nucleation.cpp nucleation.hpp par8.1d6.hpp
	$(pcompiler) $(pflags) $< -include par8.1d6.hpp -o $@ -lz

nuc8.1d7: nucleation.cpp nucleation.hpp par8.1d7.hpp
	$(pcompiler) $(pflags) $< -include par8.1d7.hpp -o $@ -lz

nuc8.1d8: nucleation.cpp nucleation.hpp par8.1d8.hpp
	$(pcompiler) $(pflags) $< -include par8.1d8.hpp -o $@ -lz

# cleanup
clean:
	rm -vf debug nucleation nuc8.1a nuc8.1b nuc8.1c \
           nuc8.1d1 nuc8.1d2 nuc8.1d3 nuc8.1d4 \
           nuc8.1d5 nuc8.1d6 nuc8.1d7 nuc8.1d8
