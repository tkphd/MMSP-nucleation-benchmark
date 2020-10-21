# Makefile
# GNU makefile for Nucleation Benchmark code
# Questions/comments to trevor.keller@nist.gov (Trevor Keller)

# includes
incdir = $(MMSP_PATH)/include

# compilers/flags
compiler = g++
flags = -O3 -fopenmp -I $(incdir)

all: mmsp2norm nuc8.4a nuc8.4b nuc8.4c nuc8.4d1 nuc8.4d2 nuc8.4d3 nuc8.4d4 nuc8.4d5 nuc8.4d6 nuc8.4d7
.PHONY: all clean

mmsp2norm: mmsp2norm.cpp nucleation.hpp
	$(compiler) $(flags) $< -o $@ -lz

debug: nucleation.cpp nucleation.hpp par8.4a.hpp
	$(compiler) $(flags) -O0 -g $< -include par8.4a.hpp -o $@ -lz

nucleation: nucleation.cpp nucleation.hpp
	$(compiler) $(flags) $< -include par8.4a.hpp -o $@ -lz

# benchmark problems
nuc8.4a: nucleation.cpp nucleation.hpp par8.4a.hpp
	$(compiler) $(flags) $< -include par8.4a.hpp -o $@ -lz

nuc8.4b: nucleation.cpp nucleation.hpp par8.4b.hpp
	$(compiler) $(flags) $< -include par8.4b.hpp -o $@ -lz

nuc8.4c: nucleation.cpp nucleation.hpp par8.4c.hpp
	$(compiler) $(flags) $< -include par8.4c.hpp -o $@ -lz

nuc8.4d1: nucleation.cpp nucleation.hpp par8.4d1.hpp
	$(compiler) $(flags) $< -include par8.4d1.hpp -o $@ -lz

nuc8.4d2: nucleation.cpp nucleation.hpp par8.4d2.hpp
	$(compiler) $(flags) $< -include par8.4d2.hpp -o $@ -lz

nuc8.4d3: nucleation.cpp nucleation.hpp par8.4d3.hpp
	$(compiler) $(flags) $< -include par8.4d3.hpp -o $@ -lz

nuc8.4d4: nucleation.cpp nucleation.hpp par8.4d4.hpp
	$(compiler) $(flags) $< -include par8.4d4.hpp -o $@ -lz

nuc8.4d5: nucleation.cpp nucleation.hpp par8.4d5.hpp
	$(compiler) $(flags) $< -include par8.4d5.hpp -o $@ -lz

nuc8.4d6: nucleation.cpp nucleation.hpp par8.4d6.hpp
	$(compiler) $(flags) $< -include par8.4d6.hpp -o $@ -lz

nuc8.4d7: nucleation.cpp nucleation.hpp par8.4d7.hpp
	$(compiler) $(flags) $< -include par8.4d7.hpp -o $@ -lz

# cleanup
clean:
	rm -vf debug mmsp2norm nucleation \
           nuc8.4a nuc8.4b nuc8.4c \
           nuc8.4d1 nuc8.4d2 nuc8.4d3 nuc8.4d4 \
           nuc8.4d5 nuc8.4d6 nuc8.4d7
