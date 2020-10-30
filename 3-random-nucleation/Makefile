# Makefile
# GNU makefile for Nucleation Benchmark code
# Questions/comments to trevor.keller@nist.gov (Trevor Keller)

# includes
incdir = $(MMSP_PATH)/include

# compilers/flags
compiler = g++
flags = -O3 -fopenmp -Dcimg_display=0 -I $(incdir)

icompiler = icc
iflags = -O3 -qopenmp -Dcimg_display=0 -I $(incdir) -std=c++11 -xCORE-AVX2

all: nucleation
.PHONY: all clean

debug: nucleation.cpp nucleation.hpp par8.3.hpp
	$(compiler) $(flags) -O0 -g $< -include par8.3.hpp -o $@ -lz

nucleation: nucleation.cpp nucleation.hpp par8.3.hpp
	$(compiler) $(flags) $< -include par8.3.hpp -o $@ -lz

# cleanup
clean:
	rm -vf debug nucleation
