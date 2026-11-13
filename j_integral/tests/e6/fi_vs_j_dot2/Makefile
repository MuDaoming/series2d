CXX ?= g++
FLINT_PREFIX ?= /usr/local
ROOT := ../../../..
reconstruct: reconstruct.cpp $(wildcard $(ROOT)/j_integral/src/IBP/kernel/*)
	$(CXX) -O3 -std=c++17 -I$(ROOT)/j_integral/src/IBP/kernel -I$(FLINT_PREFIX)/include $< -L$(FLINT_PREFIX)/lib -Wl,-rpath,$(FLINT_PREFIX)/lib -lflint -lgmp -o $@
