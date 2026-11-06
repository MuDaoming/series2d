CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra
CPPFLAGS += -I/usr/local/include
LDLIBS += -L/usr/local/lib -Wl,-rpath,/usr/local/lib -lflint -lgmp
SHARED := ../../j_integral/native
METHODS := ../../j_integral_methods/native
GLOBAL := ../../j_integral_global/native
.PHONY: all kernels
all: ../build/cut_diagonal ../build/reconstruct_partial ../build/verify_connection ../build/verify_reductions
../build/cut_diagonal: cut_diagonal.cpp $(SHARED)/prepare.hpp $(SHARED)/series.cpp $(METHODS)/direct_j.cpp $(GLOBAL)/triangular.hpp
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
kernels:
	python3 ../tests/build_kernels.py
../build/reconstruct_partial: reconstruct_partial.cpp $(SHARED)/prepare.hpp $(SHARED)/series.cpp $(METHODS)/direct_j.cpp
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
../build/verify_connection: verify_connection.cpp $(SHARED)/prepare.hpp $(SHARED)/series.cpp $(METHODS)/direct_j.cpp
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

../build/verify_reductions: verify_reductions.cpp $(SHARED)/prepare.hpp $(SHARED)/series.cpp $(METHODS)/direct_j.cpp
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
