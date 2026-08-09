CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra

all: series_basis_solver

series_basis_solver: series_basis_solver.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f series_basis_solver
