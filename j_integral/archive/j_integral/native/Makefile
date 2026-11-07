CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra
CPPFLAGS += -I/usr/local/include
LDLIBS += -L/usr/local/lib -Wl,-rpath,/usr/local/lib -lflint -lgmp
../build/j_series: series.cpp prepare.hpp
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

../build/check_series: check_series.cpp
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

../build/j_search: search.cpp
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

../build/tests: tests.cpp series.cpp prepare.hpp
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) tests.cpp -o $@ $(LDLIBS)

.PHONY: all test
all: ../build/j_series ../build/j_search ../build/check_series

test: ../build/tests
	../build/tests

../build/fi_ibp: fi_ibp.cpp search.cpp
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

../build/compare_blade: compare_blade.cpp
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

all: ../build/fi_ibp ../build/compare_blade

../build/verify_j_points: verify_j_points.cpp series.cpp prepare.hpp
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
