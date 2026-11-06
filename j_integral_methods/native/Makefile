CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra
CPPFLAGS += -I/usr/local/include
LDLIBS += -L/usr/local/lib -Wl,-rpath,/usr/local/lib -lflint -lgmp
SHARED := ../../j_integral/native
.PHONY: all test
all: ../build/direct_j ../build/verify_j_points
../build/direct_j: direct_j.cpp $(SHARED)/series.cpp $(SHARED)/prepare.hpp
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
../build/verify_j_points: $(SHARED)/verify_j_points.cpp $(SHARED)/series.cpp $(SHARED)/prepare.hpp
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
../build/tests: ../tests/tests.cpp direct_j.cpp $(SHARED)/series.cpp $(SHARED)/prepare.hpp
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
test: ../build/tests
	../build/tests
../build/symbolic_j: symbolic_j.cpp direct_j.cpp $(SHARED)/series.cpp $(SHARED)/prepare.hpp
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
all: ../build/symbolic_j
../build/baseline_series: $(SHARED)/series.cpp $(SHARED)/prepare.hpp
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
../build/baseline_search: $(SHARED)/search.cpp
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
.PHONY: baseline
baseline: ../build/baseline_series ../build/baseline_search
