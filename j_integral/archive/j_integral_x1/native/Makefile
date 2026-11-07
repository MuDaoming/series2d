CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra
CPPFLAGS += -I/usr/local/include
LDLIBS += -L/usr/local/lib -Wl,-rpath,/usr/local/lib -lflint -lgmp
SHARED := ../../j_integral/native
.PHONY: all test
all: ../build/j_series ../build/j_search ../build/verify_j_points
../build/j_series: $(SHARED)/series.cpp $(SHARED)/prepare.hpp
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
../build/j_search: $(SHARED)/search.cpp
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
../build/verify_j_points: $(SHARED)/verify_j_points.cpp $(SHARED)/series.cpp $(SHARED)/prepare.hpp
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
../build/test_coordinates: ../tests/test_coordinates.cpp $(SHARED)/series.cpp $(SHARED)/prepare.hpp
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
test: ../build/test_coordinates
	../build/test_coordinates ../examples/vac.json
