CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra
CPPFLAGS += -I/usr/local/include
LDLIBS += -L/usr/local/lib -Wl,-rpath,/usr/local/lib -lflint -lgmp
SHARED := ../../j_integral/native
METHODS := ../../j_integral_methods/native
COMMON := $(SHARED)/prepare.hpp $(SHARED)/series.cpp $(METHODS)/direct_j.cpp
TESTS := test_point_family test_triangular test_dense_coupling test_first_jet test_fbi_aliases
.PHONY: all test
all: ../build/probe ../build/global_point ../build/audit_fbi ../build/reconstruct
../build/probe: probe.cpp $(COMMON)
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
../build/global_point: global_point.cpp triangular.hpp $(COMMON)
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
../build/audit_fbi: audit_fbi.cpp triangular.hpp $(COMMON)
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
../build/reconstruct: reconstruct.cpp $(COMMON)
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
../build/test_%: ../tests/test_%.cpp triangular.hpp $(COMMON)
	mkdir -p ../build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
test: $(addprefix ../build/,$(TESTS))
	../build/test_point_family ../examples/dp.json
	../build/test_triangular
	../build/test_dense_coupling ../examples/small_coupled.json
	../build/test_first_jet ../examples/small_coupled.json
	../build/test_first_jet ../examples/minimal_unclosed.json
	../build/test_fbi_aliases ../examples/minimal_unclosed.json
	python3 ../tests/test_dispatch.py
