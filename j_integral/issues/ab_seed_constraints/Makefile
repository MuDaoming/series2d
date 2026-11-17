CXX ?= g++
CPPFLAGS += -I/usr/local/include
CXXFLAGS ?= -O3 -std=c++17
LDFLAGS += -L/usr/local/lib -Wl,-rpath,/usr/local/lib
LDLIBS += -lflint -lgmp
build/experiment: src/experiment.cpp $(wildcard ../../src/IBP/kernel/*) ../../src/IBP/main.cpp
	mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@.tmp
	mv $@.tmp $@
