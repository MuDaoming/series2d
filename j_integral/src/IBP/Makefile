CXX ?= g++
FLINT_PREFIX ?= /usr/local
CXXFLAGS ?= -O3 -std=c++17
CPPFLAGS += -I$(FLINT_PREFIX)/include
LDFLAGS += -L$(FLINT_PREFIX)/lib -Wl,-rpath,$(FLINT_PREFIX)/lib
LDLIBS += -lflint -lgmp
TARGET := build/j-ibp
SOURCES := main.cpp $(wildcard kernel/*)
.PHONY: all test clean
all: $(TARGET)
$(TARGET): $(SOURCES)
	mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) main.cpp $(LDFLAGS) $(LDLIBS) -o $@.tmp
	mv $@.tmp $@
test: $(TARGET)
	python3 tests/test_interface.py
clean:
	rm -rf build
