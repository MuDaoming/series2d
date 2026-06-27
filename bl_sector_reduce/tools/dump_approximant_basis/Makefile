CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra

ROOT := ../..
EXPAND := ../../../expand

INCLUDES := -I$(ROOT)/include -I$(EXPAND)/include
LIBS := -lflint -lgmp

TARGET := dump_approximant_basis
SRC := dump_approximant_basis.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LIBS)

clean:
	rm -f $(TARGET) *.o *.d

.PHONY: all clean
