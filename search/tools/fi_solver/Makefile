CXX      := g++
CXXFLAGS := -std=c++17 -O2 -I../../include -I../../../expand/include -I../../../reconstruct/include
LDFLAGS  := -lflint -lgmp

TARGET := fi_solver
SRC    := fi_solver.cpp

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm -f $(TARGET)
