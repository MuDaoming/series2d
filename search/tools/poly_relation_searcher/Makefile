CXX      := g++
CXXFLAGS := -std=c++17 -O2 -I../../include -I../../../expand/include -I../../../reconstruct/include
LDFLAGS  := -lflint -lgmp

TARGET := poly_relation_searcher
SRC    := poly_relation_searcher.cpp

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm -f $(TARGET)
