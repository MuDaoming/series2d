CXX := g++
CXXFLAGS := -std=c++17 -O2 -I../../include
LDFLAGS := -lginac -lcln -lflint -lgmp -lpthread

TARGET := dump_2dseries
SRC := main.cpp

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm -f $(TARGET)
