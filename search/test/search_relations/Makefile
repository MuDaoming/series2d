CXX := g++
CXXFLAGS := -std=c++17 -O2 -I../../include -I../../../expand/include -I../../../reconstruct/include
LDFLAGS := -lginac -lcln -lflint -lgmp

TARGET := search_test_runner
SRC := test_search.cpp

.PHONY: all run clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)
