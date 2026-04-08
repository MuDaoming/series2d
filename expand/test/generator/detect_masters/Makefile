CXX := g++
CXXFLAGS := -std=c++17 -O2 -I../../../../include
LDFLAGS := -lginac -lcln -lflint -lgmp

TARGET := detect_masters
SRC := *.cpp

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

run: $(TARGET)
	./$(TARGET) S config masters
clean:
	rm -f $(TARGET)
