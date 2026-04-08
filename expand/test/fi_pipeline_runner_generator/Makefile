CXX := g++
CXXFLAGS := -std=c++17 -O2 -I../../include
LDFLAGS := -lginac -lcln -lflint -lgmp

TARGET := fi_pipeline_runner
SRC := fi_pipeline_runner_generator.cpp

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

run: $(TARGET)
	./$(TARGET) S config target series

clean:
	rm -f $(TARGET)
