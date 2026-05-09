CXX := g++
CXXFLAGS := -std=c++17 -O2 -MMD -MP -I../../include
LDFLAGS := -lginac -lcln -lflint -lgmp -lpthread

TARGET := dump_masters
SRC := main.cpp
OBJ := $(SRC:.cpp=.o)
DEPS := $(OBJ:.o=.d)

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $@ $(LDFLAGS)

-include $(DEPS)

run: $(TARGET)
	./$(TARGET) S config masters

clean:
	rm -f $(TARGET) $(OBJ) $(DEPS)
