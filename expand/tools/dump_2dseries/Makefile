CXX := g++
CXXFLAGS := -std=c++17 -O2 -MMD -MP -I../../include
LDFLAGS := -lginac -lcln -lflint -lgmp -lpthread

TARGET := dump_2dseries
SRC := main.cpp
OBJ := $(SRC:.cpp=.o)
DEPS := $(OBJ:.o=.d)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $@ $(LDFLAGS)

-include $(DEPS)

clean:
	rm -f $(TARGET) $(OBJ) $(DEPS)
