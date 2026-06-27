CXX      := g++
CXXFLAGS := -std=c++17 -O2 -MMD -MP -I../../include -I../../../expand/include -I../../../reconstruct/include
LDFLAGS  := -lflint -lgmp

TARGET := absorbing_integral_solver
SRC    := absorbing_integral_solver.cpp
OBJ    := $(SRC:.cpp=.o)
DEPS   := $(OBJ:.o=.d)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $@ $(LDFLAGS)

-include $(DEPS)

clean:
	rm -f $(TARGET) $(OBJ) $(DEPS)
