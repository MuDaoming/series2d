CXX      := g++
CXXFLAGS := -std=c++17 -O2 -MMD -MP -I../../include -I../../../expand/include -I../../../reconstruct/include
LDFLAGS  := -lflint -lgmp

TARGET := absorbing_poly_searcher
SRC    := absorbing_poly_searcher.cpp
OBJ    := $(SRC:.cpp=.o)
DEPS   := $(OBJ:.o=.d)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $@ $(LDFLAGS)

-include $(DEPS)

clean:
	rm -f $(TARGET) $(OBJ) $(DEPS)
