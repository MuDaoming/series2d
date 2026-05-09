CXX      := g++
CXXFLAGS := -std=c++17 -O2 -MMD -MP -I../../include -I../../../expand/include -I../../../reconstruct/include
LDFLAGS  := -lflint -lgmp

TARGET := poly_relation_searcher
SRC    := poly_relation_searcher.cpp
OBJ    := $(SRC:.cpp=.o)
DEPS   := $(OBJ:.o=.d)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $@ $(LDFLAGS)

-include $(DEPS)

clean:
	rm -f $(TARGET) $(OBJ) $(DEPS)
