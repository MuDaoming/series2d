CXX      := g++
CXXFLAGS := -std=c++17 -O2 -MMD -MP -I../../include -I../../../expand/include -I../../../reconstruct/include
LDFLAGS  := -lflint -lgmp

TARGETS := search_de
SRC     := search_de.cpp
OBJ    := $(SRC:.cpp=.o)
DEPS   := $(OBJ:.o=.d)

.PHONY: all clean

all: $(TARGETS)

search_de: search_de.o
	$(CXX) $< -o $@ $(LDFLAGS)

-include $(DEPS)

clean:
	rm -f $(TARGETS) $(OBJ) $(DEPS)
