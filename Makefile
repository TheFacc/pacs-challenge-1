# Compiler settings
CXX = g++
CXXFLAGS = -std=c++17
LDFLAGS = -L/nix/store/46s43nm5gfmhz22q8ksbv20x19r0h0wb-muparser-2.2.3/lib -lmuparser

# objects
OBJS = main.o
EXEC = main

# targets
all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

main.o: main.cpp
	$(CXX) -c $(CXXFLAGS) $<

clean:
	rm -f $(OBJS) $(EXEC)

.PHONY: all clean
