# compiler
CXX := g++
CXXFLAGS := -std=c++17 -O3

# objects
EXEC := main
SRCS := main.cpp
OBJS := $(SRCS:.cpp=.o)

# targets
all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)

.PHONY: all clean
