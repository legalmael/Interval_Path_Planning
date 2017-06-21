CXX= g++
EXEC= Cameleon

CXXFLAGS := $(shell pkg-config --cflags ibex) -std=c++11 -Wall 
LIBS	 := $(shell pkg-config --libs  ibex)
#LIBDIR	 := $(shell pkg-config --libdir  ibex)

SRC= main.cpp box.cpp interval.cpp vibes.cpp node.cpp
OBJ= $(SRC:.cpp=.o)

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(LIBS)

.PHONY: clean mrproper

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)


	