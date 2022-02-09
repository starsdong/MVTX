OBJS = StPhysicalHelix.o StHelix.o $1.o
EXE = $1

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

INCFLAGS = -I$(ROOTSYS)/include
LDFLAGS = -L$(ROOTSYS)/lib

CXX = g++
FLAGS = -Wall -g $(INCFLAGS) $(LDFLAGS)

COMPILE = $(CXX) $(FLAGS) -c 

all: $(EXE)

$(EXE): $(OBJS)
	$(CXX) -o $(EXE) $(OBJS) $(ROOTFLAGS) $(ROOTLIBS)

%.o: %.cxx
	$(COMPILE)  $< 

