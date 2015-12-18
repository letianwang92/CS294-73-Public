HEADS = ../include
SOURCE = ../src
VPATH = . $(HEADS) $(SOURCE)/Main $(SOURCE)/Dendritic $(SOURCE)/RectArray $(SOURCE)/Timer $(SOURCE)/Writers
DIM = 2
system := $(shell uname)
CXX=g++
#CXX=clang++
CFLAGS = -O3
#CFLAGS = -g -Wall
CFLAGS += -std=c++11
ifeq ($(system),Darwin)
  CXX = clang++
endif
ifeq ($(CXX) , clang++)
	CFLAGS += -stdlib=libc++ 
endif

CFLAGS += -I. -I$(HEADS) -I$(SOURCE)/Main -I$(SOURCE)/Dendritic -I$(SOURCE)/RectArray -I$(SOURCE)/Timer -I$(SOURCE)/Writers -DDIM=$(DIM)
LIBS:= 

SRCFILES:= $(wildcard ./*.cpp $(SOURCE)/Main/*.cpp $(SOURCE)/Dendritic/*.cpp $(SOURCE)/RectArray/*.cpp $(SOURCE)/Timer/*.cpp $(SOURCE)/Writers/*.cpp)
OBJS:= $(patsubst %.cpp, %.o, $(SRCFILES))
DEPENDENTFILE:= $(patsubst %.cpp, %.d, $(SRCFILES))

%.o: %.cpp GNUmakefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

all: GNUmakefile $(OBJS)
	$(CXX) $(CFLAGS) $(OBJS) -o Dendritic2D.exe
	cd ../../documents/doxygenDocument/; sh link.sh

clean:
	rm *.o *.exe *.d *.vtk
	
realclean:
	rm -f *.o *.exe *.d *.vtk $(OBJS) $(DEPENDENTFILE) 
	cd ../../documents/doxygenDocument/; make clean

-include $(OBJS:.o=.d)
