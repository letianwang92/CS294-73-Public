HEADS = ../../../include
SOURCE = ../..
VPATH = . $(HEADS) $(SOURCE)/Dendritic $(SOURCE)/RectArray $(SOURCE)/Timer $(SOURCE)/Writers
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

CFLAGS += -I. -I$(HEADS) -I$(SOURCE)/Dendritic -I$(SOURCE)/RectArray -I$(SOURCE)/Timer -I$(SOURCE)/Writers -DDIM=$(DIM)
LIBS:= 

SRCFILES:= $(wildcard ./*.cpp $(SOURCE)/Dendritic/*.cpp $(SOURCE)/RectArray/*.cpp $(SOURCE)/Timer/*.cpp $(SOURCE)/Writers/*.cpp)
OBJS:= $(patsubst %.cpp, %.o, $(SRCFILES))
DEPENDENTFILE:= $(patsubst %.cpp, %.d, $(SRCFILES))

%.o: %.cpp GNUmakefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

all: GNUmakefile $(OBJS)
	$(CXX) $(CFLAGS) $(OBJS) -o main.exe

clean:
	rm *.o *.exe *.d *.vtk
	
realclean:
	rm *.o *.exe *.d *.vtk $(OBJS) $(DEPENDENTFILE) 
	
-include $(OBJS:.o=.d)
