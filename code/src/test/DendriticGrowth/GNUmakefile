HOME = ../..
VPATH = . $(HOME)/Dendritic $(HOME)/RectArray $(HOME)/Timer
DIM = 2
system := $(shell uname)
CXX=g++
#CXX=clang++
#CFLAGS = -O3
CFLAGS = -g -Wall
CFLAGS += -std=c++11
ifeq ($(system),Darwin)
  CXX = clang++
endif
ifeq ($(CXX) , clang++)
	CFLAGS += -stdlib=libc++ -mmacosx-version-min=10.10.5
endif

CFLAGS += -I. -I$(HOME)/Dendritic -I$(HOME)/RectArray -I$(HOME)/Timer -I$(HOME)/Writers -DDIM=$(DIM)
LIBS:= 

SRCFILES:= $(wildcard ./*.cpp $(HOME)/Dendritic/*.cpp $(HOME)/RectArray/*.cpp $(HOME)/Timer/*.cpp $(HOME)/Writers/*.cpp)
OBJS:= $(patsubst %.cpp, %.o, $(SRCFILES))
DEPENDENTFILE:= $(patsubst %.cpp, %.d, $(SRCFILES))

%.o: %.cpp GNUmakefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

testDendritic: GNUmakefile $(OBJS)
	$(CXX) $(CFLAGS) $(OBJS) -o testDendritic.exe

clean:
	rm *.o *.exe *.d
	
realclean:
	rm *.o *.exe *.d $(OBJS)
	
-include $(OBJS:.o=.d)
