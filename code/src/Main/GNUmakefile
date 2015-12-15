HOME = ..
VPATH = . $(HOME)/Dendritic $(HOME)/RectArray $(HOME)/Timer $(HOME)/Writers
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

%.o: %.cpp GNUmakefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

main: GNUmakefile $(OBJS)
	$(CXX) $(CFLAGS) $(OBJS) -o main.exe

clean:
	rm *.o *.exe *.d *.vtk
	
realclean:
	rm *.o *.exe *.d *.vtk $(OBJS) ../Dendritic/*.d ../Writers/*.d ../RectArray/*.d
	
-include $(OBJS:.o=.d)
