HOME = ../..
VPATH = . $(HOME)/Dendritic $(HOME)/RectArray $(HOME)/Timer
DIM = 2
system := $(shell uname)
CXX=g++
#CXX=clang++
CFLAGS = -O3
CFLAGS += -std=c++11
ifeq ($(CXX) , clang++)
	CFLAGS += -stdlib=libc++
endif

CFLAGS += -I. -I$(HOME)/Dendritic -I$(HOME)/RectArray -I$(HOME)/Timer -DDIM=$(DIM)
LIBS:= 

SRCFILES:= $(wildcard ./*.cpp $(HOME)/RectArray/*.cpp $(HOME)/Timer/*.cpp)
SRCFILES += $(HOME)/Dendritic/RectMDOperators.cpp
OBJS:= $(patsubst %.cpp, %.o, $(SRCFILES))

%.o: %.cpp GNUmakefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

testOperators: GNUmakefile $(OBJS)
	$(CXX) $(CFLAGS) $(OBJS) -o testOperators.exe

clean:
	rm *.o *.exe *.d
	
realclean:
	rm *.o *.exe *.d $(OBJS)
	
-include $(OBJS:.o=.d)