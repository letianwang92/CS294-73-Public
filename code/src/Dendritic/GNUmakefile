# our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.
# TODO: To be modified.
HOME = ..
VPATH= . $(HOME)/Writers $(HOME)/../RectArray $(HOME)/fftTools $(HOME)/../../resources/timer $(HOME)/Hockney
DIM=2
system := $(shell uname)
CXX=g++
#CXX=clang++
FFTWDIR = /usr/local
FFTINCFLAGS = -I$(FFTWDIR)/include
#CFLAGS = -g -Wall
CFLAGS = -O3
CFLAGS += -std=c++11
ifeq ($(CXX) , clang++)
  CFLAGS += -stdlib=libc++
endif
osuffix:=$(DIM)D.o
dsuffix:=$(DIM)D.d

CFLAGS += -I. -I$(HOME)/Writers -I$(HOME)/RectArray -I$(HOME)/fftTools -I$(HOME)/Hockney -DDIM=$(DIM) $(FFTINCFLAGS) -I $(HOME)/../../resources/timer -D DIM=$(DIM)
LIBS:=  -L$(FFTWDIR)/lib -lfftw3 -L../lib -lfft$(DIM)D -lfft1D

SRCFILES:= $(wildcard $(HOME)/../../resources/timer/*.cpp $(HOME)/RectArray/*.cpp $(HOME)/Writers/*.cpp ./*.cpp)
SRCFILES += $(HOME)/Hockney/Hockney.cpp
OBJS:=$(patsubst %.cpp,%.o, $(SRCFILES)) 

%.o: %.cpp GNUmakefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

%$(osuffix): %.cpp GNUmakefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*$(dsuffix)

../lib/libfft$(DIM)D.a:$(wildcard $(HOME)/fftTools/*.H $(HOME)/fftTools/*.cpp)
	cd $(HOME)/fftTools;rm *.o;make libfft$(DIM)D.a CXX=$(CXX) DIM=$(DIM) FFTWDIR=$(FFTWDIR) FFTINCFLAGS=$(FFTINCFLAGS)

../lib/libfft1D.a:$(wildcard $(HOME)/fftTools/*.H $(HOME)/fftTools/*.cpp)
	cd $(HOME)/fftTools;rm *.o;make libfft1D.a CXX=$(CXX) DIM=$(DIM) FFTWDIR=$(FFTWDIR) FFTINCFLAGS=$(FFTINCFLAGS)

vars:
	echo $(wildcard $(HOME)/fftTools/*.H $(HOME)/fftTools/*.cpp)
objs:
	echo $(OBJS)

vortex2D: GNUmakefile ../lib/libfft$(DIM)D.a ../lib/libfft1D.a $(OBJS) 
	$(CXX) $(CFLAGS) $(FFTINCFLAGS) $(OBJS) $(LIBS) -o vortex$(DIM)D.exe

clean:
	rm  *.o *.exe *.d

realclean:
	rm  *.exe *.d *.o $(OBJS);cd ../fftTools;make clean

-include $(OBJS:.o=.d)
