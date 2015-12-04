HOME = ../..
VPATH = . $(HOME)/Dendritic $(HOME)/RectArray $(HOME)/Timer $(HOME)/Writers
DIM = 2
system := $(shell uname)
CXX=g++
#CXX=clang++
CFLAGS = -O3
CFLAGS += -std=c++11
ifeq ($(CXX) , clang++)
  CFLAGS += -stdlib=libc++
endif

