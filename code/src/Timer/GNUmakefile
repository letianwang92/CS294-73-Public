#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

#CFLAGS = -g -Wall
CFLAGS = -O3 
CXX = g++

CH_Timer.o: CH_Timer.H CH_Timer.cpp ClockTicks.H GNUmakefile
	$(CXX) $(CFLAGS) -c -o CH_Timer.o CH_Timer.cpp
