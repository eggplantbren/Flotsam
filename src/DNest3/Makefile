# Compiler settings
CPP = g++
CFLAGS = -O2 -Wall -Wextra -ansi -pedantic -DNDEBUG -DGSL_RANGE_CHECK_OFF -DHAVE_INLINE
CLIBS = -ldnest3 -lgsl -lgslcblas -lboost_thread -lboost_system
OBJECTS = Data.o Limits.o Matrix.o TDModel.o Vector.o

default:
	g++ $(CFLAGS) -c *.cpp
	g++ -o Flotsam main.o $(OBJECTS) $(CLIBS)
	g++ -o mean_covariance mean_covariance.o $(OBJECTS) $(CLIBS)
	rm -f *.o

