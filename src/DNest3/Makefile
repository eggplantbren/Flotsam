OBJECTS = Data.o Limits.o TDModel.o

# Compiler settings
CPP = g++
# Development
# CFLAGS = -O0 -g -Wall -Wextra -ansi -pedantic
# Optimized
CFLAGS = -O2 -Wall -Wextra -ansi -pedantic -DNDEBUG -DGSL_RANGE_CHECK_OFF -DHAVE_INLINE
CLIBS = -ldnest3 -lgsl -lgslcblas -lboost_thread

# The default build rule
.cpp.o:
	$(CPP) $(CFLAGS) -o $*.o -c $*.cpp

default: main.o $(OBJECTS)
	$(CPP) -o main main.o $(OBJECTS) $(CLIBS)

clean:
	rm -f main main.o $(OBJECTS) *.pyc

