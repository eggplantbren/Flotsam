CFLAGS = -O2 -Wall -Wextra -ansi -pedantic -pthread
LIBS = -lgsl -lgslcblas -ldnest -lboost_thread -lpthread -lboost_program_options

default: main interpolate

TDModel.o: TDModel.cpp TDModel.h Data.h
	g++ $(CFLAGS) -c TDModel.cpp

Data.o: Data.cpp Data.h
	g++ $(CFLAGS) -c Data.cpp

main.o: main.cpp Data.h TDModel.h
	g++ $(CFLAGS) -c main.cpp

interpolate.o: interpolate.cpp Data.h TDModel.h
	g++ $(CFLAGS) -c interpolate.cpp

main: TDModel.o Data.o main.o
	g++ $(LFLAGS) -o main main.o TDModel.o Data.o $(LIBS)

interpolate: TDModel.o Data.o interpolate.o
	g++ $(LFLAGS) -o interpolate interpolate.o TDModel.o Data.o $(LIBS)

clean:
	rm -f *.o


