CC = g++
CFLAGS = -Wall -g
FLAGS = -std=c++11 
LDFLAGS = -L/usr/lib -lboost_program_options -lpython2.7 -lfftw3 -lm -I/usr/include/python2.7 -I/home/hermanni/windowing/eigen-eigen-323c052e1731/

main : main.cpp resample.cpp resample.h upfirdn.h
	${CC} main.cpp resample.cpp  ${CFLAGS} ${LDFLAGS} -o main

.PHONY: clean

clean:
	rm -f *.o main



#sudo g++ main.cpp -std=c++11 -I/usr/include/python2.7 -lpython2.7 -lfftw3 -lm -Wall -g -L/usr/lib -lboost_program_options
#g++ example.cpp resample.cpp  -Wall -g -L/usr/lib -lboost_program_options -o example

