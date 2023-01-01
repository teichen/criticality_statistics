SHELL = /bin/sh

# make clean

OBJS = 
CFLAGS =
CC = clang++
INCLUDES = 
LIBS =

clean:
	-rm -f *.o core *.core *.dat

.cpp.o:
	${CC} ${CFLAGS} ${INCLUDES} -c $<


