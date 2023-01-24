SHELL = /bin/sh

# make clean
# make stats_tests
# make all

OBJS = stats_tests.o ../Gofr.o ../SimSpace.o ../blocking.o ../critexp.o ../pca.o
CFLAGS =
CC = clang++
INCLUDES = 
LIBS = -L/usr/local/lib -lgsl -lgslcblas

all:stats_tests

stats_tests:${OBJS}
	${CC} ${CFLAGS} ${INCLUDES} -o $@ ${OBJS} ${LIBS}

clean:
	-rm -f *.o core *.core *.dat

.cpp.o:
	${CC} ${CFLAGS} ${INCLUDES} -c $<


