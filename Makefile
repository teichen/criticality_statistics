SHELL = /bin/sh

# make clean
# make all
# make criticality_stats

OBJS = SimSpace.o Gofr.o pca.o blocking.o
CFLAGS =
CC = clang++
INCLUDES = 
LIBS = omp

all:criticality_stats

criticality_stats:${OBJS}
	${CC} ${CFLAGS} ${INCLUDES} -o $@ ${OBJS} ${LIBS}

clean:
	-rm -f *.o core *.core *.dat

.cpp.o:
	${CC} ${CFLAGS} ${INCLUDES} -c $<


