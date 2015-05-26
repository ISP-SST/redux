CC = mpic++

GENLIB_PATH = /home/woeger/bin/lib

INCLUDE_PATH = -I${GENLIB_PATH}/gsl/include -I${GENLIB_PATH}/mpi/include

LIB_PATH = -L${GENLIB_PATH}/gsl/lib -L${GENLIB_PATH}/mpi/lib

SYS_CFLAGS = -O2 -std=c++11

#LIB_FLAGS = -lmpi -llam -lnsl -lrt -lgsl -lgslcblas -lc -lm 
LIB_FLAGS = -lgsl -lgslcblas

MPI_CFLAGS = -DUSE_STDARG -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1\
             -DHAVE_UNISTD_H=1 -DHAVE_STDARG_H=1 -DUSE_STDARG=1\
             -DMALLOC_RET_VOID=1

stf_new:  stf.o zernike.o wam.o
	${CC} ${MPI_CFLAGS} ${SYS_CFLAGS} ${INCLUDE_PATH}\
        stf.o zernike.o wam.o\
        ${LIB_PATH} ${LIB_FLAGS} -o $@

stf.o: stf.cpp
	${CC} -c ${MPI_CFLAGS} ${SYS_CFLAGS} ${INCLUDE_PATH} $*.cpp

zernike.o: zernike.cpp
	${CC} -c ${MPI_CFLAGS} ${SYS_CFLAGS} ${INCLUDE_PATH} $*.cpp

wam.o: wam.cpp
	${CC} -c ${MPI_CFLAGS} ${SYS_CFLAGS} ${INCLUDE_PATH} $*.cpp

clean:
	rm *.o stf
