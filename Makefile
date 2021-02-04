.SUFFIXES: .f90 .o

# compiler
#FC = ifort
FC = gfortran
# compile flags
#FFLAGS =  -g -check all -zmuldefs -i8
FFLAGS = -fcheck=all
# librariries
#MKLROOT = /common/intel/composer_xe_2013.2.146/mkl
#LIBS = ${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a \
#      ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a -Wl,--start-group \
#      ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a \
#      ${MKLROOT}/lib/intel64/libmkl_sequential.a \
#      ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
LIBS = 
#/usr/local/lib/libblas.a /usr/local/lib/liblapack.a \
#/usr/local/lib/libfftpack.a
# main
MAIN = soget.o
# modules
MODS = global.o balda.o dimer.o
# objects files
OBJS = secant.o quadpack_double.o dxintany.o
# program name
EXEC = soget.x

all: ${MODS} ${OBJS} ${MAIN} Makefile
	${FC} ${FFLAGS} -o ${EXEC} ${MAIN} ${MODS} ${OBJS} ${LIBS}

.f90.o: ${MODS} ${OBJS} Makefile
	${FC} ${FFLAGS} -c $<

clean:
	rm -f $(EXEC) *.o *.mod
