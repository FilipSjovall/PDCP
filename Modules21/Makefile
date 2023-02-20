# 
# Only compile the modulus 
#
 
CF77=ifort
CF90=ifort
#DEBUGFLAGS= -g -ftrapuv -check all -check bounds -check arg_temp_created -traceback -gen-interfaces -warn interfaces -check uninit -fpe0
DEBUGFLAGS= -g -traceback -check -gen-interfaces -warn interfaces #-fpe0 
#F77FLAGS= -O0 -g -traceback -fpe0 #-check -fp-model strict #-O0 -xHost ${DEBUGFLAGS}
F77FLAGS= -O3 -xHost -heap-arrays
#F77FLAGS= -O0 ${DEBUGFLAGS}
 
F90FLAGS = ${F77FLAGS} 


# modules 
# note order is important if module used in other module
SolidLib=mesh_reader.o matrix_util.o 		\
	sparse_util.o fem_util.o	memory_util.o fem_system.o	\
	matlab_util.o mater_large.o\
	elem_large_cont_2d6.o mater_hyperel.o elem_large_cont_3d.o elem_flow_3d.o	\
	mater_small.o	fem_element.o	\
	abaqus_util.o \
    ksmaxim07.o ksmaxsu07.o ksmmasu07.o 
        
         

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .x

# making objects
.f90.o:
	${CF90} ${F90FLAGS} -c  $*.f90 

.f.o:
	${CF77} ${F77FLAGS} -c $*.f 

# making module
.f90.mod:
	${CF90} ${F90FLAGS} -c $*.f90 

default: solid_mod.o

solid_mod.o: ${SolidLib}

clean:
	/bin/rm -f *.x *.o *.mod
