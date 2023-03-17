# Fortran Compiler options
CF77=ifort
CF90=ifort
F77FLAGS= -O2 -xHost -qopenmp -mkl=parallel# -03 -mkl=parallel#

#######################################################
#		 Sequential
#######################################################
# Add -check all for run-time debug checks, disables optimization, also -traceback, -g for debugging , -implicitnone -u ? -m64 / -check all
F90FLAGS = ${F77FLAGS} -w -m64 -fpp -fcheck=all -fbacktrace -Wall -wextra -g -traceback -gen-interfaces -warn interfaces  -I${SOLIDroot} 
#F90FLAGS = ${F77FLAGS} ${ADDINCFLAGS}

# Link to modules
SOLIDroot = Modules21

# Name of your executable
EXE = run_4

# Here is the name of all object files corresponding to the source
# code that you wrote in order to define the problem statement
#OBJS = opt_module.o  run_opt.o # = not needed? 
OBJS = run_4.o 

# Dependencies
DEPS = ${SOLIDroot}/matrix_util.o ${SOLIDroot}/fem_util.o	${SOLIDroot}/memory_util.o \
	${SOLIDroot}/sparse_util.o \
	${SOLIDroot}/mater_large.o \
	${SOLIDroot}/elem_large_cont_2d6.o ${SOLIDroot}/mater_hyperel.o ${SOLIDroot}/elem_large_cont_3d.o ${SOLIDroot}/elem_flow_3d.o	\
	${SOLIDroot}/mater_small.o	${SOLIDroot}/fem_element.o	\
	${SOLIDroot}/fem_system.o ${SOLIDroot}/matlab_util.o  \
	${SOLIDroot}/abaqus_util.o \
	${SOLIDroot}/Spl_geom.o ${SOLIDroot}/shape_sens.o ${SOLIDroot}/plani4axi.o ${SOLIDroot}/assembleG.o ksmaxim07.o ksmaxsu07.o ksmmasu07.o

# Additional flags for compilation (e.g., include flags)
ADDINCFLAGS = -I${SOLIDroot} -I$(MKLroot)/../.. 

#######################################################
#	        Intel MKL
#######################################################
MKLroot= /opt/intel/oneapi/mkl/2023.0.0/lib/intel64

#######################################################
#		 Sequential
#######################################################
ADDLIBS =  -L/software/matlab/R2019b/bin/glnxa64 -leng -lmat -lmx -lut -Wl,--start-group ${MKLroot}/libmkl_intel_lp64.a ${MKLroot}/libmkl_sequential.a ${MKLroot}/libmkl_core.a -Wl,--end-group -lpthread

#######################################################
#	 MATLAB
#######################################################
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/software/matlab/R2019b/bin/glnxa64:/opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/software/matlab/R2019b/bin/glnxa64:/opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin

all: subsystem $(EXE) 

.SUFFIXES: .f .f90 .mod .o .x .obj

subsystem:
	cd ${SOLIDroot} && $(MAKE)

$(EXE): $(OBJS) $(DEPS)
	$(CF90) ${F90FLAGS} $(F77LINKFLAGS) -o $@ $(OBJS) $(DEPS) $(ADDLIBS)

%.o : %.mod

.f90.o:
	${CF90} ${F90FLAGS} $(INCL) -c -o $@ $< 

.f.o:
	${CF77} ${F77FLAGS} $(INCL) -c -o $@ $< 
	
.f.obj:
	$(CF77) $(F77FLAGS) $(INCL) -c -o $@ $<
	
.f90.mod:
	${CF90} ${F90FLAGS} $(INCL) -c -o $@ $<

clean:
	rm -f $(EXE) $(OBJS) IPOPT.OUT *.mod