#dest = m1
# Compiler and Flags
FC = mpif90
#FC = ftn 
FLAGS = -O3 -ffast-math -funroll-loops -Wall -Wtabs

# user customized source subroutine in first line
OBJS = stringmod.o globals.o util.o radiations.o\
	frio.o  mpi.o arrays.o setup.o fieldio.o parameters.o  assemble.o main.o

all: PIOL-FFD

PIOL-FFD: $(OBJS)
	$(FC) $(FLAGS) $^ -o $@ 
	cp $@ ../

# independent
util.o:         util.f90
	$(FC) $(FLAGS) -c $^

frio.o:         frio.f90
	$(FC) $(FLAGS) -c $^


globals.o:      globals.f90
	$(FC) $(FLAGS) -c $^ 

#********** user customized ***********
stringmod.o: stringmod.f90
	$(FC) $(FFLAGS) -c $^
	
radiations.o: radiations.f90
	$(FC) $(FFLAGS) -c $^

assemble.o: assemble.f90
	$(FC) $(FFLAGS) -c $^
#**************************************

# dependent on m_frio
mpi.o:          mpi.f90
	$(FC) $(FLAGS) -c $^

serial.o:       serial.f90
	$(FC) $(FLAGS) -c $^
	
# dependent on globals.o
arrays.o:       arrays.f90
	$(FC) $(FLAGS) -c $^ 

# dependent on m_globals, m_collective (mpi.o), m_util
setup.o:        setup.f90
	$(FC) $(FLAGS) -c $^


# dependent on m_globals, m_collective (mpi.o), m_util, m_frio
fieldio.o:      fieldio.f90
	$(FC) $(FLAGS) -c $^

# dependent on m_globals, m_fieldio
parameters.o:   parameters.f90
	$(FC) $(FLAGS) -c $^


# depdentent on all above
main.o:         main.f90
	$(FC) $(FLAGS) -c $^


# PHONY TARGETS 
clean:
	rm -rf $(OBJS) *.mod $(bin)
