F95 = mpif90

OPTS = -O3 -fcheck=all -fbacktrace -Wall

LIBS=-L${BLASDIR} ${BLASLIB} -lpthread

OBJS = header.o coefficients.o create_matrices.o vecdot.o vecaxpy.o matmult.o bicgstab.o invit.o neutron.o save_fields.o sparsegather.o

all: neutron

header.mod: header.f90
	$(F95) $(OPTS) -c header.f90

neutron: $(OBJS)
	$(F95) $(OPTS) $(OBJS) -o neutron $(LIBS)

%.o: %.f90 header.mod
	$(F95) $(OPTS) -c $<

clean:
	rm -rf neutron *.o *.mod core solution.*

