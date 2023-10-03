F95       = mpif90 
LIB_LIST  = -L${BLASDIR} ${BLASLIB} -lpthread

OPT      = -O3 -Wall
OBJS      = vecdot.o header.o testq1.o vecaxpy.o

all: testq1

header.mod: header.f90
	$(F95) $(OPT) -c header.f90

%.o: %.f90 header.mod
	$(F95) $(OPT) -c $<

testq1: $(OBJS)
	 $(F95) $(OPT) -o testq1 $(OBJS) $(LIB_LIST)

clean: 
	rm -f *.o *.mod core testq1
