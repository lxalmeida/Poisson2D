makefile:
all: poisson

.PHONY: clean allclean sanit

include ../Makefile.common

poisson:
	cd .. ; make all

clean:
	-rm -f poisson2d poisson_* sent* recv* log* output_topology*

sanit:
	-rm -f poisson_* sent* recv* log* output_topology*

allclean:
	cd .. ; make clean

run:
	-rm -f sent*
	mpirun -np 4 ./poisson2d 200 200 1000000 5
