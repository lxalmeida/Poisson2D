makefile:

.PHONY: all clean

include Makefile.common

all: ${POISSON_PATH}/bin/poisson2d

${POISSON_PATH}/bin/poisson2d: ${POISSON_PATH}/obj/main.o ${POISSON_PATH}/obj/jacobi2d.o 
	cd ${MESH2D_PATH} ; make all
	${MPICPP} ${POISSON_PATH}/obj/main.o ${POISSON_PATH}/obj/jacobi2d.o ${OBJ_PATH}/mesh2d/mesh2d.o ${OBJ_PATH}/topology/topology.o ${OBJ_PATH}/resource/resource.o ${OBJ_PATH}/topology/topology_output.o ${OBJ_PATH}/mesh2d/ghost_zones/ghost_zones.o ${OBJ_PATH}/mesh2d/ghost_zones/halo_left.o ${OBJ_PATH}/mesh2d/ghost_zones/halo_right.o ${OBJ_PATH}/mesh2d/ghost_zones/halo_top.o ${OBJ_PATH}/mesh2d/ghost_zones/halo_bottom.o -o bin/poisson2d  ${LD_FLAGS}

${POISSON_PATH}/obj/main.o: ${POISSON_PATH}/main.cpp
	${MPICPP} ${POISSON_PATH}/main.cpp -c -o obj/main.o -I ${MESH2D_PATH}/include ${CPP_FLAGS}

${POISSON_PATH}/obj/jacobi2d.o: ${POISSON_PATH}/jacobi2d.cpp
	${MPICPP} ${POISSON_PATH}/jacobi2d.cpp -c -o obj/jacobi2d.o -I ${MESH2D_PATH}/include ${CPP_FLAGS}

clean:
	cd obj ; make clean
	cd bin ; make clean

