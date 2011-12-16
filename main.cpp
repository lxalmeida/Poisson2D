#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include "jacobi2d.h"

using namespace std;

fstream logfile;

void print_usage(void);

int main(int argc, char **argv){
    int m, n, max_iter, iter = 0, ngz = 1, numthreads = 1;
    double error;
    //double t_init, t_end;
    int rank = 0;
    
    MPI_Init(&argc, &argv);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(argc < 4) {
		if(rank == 0) {
    		print_usage();
    	}
		MPI_Finalize();
		return(255);
    }
	
    m = atoi(argv[1]);
    n = atoi(argv[2]);
    max_iter = atoi(argv[3]);

    if(argc > 4) {
    	numthreads = atoi(argv[4]);
    	if(argc == 6) {
    		ngz = atoi(argv[5]);
    	}
	}
    
    //ostringstream filename(ostringstream::out);
    
    //filename << "rank" << rank << "_ngz" << ngz << ".txt";
    
    //logfile.open(filename.str().c_str(), fstream::out);
    
    Mesh mesh(m, n, ngz, numthreads);
    
    // Inicializa mesh
    for(int i = 0; i < mesh.get_size_y(); i++){
        for(int j = 0; j < mesh.get_size_x(); j++){
            mesh.data_dest[i][j] = 1.0;
        }
    }
    // Inicializa heat zone a esquerda
    mesh.set_left_extern(15.0);
    
    // Swap nos meshes para inicializar o outro
    mesh.swap();
    // Inicializa mesh
    for(int i = 0; i < mesh.get_size_y(); i++){
        for(int j = 0; j < mesh.get_size_x(); j++){
            mesh.data_dest[i][j] = 1.0;
        }
    }
    // Inicializa heat zone a esquerda
    mesh.set_left_extern(15.0);
    
    mesh.swap();
    
    // Sobrepõe os dados dos subdomínios
    mesh.send_borders();
    mesh.sync();
    
    //mesh.print_ghost_zones();
    
    //t_init = MPI_Wtime();
    
    iter = poisson2d(mesh, max_iter, m, &error);

    mesh.gather();
    
    //t_end = MPI_Wtime();
    
    if(mesh.topology->get_rank() == 0){
        //cout << fixed << setprecision(6);
        //cout << t_end - t_init;

        /*if(iter > 0) {
            cout << "\t(" << iter << ")" << endl;
        } else {
            cout << "\t(-1)" << endl;
        }*/
    }
    
    //mesh.print_final_result();
    
    MPI_Finalize();
    
    return(0);
}

void print_usage(void) {
	cerr << "Usage: ./poisson2d m n max_iter [num_threads] [num_ghost_zones]" << endl;
	cerr << "Mandatory arguments:" << endl;
	cerr << "            m = total number of rows" << endl;
	cerr << "            n = total number of columns" << endl;
	cerr << "     max_iter = max number of iterations" << endl;
	cerr << "Optional arguments:" << endl;
	cerr << "         num_threads = number of threads (default = 1)" << endl;
	cerr << "     num_ghost_zones = number of subdomain ghost zones (default = 1)" << endl;
}
