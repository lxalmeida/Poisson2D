#ifndef JACOBI2D_H
#define JACOBI2D_H

#include <mesh2d.h>
#include <cmath>
#include <pthread.h>

struct args_thread_borders {
	struct args_t *mesh_thread_args;
	int border_index;
};

int poisson2d(Mesh &mesh, int max_iter, int iter_skip, double *error);
void *inner_jacobi_iter(void *args);
void *proc_border(void *args);
void *calc_diff_slice(void *args);
int inner_jacobi_iter(Mesh &mesh);
int outer_jacobi_iter(Mesh &mesh);
double calc_diff(Mesh &mesh);

#endif

