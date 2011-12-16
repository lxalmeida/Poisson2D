#include "jacobi2d.h"

int poisson2d(Mesh &mesh, int max_iter, int iter_skip, double *error){
    int k = 1;
    double diff = 0.0, diff_norm = 0.0;
    //extern fstream logfile;
    
    double t_ini = 0, t_fim = 0;
    //double t_prev = 100, t_curr = 0, mean;
    //int threads = 2, nsamples = 100, points = 0;
    
    if(max_iter <= 0){
        return(-255);
    }
    
    //mesh.print_mesh_info();
    //mesh.print_file_mesh();
    
    //t_ini = MPI_Wtime();
    
    /* Calculo inicial */
    inner_jacobi_iter(mesh);
    outer_jacobi_iter(mesh);
    /* Dispara a comunicacao das fronteiras */
    mesh.send_borders();
    /* Calcula os pontos do subdominio, menos as fronteiras */
    inner_jacobi_iter(mesh);
    /* Bloqueia enquanto os Isends/Irecvs nao terminarem */
    mesh.sync();
    /* Calcula os pontos das fronteiras */
    outer_jacobi_iter(mesh);
    
    //t_fim = MPI_Wtime();
    
	//t_curr = (t_fim - t_ini);
	
    /* Determina a diferenca entre as iteracoes */
    diff = calc_diff(mesh);
    /* A diferenca local de cada processo eh somada a um valor global */
    MPI_Allreduce(&diff, &diff_norm, 1, MPI_DOUBLE, MPI_SUM, mesh.topology->get_comm());
    
    //if(mesh.topology->get_rank() == 0)
    //    cout << diff_norm << endl;
    
    mesh.swap();
    
    cout << setprecision(10) << fixed;
    
    /* Itera ate que convirja, ou ate que o numero maximo de iteracoes seja
       atingido */
    while(k < max_iter && diff_norm > 1.0E-5){
    	t_ini = MPI_Wtime();
    	/* Dispara a comunicacao das fronqteiras. Envia as fronteiras de u. */
        mesh.send_borders();
        /* Calcula os pontos do subdominio, menos as fronteiras */
        inner_jacobi_iter(mesh);
        /* Bloqueia enquanto os Isends/Irecvs nao terminarem */
        mesh.sync();
        /* Calcula os pontos das fronteiras */
        outer_jacobi_iter(mesh);
        t_fim = MPI_Wtime();
        
        cout << t_fim - t_ini << endl;
        
        mesh.swap();
        
        /* Dispara a comunicacao das fronteiras. Envia as fronteiras de v. */
        mesh.send_borders();
        /* Calcula os pontos do subdominio, menos as fronteiras */
        inner_jacobi_iter(mesh);
        /* Bloqueia enquanto os Isends/Irecvs nao terminarem */
        mesh.sync();
        /* Calcula os pontos das fronteiras */
        outer_jacobi_iter(mesh);
        
        /*t_fim = MPI_Wtime();
        logfile << t_fim - t_ini << endl;*/
        
        
        /*t_fim = MPI_Wtime();
        
        t_curr += (t_fim - t_ini);
        
        if(k%nsamples == 0) {
        	mean = t_curr/nsamples;
        	// Se a diferença for maior que 10% do valor anterior
        	if((t_prev*0.1) < fabs(t_prev - mean)) {
        		cout << "t_prev = " << t_prev << "  t_curr/nsamples = " << t_curr/nsamples << endl;
				if(t_prev > mean) {
					threads += 2;
					t_prev = mean;
				} else {
					points--;
					cout << " Perdeu pontos! pontos = " << points << endl;
					if(points == 0) {
						threads -= 2;
					}
				}
				cout << "  num_threads = " << threads << endl;
        	} else {
        		points++;
        	}
			t_curr = 0.0;
        }*/
		
        mesh.swap();
        
        /* Como o algoritmo leva muitas iteracoes ate convergir, iter_skip eh
           um fator utilizado para atualizar a diferenca entre os subdominios
           a cada n iteracoes, evitando o Allreduce */        
        if((k%iter_skip) == 0){
            diff = calc_diff(mesh);
            
            MPI_Allreduce(&diff, &diff_norm, 1, MPI_DOUBLE, MPI_SUM, mesh.topology->get_comm());

            /*if(mesh.topology->get_rank() == 0){
                cout << diff_norm <<  endl;
            }*/
        }
        
        k++;
    }

    if(diff_norm > 1.0E-5){
        return(-1);
    }
    
    *error = diff_norm;
    
    return(k);
}

/* Calcula os pontos do subdominio de u utilizando v, menos as regioes de
   fronteira */
int inner_jacobi_iter(Mesh &mesh){

	/* Complexidade (n-2)*(n-2) = (n^2 - 4n - 4) */
/*	for(int i = 1; i < mesh.get_size_y()-1; i++){
		for(int j = 1; j < mesh.get_size_x()-1; j++){
			mesh.data_dest[i][j] = 0.25 * (mesh.data_source[i-1][j] +
									  mesh.data_source[i][j+1] +
									  mesh.data_source[i][j-1] +
									  mesh.data_source[i+1][j]);
		}
	}*/
	for(int i = 0; i < mesh.get_num_threads(); i++) {
		//cout << "Criando thread " << i << endl;
		pthread_create(&(mesh.thread_id[i]), &(mesh.thread_attr), inner_jacobi_iter, (void*)&(mesh.thread_args[i]));
	}
	
	for(int i = 0; i < mesh.get_num_threads(); i++) {
		//cout << "Esperando thread " << i << endl;
		pthread_join(mesh.thread_id[i], NULL);
		//cout << "  feito!" << endl;
	}
	
    return(0);
}

void *inner_jacobi_iter(void *args) {
	struct args_t *a = (struct args_t *) args;
	double **data_dest = a->mesh->data_dest, **data_source = a->mesh->data_source;
	// Cada thread tem um número de linhas igual a rows/NUM_THREADS e todas as colunas.
	// Essa política de divisão tende a aproveitar a localidade dos dados em memória.
	int lower_col = 1, upper_col = a->mesh->get_size_x()-1;
	
    for(int i = a->lower_row; i <= a->upper_row; i++){
        for(int j = lower_col; j < upper_col; j++){
            data_dest[i][j] = 0.25 * (data_source[i-1][j] +
            						  data_source[i][j+1] +
            						  data_source[i][j-1] +
            						  data_source[i+1][j]);
        }
    }
    
	pthread_exit(NULL);
}


/* Calcula as fronteiras de u utilizando v */
int outer_jacobi_iter(Mesh &mesh){
	/*int i = 0, j = 0;
	
	// Fronteira superior (linha fixa)
	i = 0;
	for(j = 0; j < mesh.get_size_x(); j++){
		mesh.data_dest[i][j] = 0.25 * (mesh.get_source_halo(i-1, j) +
								  mesh.get_source_halo(i, j+1) +
								  mesh.get_source_halo(i, j-1) +
								  mesh.get_source_halo(i+1, j));
	}

	// Fronteira da esquerda (coluna fixa)
	j = 0;    
	for(i = 0; i < mesh.get_size_y(); i++){
		mesh.data_dest[i][j] = 0.25 * (mesh.get_source_halo(i-1, j) +
								  mesh.get_source_halo(i, j+1) +
								  mesh.get_source_halo(i, j-1) +
								  mesh.get_source_halo(i+1, j));
	}
	
	// Fronteira inferior (linha fixa)
	i = mesh.get_size_y() - 1;
	for(j = 0; j < mesh.get_size_x(); j++){
		mesh.data_dest[i][j] = 0.25 * (mesh.get_source_halo(i-1, j) +
								  mesh.get_source_halo(i, j+1) +
								  mesh.get_source_halo(i, j-1) +
								  mesh.get_source_halo(i+1, j));
	}
	
	// Fronteira da direita (coluna fixa)
	j = mesh.get_size_x() - 1;
	for(i = 0; i < mesh.get_size_y(); i++){
		mesh.data_dest[i][j] = 0.25 * (mesh.get_source_halo(i-1, j) +
								  mesh.get_source_halo(i, j+1) +
								  mesh.get_source_halo(i, j-1) +
								  mesh.get_source_halo(i+1, j));
	}*/
	
	for(int i = 0; i < mesh.get_num_threads(); i++) {
		pthread_create(&(mesh.thread_id[i]), &(mesh.thread_attr), proc_border, (void*)&(mesh.thread_args[i]));
	}
	for(int i = 0; i < mesh.get_num_threads(); i++) {
		pthread_join(mesh.thread_id[i], NULL);
	}
    
    return(0);
}

void *proc_border(void *args) {
	int i = 0, j = 0;
	struct args_t *a = (struct args_t *) args;
	double **data_dest = a->mesh->data_dest;
	
	/* Fronteira superior (linha fixa) */
    for(int j = a->border_lower_col; j <= a->border_upper_col; j++){
        data_dest[0][j] = 0.25 * (a->mesh->get_source_halo(-1, j) +
                                  a->mesh->get_source_halo(0, j+1) +
                                  a->mesh->get_source_halo(0, j-1) +
                                  a->mesh->get_source_halo(1, j));
    }
    
    /* Fronteira da esquerda (coluna fixa) */
	for(i = a->border_lower_row; i <= a->border_upper_row; i++){
		data_dest[i][0] = 0.25 * (a->mesh->get_source_halo(i-1, 0) +
								a->mesh->get_source_halo(i, 1) +
								a->mesh->get_source_halo(i, -1) +
								a->mesh->get_source_halo(i+1, 0));
	}
	
	/* Fronteira inferior (linha fixa) */
	i = a->mesh->get_size_y() - 1;
	for(int j = a->border_lower_col; j <= a->border_upper_col; j++){
		data_dest[i][j] = 0.25 * (a->mesh->get_source_halo(i-1, j) +
				a->mesh->get_source_halo(i, j+1) +
				a->mesh->get_source_halo(i, j-1) +
				a->mesh->get_source_halo(i+1, j));
	}
	
	/* Fronteira da direita (coluna fixa) */
	j = a->mesh->get_size_x() - 1;
	for(i = a->border_lower_row; i <= a->border_upper_row; i++){
		data_dest[i][j] = 0.25 * (a->mesh->get_source_halo(i-1, j) +
				a->mesh->get_source_halo(i, j+1) +
				a->mesh->get_source_halo(i, j-1) +
				a->mesh->get_source_halo(i+1, j));
	}
    
	pthread_exit(NULL);
}

double calc_diff(Mesh &mesh){
    double sum = 0.0;
	/*for(int i = 0; i < mesh.get_size_y(); i++){
        for(int j = 0; j < mesh.get_size_x(); j++){
            sum += pow(mesh.data_dest[i][j] - mesh.data_source[i][j], 2);
        }
    }*/
    
	for(int i = 0; i < mesh.get_num_threads(); i++) {
		//cout << "Criando thread " << i << endl;
		pthread_create(&(mesh.thread_id[i]), &(mesh.thread_attr), calc_diff_slice, (void*)&(mesh.thread_args[i]));
	}
	
	for(int i = 0; i < mesh.get_num_threads(); i++) {
		//cout << "Esperando thread " << i << endl;
		pthread_join(mesh.thread_id[i], NULL);
		//cout << "  feito!" << endl;
		sum += mesh.thread_args[i].val;
	}
	
	return(sum);
}

void *calc_diff_slice(void *args) {
	double partial_sum = 0.0;
	struct args_t *a = (struct args_t *) args;
	double **data_dest = a->mesh->data_dest, **data_source = a->mesh->data_source;
	
	for(int i = a->border_lower_row; i < a->border_upper_row; i++){
		for(int j = 0; j < a->mesh->get_size_x(); j++){
			partial_sum += pow(data_dest[i][j] - data_source[i][j], 2);
		}
	}
	
	a->val = partial_sum;
	
	pthread_exit(NULL);
}
