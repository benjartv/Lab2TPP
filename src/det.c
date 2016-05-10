#include "detfunct.h"

int main(int argc, char *argv[]){

	//Obtener argumentos de entrada
	Options variables = get_variables(argc, argv);

	//Almacenar matriz en un arreglo
	float *buffer = read_input(variables.input_file, variables.matrix_size);

  //Calcula el que índices debe analizar cada hebra
	P_range* threads_work = calculate_work(variables.n_threads, variables.matrix_size);
	int tid;
	float determinant = 0;

  //Variables para calcular el tiempo
  double time_spent;
  struct timeval start, end;
  gettimeofday(&start, NULL);

  //Bloque paralelo, cada hebra calcula el determinante de los íncides que le corresponden y
  //luego suma los valores de cada hebra con la claúsula reduction
	#pragma omp parallel private(tid) num_threads(variables.n_threads) reduction(+:determinant)
	{
		tid = omp_get_thread_num();
    determinant = parallel_det(tid, variables.matrix_size, buffer, threads_work);
	}

  gettimeofday(&end, NULL);
  time_spent = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

  printf("Time: %f\n", time_spent);
	printf("Determinant: %f\n", determinant);

	return 0;
}
