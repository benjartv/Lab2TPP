#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <ctype.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include <time.h>

/*	struct getopt_opc{} Options
	nombre archivo matriz de entrada
	tamaño de la matriz
	numero de hebras
*/
typedef struct getopt_opc{
	char* input_file;
	int matrix_size;
	int n_threads;
}Options;

/* struct parallel_range{}P_range
	índice inicial
	índice final
*/
typedef struct parallel_range{
	int begin;
	int end;
}P_range;

void help();
float* read_input(char*, int);
Options get_variables(int, char **);
P_range* calculate_work(int, int);
float calculate_det(float*, int, int, int, int*);
float parallel_det(int, int, float*, P_range*);
