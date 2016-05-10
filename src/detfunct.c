#include "detfunct.h"


/*	void usage()
	Entrada: none
	Salida: none
	Func: imprime por stdout las diferentes opciones de uso del programa det
*/
void usage(){
	printf("usage: det \t-i input_matrix\n");
	printf("\t\t-N size of the matrix\n");
	printf("\t\t-H number of threads\n");
	printf("\t\t-h (Display this message)\n");
}


/*	int* read_input(char* input_file, matrix_size)
	Entrada: nombre archivo de matriz de entrada, tamaño de la matriz
	Salida: buffer con el contenido de la matriz de entrada
	Func: lee el contenido del archivo de entrada y almacena los valores de los elementos
		de la matriz en un buffer de flotantes
*/
float* read_input(char* input_file, int matrix_size){
	int t_elem_matrix = matrix_size*matrix_size;
	float* buffer = (float*)malloc(sizeof(float)*t_elem_matrix);
	int filedesc = open(input_file, O_RDONLY);
	read(filedesc, buffer, matrix_size*sizeof(float)*t_elem_matrix);
	close(filedesc);
	return buffer;
}


/*	P_range* calculate_work(int n_threads, int matrix_size)
	Entrada: cantidad de hebras y tamaño de la matriz
	Salida: lista de P_range, que contiene los índices que cada hebra tiene que procesar
	Funct: crea un arreglo de tamaño n_threads y distribuye el trabajo en todas las hebras, 
		considerando que si el trabajo no se puede distribuir equitativamente entre las hebras,
		por lo menos la diferencia de trabajo sea de +-1
		crea un segundo arreglo (de estrucutras) que contiene para cada hebra el primer y último
		índice que debe procesar
*/
P_range* calculate_work(int n_threads, int matrix_size){
	int i;
	int* threads_work = (int*)malloc(sizeof(int)*n_threads);
	P_range* range_work = (P_range*)malloc(sizeof(P_range)*n_threads);
	int base = matrix_size / n_threads;
	for (i = 0; i < n_threads; i++){
		threads_work[i] = base;
	}
	if (matrix_size % n_threads != 0){
		i = 0;
		int rest = matrix_size % n_threads;
		while(rest > 0){
			threads_work[i] += 1;
			rest--;
			i++;
			if (i == n_threads){
				i = 0;
			}
		}
	}
	for (i = 0; i < n_threads; i++){
		if (i == 0){
			range_work[i].begin = 0;
			range_work[i].end = threads_work[i];
		}
		else{
			range_work[i].begin=range_work[i-1].end;
			range_work[i].end=range_work[i-1].end + threads_work[i];
		}
	}
	free(threads_work);
	return range_work;
}


/*	float calculate_det(float *matrix, int mt_size, int mn_size, int n_row, int*col)
	Entrada: buffer con los elementos de la matriz, tamaño de la matriz, tamaño del minor
		fila actual, arreglo que representa que columnas no debe considerar para construir
		el minor, 1 si la columna debe considerarse, 0 en caso contrario
	Saluda: determinante de la matriz
	Func: Si el minor es igual a 2, la función retorna el determinante como la multiplicación
	cruzada de los elementos del minor.
	Sino para cada índice en la fila actual, si debe considerar la columna, llama recursivamente
	a la función calculate_det, indicando un tamaño menor del minor, la fila siguiente y el 
	arreglo de columnas(con la columna actual con valor 0, para que no sea considerada por la
	siguiente iteración). Cuando recupera el valor del determinante del minor, desmarca la columna
	para que la siguiente iteración pueda considerarla al momento de calcular el determinante del nuevo
	minor
*/
float calculate_det(float *matrix, int mt_size, int mn_size, int n_row, int*col){
	int i;
	float determinant;
	if (mn_size == 2){
		int c1 = -1, c2;
		for (i = 0; i < mt_size; i++){
			if (col[i] == 1){
				if (c1 == -1)
					c1 = i;
				else
					c2 = i;
			}
		}
		determinant = (matrix[mt_size*n_row + c1] * matrix[mt_size*(n_row+1) + c2]) - (matrix[mt_size*(n_row+1) + c1] * matrix[mt_size*n_row + c2]);
	}
	else{
		determinant = 0;
		int sub = 0;
		for (i = 0; i < mt_size; i++)
		{
			//Si la columna no ha sido marcada, calcular minor
			if (col[i] == 1){
				col[i] = 0;
				determinant += pow(-1.0, sub)*matrix[mt_size*n_row + i]*calculate_det(matrix, mt_size, mn_size-1, n_row+1, col);
				//Desmarcar columna para siguiente iteración
				col[i] = 1;
				sub++;
			}
		}
	}
	return(determinant);
}


/*	float parallel_det(int tid, int matrix_size, float*buffer, P_range* threads_work, float determinant)
	Entrada: tid de la hebra, tamaño de la matriz, arreglo con la matriz, arreglo con el trabajo
	de cada hebra(índices de la matriz)
	Salida: determinante parcial de la matriz (correspondiente a los índices de esa hebra)
	Func: crea un arreglo igual a la cantidad de columnas de la matriz que será utilizado por la función
		calculate_det para determinar el minor en cada iteración, lo inicializa en 1, que representa que cada
		columna de la matriz esta disponible para el calculo del minor.
		Luego para cada índice correspondiente a la hebra actual, marcar la columna y calcular el determinante
		del minor para multiplicarlo por el indice actual, luego desmarcar la columna y repetir el proceso
		para todos los índices.
*/
float parallel_det(int tid, int matrix_size, float*buffer, P_range* threads_work){
	float determinant = 0;
	int i,j;
    int *col_vector = (int*)malloc(sizeof(int)*matrix_size);

    //Deja todas las columnas disponibles (revisar función calculate_det)
    for (j = 0; j < matrix_size; j++){
        col_vector[j] = 1;
    }

    //Cada hebra calcula el determinante de los índices que le tocaron
    for (i = threads_work[tid].begin; i < threads_work[tid].end; i++)
    {
      col_vector[i] = 0;
      determinant += pow(-1.0, i)*buffer[i]*calculate_det(buffer, matrix_size, matrix_size-1, 1, col_vector);
      col_vector[i] = 1;
    }
    free(col_vector);
    return determinant;
}


/*	Options get_variables(int argc, char *argv[])
	Entrada: cantidad de argumentos de ejecucióm_size del programa, argumentos de ejecucióm_size
	Salida: estructura con la informacióm_size de los argumentos de ejecucióm_size
	Func: utilizando getopt() almecena los diferentes argumentos de ejecucióm_size en una
		estructura de datos
*/
Options get_variables(int argc, char *argv[]){
	if (argc < 7){
		usage();
		exit(1);
	}
	Options variables;
	int gopt;
	while((gopt = getopt(argc, argv, "hi:N:H:") ) != -1)
		switch(gopt){
			case 'h':
				usage();
				exit(1);
			case 'i':
				variables.input_file = (char*)malloc(sizeof(char)*strlen(optarg));
				strcpy(variables.input_file, optarg);
				break;
			case 'N':
				variables.matrix_size = atoi(optarg); 
				break;
			case 'H':
				variables.n_threads = atoi(optarg);
				break;
			default:
				abort();
		}
	return variables;
}