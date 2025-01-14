#ifndef AUXILARY_h
#define AUXILARY_h

void copy_vector(double *arr1, double *arr2, int size);
void copy_matrix(double **m1, double **m2, int row, int col);
void switch_matrices(double ***u, double ***unm1, double ***unm2, double ***aux);
void switch_matrices_linear(double **u, double **unm1, double **unm2, double **aux);
void print_vector(double *vector, int size);
void print_matrix(double **matrix, int row, int col);
void write_vector(double *vector, int size, FILE *path);
void write_matrix(double **matrix, int row, int col, FILE *path);
void write_dir_matrix(double **matrix, int row, int col, char *dir_path, int* index);
void write_linear_dir_matrix(double *matrix, int row, int col, char *dir_path, int* index);
void write_multiple_matrix(double **matrix, int row, int col, FILE *output);


#endif