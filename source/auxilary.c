#include <stdio.h>
#include <string.h>
#include "../include/auxilary.h"

void copy_vector(double *arr1, double *arr2, int size){
    //Copy arr1 into arr2
    
    int i;

    for(i=0; i<size; i++){
        arr1[i] = arr2[i];
    }
}

void copy_matrix(double **m1, double **m2, int row, int col){
    int i, j;

    for(j=0; j<col; j++){
        for(i=0; i<row; i++){
            m1[j][i] = m2[j][i];
        }
    }
}

void switch_matrices(double ***u, double ***unm1, double ***unm2, double ***aux){
    *aux = *u;
    *u = *unm2;
    *unm2 = *unm1;
    *unm1 = *aux;
}

void switch_matrices_linear(double **u, double **unm1, double **unm2, double **aux){
    *aux = *u;
    *u = *unm2;
    *unm2 = *unm1;
    *unm1 = *aux;
}

void print_vector(double *vector, int size){
    int i=0;
    for(i=0; i<size; i++){
        printf("%lf\t", vector[i]);
    }

    return;
}

void print_matrix(double **matrix, int row, int col){
    int i, j;
    for(i=0; i<row; i++){
        for(j=0; j<col; j++){
            printf("%lf ", matrix[i][j]);
        }
        printf("\n");
    }
}

void write_vector(double *vector, int size, FILE *path){
    int i;
    for(i=0; i<size-1; i++){
        fprintf(path, "%.8lf,", vector[i]);
    }
    fprintf(path, "%.8lf", vector[i]);
    fprintf(path, "\n");
}

void write_matrix(double **matrix, int row, int col, FILE *path){
    int i, j;
    for(i=0; i<row; i++){
        for(j=0; j<col-2; j++){
            fprintf(path, "%.4lf,", matrix[i][j]);
        }
        fprintf(path, "%.4lf", matrix[i][col-1]);
        fprintf(path, "\n");
    }
}

// Binary version
void write_dir_matrix(double **matrix, int row, int col, char *dir_path, int *index){
    int str_size = strlen(dir_path), i;
    char file[str_size+20];
    FILE *output;

    sprintf(file, "%s%d.bin", dir_path, *index);
    (*index)++;
    output = fopen(file, "wb");

    for(i=0; i<row; i++) fwrite(matrix[i], sizeof(double), col, output);
    fclose(output);
}

// Binary version
void write_linear_dir_matrix(double *matrix, int row, int col, char *dir_path, int *index){
    int str_size = strlen(dir_path), i;
    char file[str_size+20];
    FILE *output;

    sprintf(file, "%s%d.bin", dir_path, *index);
    (*index)++;
    output = fopen(file, "wb");

    fwrite(matrix, sizeof(double), row*col, output);
    fclose(output);
}

//1 Binary file version
void write_multiple_matrix(double **matrix, int row, int col, FILE *output){
    int i;
    for(i=0; i<row; i++) fwrite(matrix[i], sizeof(double), col, output);
}