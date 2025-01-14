#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include "../include/auxilary.h"

#define Nx 800
#define Ny 800
#define T 10
#define a 10.0
#define b 10.0
#define Cx 0.5
#define Cy 0.5
#define GHOST_ROWS 1 // Modular number of ghost rows
#define frames 20

#define idx(i, j, col) ((j) * (col) + (i))

double f(int i, int j, double t) {
    if (i < (Nx / 2 + Nx / 100) && i > (Nx / 2 - Nx / 100)) {
        if (j < (Ny / 2 + Ny / 100) && j > (Ny / 2 - Ny / 100)) {
            return 30 * sin(4.7123889 * t);
        }
    }
    return 0;
}

void enforce_boundary_conditions(double* U, int local_Ny, int rank, int size) {
    // Apply top boundary condition on rank 0
    if (rank == 0) {
        for (int i = 0; i < Nx; i++) {
            U[idx(i, GHOST_ROWS - 1, Nx)] = 0; // First physical row
        }
    }

    // Apply bottom boundary condition on the last rank
    if (rank == size - 1) {
        for (int i = 0; i < Nx; i++) {
            U[idx(i, local_Ny - GHOST_ROWS, Nx)] = 0; // Last physical row
        }
    }
}

double* inicializa_matriz_mpi(int rank, int size, int* n_linhas_out) {
    int n_linhas = Ny / size;
    *n_linhas_out = n_linhas + GHOST_ROWS * 2; // Add ghost rows for both top and bottom

    double* u = calloc((*n_linhas_out) * Nx, sizeof(double));
    return u;
}

void communicate_boundaries(double* u, int n_linhas, int rank, int size, MPI_Comm comm) {
    MPI_Status status;
    int cima = rank - 1;   // Rank above
    int baixo = rank + 1;  // Rank below

    // Send and receive with the bottom neighbor
    if (baixo < size) {
        MPI_Sendrecv(&u[(n_linhas - GHOST_ROWS * 2) * Nx], Nx * GHOST_ROWS, MPI_DOUBLE, baixo, 0, 
                     &u[(n_linhas - GHOST_ROWS) * Nx], Nx * GHOST_ROWS, MPI_DOUBLE, baixo, 0, MPI_COMM_WORLD, &status);
    }

    // Send and receive with the top neighbor
    if (cima >= 0) {
        MPI_Sendrecv(&u[GHOST_ROWS * Nx], Nx * GHOST_ROWS, MPI_DOUBLE, cima, 0, 
                     &u[0 * Nx], Nx * GHOST_ROWS, MPI_DOUBLE, cima, 0, MPI_COMM_WORLD, &status);
    }
}

void next_step(double* u, double* unm1, double* unm2, double Cx2, double Cy2, double dt2, int local_Ny, int rank, int size, double t, MPI_Comm comm) {
    for (int j = GHOST_ROWS; j < local_Ny - GHOST_ROWS; j++) {
        for (int i = 1; i < Nx - 1; i++) {
            u[idx(i, j, Nx)] = -unm2[idx(i, j, Nx)] + 2 * unm1[idx(i, j, Nx)]
                             + Cx2 * (unm1[idx(i + 1, j, Nx)] - 2 * unm1[idx(i, j, Nx)] + unm1[idx(i - 1, j, Nx)])
                             + Cy2 * (unm1[idx(i, j + 1, Nx)] - 2 * unm1[idx(i, j, Nx)] + unm1[idx(i, j - 1, Nx)])
                             + dt2 * f(i, rank * (local_Ny - 2 * GHOST_ROWS) + j - GHOST_ROWS, t);
        }
    }

    communicate_boundaries(u, local_Ny, rank, size, comm);
}

void gather_and_save(double *u, int rank, int size, int n_linhas, char *dir_path, int* file_index) {
    // MPI_Barrier(MPI_COMM_WORLD);
    
    double *u_linear = malloc(Nx * Ny * sizeof(double));

    int *sendcounts = NULL;
    int *displs = NULL;

    if (rank == 0) {
        sendcounts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));
    }

    int sendcount = (n_linhas - 2 * GHOST_ROWS) * Nx;
    MPI_Gather(&sendcount, 1, MPI_INT, sendcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        displs[0] = 0;
        for (int i = 1; i < size; i++) {
            displs[i] = displs[i - 1] + sendcounts[i - 1];
        }
    }

    MPI_Gatherv(&u[GHOST_ROWS * Nx], sendcount, MPI_DOUBLE,
                u_linear, sendcounts, displs, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    if (rank == 0) {
        free(sendcounts);
        free(displs);
        write_linear_dir_matrix(u_linear, Nx, Ny, dir_path, file_index);
    }
}

int main(int argc, char** argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int local_Ny = Ny / size + 2 * GHOST_ROWS; // Include ghost rows
    double dx = a / Nx, dy = b / Ny;
    double Cx2 = Cx * Cx, Cy2 = Cy * Cy, c=1.2;
    double dt = Cx * dx / c, dt2 = dt * dt;
    int Nt = floor(T / dt);
    int n_linhas;

    char dir_path[] = "output/multi_matrix/";
    int file_index = 0;

    // We have to define the number of files in a smart way, since it may not result in a integer value
    int print_interval = ceil(T/dt/frames), print_counter=print_interval;
    int n_files = floor(T/dt/print_interval)+2;

    double start_time, end_time;

    // Synchronize all processes and start the timer
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    double *u, *unm1, *unm2;

    u = inicializa_matriz_mpi(rank, size, &n_linhas);
    unm1 = inicializa_matriz_mpi(rank, size, &n_linhas);
    unm2 = inicializa_matriz_mpi(rank, size, &n_linhas);

    for (double t = 0; t < T; t += dt) {
        enforce_boundary_conditions(u, n_linhas, rank, size);
        next_step(u, unm1, unm2, Cx2, Cy2, dt2, n_linhas, rank, size, t, MPI_COMM_WORLD);

         if(print_counter >= print_interval){
            if (rank == 0)
                printf("Current time step: %.2lf/%d\n", t, T);
            gather_and_save(u, rank, size, n_linhas, dir_path, &file_index);
            print_counter=0;
        }
        print_counter++;

        double* temp = unm2;
        unm2 = unm1;
        unm1 = u;
        u = temp;
    }

    free(u);
    free(unm1);
    free(unm2);

    // Synchronize again and stop the timer
    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();

    FILE *time_file;
    time_file = fopen("output/Execution_time_mpi.txt", "w");
    fprintf(time_file, "Elapsed Time:%lf\ndx = %g \t dy = %g\nNx = %d \t Ny = %d\na=%lf\nb=%lf\nT=%d\nc=%lf\nCx=%lf \t Cy=%lf",
    end_time - start_time, dx, dy, Nx, Ny, a, b, T, c, Cx, Cy);

    MPI_Finalize();
    return 0;
}
