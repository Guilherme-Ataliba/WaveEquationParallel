#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include "auxilary.c"

/*
    This is the shared memory version
    - No halo regions are required since they can be accessed through the global memory
*/

#define Nx 800
#define Ny 800
#define T 10
#define a 10.0
#define b 10.0
#define Cx 0.5
#define Cy 0.5
#define frames 20

__device__ double f(int i, int j, double t) {
    if (i < (Nx / 2 + Nx / 100) && i > (Nx / 2 - Nx / 100)) {
        if (j < (Ny / 2 + Ny / 100) && j > (Ny / 2 - Ny / 100)) {
            return 30 * sin(4.7123889 * t);
        }
    }
    return 0.0;
}

__global__ void step_zero_kernel(double* unm2) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < Nx && j < Ny) {
        int idx = j * Nx + i;
        unm2[idx] = 0; // Initial condition
    }
}

__global__ void step_one_kernel(double* unm1, double* unm2,double Cx2, double Cy2, double dt2) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i > 0 && i < Nx - 1 && j > 0 && j < Ny - 1) {
        int idx = j * Nx + i;
        unm1[idx] = unm2[idx]
                    + 0.5 * Cx2 * (unm2[idx + 1] - 2 * unm2[idx] + unm2[idx - 1])
                    + 0.5 * Cy2 * (unm2[idx + Nx] - 2 * unm2[idx] + unm2[idx - Nx])
                    + 0.5 * dt2 * f(i, j, 1);
    }
}

__global__ void next_step_kernel(double* u, double* unm1, double* unm2, double t, double Cx2, double Cy2, double dt2) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i > 0 && i < Nx - 1 && j > 0 && j < Ny - 1) {
        int idx = j * Nx + i;
        u[idx] = -unm2[idx] + 2 * unm1[idx]
                 + Cx2 * (unm1[idx + 1] - 2 * unm1[idx] + unm1[idx - 1])
                 + Cy2 * (unm1[idx + Nx] - 2 * unm1[idx] + unm1[idx - Nx])
                 + dt2 * f(i, j, t);
    }
}

int main() {
    double dx = a / Nx, dy = b / Ny;
    double Cx2 = Cx * Cx, Cy2 = Cy * Cy, c=1.2;
    double dt = Cx * dx/c, dt2 = dt * dt;
    int Nt = floor(T / dt);

    // Create events for timing
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start, 0);

    // double *x, *y;
    double *u, *unm1, *unm2;
    double *h_u;

    h_u = (double*)malloc(Nx*Ny*sizeof(double));

    int index=0;

    // We have to define the number of files in a smart way, since it may not result in a integer value
    int print_interval = ceil(T/dt/frames), print_counter=print_interval;
    int n_files = floor(T/dt/print_interval)+2;

    // Allocate memory on the host
    // cudaMallocManaged(&x, Nx * sizeof(double));
    // cudaMallocManaged(&y, Ny * sizeof(double));
    cudaMallocManaged(&u, Nx * Ny * sizeof(double));
    cudaMallocManaged(&unm1, Nx * Ny * sizeof(double));
    cudaMallocManaged(&unm2, Nx * Ny * sizeof(double));

    // for (int i = 0; i < Nx; i++) x[i] = i * dx;
    // for (int j = 0; j < Ny; j++) y[j] = j * dy;

    dim3 blockSize(16, 16);
    dim3 gridSize((Nx + blockSize.x - 1) / blockSize.x, (Ny + blockSize.y - 1) / blockSize.y);

    step_zero_kernel<<<gridSize, blockSize>>>(unm2);

    step_one_kernel<<<gridSize, blockSize>>>(unm1, unm2, Cx2, Cy2, dt2);

    for (int t = 0; t < Nt; t++) {
        next_step_kernel<<<gridSize, blockSize>>>(u, unm1, unm2, t * dt, Cx2, Cy2, dt2);

        if(print_counter >= print_interval){
            cudaMemcpy(h_u, u, Nx*Ny*sizeof(double), cudaMemcpyDeviceToHost);
            write_linear_dir_matrix(h_u, Nx, Ny, "output/multi_matrix/", &index);
            print_counter = 0;
        }
        print_counter++;

        // Swap pointers
        double* temp = unm2;
        unm2 = unm1;
        unm1 = u;
        u = temp;
    }

    // Record the stop event
    cudaEventRecord(stop, 0);

    // Wait for the events to finish
    cudaEventSynchronize(stop);

    // Calculate elapsed time
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);

    FILE *time_file;
    time_file = fopen("output/Execution_time_cuda.txt", "w");
    fprintf(time_file, "Elapsed Time:%lf\ndx = %g \t dy = %g\nNx = %d \t Ny = %d\na=%lf\nb=%lf\nT=%d\nc=%lf\nCx=%lf \t Cy=%lf",
    milliseconds/1000, dx, dy, Nx, Ny, a, b, T, c, Cx, Cy);

    // cudaFree(x);
    // cudaFree(y);
    cudaFree(u);
    cudaFree(unm1);
    cudaFree(unm2);

    return 0;
}
