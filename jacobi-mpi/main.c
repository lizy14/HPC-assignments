#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#define TIME(a,b) (1.0*((b).tv_sec-(a).tv_sec)+0.000001*((b).tv_usec-(a).tv_usec))

extern int Init(double *data, long long L);
extern int Check(double *data, long long L);

//Can be modified
typedef struct
{
    int nx, ny, nz;
} Info;

int test(){
    return test_partition();
}
int test_partition(){
    int n = 27;
    int k = 4;
    for(int i=0; i<k; ++i){
        printf(
            "%d %d %d\n",
            i,
            number_of_elements_in_group(i, n, k),
            index_of_last_element_in_group(i, n, k)
        );

    }
    return 0;
}

// evenly distribute elements into groups
// i-th (starting from 0) group; n elements; k groups
// 27 / 4 -> 7, 7, 7, 6
int number_of_elements_in_group(int i, int n, int k){
    int q = n / k;
    int r = n % k;
    return q + (i < r? 1: 0);
}
// 27 / 4 -> 6, 13, 20, 26
int index_of_last_element_in_group(int i, int n, int k){
    int result = -1;
    for(int _=0; _<=i; ++_){
        result += number_of_elements_in_group(_, n, k);
    }
    return result;
}

//Must be modified because only single process is used now
Info setup(int NX, int NY, int NZ, int P)
{
    Info result;
    int myrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if(myrank == 0)
    {
        result.nx = NX;
        result.ny = NY;
        result.nz = NZ;
    }
    else
    {
        result.nx = 0;
        result.ny = 0;
        result.nz = 0;
    }
    return result;
}


//Must be re-written, including all the parameters
int stencil(double *A, double *B, int nx, int ny, int nz, int steps)
{
    int i, j, k, s;
#define IDX(i,j,k) ((i)*ny*nz+(j)*nz+(k))
    for(s = 0; s < steps; s ++) {
        for(i = 0; i < nx; i ++) {
            for(j = 0; j < ny; j ++) {
                for(k = 0; k < nz; k ++) {
                    double r = 0.4*A[IDX(i,j,k)];
                    if(k !=  0)   r += 0.1*A[IDX(i,j,k-1)];
                    else          r += 0.1*A[IDX(i,j,k)];
                    if(k != nz-1) r += 0.1*A[IDX(i,j,k+1)];
                    else          r += 0.1*A[IDX(i,j,k)];
                    if(j !=  0)   r += 0.1*A[IDX(i,j-1,k)];
                    else          r += 0.1*A[IDX(i,j,k)];
                    if(j != ny-1) r += 0.1*A[IDX(i,j+1,k)];
                    else          r += 0.1*A[IDX(i,j,k)];
                    if(i !=  0)   r += 0.1*A[IDX(i-1,j,k)];
                    else          r += 0.1*A[IDX(i,j,k)];
                    if(i != nx-1) r += 0.1*A[IDX(i+1,j,k)];
                    else          r += 0.1*A[IDX(i,j,k)];
                    B[IDX(i,j,k)] = r;
                }
            }
        }
        double *tmp = NULL;
        tmp = A, A = B, B = tmp;
    }
    return 0;
}

int main(int argc, char **argv) {
    return test();
/*
    double *A = NULL, *B = NULL;
    int myrank, nprocs, nx, ny, nz;
    int NX=100, NY=100, NZ=100, STEPS=10;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    NX = atoi(argv[1]);
    NY = atoi(argv[2]);
    NZ = atoi(argv[3]);
    STEPS = atoi(argv[4]);

    if(myrank == 0){
        printf("Size: %dx%dx%d, # of Steps: %d, # of procs: %d\n",
                NX, NY, NZ, STEPS, nprocs);
    }
    Info info = setup(NX, NY, NZ, nprocs); //This is a single version
    nx = info.nx, ny = info.ny, nz = info.nz;

    long long size = nx*ny*nz;

    A = (double*)malloc(size*sizeof(double));
    B = (double*)malloc(size*sizeof(double));

    Init(A, size);

    struct timeval t1, t2;
    MPI_Barrier(MPI_COMM_WORLD);
    gettimeofday(&t1, NULL);

    stencil(A, B, nx, ny, nz, STEPS);

    MPI_Barrier(MPI_COMM_WORLD);
    gettimeofday(&t2, NULL);

    if(myrank == 0){
        printf("Total time: %.6lf\n", TIME(t1,t2));
    }

    if(STEPS%2){
        Check(B, size);
    }else{
        Check(A, size);
    }

    free(A);
    free(B);

    MPI_Finalize();
    */
}
