#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#define TIME(a,b) (1.0*((b).tv_sec-(a).tv_sec)+0.000001*((b).tv_usec-(a).tv_usec))
#define debug_printf(fmt, ...) (0)
extern int Init(double *data, long long L);
extern int Check(double *data, long long L);


int test(){
    return test_partition();
}
int test_partition(){
    int n = 27;
    int k = 4;
    int i;
    for(i=0; i<k; ++i){
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
    int _;
    for(_=0; _<=i; ++_){
        result += number_of_elements_in_group(_, n, k);
    }
    return result;
}



//Must be re-written, including all the parameters
int stencil(
        double *A, double *B,
        int nx, int ny, int nz, int steps,
        int first, int last, int rank, int nprocs)
{
    #define IDX(i,j,k) ((i) + (j)*nx + (k)*nx*ny)
    int s = nx * ny;
    int i_am_head = first == 0;
    int i_am_tail = last == nz - 1;
    int i, j, k, step;
    for(step = 0; step < steps; step ++) {
        for(k = first; k <= last; k ++) {
            for(j = 0; j < ny; j ++) {
                for(i = 0; i < nx; i ++) {


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
//                    debug_printf("#%d: %d, %d, %d: r=%lf\n", rank, i,j,k, r);
                    B[IDX(i,j,k)] = r;
                }
            }
        }
        // message passing

        MPI_Request reqs[4];
        MPI_Status stats[4];

//        debug_printf("#%d: tail=%d, head=%d, s=%d, first=%d, last=%d\n", rank, i_am_tail, i_am_head, s, first, last);
        debug_printf("#%d: %lf, %lf, %lf, %lf\n", rank, B[0], B[1], B[2], B[3]);

        if(i_am_tail && i_am_head){
            // only one worker, nothing done
        }else{
            if(! i_am_head){
                MPI_Isend(B + IDX(0,0,first),
                        s, MPI_DOUBLE,
                        rank - 1,
                        0, MPI_COMM_WORLD,
                        &reqs[0]
                        );
                MPI_Irecv(B + IDX(0,0,first - 1),
                        s, MPI_DOUBLE,
                        rank - 1,
                        0, MPI_COMM_WORLD,
                        &reqs[1]
                        );
            }
            if(! i_am_tail){
                MPI_Isend(B + IDX(0,0,last),
                        s, MPI_DOUBLE,
                        rank + 1,
                        0, MPI_COMM_WORLD,
                        &reqs[2]
                        );
                MPI_Irecv(B + IDX(0,0,last + 1),
                        s, MPI_DOUBLE,
                        rank + 1,
                        0, MPI_COMM_WORLD,
                        &reqs[3]
                        );
            }
            if(! i_am_tail && ! i_am_head){
                MPI_Waitall(4, reqs, stats);
            }else if(i_am_tail){
                MPI_Waitall(2, reqs, stats);
            }else if(i_am_head){
                MPI_Waitall(2, reqs+2, stats);
            }
        }

        debug_printf("#%d synced: %lf, %lf, %lf, %lf\n", rank, B[0], B[1], B[2], B[3]);




        double *tmp = NULL;
        tmp = A, A = B, B = tmp;
    }
MPI_Barrier(MPI_COMM_WORLD);
    // send everything to master 0
    if(rank == 0){
        int _i;
        for(_i=1; _i<nprocs; ++_i){

            int _n = number_of_elements_in_group(_i, nz, nprocs);
            int _last = index_of_last_element_in_group(_i, nz, nprocs);
            int _first = _last - _n + 1;
            MPI_Status st;
            MPI_Recv(B + IDX(0,0,_first),
                    s * _n, MPI_DOUBLE,
                    _i,
                    0, MPI_COMM_WORLD,
                    &st
                    );
        }
    }else{
        MPI_Send(B + IDX(0,0,first),
                s*(last - first + 1), MPI_DOUBLE,
                0, // master
                0, MPI_COMM_WORLD
                );
    }
//return 0;

        debug_printf("#%d sent to master: %lf, %lf, %lf, %lf\n", rank, B[0], B[1], B[2], B[3]);
        MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0){
        int _i;
        for(_i=1; _i<nprocs; ++_i){
            MPI_Send(B, s * nz, MPI_DOUBLE, _i, 0, MPI_COMM_WORLD);
        }
    }else{
        MPI_Status st;
        MPI_Recv(B, s * nz, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &st);
    }

        debug_printf("#%d returning: %lf, %lf, %lf, %lf\n", rank, B[0], B[1], B[2], B[3]);
/*    MPI_Allgather(
        B,
        s * (last - first + 1),
        MPI_DOUBLE,
        A,
        s * nz,
        MPI_DOUBLE,
        MPI_COMM_WORLD
            );*/

    return 0;
}

int main(int argc, char **argv) {
    //    return test();

    double *A = NULL, *B = NULL;
    int myrank, nprocs, nx, ny, nz;
    int NX=100, NY=100, NZ=100, STEPS=10;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if(argc >= 5){
        NX = atoi(argv[1]);
        NY = atoi(argv[2]);
        NZ = atoi(argv[3]);
        STEPS = atoi(argv[4]);
    }

    if(myrank == 0){
        printf("Size: %dx%dx%d, # of Steps: %d, # of procs: %d\n",
                NX, NY, NZ, STEPS, nprocs);
    }

    nx = NX, ny = NY, nz = NZ;

    // allocate matricies
    long long size = nx * ny * nz;

    A = (double*)malloc(size*sizeof(double));
    B = (double*)malloc(size*sizeof(double));

    // initialize matrix
    Init(A, size);


    debug_printf("Initialized: %lf, %lf, %lf, %lf\n", A[0], A[1], A[2], A[3]);
    // get assignment
    int n = number_of_elements_in_group(myrank, nz, nprocs);
    int last = index_of_last_element_in_group(myrank, nz, nprocs);
    int first = last - n + 1;


    struct timeval t1, t2;
    MPI_Barrier(MPI_COMM_WORLD);
    gettimeofday(&t1, NULL);

    debug_printf("#%d: %d - %d\n", myrank, first, last);

    stencil(A, B, nx, ny, nz, STEPS, first, last, myrank, nprocs);

    MPI_Barrier(MPI_COMM_WORLD);
    gettimeofday(&t2, NULL);

    if(myrank == 0){
        printf("Total time: %.6lf\n", TIME(t1,t2));
    }
    if(STEPS % 2){
        double* tmp = B;
        B = A;
        A = tmp;
    }
    debug_printf("#%d Checking: %lf, %lf, %lf, %lf\n", myrank, B[0], B[1], B[2], B[3]);
    Check(B, size);


    free(A);
    free(B);

    MPI_Finalize();
}
