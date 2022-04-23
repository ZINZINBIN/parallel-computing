/**
 * @file main.cpp
 * @author zinzinbin
 * @brief matrix multiplication with OpenMPI + non blocking method
 * @version 0.1
 * @date 2022-04-21
 * 
 * How to execute
 * (1) mpic++ matrix_multiplication_non_blocking.cpp -o matrix_multiplication_non_blocking.out
 * (2) mpirun -np <number of processes> ./matrix_multiplication_non_blocking.out
 */

#include <iostream>
#include <mpi.h>
#include <time.h>
#include <array>
#include <vector>
#include <chrono>
#include <math.h>
#include <cmath>
#define EPS 1e-8
#define MAXPROC 16

using namespace std;

void init_random_array(float *x, int n_size, float min, float max);
void init_random_matrix(float **A, int m, int n, float min, float max);
float get_random_number(float min, float max);
float** generate_matrix(int m, int n);
float* generate_array(int m);
void check_solution_validity(float **A, float *x, float *target, int m, int n, float tolerance);

int main(int argc, char *argv[]){

    int baton = 1;
    
    int rank, size; 
    int tag;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size); // get number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get current process id

    int m = 16;
    int n = 16;
    float **A = generate_matrix(m,n);
    float *x = generate_array(n);
    float *y = generate_array(m);
    float ans = 0;

    // initiate matrix and array with each component extracted from range of (min, max) randomly
    float min = -1.0;
    float max = 1.0;
    init_random_matrix(A, m, n, min, max);
    init_random_array(x, n, min, max);
    init_random_array(y, m, min, max);

    if(rank == 0){
        cout << "number of processes : " << size << endl;
        cout << "current process id : " << rank << endl;
        cout << "matrix : (" << m << "," << n << ")" << endl;
    }

    int numsent = 0;
    int sender = 0;
    int row = 0;
    int po = 1;
    MPI_Status status;
    MPI_Request send_req[MAX_PROC], recv_req[MAX_PROC];
    float *buffer;

    if(rank == 0){
        double start = MPI_Wtime();
        numsent = 0;

        for(int i = 0; i < size - 1 && i < m; i++){
            MPI_Isend(&A[i][0], 1, MPI_FLOAT, i+1, i, MPI_COMM_WORLD, &send_req[i]); // *ptr_buffer, count, MPI_Datatype, target, tag, MPI_Comm
            MPI_Isend(&x[0], 1, MPI_FLOAT, i+1, i, MPI_COMM_WORLD, &send_req[i]);
            numsent ++;
            cout << "MPI send to rank : " << i+1 << " complete" << endl;
        }

        MPI_Waitall(size - 1, &send_req[1], &status);
        cout << "MPI rank 0 process receiving from all other processes" << endl;

        for(int i = 0; i < m; i++){
            MPI_Irecv(&ans, 1, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            sender = status.MPI_SOURCE;
            row = status.MPI_TAG;
            y[row] = ans;

            if(numsent < m){
                MPI_Send(&A[numsent][0], n, MPI_FLOAT, sender, numsent, MPI_COMM_WORLD);
                numsent ++;
            }
            else{
                MPI_Send(MPI_BOTTOM, 0, MPI_FLOAT, sender, n, MPI_COMM_WORLD);
            }
        }
        
        double finish = MPI_Wtime();
        cout << "rank : 0, process done in " << finish - start << "(s)" << endl;
    }
    else{
        MPI_Bcast(&x[0], n, MPI_FLOAT, 0, MPI_COMM_WORLD);
        buffer = (float *)malloc(sizeof(float)*n);
        po = 1;
        while(rank <= m && po){
            MPI_Recv(buffer, n, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            row = status.MPI_TAG;

            if(row < n){
                ans = 0;
                for(int j = 0; j < n; j++){
                    ans += buffer[j] * x[j];
                }
                MPI_Send(&ans, 1, MPI_FLOAT, 0, row, MPI_COMM_WORLD);
            }
            else{
                po = 0;
            }
        }

        free(buffer);
    }

    // clean up the MPI state in preparation for the processes to exit
    MPI_Finalize();

    // check parallel computing solution
    if(rank == 0){
        check_solution_validity(A,x,y,m,n,EPS);
    }

    return 0;
}

float **generate_matrix(int m, int n){
    /**
     * @brief generate m x n matrix
     * m : row
     * n : col
     */
    float **A = new float*[m];
    for (int i = 0; i < m; i++){
        A[i] = new float[n];
    }
    return A;
}

float *generate_array(int m){
    /**
     * @brief generate m x 1 array
     * m : row
     */
    float *B = new float[m];
    return B;
}

void init_random_array(float *x, int n_size, float min, float max){
    for(int i = 0; i < n_size; i++){
        x[i] = get_random_number(min, max);
    }
}

void init_random_matrix(float **A, int m, int n, float min, float max){
    // if abs(A[i][i]) > sum of A[i][j] with j != i, then jacobi method can be converged
    // check diagonally dominanted
    int abs_sum = 0;

    for (int i = 0; i < m; i++)
    {
        abs_sum = 0;
        for (int j = 0; j < n; j++)
        {   
            // initialize matrix with condition : diagonal dominant
            if(i!=j){
                A[i][j] = get_random_number(min, max);
                abs_sum += fabs(A[i][j]);
            }
        }
        A[i][i] = abs_sum + 1;
    }

}

float get_random_number(float min, float max){
    int fraction = int(max - min + 1.0);
    return float(min + (std::rand() % fraction));
}

void check_solution_validity(float **A, float *x, float *target, int m, int n, float tolerance)
{
    float *y = new float[n];
    for (int i = 0; i < m; i++)
    {
        y[i] = 0;
        for (int j = 0; j < n; j++)
        {
            y[i] += A[i][j] * x[j];
        }
    }

    for(int i = 0; i < m; i++){
        cout << "target[" << i << "] : " << target[i] << ", y[" << i << "] : " << y[i] << endl;
    }
    return;
}