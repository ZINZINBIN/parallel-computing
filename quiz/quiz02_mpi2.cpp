/**
 * How to execute
 * (1) mpic++ quiz02_mpi2.cpp -o quiz02_mpi2.out
 * (2) mpirun -np <number of processes> ./quiz02_mpi2.out
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

using namespace std;

// function for use
void init_random_array(float *x, int n_size, float min, float max);
void init_random_matrix(float **A, int m, int n, float min, float max);
float get_random_number(float min, float max);
float** generate_matrix(int m, int n);
float* generate_array(int m);
void check_solution_validity(float **A, float *x, float *target, int m, int n, float tolerance);


int main(int argc, char *argv[]){

    int baton = 1;
    int MAX_PROC = 16;
    int rank, size; 
    int tag;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size); // get number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get current process id

    MPI_Status status;


    // Ring example
    if(rank == 0){
        baton = 1;
        cout << "Ring Example" << endl;
        MPI_Send(&baton, 1, MPI_INT, 1, 999, MPI_COMM_WORLD);
        MPI_Recv(&baton, 1, MPI_INT, size - 1, 999, MPI_COMM_WORLD, &status);

        cout << "baton : " << baton << ", process : " << size - 1 << " to process : 0" << endl;
    }   
    else{
        MPI_Recv(&baton, 1, MPI_INT, rank - 1, 999, MPI_COMM_WORLD, &status);
        cout << "baton : " << baton << ", process : " << rank - 1 << " to process : " << rank << endl;

        MPI_Send(&baton, 1, MPI_INT, (rank + 1) % size, 999, MPI_COMM_WORLD);
    }

    int m = 1024;
    float **A = generate_matrix(m,m);
    float *x = generate_array(m);
    float *y = generate_array(m);
    
    int n_rows = m;
    int n_cols = m;
    int numprocs = size;
    int numsent = 0;
    int sender;
    int row;
    int po;
    float ans;
    float *buffer;

    double t_start;
    double t_end;

    // Matrix - vector multiplication
    if(rank == 0){

        t_start = MPI_Wtime();

        numsent = 0;
        MPI_Bcast(&x[0], n_cols, MPI_FLOAT, 0, MPI_COMM_WORLD);

        for(int i = 0; i < numprocs - 1 && i < n_rows; i++){
            MPI_Send(&A[i][0], n_cols, MPI_FLOAT, i + 1, i, MPI_COMM_WORLD);
            numsent ++;
        }

        for(int i = 0; i < n_rows; i++){
            // 각각의 process에서 진행한 연산을 모두 receive
            MPI_Recv(&ans, 1, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            sender = status.MPI_SOURCE;
            row = status.MPI_TAG;

            y[row] = ans;

            // numsent가 n_rows보다 작은 경우 numsent를 증가시키면서 send를 보낸다. 
            if(numsent < n_rows){
                MPI_Send(&A[numsent][0], n_cols, MPI_FLOAT, sender, numsent, MPI_COMM_WORLD);
                numsent ++;
            }
            else{
                MPI_Send(MPI_BOTTOM, 0, MPI_FLOAT, sender, n_rows, MPI_COMM_WORLD);
            }
        }

        t_end = MPI_Wtime();

        cout << "matrix multiplication time spend : " << t_end - t_start << endl;

    }
    else{
        MPI_Bcast(&x[0], n_cols, MPI_FLOAT, 0, MPI_COMM_WORLD); // rank = 0이 아닌 경우에는 따로 진행 X
        buffer = new float[n_cols];
        po = 1;

        while(rank <= n_rows && po){
            MPI_Recv(buffer, n_cols, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            row = status.MPI_TAG;

            if(row < n_rows){
                ans = 0;
                for(int j = 0; j < n_cols; j++){
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

    MPI_Finalize();
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