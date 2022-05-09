/**
 * @file homework.cpp
 * @author zinzinbin
 * @brief homework 2 
 * @version 0.1
 * @date 2022-05-03
 *
 * How to execute
 * (1) mpic++ homework.cpp -o homework.out
 * (2) mpirun -np <number of processes> ./homework.out
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

// custom function 
float **generate_matrix(int m, int n);
float *generate_array(int m);
float **convertArray2Mat(float *B, int m, int n);
float *convertMat2Array(float **A, int m, int n);
float calculate_component(int i, int j, int k, int l);
void initiate_matrix(float **A, int m, int n);
float **transpose_matrix(float **A, int m, int n);
int check_coincidency(float **A, float **B, int m, int n, float eps);
void print_matrix(float **A, int m, int n);

int main(int argc, char *argv[]){

    int size, rank;
    int N = 256;
    float eps = EPS;
    double p_time = 0; // time for execution
    double c_time = 0; // time for communication
    double exec_time = 0; // time for total execution

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // generate matrix A : m x n
    int m = N;
    int n = N;
    float **A = generate_matrix(m, n);
    initiate_matrix(A, m, n);

    if(m%size != 0 || n%size !=0){
        if(rank == 0){
            cout << "Error : dimensions of matrix should be divided by number of processes!" << endl;
        }
        MPI_Finalize();
        return -1;
    }
    else{
        if (rank == 0){
            cout << "matrix A : (" << m << "," << n << ")" << endl;
            cout << "nprocs : " << size << endl;
            if(m < 16 && n < 16){
                print_matrix(A,m,n);
            }
        }
    }

    // convert matrix A to 1D array
    float *A_1D = convertMat2Array(A, m, n);

    // transpose of A without MPI: real value
    // time check

    auto start = std::chrono::high_resolution_clock::now();

    float **A_T = transpose_matrix(A, m, n);

    auto end = std::chrono::high_resolution_clock::now();
    double exec_time_without_mpi = std::chrono::duration<double>(end - start).count();

    if(rank == 0){
        // printf("transposition without MPI(s) : %3.6f \n", exec_time_without_mpi);
        cout << endl;
        cout << "-------------------------- Matrix Transposition with MPI ---------------------------" << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // matrix a_t : transpose of A using MPI collective communication
    float **a_t = generate_matrix(n,m);
    float *a_t_1D = generate_array(m * n);

    int row_per_procs = m / size;
    int col_per_procs = n / size;

    float *A_procs = generate_array(row_per_procs * n);
    float *A_procs_t01 = generate_array(row_per_procs * n);
    float *A_procs_t02 = generate_array(col_per_procs * m);
    float *A_procs_t03 = generate_array(row_per_procs * n);

    // tranpose algorithm
    if(rank == 0){
        exec_time -= MPI_Wtime();
    }

    int row_procs_i = rank * row_per_procs;
    int row_procs_f = (rank + 1) * row_per_procs;

    if(rank == 0){
        p_time -= MPI_Wtime();
    }

    for(int i = 0; i < row_per_procs; i++){
        for(int proc = 0; proc < size; proc ++){
            for(int j = 0; j < col_per_procs; j++){
                A_procs[row_per_procs * col_per_procs * proc + col_per_procs * i + j] = A_1D[row_procs_i * n + proc * row_per_procs * col_per_procs + i * col_per_procs + j];
            }
        }
    }

    if(rank == 0){
        p_time += MPI_Wtime();
        printf("# p_time : %3.6f, process : generate A_procs\n", p_time);
        p_time = 0;
    }

    if(rank == 0){
        p_time -= MPI_Wtime();
    }

    // transpose process 1 : local transpose
    for(int proc = 0; proc < size; proc ++){
        for(int i = 0; i < row_per_procs; i++){
            for(int j = 0; j < col_per_procs; j++){
                A_procs_t01[row_per_procs * col_per_procs * proc + row_per_procs * j + i] = A_procs[col_per_procs * size * i + proc * col_per_procs + j];
            }
        }
    }

    // MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0){
        p_time += MPI_Wtime();
        printf("# p_time : %3.6f, process : transpose 01\n", p_time);
        p_time = 0;
    }

    if(rank == 0){
        c_time -= MPI_Wtime();
    }

    // transpose process 2 : AlltoAll, rearangement of block matrix
    MPI_Alltoall(A_procs_t01, row_per_procs * col_per_procs, MPI_FLOAT, A_procs_t02, row_per_procs * col_per_procs, MPI_FLOAT, MPI_COMM_WORLD);

    if(rank == 0){
        c_time += MPI_Wtime();
        printf("# c_time : %3.6f, process : MPI_Alltoall\n", c_time);
        c_time = 0;
    }

    if(rank == 0){
        p_time -= MPI_Wtime();
    }

    // transpose process 3 : rearange each block matrix with correct order
    for(int proc = 0; proc < size; proc ++){
        for(int j = 0; j < col_per_procs; j++){
            for(int i = 0; i < row_per_procs; i++){
                A_procs_t03[size * row_per_procs *j + proc * row_per_procs + i] = A_procs_t02[proc * row_per_procs * col_per_procs + j * row_per_procs + i];
            }
        }
    }

    // MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0){
        p_time += MPI_Wtime();
        printf("# p_time : %3.6f, process : transpose 03\n", p_time);
        p_time = 0;
    }

    if(rank == 0){
        c_time -= MPI_Wtime();
    }

    // transpose process 4
    MPI_Gather(A_procs_t03, row_per_procs * n, MPI_FLOAT, a_t_1D, row_per_procs * n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if(rank == 0){
        exec_time += MPI_Wtime();
        c_time += MPI_Wtime();
        printf("# c_time : %3.6f, process : MPI_Gather\n", c_time);
        c_time = 0;
    }

    // result
    a_t = convertArray2Mat(a_t_1D, n, m);

    // check the coincidence
    if(rank == 0){
        cout << "------------------------------------------------------------------------------------" << endl;
        cout << endl;

        if(check_coincidency(A_T, a_t, m, n, eps)){
            cout <<"transposition of matrix A complete, all components are same" << endl;
            printf("exec_time : %3.6f, transposition without MPI(s)\n", exec_time_without_mpi);
            printf("exec_time : %3.6f, transposition with MPI(s)\n", exec_time);
        }
        else{
            cout << "transpose of matrix A : error, wrong result" << endl;
            printf("exec_time : %3.6f, transposition without MPI(s)\n", exec_time_without_mpi);
            printf("exec_time : %3.6f, transposition with MPI(s)\n", exec_time);
        }
    }

    if (rank == 0 && m < 16 && n < 16)
    {
        print_matrix(a_t, n, m);
    }

    MPI_Finalize();

    delete A;
    delete A_T;
    delete A_1D;
    delete A_procs;
    delete A_procs_t01;
    delete A_procs_t02;
    delete A_procs_t03;
    delete a_t_1D;
    delete a_t;

    return 0;
}

float **generate_matrix(int m, int n){
    float **A = new float*[m];
    for (int i = 0; i < m; i++){
        A[i] = new float[n];
    }
    return A;
}

float *generate_array(int m){
    float *B = new float[m];
    for(int i = 0; i < m; i ++){
        B[i] = 0.0;
    }
    return B;
}

float *convertMat2Array(float **A, int m, int n){
   int array_size = m * n;
   int idx;
   float *ptr = generate_array(array_size);

    for(int i = 0; i < m; i++){
        for(int j = 0; j<n; j++){
            idx = i * n + j;
            ptr[idx] = A[i][j];
        }
    }

    return ptr;
}

float **convertArray2Mat(float *B, int m, int n){
    int idx;
    float **ptr = generate_matrix(m,n);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            idx = i * n + j;
            ptr[i][j] = B[idx];
        }
    }
    return ptr;
}

float calculate_component(int i, int j, int k, int l){
    float comp = (float) 1.0 / (sqrt((i-k) * (i-k) + (j-l) * (j-l)));
    return comp;
}

void initiate_matrix(float **A, int m, int n){

    float sum = 0;
    for(int i = 0; i < m; i++){
        for(int j = 0; j<n; j++){
            sum = 0;
            for(int k = 0; k < n; k++){
                for(int l = (-1) * 5 * n; l < 5 * n ; l++){
                    if(k != i && l != j){
                        sum += calculate_component(i, j, k, l);
                    }
                }
            }
            A[i][j] = sum;
        }
    }
}

float **transpose_matrix(float **A, int m, int n){

    float **A_t = generate_matrix(n,m);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            A_t[i][j] = A[j][i];
        }
    }
    return A_t;
}

int check_coincidency(float **A, float **B, int m, int n, float eps){
    int is_coincidence = 1;
   float loss;

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            loss = abs((A[i][j] - B[i][j]) / (A[i][j] + eps));
            
            if(loss > eps){
                is_coincidence = 0;
            }
        }
    }

    return is_coincidence;
}

void print_matrix(float **A, int m, int n){
    cout << "matrix configuration" << endl;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%3.3f \t", A[i][j]);
        }
        cout << endl;
    }
}