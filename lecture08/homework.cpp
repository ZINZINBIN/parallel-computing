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
#define MAXPROC 16

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
    int N = 8;
    float eps = EPS;
    double p_time = 0; // time for execution
    double p_time_avg = 0; // average time for execution
    double c_time = 0; // time for communication

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // generate matrix A : m x n
    int m = N;
    int n = N;
    float **A = generate_matrix(m, n);
    initiate_matrix(A, m, n);

    if (rank == 0){
        cout << "matrix A : (" << m << "," << n << ")" << endl;
        print_matrix(A,m,n);
    }

    // 1D array of A
    float *A_1D = convertMat2Array(A, m, n);

    // transpose of A : real value
    float **A_T = transpose_matrix(A, m, n);

    // matrix a_t : transpose of A using MPI collective communication
    float **a_t = generate_matrix(n,m);
    float *a_t_1D = convertMat2Array(a_t, n, m);
    
    // tranpose algorithm
    if(rank == 0){
        c_time -= MPI_Wtime();
        p_time -= MPI_Wtime();
    }

    int row_per_proc;
    int rest;

    if(m >= size){
        row_per_proc = int(m / size);
    }
    else{
        row_per_proc = 1;
    }

    rest = m % size;

    // variable used for parallel computing
    float *recv_buffer = generate_array(n * row_per_proc + n);
    float *send_buffer = generate_array(n * row_per_proc + n);
    float *shuffle_buffer = generate_array(n * row_per_proc + n);

    // define send counts and displacement
    int *send_counts = new int[size];
    int *recv_counts = new int[size];
    int *send_displs = new int[size];
    int *recv_displs = new int[size];
    int recv_sum = 0;
    int send_sum = 0;

    int sum = 0;

    for(int i = 0; i < size; i++){
        send_counts[i] = row_per_proc * n;
        if(rest >0){
            rest --;
            send_counts[i] += n;
        }
        send_displs[i] = send_sum;
        send_sum += send_counts[i];
    }

    MPI_Scatterv(A_1D, send_counts, send_displs, MPI_FLOAT, recv_buffer, row_per_proc * n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    // MPI_Scatter(A_1D, row_per_proc * n, MPI_FLOAT, recv_buffer, row_per_proc * n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < row_per_proc; i++)
    {
        for (int j = 0; j < n; j++)
        {
            send_buffer[j * row_per_proc + i] = recv_buffer[i * n + j];
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
    {
        cout << "rank 0" << endl;
        cout << "recv buffer " << endl;
        for (int i = 0; i < row_per_proc * n; i++)
        {
            printf("%3.3f \t", recv_buffer[i]);
        }
        cout << endl;
        cout << "send buffer " << endl;
        for (int i = 0; i < row_per_proc * n; i++)
        {
            printf("%3.3f \t", send_buffer[i]);
        }
        cout << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // define recv counts and displacement
    send_sum = 0;
    recv_sum = 0;
    rest = m % size;

    for(int i = 0; i < size; i++){
        send_counts[i] = row_per_proc;
        recv_counts[i] = row_per_proc * n;

        if(rest >0){
            rest --;
            send_counts[i] += 1;
            recv_counts[i] += m % size;
        }
        send_displs[i] += send_sum;
        recv_displs[i] += recv_sum;
        recv_sum += recv_counts[i];
        send_sum += send_counts[i];
    }

    MPI_Alltoallv(send_buffer, send_counts, send_displs, MPI_FLOAT, a_t_1D, recv_counts, recv_displs, MPI_FLOAT, MPI_COMM_WORLD);
    
    // MPI_Alltoall(send_buffer, row_per_proc, MPI_FLOAT, shuffle_buffer, row_per_proc, MPI_FLOAT, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);

    // if (rank == 0)
    // {
    //     cout << "shuffle buffer " << endl;
    //     for (int i = 0; i < row_per_proc * n; i++)
    //     {
    //         printf("%3.3f \t", shuffle_buffer[i]);
    //     }
    //     cout << endl;
    // }

    //MPI_Gather(shuffle_buffer, row_per_proc * n, MPI_FLOAT, a_t_1D, row_per_proc * n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // define send counts and displacement
    // int *recv_counts = new int[size];
    // int *displs = new int[size];
    // int sum = 0;

    // MPI_Gatherv(send_buffer, row_per_proc * n, MPI_FLOAT, a_t_1D, row_per_proc * n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // MPI_Barrier(MPI_COMM_WORLD);

    a_t = convertArray2Mat(a_t_1D, m, n);

    // check the coincidence
    if(rank == 0){
        int is_coincidence = check_coincidency(A_T, a_t, m, n, eps);
        c_time += MPI_Wtime();
        p_time += MPI_Wtime();

        if(is_coincidence==1){
            cout << "transpose of matrix A complete" << endl;
            cout << "p_time : " << p_time << ", c_time : " << c_time << endl; 
        }
        else{
            cout << "transpose of matrix A : error, wrong result" << endl;
            cout << "p_time : " << p_time << ", c_time : " << c_time << endl;
        }
    }

    if (rank == 0)
    {
        print_matrix(a_t, n, m);
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
    float comp = (float) 1.0 / (sqrt((i-k) * (i-k) + (j-l) * (j-l)) + EPS);
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
    int is_coincidence;

    /*
    float loss = 0;
    float total_loss = 0;

    for(int i = 0; i < m ; i++){
        for(int j = 0; j < n; j++){
            loss = abs((A[i][j] - B[i][j]) / (A[i][j] + eps));
            total_loss += loss * loss
        }
    }

    if(total_loss < eps){
        is_coincidence = 1;
    }
    else{
        is_coincidence = 0;
    }
    */

   is_coincidence = 1;
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