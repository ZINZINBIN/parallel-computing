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

float **generate_matrix(int m, int n);
float *generate_arrray(int m);
float calculate_component(int i, int j, int k, int l);
void initiate_matrix(float **A, int m, int n);

int main(int argc, char *argv[]){

    int size, rank;
    int m = 128;
    int n = 128;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    float **A = generate_matrix(m, n);
    initiate_matrix(A, m, n);

    if(rank == 0){
        cout << "matrix A : (" << m << "," << n << ")" << endl;
    }

    // tranpose algorithm


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

float calculate_component(int i, int j, int k, int l){
    float comp = 0
    comp = (float) 1.0 / (sqrt((i-k) * (i-k) + (j-l) * (j-l)) + EPS);
    return comp;
}

void initiate_matrix(float **A, int m, int n){

    float sum = 0;
    for(int i = 0; i < m; i++){
        for(int j = 0; j<n; j++){
            sum = 0;
            for(int k = 0; k < n; k++){
                for(int l = (-1) * 5 * n; l < 5 * n ; l++){
                    sum += calculate_component(i,j,k,l);
                }
            }
            A[i][j] = sum;
        }
    }
}