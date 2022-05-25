/**
 * @file main.cpp
 * @author zinzinbin
 * @brief Example Code on lecture 11 : example code for FDM
 * @version 0.1
 * @date 2022-05-20
 *
 * How to execute
 * (1) mpic++ main.cpp -o main.out
 * (2) mpirun -np <number of processes> ./main.out
 */

#include <iostream>
#include <omp.h>
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

void FivePointStencil(float *A, float *X, float *Y, int M, int N, int iters){

    int idx;
    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            
        }
    }
}

void GaussSeidel(float *A, float *X, float *Y, int M, int N, int iters){

}

float **generate_matrix(int m, int n)
{
    float **A = new float *[m];
    for (int i = 0; i < m; i++)
    {
        A[i] = new float[n];
    }
    return A;
}

float *generate_array(int m)
{
    float *B = new float[m];
    for (int i = 0; i < m; i++)
    {
        B[i] = 0.0;
    }
    return B;
}

void check_solution_validity(float **A, float *B, float *x, int m, int n, float tolerance)
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
    float diff = CalculateNormDifference(y, B, m, 0);
    if (diff < tolerance)
    {
        cout << "solution x is valid, tol : " << diff << endl;
    }
    else
    {
        cout << "solution x is not valid, tol : " << diff << endl;
    }
    return;
}

float CalculateNormDifference(float *x1, float *x2, int n)
{
    float difference = 0;
    int num_threads = 0;
    int steps = 0;

    for (int i = 0; i < n; i++)
    {
        difference += powf(x1[i] - x2[i], 2.0);
    }
    
    difference = sqrt(difference / n); // sqrt(sum of (x1-x2)^2 / n)
    return difference;
}

int main(int argc, char *argv[])
{   
    int M = 64;
    int N = 64;
    int iters = 128;

    float *A = generate_array(M * N);
    float *x = generate_array(N);
    float *b = generate_array(M);

    float tolerance = 1e-8;

    FivePointStencil(A,x,b,M,N,iters);
    




    return 0;
}