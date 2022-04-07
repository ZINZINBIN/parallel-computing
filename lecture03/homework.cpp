/**
 * @file homework.cpp
 * @author KIM JINSU(김진수, 2019-27420)
 * @brief Parallelize a program to solve a linear equation Ax = b by using the Jacobi method
 * @date 2022-04-07
 * step 1. Test your program (for checking accuracy) with a linear equation 
 * step 2. Check your program's speedup with a large matrix.
 * step 3. Describe your effort to enhance efficiency.
 */

#include <iostream>
#include <omp.h>
#include <time.h>
#include <array>
#include <vector>
#include <chrono>
#include <math.h>
#define EPS 1e-12

void init_random_array(float *x, int n_size, float min, float max);
void init_random_matrix(float **A, int m, int n, float min, float max);
float get_random_number(float min, float max);
float** generate_matrix(int m, int n);
float* generate_array(int m);
void initiate(float *x, int n_size, int method);
void copy_array(float *x, float *x_copy, int n_size, int method) ;
void check_solution_validity(float **A, float *B, float *x, int m, int n, float tolerance) ;
float CalculateNormDifference(float *x1, float *x2, int n, int is_parallel);
float *JacobiMethodSerial(float **A, float *B, int m, int n, int max_iters, float tolerance);
float *JacobiMethodParallel(float **A, float *B, int m, int n, int max_iters, float tolerance) ;
void print_solution(float *x, int n_size);

using namespace std;

int main(void){

    // Parallelize a program to solve a linear equation Ax = b
    // solve Ax = B by using Jacobi method

    cout << "### step 0. setting ###" << endl;
    cout << endl;
    int THREADS_NUM = 128;
    omp_set_num_threads(THREADS_NUM);
    cout << THREADS_NUM << "-threads set" << endl;
    cout << omp_get_num_procs() << "-procs" << endl;
    cout << omp_get_max_threads() << "-max.threads" << endl;
    cout << omp_get_num_threads() << "-threads now" << endl;
    cout << omp_get_thread_num() << "-th thread now" << endl;
    cout <<endl;
    
    /*
    * step 1. test program for checking accuarcy
    * example : A = [[4,-1,-1],[-2,6,1],[-1,1,7]], B = [3,9,-6], x = [0.75, 1.5, -0.857]
    */

    float **A_example = generate_matrix(3,3);
    float *B_example = generate_array(3);
    A_example[0][0] = 4; A_example[0][1] = -1; A_example[0][2] = -1;
    A_example[1][0] = -2; A_example[1][1] = 6; A_example[1][2] = 1;
    A_example[2][0] = -1; A_example[2][1] = 1; A_example[2][2] = 7;
    B_example[0] = 3; B_example[1] = 9; B_example[2] = -6;

    cout << "### step 1. test program for checking accuracy ######" << endl;
    cout << endl;
    float *solution_example = JacobiMethodSerial(A_example, B_example, 3, 3, 32, EPS);
    print_solution(solution_example, 3);
    cout << endl;

    /*
     * step 2. test program's speedup with a large matrix
     * matrix size : 1024 x 1024 fixed, initalized with random numbers in range of (min = 1, max = 5)
     * max iteration : 64 epoch
     */

    cout << "### step 2. test program's speedup with a large matrix ###" << endl;
    cout << endl;

    float tolerance = EPS;
    int m = 1024;
    int n = 1024;
    int max_iter = 64;
    float **A = generate_matrix(m,n);
    float *B = generate_array(m);

    // initiate matrix and array with each component extracted from range of (min, max) randomly
    float min = 1.0;
    float max = 5.0;
    init_random_matrix(A, m, n, min, max);
    init_random_array(B, m, min, max);
    
    // Jacobi method with serial computing
    cout << "------------Serial jacobi method-------------"<<endl;
    float *solution = JacobiMethodSerial(A,B,m,n,max_iter,tolerance);
    check_solution_validity(A,B,solution,m,n,tolerance);
    
    // Jacobi method with parallel computing
    cout << "------------Parallel jacobi method-----------"<<endl;
    float *solution_parallel = JacobiMethodParallel(A,B,m,n,max_iter,tolerance);
    check_solution_validity(A,B,solution_parallel, m,n,tolerance);

    /*
     * Extra(1): check speedup with varying matrix size and fixing number of threads
     * # of threads : 128 fixed
     * matrix size : 4, 16, 64, 256, ...
     * max iteration : 64 epoch
     * tolerance : 1e-8(common)
     */

    cout << endl;
    cout << "------------Check Speedup : fixed # of threads / varying matrix size-----------" << endl;
    cout << endl;

    omp_set_num_threads(128);
    for (int i = 1; i < 8; i++)
    {
        m = int(pow(4, i));
        n = m;
        A = generate_matrix(m, n);
        B = generate_array(m);

        init_random_matrix(A, m, n, min, max);
        init_random_array(B, m, min, max);
        cout << "matrix size : " << m << ", ";
        solution_parallel = JacobiMethodParallel(A, B, m, n, max_iter, tolerance);
    }

    /*
     * Extra(2): check speedup with varying number of threads and fixing matrix size
     * # of threads : 2, 4, 8, 16, 32, ...
     * matrix size : 1024 fixed
     * max iteration : 64 epoch
     * tolerance : 1e-8(common)
     */

    cout << endl;
    cout << "------------Check Speedup : fixed matrix size / varying # of threads-----------" << endl;
    cout << endl;
    m = 1024;
    n = 1024;
    A = generate_matrix(m, n);
    B = generate_array(m);
    init_random_matrix(A, m, n, min, max);
    init_random_array(B, m, min, max);

    for (int i = 1; i < 8; i++)
    {
        int thread_num = int(pow(2, i));
        cout << "# of threads : " << thread_num << ", ";
        omp_set_num_threads(thread_num);
        solution_parallel = JacobiMethodParallel(A, B, m, n, max_iter, tolerance);
    }

    return 0;
}

void print_solution(float *x, int n_size){
    cout << "[";
    for (int i = 0; i < n_size-1; i++)
    {
        cout << x[i] << ",";
    }
    cout << x[n_size-1] << "]" << endl;
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

void initiate(float *x, int n_size, int method)
{
    // serial initiate : 0
    // parallel initiate : 1 or any number
    if (method == 0)
    {
        for (int i = 0; i < n_size; i++)
        {
            x[i] = 0;
        }
    }
    else
    {
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < n_size; i++)
        {
            x[i] = 0;
        }
    }
}

void init_random_array(float *x, int n_size, float min, float max){
    for(int i = 0; i < n_size; i++){
        x[i] = get_random_number(min, max);
    }
}

void init_random_matrix(float **A, int m, int n, float min, float max){
    for (int i = 0; i < m; i++){
        for(int j = 0; j<n;j++){
            A[i][j] = get_random_number(min, max);
        }
    }
}

float get_random_number(float min, float max){
    float fraction = 1.0 / (max - min + 1.0);
    return min + (max - min + 1.0) * (std::rand() * fraction);
}

void copy_array(float *x, float *x_copy, int n_size, int method)
{
    // serial initiate : 0
    // parallel initiate : 1 or any number
    if (method == 1)
    {
        for (int i = 0; i < n_size; i++)
        {
            x_copy[i] = x[i];
        }
    }
    else
    {
    #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < n_size; i++)
        {
            x_copy[i] = x[i];
        }
    }
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
    if (CalculateNormDifference(y, B, m, 0) < tolerance)
    {
        cout << "solution x is valid" << endl;
    }
    else
    {
        cout << "solution x is not valid" << endl;
    }
}

float CalculateNormDifference(float *x1, float *x2, int n, int is_parallel)
{
    /**
     * @brief Calculate the L2 Norm of |x1 - x2|
     * x1 : float pointer(array 1)
     * x2 : float pointer(array 2)
     * n : length of array(x1 and x2 should be equal)
     */
    float difference = 0;
    if(is_parallel){
        #pragma omp parallel for schedule(dynamic) reduction(+: difference)
        for (int i = 0; i < n; i++)
        {
            difference += pow(abs(x1[n] - x2[n]), 2);
        }
    }
    else{
        for (int i = 0; i < n; i++)
        {
            difference += pow(abs(x1[n] - x2[n]), 2);
        }
    }
    difference /= n;
    difference = pow(difference, 0.5); // sqrt(sum of (x1-x2)^2 / n)
    return difference;
}

float* JacobiMethodSerial(float **A, float *B, int m, int n, int max_iters, float tolerance)
{
    /**
     * @brief Solve Ax = B with Jacobi Method without parallel computing
     * A : dim(m x n), matrix with float type
     * B : dim(m x 1), array with float type
     * m : row dim
     * n : col dim
     * max_iters : maximum iteration
     * tolerance : criteria for checking convergence
     */
    float *x = new float[n];
    float *x_new = new float[n];
    float sum = 0;
    int iteration = 0;

    // time check : start
    auto start = std::chrono::high_resolution_clock::now();

    initiate(x, n, 0);
    initiate(x_new, n, 0);

    // Jacobi method
    for (int n_iter = 0; n_iter < max_iters; n_iter++)
    {
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                sum += A[i][j] * x[j];
            }
            x_new[i] = B[i] - sum;
            x_new[i] /= A[i][i];
        }

        iteration ++;

        if (CalculateNormDifference(x_new, x, n, 0) < tolerance)
        {   
            break;
        }
        else
        {
            copy_array(x_new, x, n, 0);
        }
    }
    // time check : end
    auto end = std::chrono::high_resolution_clock::now();

    double elapsed_time_ms = std::chrono::duration<double, std::milli>(end - start).count();
    cout << "Jacobi method runtime : " << elapsed_time_ms << ", iters for convergence : " << iteration << endl;

    return x_new;
};

float* JacobiMethodParallel(float **A, float *B, int m, int n, int max_iters, float tolerance)
{
    /**
     * @brief Solve Ax = B with Jacobi Method with parallel computing
     * A : dim(m x n), matrix with float type
     * B : dim(m x 1), array with float type
     * m : row dim
     * n : col dim
     * max_iters : maximum iteration
     * tolerance : criteria for checking convergence
     * Iteration간에는 parallel은 사용하지 않는다.
     * 단, row마다 Jacobi method를 이용해 x_new를 update 하는 과정에서 parallel computing을 적용
     * row : m이라면 schuduler를 이용해 dynamic하게 thread마다 할당하여 연산을 진행
     * 이후 Join하여 새로 업데이트된 x_new 값에 대한 tolerance를 확인
     */

    float *x = new float[n];
    float *x_new = new float[n];
    float sum = 0;
    int iteration = 0;

    // time check : start
    auto start = std::chrono::high_resolution_clock::now();

    initiate(x, n, 1);
    initiate(x_new, n, 1);

    // Jacobi method
    for (int n_iter = 0; n_iter < max_iters; n_iter++)
    {
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                sum += A[i][j] * x[j];
            }
            x_new[i] = B[i] - sum;
            x_new[i] /= A[i][i];
        }

        iteration ++;

        if (CalculateNormDifference(x_new, x, n, 1) < tolerance)
        {
            break;
        }
        else
        {
            copy_array(x_new, x, n, 1);
        }
    }
    // time check : end
    auto end = std::chrono::high_resolution_clock::now();

    double elapsed_time_ms = std::chrono::duration<double, std::milli>(end - start).count();
    cout << "Jacobi method runtime : " << elapsed_time_ms << ", iters for convergence : " << iteration << endl;

    return x_new;
};
