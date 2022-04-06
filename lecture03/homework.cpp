/*
# Assignment
(1) Sections example 만들기
- 어느 thread가 어느 부분을 맡았는지 확인해보기. 여러번 실행해보기
(2) Matrix Vector Multiplication
- 곱한 값을 아는 행렬과 벡터로 시험해보기
- 큰 규모의 행렬과 벡터에 대해 schedule clause 적용해보고 어느 thread가 어느 부분을 맡았는지 확인해보기
(3) 첫번째 과제 2주 후 마감. 3월 25일 과제 오픈시 2주 후에 마감할 예정
*/

#include <iostream>
#include <omp.h>
#include <time.h>
#include <array>
#include <vector>
#include <chrono>
#include <math.h>
#define EPS 1e-8

void init_random_array(float *x, int n_size, float min, float max);
void init_random_matrix(float **A, int m, int n, float min, float max);
float get_random_number(float min, float max);
float** generate_matrix(int m, int n);
float* generate_array(int m);
void initiate(float *x, int n_size, int method);
void copy_array(float *x, float *x_copy, int n_size, int method) ;
void check_solution_validity(float **A, float *B, float *x, int m, int n, float tolerance) ;
float CalculateNormDifference(float *x1, float *x2, int n);
float *JacobiMethodSerial(float **A, float *B, int m, int n, int max_iters, float tolerance);
float *JacobiMethodParallel(float **A, float *B, int m, int n, int max_iters, float tolerance) ;
void print_solution(float *x, int n_size);

using namespace std;

int main(void){

    int THREADS_NUM = 128;
    omp_set_num_threads(THREADS_NUM);
    cout << THREADS_NUM << "-threads set" << endl;
    cout << omp_get_num_procs() << "-procs" << endl;
    cout << omp_get_max_threads() << "-max.threads" << endl;
    cout << omp_get_num_threads() << "-threads now" << endl;
    cout << omp_get_thread_num() << "-th thread now" << endl;

    // Parallelize a program to solve a linear equation Ax = b 
    // solve Ax = B by using Jacobi method
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

    // cout << "# solution for serial jacobi method"<<endl;
    // print_solution(solution, m);
    
    // Jacobi method with parallel computing
    cout << "------------Parallel jacobi method-----------"<<endl;
    float *solution_parallel = JacobiMethodParallel(A,B,m,n,max_iter,tolerance);
    check_solution_validity(A,B,solution_parallel, m,n,tolerance);

    // cout << "# solution for parallel jacobi method" << endl;
    // print_solution(solution_parallel, m);


    // Check speedup with a large matrix
    // process 1 : fixed number of threads(=128), varying matrix size(m = 4,16,64,256,1024,...)
    cout << "------------Check Speedup : fixed # of threads with varying matrix size-----------"<<endl;

    omp_set_num_threads(128);
    for(int i = 1; i < 8; i++){
        m = int(pow(4,i));
        n = m;
        A = generate_matrix(m,n);
        B = generate_array(m);

        init_random_matrix(A, m, n, min, max);
        init_random_array(B,m,min,max);
        cout << "matrix size : " << m << ", ";
        solution_parallel = JacobiMethodParallel(A,B,m,n,max_iter,tolerance);
    }
    
    // process 2 : fixed matrix size, varying number of threads
    cout << "------------Check Speedup : fixed matrix size with varying # of threads-----------"<<endl;
    m = 1024;
    n = 1024;
    A = generate_matrix(m,n);
    B = generate_array(m);
    init_random_matrix(A, m, n, min, max);
    init_random_array(B,m,min,max);

    for(int i = 1; i < 10; i++){
        int thread_num = int(pow(2,i));
        cout << "# of threads : " << thread_num << ", ";
        omp_set_num_threads(thread_num);
        solution_parallel = JacobiMethodParallel(A,B,m,n,max_iter,tolerance);
    }
    return 0;
}

void print_solution(float *x, int n_size){
    /*
    for (int i = 0; i < n_size; i++)
    {
        cout << "x[" << i << "] : " << x[i] << endl;
    }
    */
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
    if (CalculateNormDifference(y, B, m) < tolerance)
    {
        cout << "solution x is valid" << endl;
    }
    else
    {
        cout << "solution x is not valid" << endl;
    }
}

float CalculateNormDifference(float *x1, float *x2, int n)
{
    /**
     * @brief Calculate the L2 Norm of |x1 - x2|
     * x1 : float pointer(array 1)
     * x2 : float pointer(array 2)
     * n : length of array(x1 and x2 should be equal)
     */
    float difference = 0;
    #pragma omp parallel for schedule(dynamic, 1) reduction(+:difference)
    for (int i = 0; i < n; i++)
    {
        difference += pow(abs(x1[n] - x2[n]), 2);
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

        if (CalculateNormDifference(x_new, x, n) < tolerance)
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
    cout << "Jacobi method runtime : " << elapsed_time_ms << endl;

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

        if (CalculateNormDifference(x_new, x, n) < tolerance)
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
    cout << "Jacobi method runtime : " << elapsed_time_ms << endl;

    return x_new;
};
