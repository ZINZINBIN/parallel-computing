#include <iostream>
#include <omp.h>
#include <array>

using namespace std;

float **generate_matrix(int m, int n){
    float **A = new float*[m];
    for (int i = 0; i < m; i++){
        A[i] = new float[n];
    }
    return A;
}

void init_matrix(float **A, int m, int n){
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {   
            A[i][j] = 0;
        }
    }
}

struct MyType
{
    int thread_num;
};

MyType threaded_ver;
#pragma omp threadprivate(threaded_ver)
int main(void){

    int sum;
    int N = 32;
    int M = 32;
    int THREADS_NUM = 16;
    int i,k;
    int count = 256;
    array<int, 256> x = {0};
    double start, end;

    cout << THREADS_NUM << "-threads set by me" << endl;
    omp_set_num_threads(THREADS_NUM);
    cout << omp_get_num_procs() << "-procs" << endl; 
    cout << omp_get_max_threads() << "-max.threads" << endl;
    cout << omp_get_num_threads() << "-threads now" << endl;
    cout << omp_get_thread_num() << "-th thread now" << endl;
    cout << "check this region parallel : " << omp_in_parallel() << endl;
    cout << "Fork" << endl;

    for(int i = 0; i < count; i++){
        x[i] = i;
    }

    int idx = 10;

    cout << "private and first private test" << endl;

    #pragma omp parallel private(idx)
    {
        #pragma omp critical(print_first)
        {
            cout << "thread : " << omp_get_thread_num() << " - idx : " << idx << endl;
            idx = 1000 + omp_get_thread_num();
        }
    }

    cout << "idx : " << idx << endl;
    cout << "using firstprivate" << endl;

    #pragma omp parallel firstprivate(idx)
    {
        #pragma omp critical(print_first)
        {
            cout << "thread : " << omp_get_thread_num() << " - idx : " << idx << endl;
            idx = 1000 + omp_get_thread_num();
        }
    }

    cout << "idx : " << idx << endl;

    // collapse clause
    float **a = generate_matrix(M,N);
    init_matrix(a, M, N);

    #pragma omp parallel for collapse(2)
    for(int i = 0; i < M; i++){
        for(int j = 0; j< N; j++){
            a[i][j] *= 2.0;
        }
    }

    // sections
    #pragma omp parallel
    {
        #pragma omp sections
        {
            #pragma omp section
            {

            }
            #pragma omp section
            {

            }
        }
    }

    // task
    int op = 1;
    int idx_i = 1;
    int idx_j = 0;
    #pragma omp parallel
    {
        #pragma omp master
        {
            #pragma omp critical(print_task)
            {
                cout << "master thread : " << omp_get_thread_num() << endl;
            }

            while(op && idx_i < 50)
            {
                #pragma omp task firstprivate(idx_i) private(idx_j)
                {
                    for(idx_j = 0; idx_j < 128; idx_j++){
                        if(idx_j == 64 && idx_i == 32){
                            op = 0;
                        }
                    }
                }
                idx_i++;
            }
        }
    }
    
    /*
    cout << "static schedule" << endl;

    #pragma omp parallel reduction(+:sum) private(start, end)
    {
        start = omp_get_wtime();
        #pragma omp for schedule(static, 64) 
        for(int i = 0; i < count; i++)
        {
            sum += x[i];
            #pragma omp critical(printing)
            {
                cout << omp_get_thread_num() << "-thread, " << i << "-th step, "<< sum << ": sum " << endl; 
            }
        }
        end = omp_get_wtime();

        #pragma omp critical(time_check)
        {
            cout << "static schedule time spend : " << end - start << endl;
        }
    }

    cout << "fork" << endl;
    cout << "parallel sum : " << sum << " answer : " << int(count * (count - 1) / 2) << endl;
    cout << "dynamic schedule" << endl;
    
    #pragma omp parallel reduction(+:sum) private(start, end)
    {   
        start = omp_get_wtime();
        #pragma omp for schedule(dynamic, 64) 
        for(int i = 0; i < count; i++)
        {
            sum += x[i];
            #pragma omp critical(printing)
            {
                cout << omp_get_thread_num() << "-thread, " << i << "-th step, "<< sum << ": sum " << endl; 
            }
        }
        end = omp_get_wtime();

        #pragma omp critical(time_check)
        {
            cout << "static schedule time spend : " << end - start << endl;
        }
    }

    cout << "fork" << endl;
    cout << "parallel sum : " << sum << " answer : " << int(count * (count - 1) / 2) << endl;

    // section clause
    #pragma omp parallel
    {
        #pragma omp sections
        {
            #pragma omp section
            {

            }
            #pragma omp section
            {

            }
        }
    }
    */

    return 0;
}