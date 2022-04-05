#include <iostream>
#include <omp.h>
#include <array>

using namespace std;

int main(void){

    int sum;
    int N = 128;
    int M = 128;
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

    // matrix vector multiplication
    /*
    #pragma omp parallel
    {
        #pragma omp for private(i, k, sum)
        for(i = 0; i < N; i ++){
            sum = 0.0;
            for(k = 0; k < M; k++){
                sum += A[i][k] * x[k];
            }
            y[i] = sum;
        }
    }
    */
    
    return 0;
}