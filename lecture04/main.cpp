#include <iostream>
#include <omp.h>
#include <chrono>
#include <time.h>

using namespace std;

int main(void){
    int THREADS_NUM = 8;
    omp_set_num_threads(THREADS_NUM);
    cout << THREADS_NUM << "-threads set" << endl;
    cout << omp_get_num_procs() << "-procs" << endl;
    cout << omp_get_max_threads() << "-max.threads" << endl;
    cout << omp_get_num_threads() << "-threads now" << endl;
    cout << omp_get_thread_num() << "-th thread now" << endl;
    cout <<endl;

    int n = 10000000;

    int *a = new int[n];
    int *b = new int[n];
    int *c = new int[n];

    // initialize 
    for(int i = 0; i < n; i++){
        a[i] = i;
        b[i] = i;
        c[i] = i;
    }

    // serial example
    auto start = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < n; i++){
        b[i] += a[i-1];
        a[i] += c[i];
    }
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed_time_ms = std::chrono::duration<double, std::milli>(end - start).count();
    cout << "serial example runtime : " << elapsed_time_ms << endl;

    // parallel version
    start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for
    for(int i = 0; i < n; i++){
        a[i] += c[i];
        b[i+1] += a[i];
    }
    end = std::chrono::high_resolution_clock::now();
    elapsed_time_ms = std::chrono::duration<double, std::milli>(end - start).count();
    cout << "parallel example runtime : " << elapsed_time_ms << endl;


    int m = 1024;
    n = m;
    float **A = new float*[m];
    float *x = new float[n];
    float *y = new float[m];
    
    for(int i = 0; i < m; i++){
        A[i] = new float[n];
        for(int j = 0; j < n; j++){
            A[i][j] = 1;
        }
    }

    float sum;
    double start_time;
    double end_time;

    #pragma omp parllel private(start, end)
    {
        #pragma omp master
        start_time = omp_get_wtime();
        
        #pragma omp barrer

        #pragma omp for
        for(int i = 0; i < m; i++){
            sum = 0;
            for(int j = 0; j<n;j++){
                sum += A[i][j] * x[j];
            }
            y[i] = sum;
        }
        
        #pragma omp maste
        {
            end_time = omp_get_wtime();

            #pragma omp critical
            {
                cout << "better access, parallel runtime : " << end_time - start_time << endl;
            }
        }
        
    }

    #pragma omp parllel private(start, end)
    {
        #pragma omp master
        start_time = omp_get_wtime();
        
        #pragma omp barrer

        #pragma omp for
        for(int i = 0; i < m; i++){
            sum = 0;
            for(int j = 0; j<n;j++){
                sum += x[j] * A[j][i];
            }
            y[i] = sum;
        }
        
        #pragma omp maste
        {
            end_time = omp_get_wtime();

            #pragma omp critical
            {
                cout << "worse access, parallel runtime : " << end_time - start_time << endl;
            }
        }
        
    }   
    /*
    // task direcive example
    #pragma omp parallel shared(a, coord , op,
    {
        #pragma omp master
        {
            while (op && i <50){
                #pragma omp task firstprivate(i) private(j)
                {
                    for (j=0; j<10000; j++){

                    }
                }
            }
            
        {
         {
        if (a[i][j]==target){
            coord[0] = i ; coord [1] = j; op = 0;        
            }
        
        
        }
        i
        }
    }
    */

    return 0;
}

float parallel_integration(float x_min, float x_max, int n_size){
    float sum = 0;
    float x;
    int i;
    float dx = (x_max - x_min) / n_size;

    #pragma omp parallel private(i,x)
    {
        #pragma omp for reduction(+:sum)
        for(i = 0; i < n_size ; i++){
            x = x_min + float(i * dx);
            sum += 4.0 / (1.0 + x * x);
        }
    }
    return sum;
}