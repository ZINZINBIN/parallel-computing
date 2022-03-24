#include <iostream>
#include <omp.h>

using namespace std;

int main(void){

    float sum;
    int N = 128;
    int M = 128;
    int i,k;

    int x[M] = {};

    #pragma omp parallel for(int i = 0; i < count; i++)
    {
        /* code */
    }

    // matrix vector multiplication
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
    
    return 0;
}