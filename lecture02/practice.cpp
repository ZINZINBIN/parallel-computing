#include <iostream>
#include <omp.h>
#include "mimic_omp.h"
// #ifdef _OPEMMP
// id = omp_get_thread_num();
// #endif

using namespace std;

int main(void){
    int id, N;
    N = 16;
    float fraction;
    int total_count;

    // Library Routines는 serial computing에서 실행시킬 경우 에러가 발생한다
    cout << N << "-threads set by me" << endl;
    omp_set_num_threads(N);
    cout << omp_get_num_procs() << "-procs" << endl; // num_procs == max_threads?
    cout << omp_get_max_threads() << "-max.threads" << endl;
    cout << omp_get_num_threads() << "-threads now" << endl;
    cout << omp_get_thread_num() << "-th thread now" << endl;
    cout << "check this region parallel : " << omp_in_parallel() << endl;
    cout << "Fork" << endl;

    int sum = 0;
    int i;
    int a[16] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};

    #pragma omp parallel reduction(+:sum) private(i) shared(a)
    {
        i = omp_get_thread_num();
        // #pragma omp atomic
        sum += a[i] * a[i];

        #pragma omp critical(print_func)
        {
            cout << "thread : " << i << " and sum : " << sum << endl;
        }

        #pragma omp barrer
        {
            #pragma omp master
            {
                cout << "thread : " << i << " is master thread" << endl;
            }

        }
    }
    cout << "Join Process" << endl;
    cout << "sumation by parallel reduction : " << sum << endl;

    #pragma omp parallel private(id, fraction) shared(N) reduction(+:total_count)
    {
        total_count ++;
        id = omp_get_thread_num();
        fraction = (float)id /(float)N;

        #pragma omp critical(printing)
        {
            cout << "Hello, I'm thread-" << id << endl;
            cout << id << "/" << N << "=" << fraction << endl;
            cout << "check this region parallel : " << omp_in_parallel() << endl;
        }
    }

    /*
    #pragma omp parallel private(id, fraction) shared(N,a) reduction(+:total_sum)
    {

    }
    */

    cout << "Join : " << total_count << endl;

    return 0;
}