#include <iostream>
#include <omp.h>

using namespace std;

int main(void){
    int id, N;
    float fraction;
    N = 4;

    cout << N << "threads set by me" << endl;
    omp_set_num_threads(N);
    cout << omp_get_num_procs() << "procs" << endl;
    cout << omp_get_max_threads() << "max. thread" << endl;
    cout << omp_get_num_threads() << "thread now" << endl;
    cout << "Fork" << endl;

    #pragma omp parallel private(id, fraction) shared(N)
    {
        id = omp_get_thread_num();
        fraction = (float)id / (float)N;

        #pragma omp critical(printing)
        {
            cout << "Hello, I'm thread" << id << endl;
            cout << id << "/" << N << "=" << fraction << endl;
        }
    }
    cout << "Join"<< endl;
    return 0;
}