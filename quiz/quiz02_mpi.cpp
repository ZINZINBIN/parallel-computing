/**
 * How to execute
 * (1) mpic++ quiz02_mpi.cpp -o quiz02_mpi.out
 * (2) mpirun -np <number of processes> ./quiz02_mpi.out
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

using namespace std;

// function for use
void init_random_array(float *x, int n_size, float min, float max);
void init_random_matrix(float **A, int m, int n, float min, float max);
float get_random_number(float min, float max);
float** generate_matrix(int m, int n);
float* generate_array(int m);
void check_solution_validity(float **A, float *x, float *target, int m, int n, float tolerance);

// define op function
void op_func(void *a, void *b, int *len, MPI_Datatype *dt){
    int i;
    if(*dt == MPI_INT){
        for(i = 0; i<*len; i++){
            ((int *)b)[i] = 2 * ((int *)b)[i] + 2 * ((int *)a)[i];
        }
    }
}

int main(int argc, char *argv[]){

    int baton = 1;
    int MAX_PROC = 16;
    int rank, size; 
    int tag;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size); // get number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get current process id

    if(rank == 0){
        //cout << "rank : 0, example test" << endl;
    }
    else{
        //cout << "rank : " << rank << ", example test" << endl;
    }

    // define op function
    // MPI_Op newop;
    // MPI_Op_create(func_op, 1, &newop)

    int n_interval;
    int ierr;
    float dx;
    float s = 0;
    float x;
    float total;

    /**
     * @brief Numerical integration
     * integral 1/(1+x^2) from 0 to 1
     */

    if(rank == 0){
        cout << "Enter the number of inteval : ";
        cin >> n_interval;
    }

    dx = 1.0 / (float) n_interval;
    s = 0;

    ierr = MPI_Bcast(&n_interval, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(rank == 0){
        cout << "MPI Bcast ierror : " << ierr << endl;
    }

    if(ierr == MPI_ERR_INTERN){
        MPI_Abort(MPI_COMM_WORLD, ierr);
    }

    for(int i = 0; i < n_interval; i+= size){
        x = i * dx;
        s += 4 / (1 + x * x);
    }

    MPI_Reduce(&s, &total, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); // 각각의 process에서 연산된 s를 buffer인 total로 보낸다. 여기서 send_buffer는 s, recv_buffer는 total
    
    if(rank == 0){
        total += (4.0 + 2.0) / 2.0;
        total = total * dx;
        cout << "total : " << total << endl;
    }




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

void init_random_array(float *x, int n_size, float min, float max){
    for(int i = 0; i < n_size; i++){
        x[i] = get_random_number(min, max);
    }
}

void init_random_matrix(float **A, int m, int n, float min, float max){
    // if abs(A[i][i]) > sum of A[i][j] with j != i, then jacobi method can be converged
    // check diagonally dominanted
    int abs_sum = 0;

    for (int i = 0; i < m; i++)
    {
        abs_sum = 0;
        for (int j = 0; j < n; j++)
        {   
            // initialize matrix with condition : diagonal dominant
            if(i!=j){
                A[i][j] = get_random_number(min, max);
                abs_sum += fabs(A[i][j]);
            }
        }
        A[i][i] = abs_sum + 1;
    }

}

float get_random_number(float min, float max){
    int fraction = int(max - min + 1.0);
    return float(min + (std::rand() % fraction));
}

void check_solution_validity(float **A, float *x, float *target, int m, int n, float tolerance)
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

    for(int i = 0; i < m; i++){
        cout << "target[" << i << "] : " << target[i] << ", y[" << i << "] : " << y[i] << endl;
    }
    return;
}