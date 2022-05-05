/**
 * How to execute
 * (1) mpic++ quiz02_mpi2.cpp -o quiz02_mpi2.out
 * (2) mpirun -np <number of processes> ./quiz02_mpi2.out
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


int main(int argc, char *argv[]){

    int baton = 1;
    int MAX_PROC = 16;
    int rank, size; 
    int tag;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size); // get number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get current process id

    MPI_Status status;
    MPI_Request request[2];

    int value;
    int N = MAX_PROC;
    int *values = new int[N];

    MPI_Gather(&value, 1, MPI_INT, values, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Communicator
    MPI_Comm world, workers;
    MPI_Group world_group, worker_group;
    int ranks[1];

    world = MPI_COMM_WORLD;
    
    // exclude할 rank와 그를 담는 array : ranks
    int server = size - 1;
    ranks[0] = server;

    MPI_Comm_group(world, &world_group);
    MPI_Group_excl(world_group, 1, ranks, &worker_group);
    MPI_Comm_create(world, worker_group, &workers); // arg : (comm, group, newcomm, ierror : none) -> first argument comm is needed for the other information such as context

    MPI_Group_free(&worker_group);
    MPI_Group_free(&world_group);

    MPI_Comm_free(&world);
    MPI_Comm_free(&workers);

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