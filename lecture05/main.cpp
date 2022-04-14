/**
 * @file main.cpp
 * @author zinzinbin
 * @brief OpenMPI practice with some example code
 * @version 0.1
 * @date 2022-04-14
 * 
 * How to execute
 * mpicxx main.cpp : execute main.cpp with OpenMPI interface
 * maincxx main.cpp -o main.out : execute main.cpp with OpenMPI interface and save output file as main.out
 */

#include <iostream>
#include <mpi.h>

using namespace std;

int main(int argc, char *argv[]){
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0){
        cout << "Hello, World!" << endl;
    }

    MPI_Finalize();
    return 0;
}