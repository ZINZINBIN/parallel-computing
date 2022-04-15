/**
 * @file main.cpp
 * @author zinzinbin
 * @brief matrix multiplication with OpenMPI
 * @version 0.1
 * @date 2022-04-15
 * 
 * How to execute
 * (1) mpic++ matrix_multiplication.cpp -o matrix_multiplication.out
 * (2) mpirun -np <number of processes> ./matrix_multiplication.out
 */

#include <iostream>
#include <mpi.h>
using namespace std;

int main(int argc, char *argv[]){

    int baton = 1;
    
    int rank, size; 
    int tag;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size); // get number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get current process id

    if(rank == 0){
        cout << "number of processes : " << size << endl;
        cout << "current process id : " << rank << endl;
    }

    if(rank == 0){
        baton = 1;
        MPI_Send(&baton, 1, MPI_INT, 1, 999, MPI_COMM_WORLD);
        MPI_Recv(&baton, 1, MPI_INT, size -1, 999, MPI_COMM_WORLD, &status);

        cout << "baton : " << baton << ", Process : " << size - 1 << " --> Process : 0"<<endl;
    }
    else{
        MPI_Recv(&baton, 1, MPI_INT, rank - 1, 999, MPI_COMM_WORLD, &status);
        cout << "baton : " << baton << ", Process : " << rank - 1 << " --> Process : " << rank <<endl;
        MPI_Send(&baton, 1, MPI_INT, (rank + 1) % size, 999, MPI_COMM_WORLD);
    }

    // clean up the MPI state in preparation for the processes to exit
    MPI_Finalize();
    return 0;
}