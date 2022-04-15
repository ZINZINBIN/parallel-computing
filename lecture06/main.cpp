/**
 * @file main.cpp
 * @author zinzinbin
 * @brief main.cpp file for example code
 * @version 0.1
 * @date 2022-04-15
 * 
 * Ring example code with send and recive
 * How to execute
 * (1) mpic++ main.cpp -o main.out
 * (2) mpirun -np <number of processes> ./main.out
 */

#include <iostream>
#include <mpi.h>
using namespace std;

int main(int argc, char *argv[]){

    int baton = 1;
    int status = 0;
    
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