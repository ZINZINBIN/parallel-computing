#include <iostream>
#include <mpi.h>

using namespace std;

void main(){
    int rank, size;
    MPI_Init();
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0){
        cout << "Hello, World!" << endl;
    }
}