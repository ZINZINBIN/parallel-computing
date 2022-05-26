/**
 * How to execute
 * (1) mpic++ quiz03_mpi.cpp -o quiz03_mpi.out
 * (2) mpirun -np <number of processes> ./quiz03_mpi.out
 */

#include <iostream>
#include <mpi.h>

using namespace std;

int generate_random(int min, int max)
{
    int fraction = int(max - min + 1.0);
    return int(min + (std::rand() % fraction));
}

int *generate_rands_arr(int M){
    int *rands = new int[M];
    for(int i = 0; i < M; i++){
        rands[i] = generate_random(1, 10);
    }
    return rands;
}

int compute_avg(int * arr, int M){
    int avg = 0;
    for(int i = 0; i < M; i++){
        avg += arr[i];
    }
    return int(avg / M);
}

int main(int argc, char *argv[])
{
    int rank, size;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        cout << "Hello, World!" << endl;
    }

    // one to all and all to one : collective communicator
    // reference : https://mpitutorial.com/tutorials/mpi-scatter-gather-and-allgather/
    
    int value = 32;
    int element_per_procs = 32;
    int *values = new int[element_per_procs];

    if(rank == 0){
        // MPI_Bcast(value, element_per_procs, MPI_INT, 0, MPI_COMM_WORLD);
    }

    value = 16;

    int *recv_arrs = new int[element_per_procs];
    int *recv_avg_arrs = new int[size];

    if(rank == 0){
        int *random_arrs = generate_rands_arr(size * element_per_procs);
        MPI_Scatter(random_arrs, element_per_procs, MPI_INT, recv_arrs, element_per_procs, MPI_INT, 0, MPI_COMM_WORLD);
    }

    int avg_per_procs = compute_avg(recv_arrs, element_per_procs);

    MPI_Gather(&avg_per_procs, 1, MPI_INT, recv_avg_arrs, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(rank == 0){
        for(int i = 0; i < size; i++){
            cout << "avg_per_procs : " << recv_avg_arrs[i] << endl;
        }
    }

    // generate new datatype for MPI
    // datatype : contiguous, vector, indexed, struct
    // process
    // (1) construction : define, (ex : MPI_Type_contiguous)
    // (2) commission to the system : MPI_Type_commit(new_type)
    // (3) use new type 
    // (4) Destruction : cleaning (ex : MPI_Type_free(new_type))

    int col = 0;
    MPI_Datatype coltype;
    double A[5][8];
    
    // MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype old_type, MPI_Datatype *new_type) 
    // MPI_Type_commit(MPI_Datatype *new_type)  
    // matrix 내 column 데이터를 받기 위해 vecto를 이용해 저장할 예정(5 x 1)

    MPI_Type_vector(5, 1, 8, MPI_DOUBLE, &coltype);
    MPI_Type_commit(&coltype);

    int count;
    double *buffer = new double[5];

    if(rank == 0){
        MPI_Send(&A[0][col], 5, coltype, 1, 0, MPI_COMM_WORLD);
    }
    else{
        MPI_Recv(buffer, 5, coltype, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_INT, &count);
        cout << "MPI get count : " << count << endl;
    }

    // struct data type
    struct
    {
        int num;
        float x;
        double y[10];
    } data;

    int blocklength[3] = {1,1,10};
    MPI_Datatype newstruct, type3[3] = {MPI_INT, MPI_FLOAT, MPI_DOUBLE};
    MPI_Aint intext, floatext, displace[3];

    // MPI_Type_Extent : length of occupied memory에 대한 확보
    MPI_Type_extent(MPI_INT, &intext);
    MPI_Type_extent(MPI_FLOAT, &floatext);

    displace[0] = (MPI_Aint)0;
    displace[1] = intext;
    displace[2] = intext + floatext;
    //  MPI_Type_struct(int count, int *blocklength_array, MPI_Aint *displacement_array, MPI_Datatype*type_array, MPI_Datatype *newtype)
    MPI_Type_struct(3, blocklength, displace, type3, &newstruct);

    // find the size of datatype
    int sample_size;

    MPI_Type_size(MPI_DOUBLE, &sample_size);

    cout << "new struct size : " << sample_size << endl;

    MPI_Aint sample_extent;

    MPI_Type_extent(MPI_DOUBLE, &sample_extent);
    cout << "new struct extent : " << sample_extent << endl;

    MPI_Finalize();
    return 0;
}