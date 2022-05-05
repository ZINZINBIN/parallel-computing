/**
 * @file main.cpp
 * @author zinzinbin
 * @brief Example Code on lecture 08 : OpenMPI + non blocking method
 * @version 0.1
 * @date 2022-05-03
 *
 * How to execute
 * (1) mpic++ main.cpp -o main.out
 * (2) mpirun -np <number of processes> ./main.out
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
#define MAXPROC 16

struct{
    int num;
    float x;
    double y[10];
}data;

using namespace std;

int main(int argc, char *argv[]){

    int rank;
    int size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0)
    {
        cout << "Running on " << size << "-Processes" << endl;
    }

    // cout << "Greetings from process " << rank << endl;

    // create new MPI type
    int m = 16;
    int n = 16;

    // c 언어에서는 row data의 경우 주소가 인접하므로, 하나의 row에 대해 데이터를 분리해서 처리가 가능하다
    // 그러나, col data의 경우에는 주소가 연속적이지 않으므로 한번에 send 할 수 없다
    // 따라서, MPI_Type_vector를 통해 하나의 col에 대해 새로운 데이터 구조체를 정의하여 한번에 send하도록 쓸 수 있다. 
    // 참조 : http://k-atoms.ksc.re.kr/mpi/2_4_2.html

    double A[m][n]; // dynamic allocation한 array의 경우에는 데이터가 연속적이지 않을 수 있다. 그래서 vector type을 쓸 수 없을 수도 있다(주의하자)
    MPI_Datatype coltype;
    MPI_Status status;

    MPI_Type_vector(m, 1, n, MPI_DOUBLE, &coltype);
    MPI_Type_commit(&coltype);

    int dtype_size;
    int dtype_count;

    if(rank == 0){
        MPI_Type_size(coltype, &dtype_size);
        MPI_Get_count(&status, coltype, &dtype_count);
        cout << "MPI new vector type generated " << endl;
        cout << "MPI new type size : " << dtype_size << endl;
        cout << "MPI new type count : " << dtype_count << endl;
    }

    /*
    int block_length[3] = {1, 1, 10};
    MPI_Datatype newstruct, types[3] = {MPI_INT, MPI_FLOAT, MPI_DOUBLE};
    MPI_Aint intext, floatext, displace[3];

    MPI_Type_extend(MPI_INT, &intext);
    MPI_Type_extend(MPI_FLOAT, &floatext);
    
    displace[0] = (MPI_Aint)0;
    displace[1] = intext;
    displace[2] = intext + floatext;

    MPI_Type_struct(3, block_length, displace, types, &newstruct);
    */


    MPI_Type_free(&coltype); 
    MPI_Finalize();
    return 0;
}