/**
 * @file main.cpp
 * @author zinzinbin
 * @brief Using OpenMPI, numerical integration code build up
 * @version 0.1
 * @date 2022-04-14
 * 
 * How to execute
 * (1) mpic++ numerical_integration.cpp -o numerical_integration.cpp
 * (2) mpirun -np <number of processes> ./numerical_integration.out
 */

#include <iostream>
#include <mpi.h>
float func(float x, float a, float b, float c);
float integral(float xl, float xr, float a, float b, float c);

using namespace std;

float func(float x, float a, float b, float c){
    float y = a * x * x + b * x + c;
    return y;
}

float integral(float xl, float xr, float a, float b, float c){
    float yl = func(xl, a/3, b/2, c) * xl;
    float yr = func(xr, a/3, b/2, c) * xr;
    float area = yr - yl;
    return area;
}

int main(int argc, char *argv[]){
    int n;
    float xl = 0.0;
    float xr = 2.0;
    float x;
    float dx = 0;
    float sum = 0;
    float total = 0;
    float answer;
    float a = 1;
    float b = 2;
    float c = 1;
    int i;
    int rank, size; 
    int tag;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size); // get number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get current process id

    if(rank == 0){
        cout << "number of processes : " << size << endl;
        cout << "current process id : " << rank << endl;
        cout << "Enter the number of interval : ";
        cin >> n;
    }

    // Braodcast a buffer from a sending process to all other processes
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    dx = (xr - xl) / n;
    sum = 0; // 이걸 위에서 초기화했을 경우에는 해가 수렴 xx 왜 이럴까?

    for(i = rank; i < n; i+=size){
        x = xl + dx * i;
        sum += func(x,a,b,c);
    }

    if(rank == 0){
        answer = sum * dx;
        cout << "rank0 : " << answer << endl;
    }

    if(rank == 1){
        answer = sum * dx;
        cout << "rank1 : " << answer << endl;
    }

    if(rank == 2){
        answer = sum * dx;
        cout << "rank2 : " << answer << endl;
    }

    if(rank == 3){
        answer = sum * dx;
        cout << "rank3 : " << answer << endl;
    }

    // Reduction : summation of all result to rank = 0
    // question : MPI_Reduce를 쓰는 경우에는 Barrier가 default로 작동되는 구조인가 아니면 별도로 봐야 하는가? 여기서는 별도로 보더라도 적분 수행에 있어서 문제는 없음
    MPI_Reduce(&sum, &total, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

    if(rank == 0){
        answer = total * dx;
        cout << "Numerical answer : " << answer << ", Analyric answer : " << integral(xl,xr,a,b,c); 
        cout << ", relative err : " << abs(answer - integral(xl,xr,a,b,c)) / abs(integral(xl,xr,a,b,c)) << endl;
    }

    // clean up the MPI state in preparation for the processes to exit
    MPI_Finalize();
    return 0;
}