/**
 * @file main.cpp
 * @author zinzinbin
 * @brief Example Code on lecture 10 : cuda
 * @version 0.1
 * @date 2022-05-13
 *
 * How to execute
 * (1) nvcc main.cu -o main.out
 * (2) ./main.out
 */

#include <iostream>
#define nTx 4
#define nTy 4
#define BLOCK_SIZE 16

using namespace std;

__global__ void kernel(void){
    // printf("kernel function proceed\n");
}

__global__ void matmul(float *a, float *b, float *c, int M, int N){
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    int col = blockIdx.y * blockDim.y + threadIdx.y;

    float sum = 0;
    for (int i = 0; i < M; i++)
    {
        sum += a[row * M + i] * b[i * N + col];
    }
    c[row * N + col] = sum; // C[row, col] = sum of a[row,i] * b[i, col]
}

__global__ void matmul_sm(float *a, float *b, float *c, int L, int M, int N){
    __shared__ float aTile[BLOCK_SIZE][BLOCK_SIZE];
    __shared__ float bTile[BLOCK_SIZE][BLOCK_SIZE];

    int brow = blockIdx.x;
    int bcol = blockIdx.y;

    int trow = threadIdx.x;
    int tcol = threadIdx.y;

    int row = blockIdx.x * blockDim.x + threadIdx.x;
    int col = blockIdx.y * blockDim.y + threadIdx.y;

    float sum = 0;

    int nBlock = L / BLOCK_SIZE + (L % BLOCK_SIZE != 0);

    for(int k = 0; k < nBlock; k++){
        
    }

}


int main(int argc, char *argv[])
{
    int N = 8;
    size_t size = N * sizeof(float);
    float *h_a = (float *)malloc(size);
    float *h_v;
    float *d_a;
    float *d_v;

    int nBx = 2;
    int nBy = 2;

    dim3 grid(nBx,nBy);
    dim3 block(nTx,nTy);

    // example 1. asynchronous data transfer
    cout << "cuda memory allocation to host" << endl;
    //cudaHostAlloc(&d_a, size);
    cudaMalloc(&d_a, size);

    cudaMemcpyAsync(d_a, h_a, size, cudaMemcpyHostToDevice, 0); // the last argument : stream 
    kernel <<<grid,block>>> ();

    // example 1-1. use cuda stream for asynchronous data transfer
    cudaStream_t stream1, stream2;

    cudaHostAlloc(&h_v, size, cudaHostAllocDefault);
    cudaMalloc(&d_v, size);
    cudaStreamCreate(&stream1);
    cudaStreamCreate(&stream2);

    // stream 1 
    cudaMemcpyAsync(d_v, h_v, size, cudaMemcpyHostToDevice, stream1);

    // stream 2
    kernel <<<grid, block, 0, stream2>>>();

    // free memory
    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);
    cudaFreeHost(h_v);

    // example 1-3. asynchronous data transfer using openacc
    int *a = new int[N];
    int *b = new int[N];
    int *c = new int[N];
    
    #pragma acc parallel loop async(1)
    {
        for(int i = 0; i < N; i++){
            a[i] = i % 8;
        }
    }

    #pragma acc parallel loop async(1)
    {
        for (int i = 0; i < N; i++)
        {
            b[i] = 2 * a[i];
        }
    }
    
    #pragma acc wait(1) async(2)
    {
        for (int i = 0; i < N; i++)
        {
            c[i] = a[i] + b[i];
        }
    }

    for(int i = 0; i < N; i++){
        printf("a : %d, b : %d, c : %d", a[i],b[i],c[i]);
        cout<<endl;
    }

    free(a);
    free(b);
    free(c);

    // example 2. Matrix multiplication with simple example and using shared memory
    int M_row = 16;
    int M_col = 16;
    int N_row = 16;
    int N_col = 16;
    float *M = new float[M_row * M_col];
    float *N = new float[N_row * N_col];

    float *M_d;
    float *N_d;

    cudaMalloc(&M_d, M_row * M_col);
    cudaMalloc(&N_d, N_row * N_col);

    for (int i = 0; i < N_row; i++)
    {
        for (int j = 0; j < N_col; j++)
        {
            int idx = i * N_col + j;
            N[idx] = 2.0 * i + 3.0 * j;
        }
    }

    for (int i = 0; i < M_row; i++)
    {
        for (int j = 0; j < M_col; j++)
        {
            int idx = i * M_col + j;
            M[idx] = 3.0 * i + 2.0 * j;
        }
    }



    return 0;
}