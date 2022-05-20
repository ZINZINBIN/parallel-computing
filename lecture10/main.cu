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
#define nTx 16
#define nTy 16

__global__ void kernel(void){
    cout << "kernel function proceeded" << endl;
}

int main(int argc, char *argv[])
{
    int N = 32;
    size_t size = N * sizeof(float);
    float *h_a = (float *)malloc(size);
    float *h_v;
    float *d_a;
    float *d_v;

    int nBx = 4;
    int nBy = 4;

    dim3 grid(nBx,nBy);
    dim3 block(nTx,nTy);

    // example 1. asynchronous data transfer
    cout << "cuda memory allocation to host" << endl;
    cudaHostAlloc(&d_a, size);

    cudaMemcpyAsync(d_a, h_a, size, cudaMemcpyHostToDevice, 0); // the last argument : stream 
    kernel <<<grid,block>>> ();

    // example 1-1. use cuda stream for asynchronous data transfer
    cudaStream_t stream1, stream2;

    cudaHostAlloc(&h_v, size, cudaHostAllocDefault);
    cudaStreamCreate(&stream1);
    cudaStreamCreate(&stream2);

    // stream 1 
    cudaMemcpyAsync(d_v, h_v, size, cudaMemcpyHostToDevice, stream1);

    // stream 2
    kernel <<<grid, block, stream2>>>();

    // free memory
    cudaStreamDistroy(stream1);
    cudaStreamDestroy(stream2);
    cudaFreeHost(h_v);

    // example 1-3. asynchronous data transfer using openacc

    return 0;
}