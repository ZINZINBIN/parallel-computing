#include <iostream>
#define nTx 16
#define nTy 16

// 자체 실습 : Data transfer, Matrix transpose (CUDA, openacc)

// kernel definintion : subroutines executing on GPUs
__global__ void kernel(void){
    
}

__global__ void VecAdd(float *A, float *B, float *C){
    int i = threadIdx.x;
    C[i] = A[i] + B[i];
}

__global__ void simpleMultiply(float *a, float *b, float *c, int M, int N){
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    double sum = 0;
    for(int i = 0; i < M; i++){
        sum += a[row * M + i] * b[i*N + col];
    }
    c[row *N + col] = sum;
}

int main(void){
    // Data transfer
    int N = 32;
    size_t size = N * sizeof(float);
    float *h_a = (float *)malloc(size);
    float *d_a;
    cudaMalloc(&d_a, size);

    // kernel<<<1,1>>>();

    cudaMemcpy(d_a, h_a, N, cudaMemcpyHostToDevice);

    cudaMemcpy(h_a, d_a, N, cudaMemcpyDeviceToHost);

    cudaFree(d_a);
    free(h_a);

    return 0;
}