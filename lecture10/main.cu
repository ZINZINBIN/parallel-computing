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

// matrix multiplication with 1D array
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

// matrix multiplication using shared memory
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
        aTile[trow][tcol] = a[(BLOCK_SIZE*brow + trow)*M + k * BLOCK_SIZE + tcol];
        bTile[trow][tcol] = b[(k * BLOCK_SIZE + trow) * N + BLOCK_SIZE * bcol + tcol];
        __syncthreads();
        for (int i = 0; i < BLOCK_SIZE; i++){
            sum += aTile[trow][i] * bTile[i][tcol];
        }
        __syncthreads();
    }

    if(row < L && col < N){
        c[row * M + col] = sum;
    }
}

__global__ void reduction_sum(float *psum){
    int half = blockDim.x / 2;
    int id = threadIdx.x;

    while(half >= 1){
        if(id < half){
            psum[id] += psum[id + half];
        }
        half /= 2;
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
    float *N_ = new float[N_row * N_col];

    float *M_d;
    float *N_d;

    cudaMalloc(&M_d, M_row * M_col);
    cudaMalloc(&N_d, N_row * N_col);

    for (int i = 0; i < N_row; i++)
    {
        for (int j = 0; j < N_col; j++)
        {
            int idx = i * N_col + j;
            N_[idx] = 2.0 * i + 3.0 * j;
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


    // example 3. 2D array in Cuda C : memory allocation from host device with 2D array to GPU with 1D array
    int N_size = 32;
    float **example3_A = new float *[N_size];
    float *example3_x = new float [N_size];
    float *example3_y = new float [N_size];

    float *example3_d_A = new float [N_size];
    float *example3_d_x = new float [N_size];
    float *example3_d_y = new float [N_size];

    for(int i = 0; i < N_size; i++){
        example3_A[i] = new float [N_size];
    }

    cudaMalloc(&example3_d_A, N_size * N_size * sizeof(float));
    cudaMalloc(&example3_d_x, N_size * sizeof(float));
    cudaMalloc(&example3_d_y, N_size * sizeof(float));

    cudaMemcpy(example3_d_A, example3_A, N_size * N_size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(example3_d_x, example3_x, N_size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(example3_d_y, example3_y, N_size * sizeof(float), cudaMemcpyHostToDevice);

    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks(N / threadsPerBlock.x, N / threadsPerBlock.y);

    matmul<<<numBlocks,threadsPerBlock>>>(example3_d_A, example3_d_x, example3_d_y, N_size, N_size);
    
    cudaMemcpy(example3_y, example3_d_y, N_size * N_size * sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(example3_d_y);
    cudaFree(example3_d_x);
    cudaFree(example3_d_A);

    free(example3_A);
    free(example3_x);
    free(example3_y);

    // example 4. Memory allocation and transfer for 2D or 3D array to 1D array in CUDA
    int example4_N = 32;
    float **example4_A = new float *[example4_N];
    float *example4_d_A = new float[example4_N];
 
    for (int i = 0; i < example4_N; i++)
    {
        example4_A[i] = new float[example4_N];
    }

    // cudaMallocPitch(&example4_d_A, &example4_A, N_size * sizeof(float), N_size); // args : &devptr,&devpitch, N_col * sizeof(type), N_row 
    
    // example 5. multiple GPU
    int gpu_num = 0;
    cudaGetDeviceCount(&gpu_num);

    cout << "cudaGetDeviceCount : " << gpu_num << endl;

    cudaDeviceProp props;
    cudaGetDeviceProperties(&props, 1);

    cudaSetDevice(0); // using gpu : 0 only

    // example 6. Zero Copy
    // zero-copy memory는 device memory에 mapping되는 pinned-memory
    // host와 device 모두 zero copy에 접근할 수 있다. 
    // more detail info : https://junstar92.tistory.com/285
    float *a_h; 
    float *a_map;

    int n_size = 32;
    int nBytes = sizeof(float) * n_size;

    cudaGetDeviceProperties(&props, 0);

    if(!props.canMapHostMemory){
        exit(0);
    }

    cudaSetDeviceFlags(cudaDeviceMapHost);
    cudaHostAlloc(&a_h, nBytes, cudaHostAllocMapped); // cudaHostAllocMapped : mapped pinned-memory, cudaHostAllocPortable : zero-copy between GPUs
    cudaHostGetDevicePointer(&a_map, a_h, 0);
    
    // cudaFreeHost();

    // example 7. atomic functions
    // atomic function : shared memory를 이용한 계산을 통해 race conditions을 방지한다
    cout << "atomic functions example" << endl;
    int arr_a[5]= {1,2,3,4,5};
    int *arr_dev_a;
    
    cudaMalloc(&arr_dev_a, 5 * sizeof(int));
    cudaMemcpy(arr_dev_a, arr_a, 5 * sizeof(int), cudaMemcpyHostToDevice);

    int const_c = 5;
    int y;

    // y = atomicAdd(&arr_a[0], const_c);
    // cout << "y = " << y << endl;

    cudaFree(arr_dev_a);

    // example 8. Reduction example with OpenACC
    int example8_M = 8;
    int example8_N = 8;
    int example8_L = 8;

    float **example8_A = new float *[example8_M];
    float **example8_B = new float *[example8_N];
    float **example8_C = new float *[example8_M];
    for(int i = 0; i < example8_M; i++){
        example8_A[i] = new float [example8_N];
        example8_C[i] = new float [example8_L];
    }

    for (int i = 0; i < example8_N; i++)
    {
        example8_B[i] = new float[example8_L];
    }

    #pragma acc parallel loop collapse(2)
    {
        for(int i = 0; i < example8_M; i++){
            for(int j = 0; j < example8_N; j++){
                float c_ij = 0;
                #pragma acc loop reduction(+:c_ij)
                {
                    for(int k = 0; k < example8_L; k++){
                        c_ij += example8_A[i][k] * example8_B[k][j];
                    }
                }
                example8_C[i][j] = c_ij;

            }
        }
    }

    // example 9. reduction in shared memory example
    // detail : https://junstar92.tistory.com/290
    
    return 0;
}