/**
 * @file main.cpp
 * @author zinzinbin
 * @brief Example Code on lecture 09 : cuda and openacc example
 * @version 0.1
 * @date 2022-05-06
 *
 * How to execute
 * (1) nvcc main.cu -o main.out
 * (2) ./main.out
 */

#include <iostream>
#define nTx 8
#define nTy 8
using namespace std;

// 자체 실습 : Data transfer, Matrix transpose (CUDA, openacc)
// kernel definintion : subroutines executing on GPUs
__global__ void kernel(void){
    
}

__global__ void VecAdd(float *A, float *B, float *C){
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    C[idx] = A[idx] + B[idx];
    printf("c[%d] = %.3f + %.3f = %.3f\n",idx,A[idx], B[idx], C[idx]);
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

__global__ void MatrixTranspose(int *inp_mat, int *out_mat, int Nrow, int Ncol){
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int j = blockDim.y * blockIdx.y + threadIdx.y;
    out_mat[j * Nrow + i] = inp_mat[Ncol * i + j];
}

__global__ void MatrixTransposeTiled(int *inp_mat, int *out_mat, int Nrow, int Ncol){
    __shared__ int tile[nTx][nTy];
    int i = threadIdx.x;
    int j = threadIdx.y;

    int bx = blockIdx.x * blockDim.x;
    int by = blockIdx.y * blockDim.y;

    tile[i][j] = inp_mat[(bx + i) * Ncol + by + j];
    __syncthreads(); // block 내 모든 thread를 동기화
    out_mat[(by + j) * Nrow + bx + i] = tile[i][j];
}

void print_matrix(int *A, int nrow, int ncol)
{
    int idx;
    for(int i = 0; i < nrow; i++){
        for(int j = 0; j < ncol; j++){
            idx = ncol * i + j;
            printf("%3d ",A[idx]);
        }
        cout << endl;
    }
}

void print_array(int *A, int N){
    for(int i = 0; i < N; i++){
        printf("%3d ", A[i]);
    }
    cout << endl;
}

int main(void){
    // Data transfer
    int N = 32;
    size_t size = N * sizeof(float);
    float *h_a = (float *)malloc(size);
    float *d_a;

    // example 1. cuda memory allocation, copy and deallocation
    cout << "cuda memory allocation" << endl;
    cudaMalloc(&d_a, size);

    kernel <<<1,1>>>();

    cout << "cuda memory copy" << endl;
    cudaMemcpy(d_a, h_a, N, cudaMemcpyHostToDevice);

    cudaMemcpy(h_a, d_a, N, cudaMemcpyDeviceToHost);

    cout << "cuda memory deallocation" << endl;
    cudaFree(d_a);
    free(h_a);

    // example 2. vector addition
    size_t vec_size = N * sizeof(float);
    float *a = (float *)malloc(vec_size);
    float *b = (float *)malloc(vec_size);
    float *c = (float *)malloc(vec_size);

    // float *d_a;
    float *d_b;
    float *d_c;

    for(int i = 0; i < N; i++){
        a[i] = i;
        b[i] = 2 * i;
        c[i] = 0;
    }

    cudaMalloc(&d_a, size);
    cudaMalloc(&d_b, size);
    cudaMalloc(&d_c, size);

    cudaMemcpy(d_a, a, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_c, c, size, cudaMemcpyHostToDevice);

    VecAdd <<<2,N>>> (d_a, d_b, d_c);

    cudaMemcpy(a, d_a, size, cudaMemcpyDeviceToHost);
    cudaMemcpy(b, d_b, size, cudaMemcpyDeviceToHost);
    cudaMemcpy(c, d_c, size, cudaMemcpyDeviceToHost);

    for(int i = 0; i < N; i++){
        cout << "c[" << i << "] = " << c[i] << endl;
    }

    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);
    free(a);
    free(b);
    free(c);

    // example 3. Matrix Transpose
    int Nrow = 16;
    int Ncol = 16;
    size_t mat_size = Nrow * Ncol * sizeof(int);
    int *inp_mat_h = (int *)malloc(mat_size);
    int *out_mat_h = (int *)malloc(mat_size);
    int *inp_mat_d;
    int *out_mat_d;

    int nBx = Nrow / nTx + (Nrow % nTx != 0);
    int nBy = Ncol / nTy + (Ncol % nTy != 0);
    dim3 grid(nBx, nBy); // nBx x nBy 형태의 block들로 구성된 grid
    dim3 threads(nTx,nTy);

    for(int i = 0; i < Nrow * Ncol; i++){
        inp_mat_h[i] = i;
    }

    cudaMalloc(&inp_mat_d, mat_size);
    cudaMalloc(&out_mat_d, mat_size);

    cudaMemcpy(inp_mat_d, inp_mat_h, mat_size, cudaMemcpyHostToDevice);

    MatrixTranspose <<<grid, threads>>> (inp_mat_d, out_mat_d, Nrow, Ncol);

    cudaMemcpy(out_mat_h, out_mat_d, mat_size, cudaMemcpyDeviceToHost);
    
    cout << "# input matrix" << endl;
    print_matrix(inp_mat_h, Nrow, Ncol);

    cout <<"# output matrix" << endl;
    print_matrix(out_mat_h, Nrow, Ncol);

    // example 4. matrix transpose with shared variable tile
    MatrixTransposeTiled <<<grid, threads>>>(inp_mat_d, out_mat_d, Nrow, Ncol);

    cudaMemcpy(out_mat_h, out_mat_d, mat_size, cudaMemcpyDeviceToHost);

    cout << "# output matrix for example 4" << endl;
    print_matrix(out_mat_h, Nrow, Ncol);

    cudaFree(inp_mat_d);
    cudaFree(out_mat_d);
    free(inp_mat_h);
    free(out_mat_h);

    // example 5. openacc
    N = 8;
    size_t a_size = N * sizeof(int);
    int *x = (int *)malloc(a_size);
    int *y = (int *)malloc(a_size);

    #pragma acc kernels
    {
        for(int i = 0; i < N; i++){
            x[i] = i;
            y[i] = 2 * i;
        }
        for(int i = 0; i < N; i++){
            y[i] += 3.0 * x[i];
        }
    }

    cout << "example 5 code" << endl;
    print_array(y,N);

    free(x);
    free(y);

    // matrix transpose with openacc
    int a_row = 16;
    int a_col = 16;
    int *A = new int [a_row * a_col];
    int *B = new int [a_row * a_col];

    for(int i = 0; i < a_row; i++){
        for(int j = 0; j < a_col; j++){
            int idx = i * a_col + j;
            A[idx] = i * j - i + j;
        }
    }

    #pragma acc parallel loop
    {
        for(int j = 0; j < a_col; j++){
            #pragma acc loop
            {
                for(int i = 0; i < a_row; i++){
                    B[j * a_col + i] = A[i * a_col + j];
                }
            }
        }
    }

    cout << "matrix transpose using openacc" << endl;
    cout << "original" << endl;
    print_matrix(A, a_row, a_col);

    cout << "transpose" << endl;
    print_matrix(B, a_row, a_col);

    free(A);
    free(B);

    // example 6. data directive of openacc : explicit data transfer
    float err = 1;
    float tol = 0.0001;
    int iter = 0;
    int iter_max = 256;

    a_row = 16;
    a_col = 16;
    float *M = new float[a_row * a_col];
    float *M_new = new float[a_row * a_col];

    for (int i = 0; i < a_row; i++)
    {
        for (int j = 0; j < a_col; j++)
        {
            int idx = i * a_col + j;
            M[idx] = 2.0 * i + 3.0 * j;
        }
    }

    #pragma acc data copy(M), create(M_new)
    {
        while(err > tol && iter < iter_max){
            iter += 1;
            err = 0;
            #pragma acc parallel loop reduction(max:err)
            {
                for(int i = 1; i < a_row - 1; i++){
                    for(int j = 1; j < a_col - 1; j++){
                        M_new[i * a_col + j] = 0.25 * (
                            M[i * a_col + j - 1] +
                            M[i * a_col + j + 1] +
                            M[i * a_col + j - a_col] +
                            M[i * a_col + j + a_col]
                        );
                        err = max(err, abs(M_new[i * a_col + j] - M[i * a_col + j]));
                    }
                }
            }
        }
    }

    cout << "error : " << err << ", iter : " << iter << endl;
    free(M_new);
    free(M);

    return 0;
}