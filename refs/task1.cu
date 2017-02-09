#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

// RG*RG*MAXN must fit within mytype

#define MAXN 100000
#define RG 10
#define USECPSEC 1000000ULL
#define nTPB 256
#define DSIZE 8192


//cuda error checking macros
#ifdef DEBUG
#define CUDA_CALL(F)  if( (F) != cudaSuccess ) \
  {printf("Error %s at %s:%d\n", cudaGetErrorString(cudaGetLastError()), \
   __FILE__,__LINE__); exit(-1);} 
#define CUDA_CHECK()  if( (cudaPeekAtLastError()) != cudaSuccess ) \
  {printf("Error %s at %s:%d\n", cudaGetErrorString(cudaGetLastError()), \
   __FILE__,__LINE__-1); exit(-1);} 
#else
#define CUDA_CALL(F) (F)
#define CUDA_CHECK() 
#endif


typedef float mytype;
// host function to compute convolution reference results
void conv(const mytype *A, const mytype *B, mytype* C, int N) {

    for (int j = 0; j < ((2*N)-1); ++j){ // iterate over columns of result
        mytype my_sum = 0;
        for (int i = 0; i < N; ++i)  // iterate down each column
          if (((j < N) && (i <= j)) || ((j >= N) && (i > (j-N)))) my_sum += A[i]*B[j-i];
        C[j] = my_sum;}
}
// host function - alternate realization
void conv2(const mytype *A, const mytype *B, mytype* C, int N) {

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            C[i + j] += A[i] * B[j];
}
// timing measurement function
unsigned long long dtime_usec(unsigned long long prev){
  timeval tv1;
  gettimeofday(&tv1,0);
  return ((tv1.tv_sec * USECPSEC)+tv1.tv_usec) - prev;
}
// convolution GPU kernel
// Task 1
__global__ void conv_Kernel(const mytype * __restrict__ A, const mytype * __restrict__ B, mytype *C, const int N){
    int idx = threadIdx.x+blockDim.x*blockIdx.x;
    if (idx < (2*N)-1){
      mytype my_sum = 0;
// create a for-loop below to sum the products in a single column in Diagram 1
      for (int i = 0; i < N; i++)
        if (((idx < N) && (i <= idx)) || ((idx >= N) && (i > (idx-N)))) my_sum += A[i]*B[idx-i];
      C[idx] = my_sum;
    }
}


int main(int argc, char *argv[]){

  mytype *d_A, *A, *d_B, *B, *d_C, *C, *h_C;
  int my_N = DSIZE;
  if ((my_N < 1) || (my_N > MAXN)) {printf("N out of range\n"); return 1;}
// allocate host data
  A   = (mytype *)malloc(my_N*sizeof(mytype));
  B   = (mytype *)malloc(my_N*sizeof(mytype));
  C   = (mytype *)malloc(((2*my_N)-1)*sizeof(mytype));
  h_C = (mytype *)malloc(((2*my_N)-1)*sizeof(mytype));
// allocate device data
  CUDA_CALL(cudaMalloc(&d_A, my_N*sizeof(mytype)));
  CUDA_CALL(cudaMalloc(&d_B, my_N*sizeof(mytype)));
  CUDA_CALL(cudaMalloc(&d_C, ((2*my_N)-1)*sizeof(mytype)));
//initialize host input data
  for (int i=0; i < my_N; i++){
    A[i] = rand()%RG;
    B[i] = rand()%RG;}
//zero out host result data
  for (int i=0; i < (2*my_N)-1; i++){
    C[i]   = 0;
    h_C[i] = 0;}
//begin timing for host reference function
  unsigned long cpu_time = dtime_usec(0);
  conv(A, B, C, my_N);
  cpu_time = dtime_usec(cpu_time);
//initialize device result data
  CUDA_CALL(cudaMemset(d_C, 0, ((2*my_N)-1)*sizeof(mytype)));
//begin timing for host reference function
  unsigned long gpu_time = dtime_usec(0);
//copy host input data to device
  CUDA_CALL(cudaMemcpy(d_A, A, my_N*sizeof(mytype), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(d_B, B, my_N*sizeof(mytype), cudaMemcpyHostToDevice));
//run convolution kernel on GPU
  conv_Kernel<<<(((2*my_N)-1)+nTPB-1)/nTPB,nTPB>>>(d_A, d_B, d_C, my_N);
  CUDA_CHECK();
//copy results from device to host
  CUDA_CALL(cudaMemcpy(h_C, d_C, ((2*my_N)-1)*sizeof(mytype), cudaMemcpyDeviceToHost));
  gpu_time = dtime_usec(gpu_time);
//check validity of results
  for (int i = 0; i < ((2*my_N)-1); i++) if (C[i] != h_C[i]) {printf("FAIL at %d, cpu: %f, gpu %f\n", i, C[i], h_C[i]); return 1;}
//print timing and speed comparison
  printf("PASS.  cpu time: %lu us, gpu time: %lu us\n", cpu_time, gpu_time);
  printf("Speedup: cpu/gpu = %f\n", cpu_time/(float)gpu_time);
//all host and device allocated data will be implicitly freed at program termination
  return 0;
}
