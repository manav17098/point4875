#include<stdio.h>
#include<math.h>
#include<complex.h>
#include<cufft.h>
#define Nx 10
#define BATCH 1
typedef double complex cplx;

int main()
{
	cufftHandle plan;
	cuDoubleComplex *data;
	cuDoubleComplex dataH[Nx];
	for(int i=0;i<Nx;i++)
	{
		dataH[i].x=i;
		dataH[i].y=0.0;
	}
	cudaMalloc((void**)&data,sizeof(cuDoubleComplex)*Nx*BATCH);
	if (cudaGetLastError() != cudaSuccess)
	{
		fprintf(stderr, "Cuda error: Failed to allocate\n");
		return 0;	
	}
	if (cufftPlan1d(&plan, Nx, CUFFT_Z2Z,BATCH) != CUFFT_SUCCESS)
	{
		fprintf(stderr, "CUFFT error: Plan creation failed");
		return 0;	
	}	
	size_t sizeData= Nx*sizeof(cuDoubleComplex);
	cudaMemcpy(data,dataH,sizeData,cudaMemcpyHostToDevice);
	if (cufftExecZ2Z(plan, data, data, CUFFT_FORWARD) != CUFFT_SUCCESS)
	{
		fprintf(stderr, "CUFFT error: ExecC2C Forward failed");
		return 0;
	}
	cudaMemcpy(dataH,data,sizeData,cudaMemcpyDeviceToHost);
	printf("\n");
	for(int i=0;i<Nx;i++)
		printf("%f:%f ",dataH[i].x,dataH[i].y);
	return 0;
}
