#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <math.h>
#include <cuda_runtime.h>
#include <ctime>

//~ #include <thrust/reduce.h>
//~ #include <reduction.h>


extern "C" void apply_bc_cuda_(double* p_2);
extern "C" void catch_divergence_cuda_(double res2,int ierr,int it);
extern "C" void collect_res2_cuda_(double* res2, double* tres2, int* ierr, int* it);

static void HandleError( cudaError_t err, const char *file,int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))


void printGPUprops(void){
int count;
cudaDeviceProp prop;
HANDLE_ERROR( cudaGetDeviceCount( &count ) );
	for (int i=0; i< count; i++) {
	HANDLE_ERROR( cudaGetDeviceProperties( &prop, i ) );
	printf( "-- General Information for device %d ---\n", i );
	printf( "Name:%s\n", prop.name );
	printf( "Compute capability: %d.%d\n", prop.major, prop.minor );
	printf( "Clock rate: %d\n", prop.clockRate );
	printf( "Device copy overlap:" );
	if (prop.deviceOverlap)
		printf( "Enabled\n" );
	else
		printf( "Disabled\n" );
		printf( "Kernel execition timeout :" );
	if (prop.kernelExecTimeoutEnabled)
		printf( "Enabled\n" );
	else
		printf( "Disabled\n" );
		printf( "--- Memory Information for device %d ---\n", i );
		printf( "Total global mem:%ld\n", prop.totalGlobalMem );
		printf( "Total constant Mem: %ld\n", prop.totalConstMem );
		printf( "Max mem pitch: %ld\n", prop.memPitch );
		printf( "Texture Alignment:%ld\n", prop.textureAlignment );
		printf( "--- MP Information for device %d ---\n", i );
		printf( "Multiprocessor count:%d\n",prop.multiProcessorCount );
		printf( "Shared mem per mp:%ld\n", prop.sharedMemPerBlock );
		printf( "Registers per mp:%d\n", prop.regsPerBlock );
		printf( "Threads in warp:%d\n", prop.warpSize );
		printf( "Max threads per block: %d\n", prop.maxThreadsPerBlock );
		printf( "Max thread dimensions: (%d, %d, %d)\n", prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2] );
		printf( "Max grid dimensions: (%d, %d, %d)\n", prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2] );
		printf( "\n" );
	}

}

//////////// CUDA functions ///////////////

__device__ __host__
int convert_3Dto1Dsingle(int *coord, int *dim_d)
{
	int index = coord[0] + coord[1]*(dim_d[0]-1) + coord[2]*((dim_d[0]-1)*(dim_d[1]-1));
	printf("%d\n", index);
	return index;
}
__device__ __host__
void convert_3Dto1D(int *coord, int *dim_d, int* vec, int i)
{
	vec[i] = coord[0] + 
			 coord[1]*(dim_d[0]) + 
			 coord[2]*((dim_d[0])*(dim_d[1]));
}

__device__ __host__
void convert_4Dto1D(int *coord, int *dim_d,  int* vec, int i)
{
	vec[i] = coord[0] + coord[1]*dim_d[0] + 
						coord[2]*(dim_d[0]*dim_d[1])+
						coord[3]*(dim_d[0]*dim_d[1]*dim_d[2]);
}

__device__ __host__
void convert_1Dto3D(int i, int *dim_d, int *coord)
{
	coord[0] = i%dim_d[0];
	coord[1] = (i/dim_d[0])%dim_d[1];
	coord[2] = (i/dim_d[0]/dim_d[1]);
}


__device__ __host__
void indexP(int *iP, int *dim_d, int *coord)
{	
	// original
	int tempC[3]={coord[0],coord[1],coord[2]};
	convert_3Dto1D(tempC,dim_d,iP,0);
	// +X
	tempC[0]=coord[0]+1;
	tempC[1]=coord[1];
	tempC[2]=coord[2];
	convert_3Dto1D(tempC,dim_d,iP,1);
	// -X
	tempC[0]=coord[0]-1;
	tempC[1]=coord[1];
	tempC[2]=coord[2];
	convert_3Dto1D(tempC,dim_d,iP,2);
	// +Y
	tempC[0]=coord[0];
	tempC[1]=coord[1]+1;
	tempC[2]=coord[2];
	convert_3Dto1D(tempC,dim_d,iP,3);
	// -Y
	tempC[0]=coord[0];
	tempC[1]=coord[1]-1;
	tempC[2]=coord[2];
	convert_3Dto1D(tempC,dim_d,iP,4);
	// +Z
	tempC[0]=coord[0];
	tempC[1]=coord[1];
	tempC[2]=coord[2]+1;
	convert_3Dto1D(tempC,dim_d,iP,5);
	// -Z
	tempC[0]=coord[0];
	tempC[1]=coord[1];
	tempC[2]=coord[2]-1;
	convert_3Dto1D(tempC,dim_d,iP,6);
}

__device__ __host__
void indexA(int *iA, int *dim_d, int *coord)
{
	int tempC[4] = {coord[0],coord[1],coord[2],0};
	convert_4Dto1D(tempC,dim_d,iA,0);
	tempC[3] = 1;
	convert_4Dto1D(tempC,dim_d,iA,1);
	tempC[3] = 2;
	convert_4Dto1D(tempC,dim_d,iA,2);
	tempC[3] = 3;
	convert_4Dto1D(tempC,dim_d,iA,3);
	tempC[3] = 4;
	convert_4Dto1D(tempC,dim_d,iA,4);
	tempC[3] = 5;
	convert_4Dto1D(tempC,dim_d,iA,5);
	tempC[3] = 6;
	convert_4Dto1D(tempC,dim_d,iA,6);
	tempC[3] = 7;
	convert_4Dto1D(tempC,dim_d,iA,7);
	}




__global__
void p_iter(int *dim_d, int dimp, int dimStencil, int Ng, double *A_d, double *p_1, double *p_2, double beta)
{	
	int iA[8];
	int iP[7];
	int coordP[3];
	int coordA[3];
	
	int dimS[3];
	dimS[0] =*(dim_d) -2*Ng;
	dimS[1] =*(dim_d+1) -2*Ng;
	dimS[2] =*(dim_d+2) -2*Ng;
	
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	for (int i = index; i < dimStencil; i += stride)
	{
		convert_1Dto3D(i,dimS,coordA);
		coordP[0]=coordA[0]+Ng;
		coordP[1]=coordA[1]+Ng;
		coordP[2]=coordA[2]+Ng;

		indexP(iP, dim_d, coordP);
		indexA(iA, dimS, coordA);
		p_2[iP[0]]=(1.0-beta)*p_1[iP[0]] + (beta/A_d[iA[6]])*(   
		        A_d[iA[0]] * p_1[iP[2]] + A_d[iA[1]] * p_1[iP[1]] +  
		        A_d[iA[2]] * p_1[iP[4]] + A_d[iA[3]] * p_1[iP[3]] +  
		        A_d[iA[4]] * p_1[iP[6]] + A_d[iA[5]] * p_1[iP[5]] + A_d[iA[7]]);
		        
	}
}

__global__
void L1L2_norm(int *dim_d, int dimStencil, int Ng, double *A_d, double *p_1, double *res2, int norm)
{	
	int iA[8];
	int iP[7];
	int coordP[3];
	int coordA[3];
	
	int dimS[3];
	dimS[0] =*(dim_d)   -2*Ng;
	dimS[1] =*(dim_d+1) -2*Ng;
	dimS[2] =*(dim_d+2) -2*Ng;
	
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	for (int i = index; i < dimStencil; i += stride)
	{
		convert_1Dto3D(i,dimS,coordA);
		coordP[0]=coordA[0]+Ng;
		coordP[1]=coordA[1]+Ng;
		coordP[2]=coordA[2]+Ng;

		indexP(iP, dim_d, coordP);
		indexA(iA, dimS, coordA);

		res2[i]=powf(abs(-p_1[iP[0]]*A_d[iA[6]] +    
		        A_d[iA[0]] * p_1[iP[2]] + A_d[iA[1]] * p_1[iP[1]] +  
		        A_d[iA[2]] * p_1[iP[4]] + A_d[iA[3]] * p_1[iP[3]] +  
		        A_d[iA[4]] * p_1[iP[6]] + A_d[iA[5]] * p_1[iP[5]] + A_d[iA[7]]), __int2float_rd(norm));
		
		        
	}
}

__global__
void catchDivergenceA(int dimA, double *A_d, int catchDiv)
{	
	
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	for (int i = index; i < dimA; i += stride)
	{
		if (isnan(A_d[i])){catchDiv = 1;}
		        
	}
}

__global__
void catchDivergenceP(int dimp, double *p_d, int catchDiv)
{	
	
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	for (int i = index; i < dimp; i += stride)
	{
		if (isnan(p_d[i])){catchDiv = 1;}
		        
	}
}


/////////// C-CUDA functions ////////////
extern "C" void linearCUDA_(double *A, double *p, int* dimM, int *Ng, double* maxError, double* beta, int *maxit, int* it, int* ierr, int* norm, double* tres2)
{
	static double *A_d=NULL, *p_1=NULL, *p_2=NULL, *res2=NULL;
	static int *dim_d=NULL;
	int blockSize = 512;
	int dimMA[3] = {dimM[0]-2*(*Ng),dimM[1]-2*(*Ng),dimM[2]-2*(*Ng)};
	int dimA = (dimMA[0]*dimMA[1]*dimMA[2]*8);
	int dimp = dimM[0]*dimM[1]*dimM[2];
	int dimStencil = (dimMA[0]*dimMA[1]*dimMA[2]);
	double *res2_reduced;
	double init_res2;
	int *CoiA;
	//~ bool catchDiv = false;
	int catchDiv=0;
	//~ printGPUprops();
	clock_t t1=0.0, t2=0.0, t3=0.0, t4=0.0, t5=0.0, t0=0.0, t6=0.0;
		
	int numBlocks = (dimStencil + blockSize - 1) / blockSize;
	if (numBlocks>65535){numBlocks=65535;}
	
	t0 = clock();
//  // ----------------------Allocate Unified Memory – accessible from CPU or GPU
	if (p_1==NULL){
	cudaMallocManaged(&A_d, dimA*sizeof(double));
	cudaMallocManaged(&p_1, dimp*sizeof(double));
	cudaMallocManaged(&p_2, dimp*sizeof(double));
	cudaMallocManaged(&res2, dimStencil*sizeof(double));
	cudaMallocManaged(&dim_d, 3*sizeof(int));
	cudaMallocManaged(&CoiA, dimStencil*13*sizeof(int));
	}
	
	cudaMemcpy(A_d, A, dimA*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(p_2, p, dimp*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(p_1, p, dimp*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dim_d, dimM, 3*sizeof(int), cudaMemcpyHostToDevice);
	t1=t1 + clock() - t0;

//  // ------------------------------Run kernel on the GPU
	for ((*it)=0;(*it)<(*maxit);(*it)++){
		// Calls one iteration step alternating the matrixes
		if ((*it)%2==0){
			// Calls one iteration step output p_2
			t0 = clock();
			p_iter<<<numBlocks, blockSize>>>(dim_d, dimp, dimStencil, (*Ng), A_d, p_1, p_2, (*beta));
			cudaDeviceSynchronize();
			t2 = t2 +clock() -t0;
			t0 = clock();
			// impose the BC through Fortran on p_2
			apply_bc_cuda_(p_2);
			cudaDeviceSynchronize();
			t3 = t3 +clock() -t0;
			t0 = clock();
			// Compute the res2 on p_2 for each element and create the vector res_2
			L1L2_norm<<<numBlocks, blockSize>>>(dim_d, dimStencil, (*Ng), A_d, p_2, res2, (*norm));
			cudaDeviceSynchronize();
			t4 = t4 +clock() -t0;

		}
		else if((*it)%2==1){
			t0 = clock();
			// Calls one iteration step output p_1
			p_iter<<<numBlocks, blockSize>>>(dim_d, dimp, dimStencil, (*Ng), A_d, p_2, p_1, (*beta));
			cudaDeviceSynchronize();
			t2 = t2 +clock() -t0;
			t0 = clock();
			// impose the BC through Fortran on p_1
			apply_bc_cuda_(p_1);
			cudaDeviceSynchronize();
			t3 = t3 +clock() -t0;
			t0 = clock();
			// Compute the res2 on p_2 for each element and create the vector res_2
			L1L2_norm<<<numBlocks, blockSize>>>(dim_d, dimStencil, (*Ng), A_d, p_1, res2, (*norm));
			cudaDeviceSynchronize();
			t4 = t4 +clock() -t0;
			t0 = clock();
		}
		// Reduces the res_2 vector to res2_reduced and divedes it for the number of elements TO IMPROVE
		init_res2 = 0.0;
		for (int i=0 ; i < dimStencil ; i++){
			init_res2=init_res2 + res2[i];
		}
		t5 = t5 +clock() -t0;
		init_res2=init_res2/dimStencil;
		res2_reduced=&init_res2;
		catchDivergenceA<<<numBlocks, blockSize>>>(dimA,A_d,catchDiv);
		cudaDeviceSynchronize();
		if     ((*it)%2==0){
			catchDivergenceP<<<numBlocks, blockSize>>>(dimp,p_2,catchDiv);
			cudaDeviceSynchronize();
			} 
		else if((*it)%2==1){
			catchDivergenceP<<<numBlocks, blockSize>>>(dimp,p_1,catchDiv);
			cudaDeviceSynchronize();
			} 
		
		if (catchDiv==1){catch_divergence_cuda_((*res2_reduced),(*ierr),(*it));}
		t0 = clock();
			
		cudaDeviceSynchronize();
		collect_res2_cuda_(res2_reduced,tres2,ierr,it);
		t6 = t6 +clock() -t0;
		if ((*norm)==2) {(*tres2)=sqrt((*tres2));}
		cudaDeviceSynchronize();
		
		if (*tres2<*maxError){break;}
		
	}
	
//  // -----------------------Wait for GPU to finish before accessing on host
	cudaDeviceSynchronize();
	

t0 = clock();
//  // -----------------------Synchronize the memory
	if     ((*it)%2==0){
		cudaMemcpy(p, p_2, dimp*sizeof(double), cudaMemcpyDeviceToHost);}
	else if((*it)%2==1){
		cudaMemcpy(p, p_1, dimp*sizeof(double), cudaMemcpyDeviceToHost);}
	cudaMemcpy(A, A_d, dimA*sizeof(double), cudaMemcpyDeviceToHost);
	t1 = t1 +clock() -t0;
	
//	printf("%f %f %f %f %f %f \n time elapsed",(double) t1/CLOCKS_PER_SEC,(double)t2/CLOCKS_PER_SEC,(double)t3/ CLOCKS_PER_SEC,(double)t4/ CLOCKS_PER_SEC,(double)t5/ CLOCKS_PER_SEC,(double)t6/ CLOCKS_PER_SEC);
//~ //  // --------------------------------Free memory
  //~ cudaFree(A_d);
  //~ cudaFree(p_1);
  //~ cudaFree(p_2);
  //~ cudaFree(res2);
  //~ cudaFree(dim_d);
  
}
