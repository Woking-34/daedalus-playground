///////////////////////////////////////////////////////////////////////////////////////////
//  Copyright (C) 2002 - 2015, Huamin Wang
//  All rights reserved.
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions
//  are met:
//     1. Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//     2. Redistributions in binary form must reproduce the above copyright
//        notice, this list of conditions and the following disclaimer in the
//        documentation and/or other materials provided with the distribution.
//     3. The names of its contributors may not be used to endorse or promote
//        products derived from this software without specific prior written
//        permission.
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//	NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////////////////
//  Class CUDA_PROJECTIVE_MESH
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef __WHMIN_CUDA_PROJECTIVE_MESH_H__
#define __WHMIN_CUDA_PROJECTIVE_MESH_H__
#include "../lib/TIMER.h"
#include "../lib/DYNAMIC_MESH.h"

#define GRAVITY			-9.8
#define RADIUS_SQUARED	0.00025


///////////////////////////////////////////////////////////////////////////////////////////
//  Control kernel
///////////////////////////////////////////////////////////////////////////////////////////
__global__ void Control_Kernel(float* X, float *more_fixed, float control_mag, const int number, const int select_v)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i>=number)	return;

	more_fixed[i]=0;
	if(select_v!=-1)
	{
		float dist2=0;
		dist2+=(X[i*3+0]-X[select_v*3+0])*(X[i*3+0]-X[select_v*3+0]);
		dist2+=(X[i*3+1]-X[select_v*3+1])*(X[i*3+1]-X[select_v*3+1]);
		dist2+=(X[i*3+2]-X[select_v*3+2])*(X[i*3+2]-X[select_v*3+2]);
		if(dist2<RADIUS_SQUARED)	more_fixed[i]=control_mag;		
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
//  Basic update kernel
///////////////////////////////////////////////////////////////////////////////////////////
__global__ void Update_Kernel(float* X, float* V, float* fixed, const float* more_fixed, const float damping, const float t, const int number, const float dir_x, const float dir_y, const float dir_z)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i>=number)	return;

	if(more_fixed[i]!=0)
	{
		X[i*3+0]+=dir_x;
		X[i*3+1]+=dir_y;
		X[i*3+2]+=dir_z;		
		V[i*3+0] =0;
		V[i*3+1] =0;
		V[i*3+2] =0;
		return;
	}

	if(fixed[i]!=0)
	{
		V[i*3+0] =0;
		V[i*3+1] =0;
		V[i*3+2] =0;
		return;
	}

	//Apply damping
	V[i*3+0]*=damping;
	V[i*3+1]*=damping;
	V[i*3+2]*=damping;
	//Apply gravity
	V[i*3+1]+=GRAVITY*t;

	//Position update
	X[i*3+0]=X[i*3+0]+V[i*3+0]*t;
	X[i*3+1]=X[i*3+1]+V[i*3+1]*t;
	X[i*3+2]=X[i*3+2]+V[i*3+2]*t;
}

///////////////////////////////////////////////////////////////////////////////////////////
//  Laplacian Damping
///////////////////////////////////////////////////////////////////////////////////////////
__global__ void Laplacian_Damping_Kernel(const float* V, float* next_V, const float* fixed, const float* more_fixed, const int* all_VV, const int* all_vv_num, const int number, const float rate)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i>=number)	return;

	next_V[i*3+0]=0;
	next_V[i*3+1]=0;
	next_V[i*3+2]=0;
	if(more_fixed[i]!=0)	return;
	if(fixed[i]!=0)			return;

	for(int index=all_vv_num[i]; index<all_vv_num[i+1]; index++)
	{
		int j=all_VV[index];
		next_V[i*3+0]+=V[j*3+0]-V[i*3+0];
		next_V[i*3+1]+=V[j*3+1]-V[i*3+1];
		next_V[i*3+2]+=V[j*3+2]-V[i*3+2];
	}
	next_V[i*3+0]=V[i*3+0]+next_V[i*3+0]*rate;
	next_V[i*3+1]=V[i*3+1]+next_V[i*3+1]*rate;
	next_V[i*3+2]=V[i*3+2]+next_V[i*3+2]*rate;
}

///////////////////////////////////////////////////////////////////////////////////////////
//  Constraint Kernel 0
///////////////////////////////////////////////////////////////////////////////////////////
__global__ void Constraint_0_Kernel(const float* X, float* init_B, float* new_VC, const float *fixed, const float* more_fixed, const float inv_t, const float spring_k, const int number)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i>=number)	return;

	float c=(1+fixed[i]+more_fixed[i])*inv_t*inv_t;
	init_B[i*3+0]=c*X[i*3+0];
	init_B[i*3+1]=c*X[i*3+1];
	init_B[i*3+2]=c*X[i*3+2];
	new_VC[i]+=c;
}

///////////////////////////////////////////////////////////////////////////////////////////
//  Constraint Kernel 1
///////////////////////////////////////////////////////////////////////////////////////////
__global__ void Constraint_1_Kernel(const float* X, const float* init_B, float* next_X, const int* all_VV, const float* all_VL, const float* all_VW, const float* new_VC, const int* all_vv_num, const float spring_k, const int number)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i>=number)	return;

	float b[3];
	b[0]=init_B[i*3+0];
	b[1]=init_B[i*3+1];
	b[2]=init_B[i*3+2];
	
	for(int index=all_vv_num[i]; index<all_vv_num[i+1]; index++)
	{
		int j=all_VV[index];

		// Remove the off-diagonal (Jacobi method)
		b[0]-=all_VW[index]*X[j*3+0];
		b[1]-=all_VW[index]*X[j*3+1];
		b[2]-=all_VW[index]*X[j*3+2];

		// Add the other part of b
		if(all_VL[index]==-1)	continue;
		float d[3];
		d[0]=X[i*3+0]-X[j*3+0];
		d[1]=X[i*3+1]-X[j*3+1];
		d[2]=X[i*3+2]-X[j*3+2];
		float new_L=spring_k*all_VL[index]/sqrtf(DOT(d, d));
		b[0]+=d[0]*new_L;
		b[1]+=d[1]*new_L;
		b[2]+=d[2]*new_L;
	}

	float c=1/new_VC[i];
	next_X[i*3+0]=b[0]*c;
	next_X[i*3+1]=b[1]*c;
	next_X[i*3+2]=b[2]*c;
}

///////////////////////////////////////////////////////////////////////////////////////////
//  Constraint Kernel 2
///////////////////////////////////////////////////////////////////////////////////////////
__global__ void Constraint_2_Kernel(float* prev_X, float* X, float* next_X, float omega, int number, float under_relax)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i>=number)	return;
	
	float displace[3];
	displace[0]=next_X[i*3+0]-X[i*3+0];
	displace[1]=next_X[i*3+1]-X[i*3+1];
	displace[2]=next_X[i*3+2]-X[i*3+2];

	next_X[i*3+0]=displace[0]*under_relax+X[i*3+0];
	next_X[i*3+1]=displace[1]*under_relax+X[i*3+1];
	next_X[i*3+2]=displace[2]*under_relax+X[i*3+2];

	next_X[i*3+0]=omega*(next_X[i*3+0]-prev_X[i*3+0])+prev_X[i*3+0];
	next_X[i*3+1]=omega*(next_X[i*3+1]-prev_X[i*3+1])+prev_X[i*3+1];
	next_X[i*3+2]=omega*(next_X[i*3+2]-prev_X[i*3+2])+prev_X[i*3+2];
}

///////////////////////////////////////////////////////////////////////////////////////////
//  Constraint Kernel 3
///////////////////////////////////////////////////////////////////////////////////////////
__global__ void Constraint_3_Kernel(float* X, float* init_B, float* V, const float* fixed, const float* more_fixed, float inv_t, int number)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i>=number)	return;

	float c=(1+fixed[i]+more_fixed[i])*inv_t*inv_t;
	V[i*3+0]+=(X[i*3+0]-init_B[i*3+0]/c)*inv_t;
	V[i*3+1]+=(X[i*3+1]-init_B[i*3+1]/c)*inv_t;
	V[i*3+2]+=(X[i*3+2]-init_B[i*3+2]/c)*inv_t;
}


///////////////////////////////////////////////////////////////////////////////////////////
//  Class
///////////////////////////////////////////////////////////////////////////////////////////
template <class TYPE>
class CUDA_PROJECTIVE_MESH: public DYNAMIC_MESH<TYPE>
{
public:
    using   BASE_MESH<TYPE>::max_number;
    using   BASE_MESH<TYPE>::number;
    using   BASE_MESH<TYPE>::t_number;
    using   BASE_MESH<TYPE>::X;
    using   BASE_MESH<TYPE>::M;
    using   BASE_MESH<TYPE>::T;
    using   BASE_MESH<TYPE>::VN;    
    using   MESH<TYPE>::e_number;
    using   MESH<TYPE>::E;
    using   MESH<TYPE>::ET;
    using   DYNAMIC_MESH<TYPE>::V;
    
	TYPE	cost[32];
	int		cost_ptr;
	TYPE	fps;

	TYPE	rho;
	TYPE	control_mag;
	TYPE	air_damping;
	int		lap_damping;
	TYPE	under_relax;

	// Edge length
	TYPE*	L;

	TYPE*	rest_X;	// rest X
	TYPE*	b;
	
	//The all neighborhood data structure needed to form the matrix
	int*	all_VV;
	TYPE*	all_VL;
	TYPE*	all_VW;
	TYPE*	all_VC;
	int*	all_vv_num;
	
	TYPE*	fixed;

	TYPE	spring_k;
	TYPE	bending_k;
	
	//CUDA data
	TYPE*	dev_X;
	TYPE*	dev_V;
	TYPE*	dev_next_X;		// next X		(for temporary storage)
	TYPE*	dev_prev_X;		// previous X	(for Chebyshev)
	TYPE*	dev_fixed;
	TYPE*	dev_more_fixed;
	TYPE*	dev_init_B;		// Initialized momentum condition in B

	int*	dev_T;

	int*	dev_all_VV;
	TYPE*	dev_all_VL;
	TYPE*	dev_all_VW;
	TYPE*	dev_all_VC;
	TYPE*	dev_new_VC;
	int*	dev_all_vv_num;
	

///////////////////////////////////////////////////////////////////////////////////////////
//  Constructor and Deconstructor
///////////////////////////////////////////////////////////////////////////////////////////
	CUDA_PROJECTIVE_MESH()
	{	
		rest_X		= new TYPE	[max_number*3];
		b			= new TYPE	[max_number*3];
		L			= new TYPE	[max_number*3];
		fixed		= new TYPE	[max_number  ];

		all_VV		= new int	[max_number*12];
		all_VL		= new TYPE	[max_number*12];
		all_VW		= new TYPE	[max_number*12];
		all_VC		= new TYPE	[max_number   ];	
		all_vv_num	= new int	[max_number   ];

		cost_ptr	= 0;
		fps			= 0;
		rho			= 0.9992;
		control_mag	= 400;
		air_damping	= 0.999;
		lap_damping = 4;
		under_relax	= 1;

		spring_k	= 3000000;	//5000000
		bending_k	= 1;

		dev_X		= 0;
		dev_V		= 0;
		dev_init_B	= 0;
		dev_next_X	= 0;
		dev_prev_X	= 0;

		dev_T		= 0;

		dev_fixed	= 0;
		dev_more_fixed = 0;
		
		dev_all_VV	= 0;
		dev_all_VL	= 0;
		dev_all_VW	= 0;
		dev_all_VC	= 0;
		dev_new_VC	= 0;
		dev_all_vv_num	= 0;
	}
	
	~CUDA_PROJECTIVE_MESH()
	{
		if(rest_X)		delete[] rest_X;
		if(b)			delete[] b;
		if(L)			delete[] L;
		if(fixed)		delete[] fixed;

		if(all_VV)		delete[] all_VV;
		if(all_VL)		delete[] all_VL;
		if(all_VW)		delete[] all_VW;
		if(all_VC)		delete[] all_VC;
		if(all_vv_num)	delete[] all_vv_num;

		if(dev_X)			cudaFree(dev_X);
		if(dev_V)			cudaFree(dev_V);
		if(dev_next_X)		cudaFree(dev_next_X);
		if(dev_prev_X)		cudaFree(dev_prev_X);

		if(dev_T)			cudaFree(dev_T);

		if(dev_fixed)		cudaFree(dev_fixed);
		if(dev_more_fixed)	cudaFree(dev_more_fixed);
		if(dev_init_B)		cudaFree(dev_init_B);

		if(dev_all_VV)		cudaFree(dev_all_VV);
		if(dev_all_VL)		cudaFree(dev_all_VL);
		if(dev_all_VW)		cudaFree(dev_all_VW);
		if(dev_all_VC)		cudaFree(dev_all_VC);
		if(dev_new_VC)		cudaFree(dev_new_VC);
		if(dev_all_vv_num)	cudaFree(dev_all_vv_num);
	}

///////////////////////////////////////////////////////////////////////////////////////////
//  Initialization functions
///////////////////////////////////////////////////////////////////////////////////////////
	void Initialize(TYPE t)
	{
		memcpy(rest_X, X, sizeof(TYPE)*number*3);

		//Initialize connectivity
		this->Build_Connectivity();
		
		//Initialize edge lengths
		for(int e=0; e<e_number; e++)
		{
			int i=E[e*2+0];
			int j=E[e*2+1];
			L[e]=Distance(&X[i*3], &X[j*3]);
		}

		Build_All_Neighborhood(t);
		Allocate_GPU_Memory();
		Build_Matrix(rest_X);
	}

	void Allocate_GPU_Memory()
	{
		//Allocate CUDA memory
		cudaMalloc((void**)&dev_X,				sizeof(TYPE)*3*number);
		cudaMalloc((void**)&dev_V,				sizeof(TYPE)*3*number);
		cudaMalloc((void**)&dev_next_X,			sizeof(TYPE)*3*number);
		cudaMalloc((void**)&dev_prev_X,			sizeof(TYPE)*3*number);
		cudaMalloc((void**)&dev_fixed,			sizeof(TYPE)*number);
		cudaMalloc((void**)&dev_more_fixed,		sizeof(TYPE)*number);
		cudaMalloc((void**)&dev_init_B,			sizeof(TYPE)*3*number);
		cudaMalloc((void**)&dev_T,				sizeof(int )*3*t_number);

		cudaMalloc((void**)&dev_all_VV,			sizeof(int )*all_vv_num[number]);
		cudaMalloc((void**)&dev_all_VL,			sizeof(TYPE)*all_vv_num[number]);
		cudaMalloc((void**)&dev_all_VW,			sizeof(TYPE)*all_vv_num[number]);
		cudaMalloc((void**)&dev_all_VC,			sizeof(TYPE)*number);
		cudaMalloc((void**)&dev_new_VC,			sizeof(TYPE)*number);
		cudaMalloc((void**)&dev_all_vv_num,		sizeof(int)*(number+1));

		//Copy data into CUDA memory
		cudaMemcpy(dev_X,			X,			sizeof(TYPE)*3*number,				cudaMemcpyHostToDevice);
		cudaMemcpy(dev_V,			V,			sizeof(TYPE)*3*number,				cudaMemcpyHostToDevice);
		cudaMemcpy(dev_prev_X,		X,			sizeof(TYPE)*3*number,				cudaMemcpyHostToDevice);
		cudaMemcpy(dev_next_X,		X,			sizeof(TYPE)*3*number,				cudaMemcpyHostToDevice);
		cudaMemcpy(dev_fixed,		fixed,		sizeof(TYPE)*number,				cudaMemcpyHostToDevice);		
		cudaMemcpy(dev_T,			T,			sizeof(int )*3*t_number,			cudaMemcpyHostToDevice);
		
		cudaMemset(dev_more_fixed,  0,			sizeof(TYPE)*number);

		cudaMemcpy(dev_all_VV,		all_VV,		sizeof(int )*all_vv_num[number],	cudaMemcpyHostToDevice);
		cudaMemcpy(dev_all_VL,		all_VL,		sizeof(TYPE)*all_vv_num[number],	cudaMemcpyHostToDevice);
		cudaMemcpy(dev_all_VW,		all_VW,		sizeof(TYPE)*all_vv_num[number],	cudaMemcpyHostToDevice);
		cudaMemcpy(dev_all_VC,		all_VC,		sizeof(TYPE)*number,				cudaMemcpyHostToDevice);
		cudaMemcpy(dev_all_vv_num,  all_vv_num, sizeof(int)*(number+1),				cudaMemcpyHostToDevice);
	}
	
	void Build_All_Neighborhood(const float t)
	{
		//Step 1: create all edges, including original and bending edges
		int* all_E=new int[e_number*8];
		int  all_e_number=0;
		for(int e=0; e<e_number; e++)
		{
			//Add original edges
			all_E[all_e_number*2+0]=E[e*2+0];
			all_E[all_e_number*2+1]=E[e*2+1];
			all_e_number++;
			all_E[all_e_number*2+0]=E[e*2+1];
			all_E[all_e_number*2+1]=E[e*2+0];
			all_e_number++;

			//Add bending edges
			int t0=ET[e*2+0];
			int t1=ET[e*2+1];
			if(t0==-1 || t1==-1)	continue;
			int v2=T[t0*3+0]+T[t0*3+1]+T[t0*3+2]-E[e*2+0]-E[e*2+1];
			int v3=T[t1*3+0]+T[t1*3+1]+T[t1*3+2]-E[e*2+0]-E[e*2+1];			
			all_E[all_e_number*2+0]=v2;
			all_E[all_e_number*2+1]=v3;
			all_e_number++;
			all_E[all_e_number*2+0]=v3;
			all_E[all_e_number*2+1]=v2;
			all_e_number++;
		}
		Quick_Sort_BE(all_E, 0, all_e_number-1);

		//Step 2: Set all_vv_num and all_VV
		int e=0;
		int all_vv_ptr=0;
		for(int i=0; i<number; i++)
		{
			all_vv_num[i]=all_vv_ptr;
			for(; e<all_e_number; e++)
			{			
				if(all_E[e*2]!=i)						break;		// not in the right vertex
				if(e!=0 && all_E[e*2+1]==all_E[e*2-1])	continue;	// duplicate
				all_VV[all_vv_ptr++]=all_E[e*2+1];
			}
		}
		all_vv_num[number]=all_vv_ptr;
		delete[] all_E;
	}

	void Build_Matrix(float* X)
	{
		for(int i=0; i<all_vv_num[number]; i++)		all_VL[i]=-1;
		for(int i=0; i<all_vv_num[number]; i++)		all_VW[i]= 0;
		for(int i=0; i<number; i++)					all_VC[i]= 0;
		
		for(int e=0; e<e_number; e++)
		{
			int v[4];
			v[0]=E[e*2+0];
			v[1]=E[e*2+1];
			
			//First, handle spring length
			TYPE l=Distance(&X[E[e*2+0]*3], &X[E[e*2+1]*3]);
			all_VL[Find_Neighbor(v[0], v[1])]=l;
			all_VL[Find_Neighbor(v[1], v[0])]=l;
			all_VC[v[0]]+=spring_k;
			all_VC[v[1]]+=spring_k;
			all_VW[Find_Neighbor(v[0], v[1])]-=spring_k;
			all_VW[Find_Neighbor(v[1], v[0])]-=spring_k;

			//Next, handle bending weights
			int t0=ET[e*2+0];
			int t1=ET[e*2+1];
			if(t0==-1 || t1==-1)	continue;
			v[2]=T[t0*3+0]+T[t0*3+1]+T[t0*3+2]-v[0]-v[1];
			v[3]=T[t1*3+0]+T[t1*3+1]+T[t1*3+2]-v[0]-v[1];
			TYPE c01=Cotangent(&X[v[0]*3], &X[v[1]*3], &X[v[2]*3]);
			TYPE c02=Cotangent(&X[v[0]*3], &X[v[1]*3], &X[v[3]*3]);
			TYPE c03=Cotangent(&X[v[1]*3], &X[v[0]*3], &X[v[2]*3]);
			TYPE c04=Cotangent(&X[v[1]*3], &X[v[0]*3], &X[v[3]*3]);			
			TYPE area0=sqrt(Area_Squared(&X[v[0]*3], &X[v[1]*3], &X[v[2]*3]));
			TYPE area1=sqrt(Area_Squared(&X[v[0]*3], &X[v[1]*3], &X[v[3]*3]));
			TYPE weight=1/(area0+area1);
			TYPE k[4];
			k[0]= c03+c04;
			k[1]= c01+c02;
			k[2]=-c01-c03;
			k[3]=-c02-c04;
						
			for(int i=0; i<4; i++)
			for(int j=0; j<4; j++)
			{
				if(i==j)	all_VC[v[i]]+=k[i]*k[j]*bending_k*weight;
				else		all_VW[Find_Neighbor(v[i], v[j])]+=k[i]*k[j]*bending_k*weight;
			}
		}

		cudaMemcpy(dev_all_VL,		all_VL,		sizeof(TYPE)*all_vv_num[number],	cudaMemcpyHostToDevice);
		cudaMemcpy(dev_all_VW,		all_VW,		sizeof(TYPE)*all_vv_num[number],	cudaMemcpyHostToDevice);
		cudaMemcpy(dev_all_VC,		all_VC,		sizeof(TYPE)*number,				cudaMemcpyHostToDevice);
	}

	int Find_Neighbor(int i, int j)
	{
		for(int index=all_vv_num[i]; index<all_vv_num[i+1]; index++)
			if(all_VV[index]==j)	return index;
		printf("ERROR: failed to find the neighbor in all_VV.\n"); getchar();
		return -1;
	}

	void Quick_Sort_BE(int a[], int l, int r)
	{				
		if(l>=r)	return;
		int j=Quick_Sort_Partition_BE(a, l, r);
		Quick_Sort_BE(a, l, j-1);
		Quick_Sort_BE(a, j+1, r);		
	}
	
	int Quick_Sort_Partition_BE(int a[], int l, int r) 
	{
		int pivot[2], i, j;
		pivot[0] = a[l*2+0];
		pivot[1] = a[l*2+1];
		i = l; j = r+1;		
		while( 1)
		{
			do ++i; while( (a[i*2]<pivot[0] || a[i*2]==pivot[0] && a[i*2+1]<=pivot[1]) && i <= r );
			do --j; while(  a[j*2]>pivot[0] || a[j*2]==pivot[0] && a[j*2+1]> pivot[1] );
			if( i >= j ) break;
			//Swap i and j
			Swap(a[i*2+0], a[j*2+0]);
			Swap(a[i*2+1], a[j*2+1]);
		}
		//Swap l and j
		Swap(a[l*2+0], a[j*2+0]);
		Swap(a[l*2+1], a[j*2+1]);
		return j;
	} 

///////////////////////////////////////////////////////////////////////////////////////////
//  Control function
///////////////////////////////////////////////////////////////////////////////////////////
	void Reset_More_Fixed(int select_v)
	{
		int threadsPerBlock = 64;
		int blocksPerGrid = (number + threadsPerBlock - 1) / threadsPerBlock;
		//Need it here
		Control_Kernel << <blocksPerGrid, threadsPerBlock>> >(dev_X, dev_more_fixed, control_mag, number, select_v);
	}

///////////////////////////////////////////////////////////////////////////////////////////
//  Simulation functions
///////////////////////////////////////////////////////////////////////////////////////////
	void Apply_Constraints(TYPE t)
	{}
	
	void Update(TYPE t)
	{}

	void Update(TYPE t, int iterations, TYPE dir[])
	{	
		TIMER timer;
		int threadsPerBlock = 256;
		int blocksPerGrid = (number + threadsPerBlock - 1) / threadsPerBlock;
		
		// Damping + Basic update
		for(int l=0; l<lap_damping; l++)
		{
			Laplacian_Damping_Kernel<< <blocksPerGrid, threadsPerBlock>> >(dev_V, dev_next_X, dev_fixed, dev_more_fixed, dev_all_VV, dev_all_vv_num, number, 0.1);		
			cudaMemcpy(dev_V, dev_next_X, sizeof(TYPE)*number*3, cudaMemcpyDeviceToDevice);
		}
		Update_Kernel << <blocksPerGrid, threadsPerBlock>> >(dev_X, dev_V, dev_fixed, dev_more_fixed, air_damping, t, number, dir[0], dir[1], dir[2]);


		//Set up data
		cudaMemcpy(dev_new_VC, dev_all_VC, sizeof(TYPE)*number, cudaMemcpyDeviceToDevice);
		Constraint_0_Kernel << <blocksPerGrid, threadsPerBlock>> >(dev_X, dev_init_B, dev_new_VC, dev_fixed, dev_more_fixed, 1/t, spring_k, number);
		cudaMemcpy(dev_prev_X, dev_X, sizeof(TYPE)*3*number, cudaMemcpyDeviceToDevice);
		
		TYPE omega;
		for(int l=0; l<iterations; l++)
		{
			Constraint_1_Kernel << <blocksPerGrid, threadsPerBlock>> >(dev_X, dev_init_B, dev_next_X, dev_all_VV, dev_all_VL, dev_all_VW, dev_new_VC, dev_all_vv_num, spring_k, number);
									
			if(l<=10)		omega=1;
			else if(l==11)	omega=2/(2-rho*rho);
			else			omega=4/(4-rho*rho*omega);
			//omega=1;

			Constraint_2_Kernel<< <blocksPerGrid, threadsPerBlock>> >(dev_prev_X, dev_X, dev_next_X, omega, number, under_relax);
			Swap(dev_X, dev_prev_X);
			Swap(dev_X, dev_next_X);
		}
		Constraint_3_Kernel<< <blocksPerGrid, threadsPerBlock>> >(dev_X, dev_init_B, dev_V, dev_fixed, dev_more_fixed, 1/t, number);
		cudaMemcpy(X, dev_X, sizeof(TYPE)*3*number, cudaMemcpyDeviceToHost);


		cost[cost_ptr]=timer.Get_Time();
		cost_ptr=(cost_ptr+1)%8;
		fps=8/(cost[0]+cost[1]+cost[2]+cost[3]+cost[4]+cost[5]+cost[6]+cost[7]);
	}

};


#endif
