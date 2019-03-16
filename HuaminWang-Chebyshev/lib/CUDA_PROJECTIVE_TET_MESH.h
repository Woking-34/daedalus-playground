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
//  Class CUDA_PROJECTIVE_TET_MESH
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef	__WHMIN_CUDA_PROJECTIVE_TET_MESH_H__
#define __WHMIN_CUDA_PROJECTIVE_TET_MESH_H__
#include "TET_MESH.h"
#include "TIMER.h"

#define GRAVITY			-9.8
#define RADIUS_SQUARED	0.002//0.01


///////////////////////////////////////////////////////////////////////////////////////////
//  math kernels
///////////////////////////////////////////////////////////////////////////////////////////
__device__ void dev_Matrix_Product_3(const float *A, const float *B, float *R)				//R=A*B
{
	R[0]=A[0]*B[0]+A[1]*B[3]+A[2]*B[6];
	R[1]=A[0]*B[1]+A[1]*B[4]+A[2]*B[7];
	R[2]=A[0]*B[2]+A[1]*B[5]+A[2]*B[8];
	R[3]=A[3]*B[0]+A[4]*B[3]+A[5]*B[6];
	R[4]=A[3]*B[1]+A[4]*B[4]+A[5]*B[7];
	R[5]=A[3]*B[2]+A[4]*B[5]+A[5]*B[8];
	R[6]=A[6]*B[0]+A[7]*B[3]+A[8]*B[6];
	R[7]=A[6]*B[1]+A[7]*B[4]+A[8]*B[7];
	R[8]=A[6]*B[2]+A[7]*B[5]+A[8]*B[8];
}

__device__ void dev_Matrix_Substract_3(float *A, float *B, float *R)						//R=A-B
{
	for(int i=0; i<9; i++)	R[i]=A[i]-B[i];
}

__device__ void dev_Matrix_Product(float *A, float *B, float *R, int nx, int ny, int nz)	//R=A*B
{
	memset(R, 0, sizeof(float)*nx*nz);
	for(int i=0; i<nx; i++)
	for(int j=0; j<nz; j++)
	for(int k=0; k<ny; k++)
		R[i*nz+j]+=A[i*ny+k]*B[k*nz+j];
}

__device__ void Get_Rotation(float F[3][3], float R[3][3])
{
    float C[3][3];
    memset(&C[0][0], 0, sizeof(float)*9);
    for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    for(int k=0; k<3; k++)
        C[i][j]+=F[k][i]*F[k][j];
    
    float C2[3][3];
    memset(&C2[0][0], 0, sizeof(float)*9);
    for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    for(int k=0; k<3; k++)
        C2[i][j]+=C[i][k]*C[j][k];
    
    float det    =   F[0][0]*F[1][1]*F[2][2]+
                    F[0][1]*F[1][2]*F[2][0]+
                    F[1][0]*F[2][1]*F[0][2]-
                    F[0][2]*F[1][1]*F[2][0]-
                    F[0][1]*F[1][0]*F[2][2]-
                    F[0][0]*F[1][2]*F[2][1];
    
    float I_c    =   C[0][0]+C[1][1]+C[2][2];
    float I_c2   =   I_c*I_c;
    float II_c   =   0.5*(I_c2-C2[0][0]-C2[1][1]-C2[2][2]);
    float III_c  =   det*det;
    float k      =   I_c2-3*II_c;
    
    float inv_U[3][3];
    if(k<1e-10f)
    {
        float inv_lambda=1/sqrt(I_c/3);
        memset(inv_U, 0, sizeof(float)*9);
        inv_U[0][0]=inv_lambda;
        inv_U[1][1]=inv_lambda;
        inv_U[2][2]=inv_lambda;
    }
    else
    {
        float l = I_c*(I_c*I_c-4.5*II_c)+13.5*III_c;
        float k_root = sqrt(k);
        float value=l/(k*k_root);
        if(value<-1.0) value=-1.0;
        if(value> 1.0) value= 1.0;
        float phi = acos(value);
        float lambda2=(I_c+2*k_root*cos(phi/3))/3.0;
        float lambda=sqrt(lambda2);
        
        float III_u = sqrt(III_c);
        if(det<0)   III_u=-III_u;
        float I_u = lambda + sqrt(-lambda2 + I_c + 2*III_u/lambda);
        float II_u=(I_u*I_u-I_c)*0.5;
        
        float U[3][3];
        float inv_rate, factor;
        
        inv_rate=1/(I_u*II_u-III_u);
        factor=I_u*III_u*inv_rate;
        
        memset(U, 0, sizeof(float)*9);
        U[0][0]=factor;
        U[1][1]=factor;
        U[2][2]=factor;
        
        factor=(I_u*I_u-II_u)*inv_rate;
        for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            U[i][j]+=factor*C[i][j]-inv_rate*C2[i][j];
        
        inv_rate=1/III_u;
        factor=II_u*inv_rate;
        memset(inv_U, 0, sizeof(float)*9);
        inv_U[0][0]=factor;
        inv_U[1][1]=factor;
        inv_U[2][2]=factor;
        
        factor=-I_u*inv_rate;
        for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            inv_U[i][j]+=factor*U[i][j]+inv_rate*C[i][j];
    }
    
    memset(&R[0][0], 0, sizeof(float)*9);
    for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    for(int k=0; k<3; k++)
        R[i][j]+=F[i][k]*inv_U[k][j];    
}

///////////////////////////////////////////////////////////////////////////////////////////
//  Control kernel
///////////////////////////////////////////////////////////////////////////////////////////
__global__ void Control_Kernel(float* X, float *more_fixed, const float control_mag, const int number, const int select_v)
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
__global__ void Update_Kernel(float* X, float* V, const float *fixed, const float *more_fixed, const float damping, const float t, const int number, const float dir_x, const float dir_y, const float dir_z)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i>=number)	return;

	if(fixed[i]!=0)	return;

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

	//Apply damping
	V[i*3+0]*=damping;
	V[i*3+1]*=damping;
	V[i*3+2]*=damping;
	//Apply gravity
	V[i*3+1]+=GRAVITY*t;
	//Position update
	X[i*3+0]+=V[i*3+0]*t;
	X[i*3+1]+=V[i*3+1]*t;
	X[i*3+2]+=V[i*3+2]*t;	
}

///////////////////////////////////////////////////////////////////////////////////////////
//  Tet Constraint Kernel
///////////////////////////////////////////////////////////////////////////////////////////
__global__ void Tet_Constraint_Kernel(const float* X, const int* Tet, const float* inv_Dm, const float* Vol, float* Tet_Temp, const float elasticity, const int tet_number, const int l)
{
	int t = blockDim.x * blockIdx.x + threadIdx.x;
	if(t>=tet_number)	return;
	
	int p0=Tet[t*4+0]*3;
	int p1=Tet[t*4+1]*3;
	int p2=Tet[t*4+2]*3;
	int p3=Tet[t*4+3]*3;

	const float* idm=&inv_Dm[t*9];

	float Ds[9];
	Ds[0]=X[p1+0]-X[p0+0];
	Ds[3]=X[p1+1]-X[p0+1];
	Ds[6]=X[p1+2]-X[p0+2];
	Ds[1]=X[p2+0]-X[p0+0];
	Ds[4]=X[p2+1]-X[p0+1];
	Ds[7]=X[p2+2]-X[p0+2];
	Ds[2]=X[p3+0]-X[p0+0];
	Ds[5]=X[p3+1]-X[p0+1];
	Ds[8]=X[p3+2]-X[p0+2];

	float F[9], R[9], B[3], C[9];
	float new_R[9];
	dev_Matrix_Product_3(Ds, idm, F);
	
	Get_Rotation((float (*)[3])F, (float (*)[3])new_R);
	
	float half_matrix[3][4], result_matrix[3][4];
	half_matrix[0][0]=-idm[0]-idm[3]-idm[6];
	half_matrix[0][1]= idm[0];
	half_matrix[0][2]= idm[3];
	half_matrix[0][3]= idm[6];
	half_matrix[1][0]=-idm[1]-idm[4]-idm[7];
	half_matrix[1][1]= idm[1];
	half_matrix[1][2]= idm[4];
	half_matrix[1][3]= idm[7];
	half_matrix[2][0]=-idm[2]-idm[5]-idm[8];
	half_matrix[2][1]= idm[2];
	half_matrix[2][2]= idm[5];
	half_matrix[2][3]= idm[8];

	dev_Matrix_Substract_3(new_R, F, new_R);
	dev_Matrix_Product(new_R, &half_matrix[0][0], &result_matrix[0][0], 3, 3, 4);
			
	float rate=Vol[t]*elasticity;
	Tet_Temp[t*12+ 0]=result_matrix[0][0]*rate;
	Tet_Temp[t*12+ 1]=result_matrix[1][0]*rate;
	Tet_Temp[t*12+ 2]=result_matrix[2][0]*rate;
	Tet_Temp[t*12+ 3]=result_matrix[0][1]*rate;
	Tet_Temp[t*12+ 4]=result_matrix[1][1]*rate;
	Tet_Temp[t*12+ 5]=result_matrix[2][1]*rate;
	Tet_Temp[t*12+ 6]=result_matrix[0][2]*rate;
	Tet_Temp[t*12+ 7]=result_matrix[1][2]*rate;
	Tet_Temp[t*12+ 8]=result_matrix[2][2]*rate;
	Tet_Temp[t*12+ 9]=result_matrix[0][3]*rate;
	Tet_Temp[t*12+10]=result_matrix[1][3]*rate;
	Tet_Temp[t*12+11]=result_matrix[2][3]*rate;
}

///////////////////////////////////////////////////////////////////////////////////////////
//  Constraint Kernel 0
///////////////////////////////////////////////////////////////////////////////////////////
__global__ void Constraint_0_Kernel(const float* X, float* init_B, float* VC, const float* fixed, const float* more_fixed, const float inv_t, const int number)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i>=number)	return;

	VC[i]=(1+fixed[i]+more_fixed[i])*inv_t*inv_t;
	init_B[i*3+0]=VC[i]*X[i*3+0];
	init_B[i*3+1]=VC[i]*X[i*3+1];
	init_B[i*3+2]=VC[i]*X[i*3+2];
}

///////////////////////////////////////////////////////////////////////////////////////////
//  Constraint Kernel 1
///////////////////////////////////////////////////////////////////////////////////////////
__global__ void Constraint_1_Kernel(const float* X, const float* init_B, const float* VC, float* next_X, const float* Tet_Temp, const float* MD, const int* VTT, const int* vtt_num, const int number)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i>=number)	return;

	double b[3];
	b[0]=init_B[i*3+0]+MD[i]*X[i*3+0];
	b[1]=init_B[i*3+1]+MD[i]*X[i*3+1];
	b[2]=init_B[i*3+2]+MD[i]*X[i*3+2];
			
	for(int index=vtt_num[i]; index<vtt_num[i+1]; index++)
	{
		b[0]+=Tet_Temp[VTT[index]*3+0];
		b[1]+=Tet_Temp[VTT[index]*3+1];
		b[2]+=Tet_Temp[VTT[index]*3+2];
	}

	next_X[i*3+0]=b[0]/(VC[i]+MD[i]);
	next_X[i*3+1]=b[1]/(VC[i]+MD[i]);
	next_X[i*3+2]=b[2]/(VC[i]+MD[i]);
}

///////////////////////////////////////////////////////////////////////////////////////////
//  Constraint Kernel 2
///////////////////////////////////////////////////////////////////////////////////////////
__global__ void Constraint_2_Kernel(float* prev_X, float* X, float* next_X, float omega, int number)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i>=number)	return;	

	next_X[i*3+0]=(next_X[i*3+0]-X[i*3+0])*0.666+X[i*3+0];
	next_X[i*3+1]=(next_X[i*3+1]-X[i*3+1])*0.666+X[i*3+1];
	next_X[i*3+2]=(next_X[i*3+2]-X[i*3+2])*0.666+X[i*3+2];

	next_X[i*3+0]=omega*(next_X[i*3+0]-prev_X[i*3+0])+prev_X[i*3+0];
	next_X[i*3+1]=omega*(next_X[i*3+1]-prev_X[i*3+1])+prev_X[i*3+1];
	next_X[i*3+2]=omega*(next_X[i*3+2]-prev_X[i*3+2])+prev_X[i*3+2];
}

///////////////////////////////////////////////////////////////////////////////////////////
//  Constraint Kernel 3
///////////////////////////////////////////////////////////////////////////////////////////
__global__ void Constraint_3_Kernel(float* X, float* init_B, float* V, const float *fixed, const float *more_fixed, float inv_t, int number)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i>=number)	return;

	float c=(1+fixed[i]+more_fixed[i])*inv_t*inv_t;
	V[i*3+0]+=(X[i*3+0]-init_B[i*3+0]/c)*inv_t;
	V[i*3+1]+=(X[i*3+1]-init_B[i*3+1]/c)*inv_t;
	V[i*3+2]+=(X[i*3+2]-init_B[i*3+2]/c)*inv_t;
}



///////////////////////////////////////////////////////////////////////////////////////////
//  class CUDA_PROJECTIVE_TET_MESH
///////////////////////////////////////////////////////////////////////////////////////////
template <class TYPE>
class CUDA_PROJECTIVE_TET_MESH: public TET_MESH<TYPE> 
{
public:
	TYPE	cost[8];
	int		cost_ptr;
	TYPE	fps;
	TYPE*	old_X;
	TYPE*	V;
	TYPE*	fixed;

	TYPE	rho;
	TYPE	elasticity;
	TYPE	control_mag;
	TYPE	damping;

	TYPE*	MD;			//matrix diagonal
	TYPE*	TQ;
	TYPE*	Tet_Temp;
	int*	VTT;		//The index list mapping to Tet_Temp
	int*	vtt_num;

	//CUDA data
	TYPE*	dev_X;
	TYPE*	dev_E;
	TYPE*	dev_V;
	TYPE*	dev_next_X;		// next X		(for temporary storage)
	TYPE*	dev_prev_X;		// previous X	(for Chebyshev)
	TYPE*	dev_fixed;
	TYPE*	dev_more_fixed;
	TYPE*	dev_init_B;		// Initialized momentum condition in B

	TYPE*	dev_Dm;
	TYPE*	dev_inv_Dm;
	TYPE*	dev_Vol;
	int*	dev_Tet;
	TYPE*	dev_TQ;
	TYPE*	dev_Tet_Temp;
	TYPE*	dev_VC;
	TYPE*	dev_MD;
	int*	dev_VTT;
	int*	dev_vtt_num;

	TYPE*	error;
	TYPE*	dev_error;

	CUDA_PROJECTIVE_TET_MESH()
	{
		cost_ptr= 0;

		old_X	= new TYPE	[max_number*3];
		V		= new TYPE	[max_number*3];
		fixed	= new TYPE	[max_number  ];

		MD		= new TYPE	[max_number  ];
		TQ		= new TYPE	[max_number*4];
		Tet_Temp= new TYPE	[max_number*24];

		VTT		= new int	[max_number*4];
		vtt_num	= new int	[max_number  ];

		error	= new TYPE	[max_number*3];

		fps			= 0;
		elasticity	= 3000000; //5000000
		control_mag	= 10;
		rho			= 0.9992;
		damping		= 0.9995;

		memset(		V, 0, sizeof(TYPE)*max_number*3);
		memset(	fixed, 0, sizeof(int )*max_number  );

		// GPU data
		dev_X			= 0;
		dev_E			= 0;
		dev_V			= 0;
		dev_next_X		= 0;
		dev_prev_X		= 0;
		dev_fixed		= 0;
		dev_more_fixed	= 0;
		dev_init_B		= 0;

		dev_Dm			= 0;
		dev_inv_Dm		= 0;
		dev_Vol			= 0;
		dev_Tet			= 0;
		dev_TQ			= 0;
		dev_Tet_Temp	= 0;
		dev_VC			= 0;
		dev_MD			= 0;
		dev_VTT			= 0;
		dev_vtt_num		= 0;

		dev_error		= 0;
	}
	
	~CUDA_PROJECTIVE_TET_MESH()
	{
		if(old_X)			delete[] old_X;
		if(V)				delete[] V;
		if(fixed)			delete[] fixed;
		if(MD)				delete[] MD;
		if(TQ)				delete[] TQ;
		if(Tet_Temp)		delete[] Tet_Temp;
		if(VTT)				delete[] VTT;
		if(vtt_num)			delete[] vtt_num;
		if(error)			delete[] error;

		//GPU Data
		if(dev_X)			cudaFree(dev_X);
		if(dev_E)			cudaFree(dev_E);
		if(dev_V)			cudaFree(dev_V);
		if(dev_next_X)		cudaFree(dev_next_X);
		if(dev_prev_X)		cudaFree(dev_prev_X);
		if(dev_fixed)		cudaFree(dev_fixed);
		if(dev_more_fixed)	cudaFree(dev_more_fixed);
		if(dev_init_B)		cudaFree(dev_init_B);

		if(dev_Dm)			cudaFree(dev_Dm);
		if(dev_inv_Dm)		cudaFree(dev_inv_Dm);
		if(dev_Vol)			cudaFree(dev_Vol);
		if(dev_Tet)			cudaFree(dev_Tet);
		if(dev_TQ)			cudaFree(dev_TQ);
		if(dev_Tet_Temp)	cudaFree(dev_Tet_Temp);
		if(dev_VC)			cudaFree(dev_VC);
		if(dev_MD)			cudaFree(dev_MD);
		if(dev_VTT)			cudaFree(dev_VTT);
		if(dev_vtt_num)		cudaFree(dev_vtt_num);

		if(dev_error)		cudaFree(dev_error);
	}

///////////////////////////////////////////////////////////////////////////////////////////
//  Initialize functions
///////////////////////////////////////////////////////////////////////////////////////////
	void Initialize(TYPE t)
	{
		TET_MESH<TYPE>::Initialize();
		Initialize_MD();
		Build_VTT();

		for(int t=0; t<tet_number; t++)
		{
			TQ[t*4+0]=0;
			TQ[t*4+1]=0;
			TQ[t*4+2]=0;
			TQ[t*4+3]=1;
		}
		Allocate_GPU_Memory();
	}
	
	void Initialize_MD()
	{
		memset(MD, 0, sizeof(TYPE)*number);
		for(int t=0; t<tet_number; t++)
		{
			int*	v=&Tet[t*4];

			TYPE	idm[12];			
			memcpy(idm, &inv_Dm[t*9], sizeof(TYPE)*9);
			idm[ 9]=-(idm[0]+idm[3]+idm[6]);
			idm[10]=-(idm[1]+idm[4]+idm[7]);
			idm[11]=-(idm[2]+idm[5]+idm[8]);

			TYPE	M[12][12];
			for(int i=0; i<12; i++)
			{
				TYPE sum=0;
				for(int j=0; j<12; j++)
				{
					M[i][j]=idm[i]*idm[j];
					if(i!=j)	sum+=fabs(M[i][j]);
				}				
			}

			MD[v[0]]+=(idm[0]+idm[3]+idm[6])*(idm[0]+idm[3]+idm[6])*Vol[t]*elasticity;
			MD[v[0]]+=(idm[1]+idm[4]+idm[7])*(idm[1]+idm[4]+idm[7])*Vol[t]*elasticity;
			MD[v[0]]+=(idm[2]+idm[5]+idm[8])*(idm[2]+idm[5]+idm[8])*Vol[t]*elasticity;

			MD[v[1]]+=idm[0]*idm[0]*Vol[t]*elasticity;
			MD[v[1]]+=idm[1]*idm[1]*Vol[t]*elasticity;
			MD[v[1]]+=idm[2]*idm[2]*Vol[t]*elasticity;

			MD[v[2]]+=idm[3]*idm[3]*Vol[t]*elasticity;
			MD[v[2]]+=idm[4]*idm[4]*Vol[t]*elasticity;
			MD[v[2]]+=idm[5]*idm[5]*Vol[t]*elasticity;

			MD[v[3]]+=idm[6]*idm[6]*Vol[t]*elasticity;
			MD[v[3]]+=idm[7]*idm[7]*Vol[t]*elasticity;
			MD[v[3]]+=idm[8]*idm[8]*Vol[t]*elasticity;
		}		
	}

	void Build_VTT()
	{
		int* _VTT=new int[tet_number*8];
		for(int t=0; t<tet_number*4; t++)
		{
			_VTT[t*2+0]=Tet[t];
			_VTT[t*2+1]=t;
		}
		Quick_Sort_VTT(_VTT, 0, tet_number*4-1);
		
		for(int i=0, v=-1; i<tet_number*4; i++)
		{
			if(_VTT[i*2+0]!=v)	//start a new vertex
			{
				v=_VTT[i*2+0];
				vtt_num[v]=i;
			}
			VTT[i]=_VTT[i*2+1];
		}
		vtt_num[number]=tet_number*4;		
		delete[] _VTT;
	}	

	void Quick_Sort_VTT(int a[], int l, int r)
	{				
		if(l>=r)	return;
		int j=Quick_Sort_Partition_VTT(a, l, r);
		Quick_Sort_VTT(a, l, j-1);
		Quick_Sort_VTT(a, j+1, r);		
	}
	
	int Quick_Sort_Partition_VTT(int a[], int l, int r) 
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

	void Allocate_GPU_Memory()
	{
		//Allocate CUDA memory
		cudaMalloc((void**)&dev_X,			sizeof(int )*3*number);
		cudaMalloc((void**)&dev_E,			sizeof(int )*3*number);
		cudaMalloc((void**)&dev_V,			sizeof(TYPE)*3*number);
		cudaMalloc((void**)&dev_next_X,		sizeof(TYPE)*3*number);
		cudaMalloc((void**)&dev_prev_X,		sizeof(TYPE)*3*number);
		cudaMalloc((void**)&dev_fixed,		sizeof(TYPE)*  number);
		cudaMalloc((void**)&dev_more_fixed, sizeof(TYPE)*  number);
		cudaMalloc((void**)&dev_init_B,		sizeof(TYPE)*3*number);

		cudaMalloc((void**)&dev_Dm,			sizeof(int )*tet_number*9);
		cudaMalloc((void**)&dev_inv_Dm,		sizeof(int )*tet_number*9);
		cudaMalloc((void**)&dev_Vol,		sizeof(int )*tet_number);
		cudaMalloc((void**)&dev_Tet,		sizeof(int )*tet_number*4);
		cudaMalloc((void**)&dev_TQ,			sizeof(TYPE)*tet_number*4);
		cudaMalloc((void**)&dev_Tet_Temp,	sizeof(TYPE)*tet_number*12);
		
		cudaMalloc((void**)&dev_VC,			sizeof(TYPE)*number);
		cudaMalloc((void**)&dev_MD,			sizeof(TYPE)*number);
		cudaMalloc((void**)&dev_VTT,		sizeof(int )*tet_number*4);
		cudaMalloc((void**)&dev_vtt_num,	sizeof(int )*(number+1));

		cudaMalloc((void**)&dev_error,		sizeof(TYPE)*3*number);

		//Copy data into CUDA memory
		cudaMemcpy(dev_X,			X,			sizeof(TYPE)*3*number,		cudaMemcpyHostToDevice);
		cudaMemcpy(dev_V,			V,			sizeof(TYPE)*3*number,		cudaMemcpyHostToDevice);
		cudaMemcpy(dev_prev_X,		X,			sizeof(TYPE)*3*number,		cudaMemcpyHostToDevice);
		cudaMemcpy(dev_next_X,		X,			sizeof(TYPE)*3*number,		cudaMemcpyHostToDevice);
		cudaMemcpy(dev_fixed,		fixed,		sizeof(TYPE)*number,		cudaMemcpyHostToDevice);
		cudaMemset(dev_more_fixed,  0,			sizeof(TYPE)*number);	

		cudaMemcpy(dev_Dm,			Dm,			sizeof(int)*tet_number*9,	cudaMemcpyHostToDevice);
		cudaMemcpy(dev_inv_Dm,		inv_Dm,		sizeof(int)*tet_number*9,	cudaMemcpyHostToDevice);
		cudaMemcpy(dev_Vol,			Vol,		sizeof(int)*tet_number,		cudaMemcpyHostToDevice);
		cudaMemcpy(dev_Tet,			Tet,		sizeof(int)*tet_number*4,	cudaMemcpyHostToDevice);

		cudaMemcpy(dev_MD,			MD,			sizeof(TYPE)*number,		cudaMemcpyHostToDevice);
		cudaMemcpy(dev_VTT,			VTT,		sizeof(int )*tet_number*4,	cudaMemcpyHostToDevice);
		cudaMemcpy(dev_TQ,			TQ,			sizeof(TYPE)*tet_number*4,	cudaMemcpyHostToDevice);
		cudaMemcpy(dev_vtt_num,		vtt_num,	sizeof(int )*(number+1),	cudaMemcpyHostToDevice);
	}


///////////////////////////////////////////////////////////////////////////////////////////
//  Control functions
///////////////////////////////////////////////////////////////////////////////////////////
 	void Reset_More_Fixed(int select_v)
	{
		int threadsPerBlock = 64;
		int blocksPerGrid = (number + threadsPerBlock - 1) / threadsPerBlock;
		Control_Kernel << <blocksPerGrid, threadsPerBlock>> >(dev_X, dev_more_fixed, control_mag, number, select_v);
	}

///////////////////////////////////////////////////////////////////////////////////////////
//  Update functions
///////////////////////////////////////////////////////////////////////////////////////////
	void Update(TYPE t, int iterations, TYPE dir[])
	{
		int threadsPerBlock = 64;
		int blocksPerGrid = (number + threadsPerBlock - 1) / threadsPerBlock;
		int tet_threadsPerBlock = 64;
		int tet_blocksPerGrid = (tet_number + tet_threadsPerBlock - 1) / tet_threadsPerBlock;

		TIMER timer;

		// Step 1: Basic update
		Update_Kernel << <blocksPerGrid, threadsPerBlock>> >(dev_X, dev_V, dev_fixed, dev_more_fixed, damping, t, number, dir[0], dir[1], dir[2]);

		// Step 2: Set up X data
		Constraint_0_Kernel << <blocksPerGrid, threadsPerBlock>> >(dev_X, dev_init_B, dev_VC, dev_fixed, dev_more_fixed, 1/t, number);
		
		// Step 3: Running iterations
		TYPE omega;
		for(int l=0; l<iterations; l++)
		{	
			Tet_Constraint_Kernel << <tet_blocksPerGrid, tet_threadsPerBlock>> >(dev_X, dev_Tet, dev_inv_Dm, dev_Vol, dev_Tet_Temp, elasticity, tet_number, l);
			Constraint_1_Kernel << <blocksPerGrid, threadsPerBlock>> >(dev_X, dev_init_B, dev_VC, dev_next_X, dev_Tet_Temp, dev_MD, dev_VTT, dev_vtt_num, number);

			if(l<=10)		omega=1;
			else if(l==11)	omega=2/(2-rho*rho);
			else			omega=4/(4-rho*rho*omega);
			//omega=1;
			
			Constraint_2_Kernel<< <blocksPerGrid, threadsPerBlock>> >(dev_prev_X, dev_X, dev_next_X, omega, number);
			Swap(dev_X, dev_prev_X);
			Swap(dev_X, dev_next_X);
		}
		
		// Step 4: Finalizing update
		Constraint_3_Kernel<< <blocksPerGrid, threadsPerBlock>> >(dev_X, dev_init_B, dev_V, dev_fixed, dev_more_fixed, 1/t, number);
		
		//Output to main memory for rendering
		cudaMemcpy(X, dev_X, sizeof(TYPE)*3*number, cudaMemcpyDeviceToHost);

		cost[cost_ptr]=timer.Get_Time();
		cost_ptr=(cost_ptr+1)%8;
		fps=8/(cost[0]+cost[1]+cost[2]+cost[3]+cost[4]+cost[5]+cost[6]+cost[7]);
	}

///////////////////////////////////////////////////////////////////////////////////////////
//  IO functions
///////////////////////////////////////////////////////////////////////////////////////////
	void Write(std::fstream &output)
	{
		Write_Binaries(output, X, number*3);
		Write_Binaries(output, V, number*3);
	}

	bool Write_File(const char *file_name)
	{
		std::fstream output; 
		output.open(file_name,std::ios::out|std::ios::binary);
		if(!output.is_open())	{printf("Error, file not open.\n"); return false;}
		Write(output);
		output.close();
		return true;
	}
	
	void Read(std::fstream &input)
	{
		Read_Binaries(input, X, number*3);
		Read_Binaries(input, V, number*3);
	}

	bool Read_File(const char *file_name)
	{
		std::fstream input; 
		input.open(file_name,std::ios::in|std::ios::binary);
		if(!input.is_open())	{printf("Error, file not open.\n");	return false;}
		Read(input);
		input.close();
		return true;
	}
};


#endif