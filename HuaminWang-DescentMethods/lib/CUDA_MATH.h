///////////////////////////////////////////////////////////////////////////////////////////
//  Copyright (C) 2016 - 2016, Huamin Wang
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
//  CUDA_MATH library.
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef __WHMIN_CUDA_MATH_H__
#define __WHMIN_CUDA_MATH_H__
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include "svd3_cuda.h"


///////////////////////////////////////////////////////////////////////////////////////////
//  R = A - B
///////////////////////////////////////////////////////////////////////////////////////////
__host__ __device__ __forceinline__
 void dev_Matrix_Substract_3(float *A, float *B, float *R)								
{
	for (int i = 0; i < 9; i++)
		R[i]=A[i]-B[i];
}

///////////////////////////////////////////////////////////////////////////////////////////
//  R = A * B
///////////////////////////////////////////////////////////////////////////////////////////
__host__ __device__ __forceinline__
 void dev_Matrix_Product(float *A, float *B, float *R, int nx, int ny, int nz)			
{
	memset(R, 0, sizeof(float)*nx*nz);
	for (int i = 0; i < nx; i++)
	for (int j = 0; j < nz; j++)
	for (int k = 0; k < ny; k++)
		R[i*nz + j] += A[i*ny + k] * B[k*nz + j];
}

__host__ __device__ __forceinline__
 void dev_Matrix_Product_3(const float *A, const float *B, float *R)					
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

///////////////////////////////////////////////////////////////////////////////////////////
//  R = A * B'
///////////////////////////////////////////////////////////////////////////////////////////
__host__ __device__ __forceinline__
 void dev_Matrix_T_Product_3(const float *A, const float *B, float *R)
{
	R[0]=A[0]*B[0]+A[3]*B[3]+A[6]*B[6];
	R[1]=A[0]*B[1]+A[3]*B[4]+A[6]*B[7];
	R[2]=A[0]*B[2]+A[3]*B[5]+A[6]*B[8];
	R[3]=A[1]*B[0]+A[4]*B[3]+A[7]*B[6];
	R[4]=A[1]*B[1]+A[4]*B[4]+A[7]*B[7];
	R[5]=A[1]*B[2]+A[4]*B[5]+A[7]*B[8];
	R[6]=A[2]*B[0]+A[5]*B[3]+A[8]*B[6];
	R[7]=A[2]*B[1]+A[5]*B[4]+A[8]*B[7];
	R[8]=A[2]*B[2]+A[5]*B[5]+A[8]*B[8];
}

__host__ __device__ __forceinline__
 void dev_Matrix_Product_T_3(const float *A, const float *B, float *R)
{	
	R[0]=A[0]*B[0]+A[1]*B[1]+A[2]*B[2];
	R[1]=A[0]*B[3]+A[1]*B[4]+A[2]*B[5];
	R[2]=A[0]*B[6]+A[1]*B[7]+A[2]*B[8];
	R[3]=A[3]*B[0]+A[4]*B[1]+A[5]*B[2];
	R[4]=A[3]*B[3]+A[4]*B[4]+A[5]*B[5];
	R[5]=A[3]*B[6]+A[4]*B[7]+A[5]*B[8];
	R[6]=A[6]*B[0]+A[7]*B[1]+A[8]*B[2];
	R[7]=A[6]*B[3]+A[7]*B[4]+A[8]*B[5];
	R[8]=A[6]*B[6]+A[7]*B[7]+A[8]*B[8];
}

///////////////////////////////////////////////////////////////////////////////////////////
// r = A * x
///////////////////////////////////////////////////////////////////////////////////////////
__host__ __device__ __forceinline__
void dev_Matrix_Vector_Product_3(float *A, float *x, float *r)
{
	r[0]=A[0]*x[0]+A[1]*x[1]+A[2]*x[2];
	r[1]=A[3]*x[0]+A[4]*x[1]+A[5]*x[2];
	r[2]=A[6]*x[0]+A[7]*x[1]+A[8]*x[2];
}

///////////////////////////////////////////////////////////////////////////////////////////
//  return (A * B')[i, j]
///////////////////////////////////////////////////////////////////////////////////////////
__host__ __device__ __forceinline__
float dev_Matrix_Product_T(float *A, float *B, int nx, int ny, int nz, int i, int j)
{
    float ret=0;
    for(int k=0; k<ny; k++)
        ret+=A[i*ny+k]*B[j*ny+k];
    return ret;
}


#endif
