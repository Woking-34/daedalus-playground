///////////////////////////////////////////////////////////////////////////////////////////
//  Copyright (C) 2002 - 2016, Huamin Wang
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
//  Class CUDA_HYPER_TET_MESH
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef	__WHMIN_CUDA_HYPER_TET_MESH_H__
#define __WHMIN_CUDA_HYPER_TET_MESH_H__

#include "PROGRESSING_BAR.h"
#include "TET_MESH.h"
#include "CUDA_MATH.h"

#define	StVK_MODEL			0
#define NH_MODEL			1
#define MR_MODEL			2
#define	FUNG_MODEL			3
#define RADIUS_SQUARED		0.01


///////////////////////////////////////////////////////////////////////////////////////////
//  Control kernel
///////////////////////////////////////////////////////////////////////////////////////////
__global__ void Control_Kernel(float* X, float* fixed, float *more_fixed, float* offset_X, const float control_mag, const int number, const int select_v)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i>=number)	return;

	more_fixed[i]=0;
	if(fixed[i]==0 && select_v!=-1)
	{
		offset_X[i*3+0]=X[i*3+0]-X[select_v*3+0];
		offset_X[i*3+1]=X[i*3+1]-X[select_v*3+1];
		offset_X[i*3+2]=X[i*3+2]-X[select_v*3+2];

		float dist2=offset_X[i*3+0]*offset_X[i*3+0]+offset_X[i*3+1]*offset_X[i*3+1]+offset_X[i*3+2]*offset_X[i*3+2];
		//dist2+=(X[i*3+0]-X[select_v*3+0])*(X[i*3+0]-X[select_v*3+0]);
		//dist2+=(X[i*3+1]-X[select_v*3+1])*(X[i*3+1]-X[select_v*3+1]);
		//dist2+=(X[i*3+2]-X[select_v*3+2])*(X[i*3+2]-X[select_v*3+2]);
		if(dist2<RADIUS_SQUARED)	more_fixed[i]=control_mag*(1-sqrt(dist2/RADIUS_SQUARED));
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
//  Basic update kernel
///////////////////////////////////////////////////////////////////////////////////////////
__global__ void Update_Kernel(float* X, float* V, const float* prev_V, float* S, const float *fixed, const float *more_fixed, const float *offset_X, float* fixed_X, const float t, const int number, const float dir_x, const float dir_y, const float dir_z)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i>=number)	return;

	if(more_fixed[i])
	{
		//fixed_X[i*3+0]=X[i*3+0]+dir_x;
		//fixed_X[i*3+1]=X[i*3+1]+dir_y;
		//fixed_X[i*3+2]=X[i*3+2]+dir_z;

		fixed_X[i*3+0]=offset_X[i*3+0]+dir_x;
		fixed_X[i*3+1]=offset_X[i*3+1]+dir_y;
		fixed_X[i*3+2]=offset_X[i*3+2]+dir_z;
	}

	//V[i*3+1]-=2*t;

	//Calculate S
    S[i*3+0]=X[i*3+0]+V[i*3+0]*t;
	S[i*3+1]=X[i*3+1]+V[i*3+1]*t;
    S[i*3+2]=X[i*3+2]+V[i*3+2]*t;
	//Initialize position
	X[i * 3 + 0] += (V[i * 3 + 0] + (V[i * 3 + 0] - prev_V[i * 3 + 0])*0.2)*t;
	X[i * 3 + 1] += (V[i * 3 + 1] + (V[i * 3 + 1] - prev_V[i * 3 + 1])*0.2)*t;
	X[i * 3 + 2] += (V[i * 3 + 2] + (V[i * 3 + 2] - prev_V[i * 3 + 2])*0.2)*t;
}

///////////////////////////////////////////////////////////////////////////////////////////
//  Tet Constraint Kernel
///////////////////////////////////////////////////////////////////////////////////////////

__global__ void Compute_FM_Kernel(const float* X, const int* Tet, const float* inv_Dm, const float* Vol, float* lambda, float* Force, float* C, float* ext_C, float *E,
	const int model, float stiffness_0, float stiffness_1, float stiffness_2, float stiffness_3, float stiffness_p, const int tet_number, const float lower_bound, const float upper_bound, const bool update_C=true)
{
	int t = blockDim.x * blockIdx.x + threadIdx.x;
	if(t>=tet_number)	return;

	stiffness_0=-Vol[t]*stiffness_0;
    stiffness_1=-Vol[t]*stiffness_1;
    stiffness_2=-Vol[t]*stiffness_2;
    //No velocity access in this function
    int p0=Tet[t*4+0]*3;
    int p1=Tet[t*4+1]*3;
    int p2=Tet[t*4+2]*3;
    int p3=Tet[t*4+3]*3;


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
        
    float F[9], U[9], Sigma[3], V[9];
	dev_Matrix_Product_3(Ds, &inv_Dm[t*9], F);   
	
	float W[9];
	svd(F[0], F[1], F[2], F[3], F[4], F[5], F[6], F[7], F[8],
		U[0], U[1], U[2], U[3], U[4], U[5], U[6], U[7], U[8], 
		Sigma[0], W[1], W[2], W[3], Sigma[1], W[5], W[6], W[7], Sigma[2],
		V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], V[8]);
	

	

	float interpolator0, interpolator1, interpolator2;
	if (Sigma[0] < 1)	interpolator0 = (lower_bound - Sigma[0])*10.0;	//10
	//else				interpolator0 = (Sigma[0] - upper_bound)*0.5;
	if (Sigma[1] < 1)	interpolator1 = (lower_bound - Sigma[1])*10.0;	//10
	//else				interpolator1 = (Sigma[1] - upper_bound)*0.5;
	if (Sigma[2] < 1)	interpolator2 = (lower_bound - Sigma[2])*10.0;	//10
	//else				interpolator2 = (Sigma[2] - upper_bound)*0.5;

	float interpolator = MAX(MAX(interpolator0, interpolator1), interpolator2);
	if(interpolator<0)		interpolator=0;
	if(interpolator>1)		interpolator=1;

	lambda[t] = MAX(lambda[t], interpolator);
	//lambda[t] = 1;
	interpolator = lambda[t];
	stiffness_0 *= (1 - interpolator);
	stiffness_1 *= (1 - interpolator);
	stiffness_2 *= (1 - interpolator);

   
	float I			= Sigma[0]*Sigma[0]+Sigma[1]*Sigma[1]+Sigma[2]*Sigma[2];
    float J			= Sigma[0]*Sigma[1]*Sigma[2];
    float II		= Sigma[0]*Sigma[0]*Sigma[0]*Sigma[0]+Sigma[1]*Sigma[1]*Sigma[1]*Sigma[1]+Sigma[2]*Sigma[2]*Sigma[2]*Sigma[2];
	float III		= J*J;
	float rcbrt_III	= 1.0/cbrt(III);
	float factor_1	= ONE_THIRD*rcbrt_III/III;
	float factor_2	= ONE_THIRD* factor_1/III;
	float dEdI		= 0;
	float dEdII		= 0;
	float dEdIII	= 0;
	float H[3][3];
	memset(&H[0][0], 0, sizeof(float) * 9);


	float energy = 0;

	if(model==StVK_MODEL)
	{
		dEdI	+= stiffness_0*(I-3)*0.25-stiffness_1*0.5;
		dEdII	+= stiffness_1*0.25;
		dEdIII	+= 0;
		H[0][0]	+= stiffness_0*0.25;
		energy-=stiffness_0*(I-3)*(I-3)*0.125+stiffness_1*(II-2*I+3)*0.25;
	}
	if(model==NH_MODEL || model==MR_MODEL || model==FUNG_MODEL)
	{
		dEdI	+= rcbrt_III*stiffness_0;
		dEdII	+= 0;
		dEdIII	+= -factor_1*stiffness_0*I;
		H[0][2]	+= -factor_1*stiffness_0;
		H[2][2]	+=  factor_2*stiffness_0*I*4;
		// Volume correction
		dEdIII	+=stiffness_1*(J-1)/J;
		H[2][2]	+=stiffness_1/(2*III*J);

		//Just do NH energy for now...
		energy-=(stiffness_0*(I *rcbrt_III-3) + stiffness_1*(J-1)*(J-1));
	}
	if(model==MR_MODEL)
	{
		float two_term_a	= stiffness_2*rcbrt_III*I;
		float two_term_b	= stiffness_2*rcbrt_III*(I*I-II);
		dEdI	+= rcbrt_III*two_term_a;
		dEdII	+= -0.5*stiffness_2*rcbrt_III*rcbrt_III;
		dEdIII	+= -factor_1*two_term_b;
		H[0][0]	+= stiffness_2*rcbrt_III*rcbrt_III;
		H[0][2]	+= -factor_1*two_term_a*2;
		H[1][2]	+= factor_1*stiffness_2*rcbrt_III;
		H[2][2]	+= factor_2*two_term_b*5;

		energy-=stiffness_2*(0.5*rcbrt_III*rcbrt_III*(I*I-II)-3);
	}		
	if(model==FUNG_MODEL)
	{
		float value=stiffness_3*(rcbrt_III*I-3);

		if(value<4.0)
		{
			float exp_term	= expf(stiffness_3*(rcbrt_III*I-3))*stiffness_2*stiffness_3;
			dEdI	+= exp_term*rcbrt_III;
			dEdIII	+= -factor_1*I*exp_term;
			H[0][0]	+= rcbrt_III*stiffness_3*rcbrt_III*exp_term;
			H[0][2]	+= -factor_1*exp_term*(1+stiffness_3*rcbrt_III*I);
			H[2][2]	+= factor_2*I*exp_term* (4 + stiffness_3*I*rcbrt_III);
		}
		energy-=stiffness_2*(exp(MIN(value, 4.0))-1);
	}
	

	// Make H symmetric
	H[1][0]	= H[0][1];
	H[2][0]	= H[0][2];
	H[2][1]	= H[1][2];

	float P0 = 2 * dEdI*Sigma[0] + 4 * dEdII*Sigma[0] * Sigma[0] * Sigma[0] + 2 * dEdIII*III / Sigma[0];
	float P1 = 2 * dEdI*Sigma[1] + 4 * dEdII*Sigma[1] * Sigma[1] * Sigma[1] + 2 * dEdIII*III / Sigma[1];
	float P2 = 2 * dEdI*Sigma[2] + 4 * dEdII*Sigma[2] * Sigma[2] * Sigma[2] + 2 * dEdIII*III / Sigma[2];
	
    float PV_transpose[9], P[9], G[9];
    PV_transpose[0]=P0*V[0];
    PV_transpose[1]=P0*V[3];
    PV_transpose[2]=P0*V[6];
    PV_transpose[3]=P1*V[1];
    PV_transpose[4]=P1*V[4];
    PV_transpose[5]=P1*V[7];
    PV_transpose[6]=P2*V[2];
    PV_transpose[7]=P2*V[5];
    PV_transpose[8]=P2*V[8];
    dev_Matrix_Product_3(U, PV_transpose, P);
    dev_Matrix_Product_T_3(P, &inv_Dm[t*9], G);

	float force[12];
    force[ 0]=-(G[0]+G[1]+G[2]);
    force[ 1]=-(G[3]+G[4]+G[5]);
    force[ 2]=-(G[6]+G[7]+G[8]);
    force[ 3]=G[0];
    force[ 4]=G[3];
    force[ 5]=G[6];
    force[ 6]=G[1];
    force[ 7]=G[4];
    force[ 8]=G[7];
    force[ 9]=G[2];
    force[10]=G[5];
    force[11]=G[8];

	//Strain limiting part
	float new_R[9];
	dev_Matrix_Product_T_3(U, V, new_R);
	const float* idm = &inv_Dm[t * 9];
	float half_matrix[3][4], result_matrix[3][4];
	half_matrix[0][0] = -idm[0] - idm[3] - idm[6];
	half_matrix[0][1] =  idm[0];
	half_matrix[0][2] =  idm[3];
	half_matrix[0][3] =  idm[6];
	half_matrix[1][0] = -idm[1] - idm[4] - idm[7];
	half_matrix[1][1] =  idm[1];
	half_matrix[1][2] =  idm[4];
	half_matrix[1][3] =  idm[7];
	half_matrix[2][0] = -idm[2] - idm[5] - idm[8];
	half_matrix[2][1] =  idm[2];
	half_matrix[2][2] =  idm[5];
	half_matrix[2][3] =  idm[8];

	dev_Matrix_Substract_3(new_R, F, new_R);
	dev_Matrix_Product(new_R, &half_matrix[0][0], &result_matrix[0][0], 3, 3, 4);

	float pd_stiffness = Vol[t] * stiffness_p * interpolator;
	force[ 0] += result_matrix[0][0] * pd_stiffness;
	force[ 1] += result_matrix[1][0] * pd_stiffness;
	force[ 2] += result_matrix[2][0] * pd_stiffness;
	force[ 3] += result_matrix[0][1] * pd_stiffness;
	force[ 4] += result_matrix[1][1] * pd_stiffness;
	force[ 5] += result_matrix[2][1] * pd_stiffness;
	force[ 6] += result_matrix[0][2] * pd_stiffness;
	force[ 7] += result_matrix[1][2] * pd_stiffness;
	force[ 8] += result_matrix[2][2] * pd_stiffness;
	force[ 9] += result_matrix[0][3] * pd_stiffness;
	force[10] += result_matrix[1][3] * pd_stiffness;
	force[11] += result_matrix[2][3] * pd_stiffness;
	
	//PD energy
	float pd_energy=0;
	pd_energy+=new_R[0]*new_R[0];
	pd_energy+=new_R[1]*new_R[1];
	pd_energy+=new_R[2]*new_R[2];
	pd_energy+=new_R[3]*new_R[3];
	pd_energy+=new_R[4]*new_R[4];
	pd_energy+=new_R[5]*new_R[5];
	pd_energy+=new_R[6]*new_R[6];
	pd_energy+=new_R[7]*new_R[7];
	pd_energy+=new_R[8]*new_R[8];
	energy+=pd_energy*pd_stiffness*0.5;


	atomicAdd(&Force[p0 + 0], force[ 0]);
	atomicAdd(&Force[p0 + 1], force[ 1]);
	atomicAdd(&Force[p0 + 2], force[ 2]);
	atomicAdd(&Force[p1 + 0], force[ 3]);
	atomicAdd(&Force[p1 + 1], force[ 4]);
	atomicAdd(&Force[p1 + 2], force[ 5]);
	atomicAdd(&Force[p2 + 0], force[ 6]);
	atomicAdd(&Force[p2 + 1], force[ 7]);
	atomicAdd(&Force[p2 + 2], force[ 8]);
	atomicAdd(&Force[p3 + 0], force[ 9]);
	atomicAdd(&Force[p3 + 1], force[10]);
	atomicAdd(&Force[p3 + 2], force[11]);
	
	float value0 = half_matrix[0][0]*half_matrix[0][0]+half_matrix[1][0]*half_matrix[1][0]+half_matrix[2][0]*half_matrix[2][0];
	float value1 = half_matrix[0][1]*half_matrix[0][1]+half_matrix[1][1]*half_matrix[1][1]+half_matrix[2][1]*half_matrix[2][1];
	float value2 = half_matrix[0][2]*half_matrix[0][2]+half_matrix[1][2]*half_matrix[1][2]+half_matrix[2][2]*half_matrix[2][2];
	float value3 = half_matrix[0][3]*half_matrix[0][3]+half_matrix[1][3]*half_matrix[1][3]+half_matrix[2][3]*half_matrix[2][3];
	atomicAdd(&ext_C[Tet[t * 4 + 0]], value0*pd_stiffness);
	atomicAdd(&ext_C[Tet[t * 4 + 1]], value1*pd_stiffness);
	atomicAdd(&ext_C[Tet[t * 4 + 2]], value2*pd_stiffness);
	atomicAdd(&ext_C[Tet[t * 4 + 3]], value3*pd_stiffness);

	//Distribute the energy to the vertices
	atomicAdd(&E[Tet[t * 4 + 0]], energy);


	
	if(update_C==false)	return;
    //Now compute the stiffness matrix
	float alpha[3][3], beta[3][3], gamma[3][3];
	float vi[3], vj[3], r[3];
    for(int i=0; i<3; i++)
    for(int j=i; j<3; j++)
    {
        alpha[i][j]=2*dEdI+4*(Sigma[i]*Sigma[i]+Sigma[j]*Sigma[j])*dEdII;
        beta[i][j]=4*Sigma[i]*Sigma[j]*dEdII-2*III*dEdIII/(Sigma[i]*Sigma[j]);

		vi[0]=2*Sigma[i];
		vi[1]=4*Sigma[i]*Sigma[i]*Sigma[i];
		vi[2]=2*III/Sigma[i];

		vj[0]=2*Sigma[j];
		vj[1]=4*Sigma[j]*Sigma[j]*Sigma[j];
		vj[2]=2*III/Sigma[j];

        dev_Matrix_Vector_Product_3(&H[0][0], vj, r);
        gamma[i][j]=DOT(vi, r)+4*III*dEdIII/(Sigma[i]*Sigma[j]);
    }
        
    float dGdF[12][9];    
	//G is related to force, according to [TSIF05], (g0, g3, g6), (g1, g4, g7), (g2, g5, g8), (g9, g10, g11)
	float dF[9], temp0[9], temp1[9];
    for(int i=0; i<9; i++)
    {
        memset(&dF, 0, sizeof(float)*9);
        dF[i]=1;
            
        dev_Matrix_Product_3(dF, V, temp0);
        dev_Matrix_T_Product_3(U, temp0, temp1);
            
        //diagonal
        temp0[0]=(alpha[0][0]+beta[0][0]+gamma[0][0])*temp1[0]+gamma[0][1]*temp1[4]+gamma[0][2]*temp1[8];
        temp0[4]=gamma[0][1]*temp1[0]+(alpha[1][1]+beta[1][1]+gamma[1][1])*temp1[4]+gamma[1][2]*temp1[8];
        temp0[8]=gamma[0][2]*temp1[0]+gamma[1][2]*temp1[4]+(alpha[2][2]+beta[2][2]+gamma[2][2])*temp1[8];
        //off-diagonal
        temp0[1]=alpha[0][1]*temp1[1]+beta[0][1]*temp1[3];
        temp0[3]=alpha[0][1]*temp1[3]+beta[0][1]*temp1[1];
        temp0[2]=alpha[0][2]*temp1[2]+beta[0][2]*temp1[6];
        temp0[6]=alpha[0][2]*temp1[6]+beta[0][2]*temp1[2];
        temp0[5]=alpha[1][2]*temp1[5]+beta[1][2]*temp1[7];
        temp0[7]=alpha[1][2]*temp1[7]+beta[1][2]*temp1[5];
            
        dev_Matrix_Product_T_3(temp0, V, temp1);
        dev_Matrix_Product_3(U, temp1, temp0);
        dev_Matrix_Product_T_3(temp0, &inv_Dm[t*9], &dGdF[i][0]);
    }
        
    //Transpose dGdF
    float temp;
    for(int i=0; i<9; i++) for(int j=i+1; j<9; j++)
        SWAP(dGdF[i][j], dGdF[j][i]);
        
    for(int j=0; j< 9; j++)
    {
        dGdF[ 9][j]=-dGdF[0][j]-dGdF[1][j]-dGdF[2][j];
        dGdF[10][j]=-dGdF[3][j]-dGdF[4][j]-dGdF[5][j];
        dGdF[11][j]=-dGdF[6][j]-dGdF[7][j]-dGdF[8][j];
    }
        
    float new_idm[4][3];
    new_idm[0][0]=-inv_Dm[t*9+0]-inv_Dm[t*9+3]-inv_Dm[t*9+6];
    new_idm[0][1]=-inv_Dm[t*9+1]-inv_Dm[t*9+4]-inv_Dm[t*9+7];
	new_idm[0][2]=-inv_Dm[t*9+2]-inv_Dm[t*9+5]-inv_Dm[t*9+8];
	new_idm[1][0]= inv_Dm[t*9+0];
	new_idm[1][1]= inv_Dm[t*9+1];
	new_idm[1][2]= inv_Dm[t*9+2];
	new_idm[2][0]= inv_Dm[t*9+3];
	new_idm[2][1]= inv_Dm[t*9+4];
	new_idm[2][2]= inv_Dm[t*9+5];
	new_idm[3][0]= inv_Dm[t*9+6];
	new_idm[3][1]= inv_Dm[t*9+7];
	new_idm[3][2]= inv_Dm[t*9+8];

	atomicAdd(&C[p0*3+0], -dev_Matrix_Product_T(&new_idm[0][0], &dGdF[ 9][0], 4, 3, 3, 0, 0));
	atomicAdd(&C[p0*3+4], -dev_Matrix_Product_T(&new_idm[0][0], &dGdF[10][0], 4, 3, 3, 0, 1));
	atomicAdd(&C[p0*3+8], -dev_Matrix_Product_T(&new_idm[0][0], &dGdF[11][0], 4, 3, 3, 0, 2));
	atomicAdd(&C[p1*3+0], -dev_Matrix_Product_T(&new_idm[0][0], &dGdF[ 0][0], 4, 3, 3, 1, 0));
	atomicAdd(&C[p1*3+4], -dev_Matrix_Product_T(&new_idm[0][0], &dGdF[ 3][0], 4, 3, 3, 1, 1));
	atomicAdd(&C[p1*3+8], -dev_Matrix_Product_T(&new_idm[0][0], &dGdF[ 6][0], 4, 3, 3, 1, 2));
	atomicAdd(&C[p2*3+0], -dev_Matrix_Product_T(&new_idm[0][0], &dGdF[ 1][0], 4, 3, 3, 2, 0));
	atomicAdd(&C[p2*3+4], -dev_Matrix_Product_T(&new_idm[0][0], &dGdF[ 4][0], 4, 3, 3, 2, 1));
	atomicAdd(&C[p2*3+8], -dev_Matrix_Product_T(&new_idm[0][0], &dGdF[ 7][0], 4, 3, 3, 2, 2));
	atomicAdd(&C[p3*3+0], -dev_Matrix_Product_T(&new_idm[0][0], &dGdF[ 2][0], 4, 3, 3, 3, 0));
	atomicAdd(&C[p3*3+4], -dev_Matrix_Product_T(&new_idm[0][0], &dGdF[ 5][0], 4, 3, 3, 3, 1));
	atomicAdd(&C[p3*3+8], -dev_Matrix_Product_T(&new_idm[0][0], &dGdF[ 8][0], 4, 3, 3, 3, 2));
}


///////////////////////////////////////////////////////////////////////////////////////////
//  Constraint Kernel 1
///////////////////////////////////////////////////////////////////////////////////////////
__global__ void Constraint_1_Kernel(const float* M, const float* X, const float* prev_X, const float* V, float* E, float* G, float* P, const float* S, float* next_X, 
	const float* fixed, const float* more_fixed, const float* fixed_X, float *F, const float* C, const float* ext_C, const float stepping, const int number, const float t, const float inv_t, const float gravity)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i>=number)	return;

	float oc = M[i]*inv_t*inv_t;
	float c  = M[i]*inv_t*inv_t+fixed[i]+more_fixed[i];

	// Get Force
	float error[3];
	error[0]=oc*(S[i*3+0]-X[i*3+0])+(c-oc)*(fixed_X[i*3+0]-X[i*3+0])+F[i*3+0];
	error[1]=oc*(S[i*3+1]-X[i*3+1])+(c-oc)*(fixed_X[i*3+1]-X[i*3+1])+F[i*3+1]+gravity*M[i];
	error[2]=oc*(S[i*3+2]-X[i*3+2])+(c-oc)*(fixed_X[i*3+2]-X[i*3+2])+F[i*3+2];
	
	// Update Energy
	float energy=0;
	energy+=oc*(S[i*3+0]-X[i*3+0])*(S[i*3+0]-X[i*3+0]);					// kinetic energy
	energy+=oc*(S[i*3+1]-X[i*3+1])*(S[i*3+1]-X[i*3+1]);					// kinetic energy
	energy+=oc*(S[i*3+2]-X[i*3+2])*(S[i*3+2]-X[i*3+2]);					// kinetic energy
	energy+=(c-oc)*(fixed_X[i*3+0]-X[i*3+0])*(fixed_X[i*3+0]-X[i*3+0]);	// fixed energy
	energy+=(c-oc)*(fixed_X[i*3+1]-X[i*3+1])*(fixed_X[i*3+1]-X[i*3+1]);	// fixed energy
	energy+=(c-oc)*(fixed_X[i*3+2]-X[i*3+2])*(fixed_X[i*3+2]-X[i*3+2]);	// fixed energy
	energy*=0.5f;
	energy+=-gravity*M[i]*X[i*3+1];										// gravity energy
	E[i]+=energy;

	float cx=C[i*9+0]+c+ext_C[i];
	float cy=C[i*9+4]+c+ext_C[i];
	float cz=C[i*9+8]+c+ext_C[i];
	
	// Update Gradient
	G[i]=error[0]*error[0]+error[1]*error[1]+error[2]*error[2];
	// Update Force
	F[i*3+0]=error[0];
	F[i*3+1]=error[1];
	F[i*3+2]=error[2];
	// Update position
	next_X[i*3+0]=X[i*3+0]+error[0]*stepping/cx;
	next_X[i*3+1]=X[i*3+1]+error[1]*stepping/cy;
	next_X[i*3+2]=X[i*3+2]+error[2]*stepping/cz;
}

///////////////////////////////////////////////////////////////////////////////////////////
//  Constraint Kernel 2
///////////////////////////////////////////////////////////////////////////////////////////
__global__ void Constraint_2_Kernel(float* prev_X, float* next_X, float omega, int number)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i>=number)	return;	

	next_X[i*3+0]=omega*(next_X[i*3+0]-prev_X[i*3+0])+prev_X[i*3+0];
	next_X[i*3+1]=omega*(next_X[i*3+1]-prev_X[i*3+1])+prev_X[i*3+1];
	next_X[i*3+2]=omega*(next_X[i*3+2]-prev_X[i*3+2])+prev_X[i*3+2];
}

///////////////////////////////////////////////////////////////////////////////////////////
//  Constraint Kernel 3
///////////////////////////////////////////////////////////////////////////////////////////
__global__ void Constraint_3_Kernel(float* X, float* S, float* V, float* prev_V, const float *fixed, const float *more_fixed, float t, float inv_t, int number)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i>=number)	return;
	
	prev_V[i * 3 + 0] = V[i * 3 + 0];
	prev_V[i * 3 + 1] = V[i * 3 + 1];
	prev_V[i * 3 + 2] = V[i * 3 + 2];

	V[i*3+0]+=(X[i*3+0]-S[i*3+0])*inv_t;
	V[i*3+1]+=(X[i*3+1]-S[i*3+1])*inv_t;
	V[i*3+2]+=(X[i*3+2]-S[i*3+2])*inv_t;
}


template <class TYPE>
class CUDA_HYPER_TET_MESH: public TET_MESH<TYPE> 
{
public:
	TYPE	cost[64];
	int		cost_ptr;
	TYPE	fps;
	TYPE	stepping;

	TYPE*	V;
	TYPE*	fixed;
	TYPE*	more_fixed;
	TYPE*	fixed_X;

	TYPE	rho;
	TYPE	model;
	TYPE	stiffness_0;
	TYPE	stiffness_1;
	TYPE	stiffness_2;
	TYPE	stiffness_3;
	TYPE	stiffness_p;
	TYPE	lower_bound;
	TYPE	upper_bound;
	TYPE	control_mag;
	TYPE	gravity;

	//CUDA data
	TYPE*	dev_M;
	TYPE*	dev_X;
	TYPE*	dev_V;
	TYPE*	dev_prev_V;
	TYPE*	dev_F;
	TYPE*	dev_start_F;	// only contain part of viscosity
	TYPE*	dev_S;			// Initialized S
	TYPE*	dev_C;			// diagonal part
	TYPE*	dev_ext_C;		// extra diagonal part
	TYPE*	dev_next_X;		// next X		(for temporary storage)
	TYPE*	dev_prev_X;		// previous X	(for acceleration)
	TYPE*	dev_last_X;		// backup X
	
	TYPE*	dev_fixed;
	TYPE*	dev_more_fixed;
	TYPE*	dev_fixed_X;
	TYPE*	dev_offset_X;

	TYPE*	dev_inv_Dm;
	TYPE*	dev_Vol;
	int*	dev_Tet;
		
	TYPE*	dev_lambda;
	TYPE*	dev_last_lambda;

	TYPE*	dev_G;
	TYPE*	dev_E;			// energy
	TYPE*	dev_P;			// the product of descent directions

	int		profile_l[64];
	TYPE	profile_v[64];
	int		profile_length;

	int		threadsPerBlock;
	int		blocksPerGrid;
	int		tet_threadsPerBlock;
	int		tet_blocksPerGrid;



	CUDA_HYPER_TET_MESH()
	{
		cost_ptr= 0;

		V			= new TYPE	[max_number*3 ];
		fixed		= new TYPE	[max_number   ];
		more_fixed	= new TYPE	[max_number   ];	
		fixed_X		= new TYPE	[max_number*3 ];

		// Default parameters
		fps			= 0;
		control_mag	= 10;
		rho			= 0.9992;
		stepping	= 0.3;

		model		= NH_MODEL;
		stiffness_0	= 2000000;	//2000000
        stiffness_1	= 2000000;	//2000000
        stiffness_2	= 2000000;	//2000000
		stiffness_3 = 1.0;
		stiffness_p	= 2000000;
		lower_bound	= 0.4;
		upper_bound	= 3.0;
		gravity		= -9.8;

		memset(		V, 0, sizeof(TYPE)*max_number*3);
		memset(	fixed, 0, sizeof(int )*max_number  );

		// GPU data
		dev_M			= 0;
		dev_X			= 0;
		dev_V			= 0;
		dev_prev_V		= 0;

		dev_F			= 0;
		dev_start_F		= 0;
		dev_S			= 0;
		dev_C			= 0;
		dev_ext_C		= 0;
		dev_next_X		= 0;
		dev_prev_X		= 0;
		dev_last_X		= 0;

		dev_fixed		= 0;
		dev_more_fixed	= 0;
		dev_fixed_X		= 0;
		dev_offset_X	= 0;

		dev_inv_Dm		= 0;
		dev_Vol			= 0;
		dev_Tet			= 0;
		
		dev_lambda		= 0;
		dev_last_lambda	= 0;

		dev_E			= 0;
		dev_G			= 0;
		dev_P			= 0;

		profile_l[0] = 1;
		profile_l[1] = 7;
		profile_l[2] = 10;
		profile_v[0] = 0.96;
		profile_v[1] = 0.991;
		profile_v[2] = 0.999;
		profile_length = 3;
	}
	
	~CUDA_HYPER_TET_MESH()
	{
		if(V)				delete[] V;
		if(fixed)			delete[] fixed;
		if(more_fixed)		delete[] more_fixed;
		
		//GPU Data
		if(dev_M)			cudaFree(dev_M);
		if(dev_X)			cudaFree(dev_X);
		if(dev_V)			cudaFree(dev_V);
		if(dev_prev_V)		cudaFree(dev_prev_V);
		if(dev_F)			cudaFree(dev_F);
		if(dev_start_F)		cudaFree(dev_start_F);
		if(dev_S)			cudaFree(dev_S);
		if(dev_C)			cudaFree(dev_C);
		if(dev_ext_C)		cudaFree(dev_ext_C);
		if(dev_next_X)		cudaFree(dev_next_X);
		if(dev_prev_X)		cudaFree(dev_prev_X);
		if(dev_last_X)		cudaFree(dev_last_X);
		
		if(dev_fixed)		cudaFree(dev_fixed);
		if(dev_more_fixed)	cudaFree(dev_more_fixed);
		if(dev_fixed_X)		cudaFree(dev_fixed_X);
		if(dev_offset_X)	cudaFree(dev_offset_X);

		if(dev_inv_Dm)		cudaFree(dev_inv_Dm);
		if(dev_Vol)			cudaFree(dev_Vol);
		if(dev_Tet)			cudaFree(dev_Tet);

		if(dev_lambda)		cudaFree(dev_lambda);
		if(dev_last_lambda)	cudaFree(dev_last_lambda);

		if(dev_E)			cudaFree(dev_E);
		if(dev_G)			cudaFree(dev_G);
		if(dev_P)			cudaFree(dev_P);
	}


///////////////////////////////////////////////////////////////////////////////////////////
//  Initialize functions
///////////////////////////////////////////////////////////////////////////////////////////
	
	void Initialize(TYPE t)
	{
		TET_MESH<TYPE>::Initialize();

		threadsPerBlock = 64;
		blocksPerGrid = (number + threadsPerBlock - 1) / threadsPerBlock;
		tet_threadsPerBlock = 64;
		tet_blocksPerGrid = (tet_number + tet_threadsPerBlock - 1) / tet_threadsPerBlock;

		Allocate_GPU_Memory();
	}

	void Allocate_GPU_Memory()
	{
		//Allocate CUDA memory
		cudaMalloc((void**)&dev_M,			sizeof(int )*number  );
		cudaMalloc((void**)&dev_X,			sizeof(int )*number*3);
		cudaMalloc((void**)&dev_V,			sizeof(TYPE)*number*3);
		cudaMalloc((void**)&dev_prev_V,		sizeof(TYPE)*number*3);
		cudaMalloc((void**)&dev_F,			sizeof(TYPE)*number*3);
		cudaMalloc((void**)&dev_start_F,	sizeof(TYPE)*number*3);
		cudaMalloc((void**)&dev_S,			sizeof(TYPE)*number*3);
		cudaMalloc((void**)&dev_C,			sizeof(TYPE)*number*9);
		cudaMalloc((void**)&dev_ext_C,		sizeof(TYPE)*number  );
		cudaMalloc((void**)&dev_next_X,		sizeof(TYPE)*number*3);
		cudaMalloc((void**)&dev_prev_X,		sizeof(TYPE)*number*3);
		cudaMalloc((void**)&dev_last_X,		sizeof(TYPE)*number*3);

		cudaMalloc((void**)&dev_fixed,		sizeof(TYPE)*number  );
		cudaMalloc((void**)&dev_more_fixed, sizeof(TYPE)*number  );
		cudaMalloc((void**)&dev_fixed_X,	sizeof(TYPE)*number*3);
		cudaMalloc((void**)&dev_offset_X,	sizeof(TYPE)*number*3);

		cudaMalloc((void**)&dev_inv_Dm,		sizeof(TYPE)*tet_number*9);
		cudaMalloc((void**)&dev_Vol,		sizeof(int )*tet_number  );
		cudaMalloc((void**)&dev_Tet,		sizeof(int )*tet_number*4);
		
		cudaMalloc((void**)&dev_lambda,		sizeof(TYPE)*tet_number  );
		cudaMalloc((void**)&dev_last_lambda,sizeof(TYPE)*tet_number	 );

		cudaMalloc((void**)&dev_G,			sizeof(TYPE)*number);
		cudaMalloc((void**)&dev_E,			sizeof(TYPE)*number);
		cudaMalloc((void**)&dev_P,			sizeof(TYPE)*number);		

		//Copy data into CUDA memory
		cudaMemcpy(dev_M,			M,			sizeof(TYPE)*number,		cudaMemcpyHostToDevice);
		cudaMemcpy(dev_X,			X,			sizeof(TYPE)*3*number,		cudaMemcpyHostToDevice);
		cudaMemcpy(dev_V,			V,			sizeof(TYPE)*3*number,		cudaMemcpyHostToDevice);
		cudaMemcpy(dev_prev_V,		V,			sizeof(TYPE)*3*number,		cudaMemcpyHostToDevice);		
		cudaMemcpy(dev_prev_X,		X,			sizeof(TYPE)*3*number,		cudaMemcpyHostToDevice);
		cudaMemcpy(dev_next_X,		X,			sizeof(TYPE)*3*number,		cudaMemcpyHostToDevice);
		cudaMemcpy(dev_fixed,		fixed,		sizeof(TYPE)*number,		cudaMemcpyHostToDevice);
		cudaMemcpy(dev_fixed_X,		fixed_X,	sizeof(TYPE)*3*number,		cudaMemcpyHostToDevice);
		cudaMemset(dev_more_fixed,  0,			sizeof(TYPE)*number);	

		cudaMemcpy(dev_fixed_X,		X,			sizeof(TYPE)*3*number,		cudaMemcpyHostToDevice);
		
		cudaMemcpy(dev_inv_Dm,		inv_Dm,		sizeof(int)*tet_number*9,	cudaMemcpyHostToDevice);
		cudaMemcpy(dev_Vol,			Vol,		sizeof(int)*tet_number,		cudaMemcpyHostToDevice);
		cudaMemcpy(dev_Tet,			Tet,		sizeof(int)*tet_number*4,	cudaMemcpyHostToDevice);

	}

///////////////////////////////////////////////////////////////////////////////////////////
//  Constraint functions
///////////////////////////////////////////////////////////////////////////////////////////

	void Reset_More_Fixed(int select_v)
	{
		int threadsPerBlock = 64;
		int blocksPerGrid = (number + threadsPerBlock - 1) / threadsPerBlock;
		Control_Kernel << <blocksPerGrid, threadsPerBlock>> >(dev_X, dev_fixed, dev_more_fixed, dev_offset_X, control_mag, number, select_v);
		cudaMemcpy(more_fixed, dev_more_fixed, sizeof(TYPE)*number, cudaMemcpyDeviceToHost);
		printf("select: %d\n", select_v);
	}

	void Clear_Velocity()
	{ 
		cudaMemset(dev_V, 0, sizeof(TYPE)*3*number);
	}

	void Evaluate(TYPE stepping, TYPE t, bool update_C)
	{
		cudaMemset(dev_F,     0, sizeof(TYPE)*number*3);
		cudaMemset(dev_ext_C, 0, sizeof(TYPE)*number  );
		cudaMemset(dev_G,     0, sizeof(TYPE)*number  );
		cudaMemset(dev_E,     0, sizeof(TYPE)*number  );
		cudaMemset(dev_P,     0, sizeof(TYPE)*number  );
		Compute_FM_Kernel << <tet_blocksPerGrid, tet_threadsPerBlock>> >(dev_X, dev_Tet, dev_inv_Dm, dev_Vol, dev_lambda, dev_F, dev_C, dev_ext_C, dev_E,
			model, stiffness_0, stiffness_1, stiffness_2, stiffness_3, stiffness_p, tet_number, lower_bound, upper_bound, update_C);		
		Constraint_1_Kernel << <blocksPerGrid, threadsPerBlock>> >(dev_M, dev_X, dev_prev_X, dev_V, dev_E, dev_G, dev_P, dev_S, dev_next_X, dev_fixed, dev_more_fixed, dev_fixed_X, dev_F, dev_C, dev_ext_C, stepping, number, t, 1/t, gravity);
	}

	float Get_Gradient_Magnitude()
	{
		thrust::device_ptr<TYPE> dev_g_ptr(dev_G);
		return thrust::reduce(dev_g_ptr, dev_g_ptr+number);
	}

	float Get_Energy_Magnitude()
	{
		thrust::device_ptr<TYPE> dev_e_ptr(dev_E);
		return thrust::reduce(dev_e_ptr, dev_e_ptr+number);
	}


	float Update(TYPE t, int iterations, TYPE dir[])
	{		
		TIMER timer;

		TYPE rho		= 0.9992;
		TYPE theta		= 1;
		TYPE omega;

		// Update S and initialize X
		cudaMemcpy(dev_S, dev_X, sizeof(TYPE)*3*number, cudaMemcpyDeviceToDevice);
		Update_Kernel << <blocksPerGrid, threadsPerBlock>> >(dev_X, dev_V, dev_prev_V, dev_S, dev_fixed, dev_more_fixed, dev_offset_X, dev_fixed_X, t, number, dir[0], dir[1], dir[2]);

		cudaMemset(dev_lambda, 0, sizeof(TYPE)*tet_number);
		
		//stepping=0.001;
		if(stepping<0.1)	stepping=0.1;
		stepping/=0.7;


		TYPE last_energy;		
		bool first=true;

		for(int l=0; l<iterations; l++)
		{
			if(stepping<0.01)	break;
			bool update_C=l%32==0;


			//always update
			if(update_C)	cudaMemset(dev_C,	  0, sizeof(TYPE)*number*9);
			Evaluate(stepping, t, update_C);


			if(l%8==0)
			{
				float energy=Get_Energy_Magnitude();
				//	thrust::device_ptr<TYPE> dev_c_ptr(dev_C);
				//	printf("C sum: %ef\n", thrust::reduce(dev_c_ptr+0, dev_c_ptr+number*3));
				//	printf("C sum: %ef\n", thrust::reduce(dev_c_ptr+number*3, dev_c_ptr+number*6));
				//	printf("C sum: %ef\n", thrust::reduce(dev_c_ptr+number*6, dev_c_ptr+number*9));
				//	printf("C sum: %ef\n", thrust::reduce(dev_c_ptr+0, dev_c_ptr+number*9));

				if(l==0 && first==false && last_energy>0 && energy>100*last_energy)
				{
					printf("what??? %f, %f\n", energy, last_energy);
					getchar();
				}

				if(l!=0 && (isnan(energy) || energy>last_energy))	//need to go back
				{
					cudaMemcpy(dev_lambda, dev_last_lambda, sizeof(TYPE)*tet_number, cudaMemcpyDeviceToDevice);
					cudaMemcpy(dev_X, dev_last_X, sizeof(TYPE)*number*3, cudaMemcpyDeviceToDevice);					

					iterations-=l-8;
					l=-1;
					stepping*=0.7;
					theta=1;

					first=false;
					continue;
				}
				else												//back it up
				{
					cudaMemcpy(dev_last_lambda, dev_lambda, sizeof(TYPE)*tet_number, cudaMemcpyDeviceToDevice);
					cudaMemcpy(dev_last_X, dev_X, sizeof(TYPE)*number*3, cudaMemcpyDeviceToDevice);
					last_energy=energy;
				}
			}
				

			if (l == 0)				{ rho = 0;	omega = 1; }
			omega = 4 / (4 - rho*rho*omega);
			for(int i=0; i<profile_length; i++)
			{
				if(l==profile_l[i]-1)
				{
					rho=0;
					omega=1;
				}
				if(l==profile_l[i])
				{
					rho=profile_v[i];
					omega = 2 / (2 - rho*rho);
					break;
				}
			}
	
			if(omega!=1)			
				Constraint_2_Kernel<< <blocksPerGrid, threadsPerBlock>> >(dev_prev_X, dev_next_X, omega, number);			

			Swap(dev_X, dev_prev_X);
			Swap(dev_X, dev_next_X);
		}

		// End the update
		Constraint_3_Kernel<< <blocksPerGrid, threadsPerBlock>> >(dev_X, dev_S, dev_V, dev_prev_V, dev_fixed, dev_more_fixed, t, 1/t, number);		
		cudaMemcpy(X, dev_X, sizeof(TYPE)*3*number, cudaMemcpyDeviceToHost);


		// Perform timing statistics
		cost[cost_ptr]=timer.Get_Time();
		cost_ptr=(cost_ptr+1)%64;
		fps=0;
		for(int i=0; i<64; i++)
			fps+=cost[i];
		fps=64.0/fps;
		
		return 0;
	}


	void Write(std::fstream &output)
	{
		cudaMemcpy(X, dev_X, sizeof(TYPE)*3*number, cudaMemcpyDeviceToHost);
		cudaMemcpy(V, dev_V, sizeof(TYPE)*3*number, cudaMemcpyDeviceToHost);

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