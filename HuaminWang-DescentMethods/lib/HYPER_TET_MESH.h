///////////////////////////////////////////////////////////////////////////////////////////
//  Copyright (C) 2002 - 2014, Huamin Wang
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
//  Class HYPER_TET_MESH
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef	__WHMIN_HYPER_TET_MESH_H__
#define __WHMIN_HYPER_TET_MESH_H__
#include "../../eigen/Eigen/Sparse"
#include "TET_MESH.h"
#include "LBFGS.h"
#include "PROGRESSING_BAR.h"
#include "PCG_SOLVER.h"
#include <Eigen/PardisoSupport>

template <class TYPE>
class HYPER_TET_MESH: public TET_MESH<TYPE>
{
public:
    using   TET_MESH<TYPE>::max_number;
    using   TET_MESH<TYPE>::number;
    using   TET_MESH<TYPE>::tet_number;
    using   TET_MESH<TYPE>::X;
    using   TET_MESH<TYPE>::Tet;
    using   TET_MESH<TYPE>::Vol;
    using   TET_MESH<TYPE>::inv_Dm;
    using   TET_MESH<TYPE>::M;
    
	TYPE*	lambda;
	TYPE*	ext_C;

    TYPE*   next_X;
    TYPE*   prev_X;
    TYPE*   Y;
	TYPE*	prev_Y;
	TYPE*	V;
	TYPE*	prev_V;
    TYPE*   S;
    TYPE*	B;
	TYPE*	D;
    TYPE*   Force;
	int*	fixed;
    TYPE*   fixed_X;
    
    TYPE*   Q;          //quaternion array
    
    TYPE*	C;

	TYPE	gravity;
	TYPE	stiffness_0, stiffness_1, stiffness_2, stiffness_3;
	TYPE	stiffness_p;

	TYPE*	F_Temp;
    TYPE*   C_Temp;
    TYPE*   TK;         // K per tetrahedron
    TYPE*   TU;         // U per tetrahedron
    TYPE*   TV;         // V per tetrahedron
    TYPE*   T_alpha;    // alpha per tetrahedron
    TYPE*   T_beta;     // beta  per tetrahedron
    TYPE*   T_gamma;    // gamma per tetrahedron
    TYPE*   T_A;
    TYPE*   T_B01;
    TYPE*   T_B02;
    TYPE*   T_B12;
    
	int*	VTT;		// The index list mapping to Tet_Temp
	int*	vtt_num;

	int*	E;
	int*	TE;
	int		e_number;
	int*	VV;
	int*	VE;
	int*	vv_num;

    TYPE*   cg_r;
    TYPE*   cg_z;
    TYPE*   cg_p;
    TYPE*   cg_Ap;
    
	LBFGS<TYPE>*			lbfgs_solver;
	PCG_SOLVER<TYPE>		pcg_solver;
    
	HYPER_TET_MESH()
	{
		lambda	= new TYPE	[max_number* 3];
		ext_C	= new TYPE	[max_number* 3];

        next_X  = new TYPE  [max_number* 3];
        prev_X  = new TYPE  [max_number* 3];
        Y       = new TYPE  [max_number* 3];
		prev_Y	= new TYPE	[max_number* 3];
		B       = new TYPE	[max_number* 3];
		V		= new TYPE	[max_number* 3];
		prev_V	= new TYPE	[max_number* 3];
        S       = new TYPE  [max_number* 3];
		C		= new TYPE	[max_number* 3];
		D		= new TYPE	[max_number* 3];
        Q       = new TYPE  [max_number*12];
        Force   = new TYPE  [max_number* 3];
		fixed	= new int	[max_number   ];
        fixed_X = new TYPE  [max_number* 3];
        
		E			= new int	[max_number*2];
		TE			= new int	[max_number*6];
		VV			= new int	[max_number*6];
		VE			= new int	[max_number*6];
		vv_num		= new int	[max_number];

        
		F_Temp  = new TYPE	[max_number*24];
        C_Temp  = new TYPE  [max_number*24];
        TK      = new TYPE  [max_number*144];
        TU      = new TYPE  [max_number* 9];
        TV      = new TYPE  [max_number* 9];
        T_alpha = new TYPE  [max_number* 9];
        T_beta  = new TYPE  [max_number* 9];
        T_gamma = new TYPE  [max_number* 9];
        T_A     = new TYPE  [max_number* 9];
        T_B01   = new TYPE  [max_number* 4];
        T_B02   = new TYPE  [max_number* 4];
        T_B12   = new TYPE  [max_number* 4];

		VTT		= new int	[max_number* 4];
		vtt_num	= new int	[max_number   ];
        
        cg_r    = new TYPE  [max_number* 3];
        cg_z    = new TYPE  [max_number* 3];
        cg_p    = new TYPE  [max_number* 3];
        cg_Ap    = new TYPE  [max_number* 3];

		gravity		= -9.8;
        
		stiffness_0	= 2000000;//2000000
        stiffness_1	= 1000000;//2000000
		stiffness_2 = 4000000;
		stiffness_3 = 3.5;

		memset(		V, 0, sizeof(TYPE)*max_number*3);
		memset(prev_V, 0, sizeof(TYPE)*max_number*3);
		memset(	fixed, 0, sizeof(int )*max_number  );

		lbfgs_solver	= 0;
	}
	
    
	~HYPER_TET_MESH()
	{
		if(lambda)		delete[] lambda;
		if(ext_C)		delete[] ext_C;

        if(next_X)      delete[] next_X;
        if(prev_X)      delete[] prev_X;
        if(Y)           delete[] Y;
		if(prev_Y)		delete[] prev_Y;
		if(B)           delete[] B;
        if(Force)       delete[] Force;
		if(V)			delete[] V;
		if(prev_V)		delete[] prev_V;
        if(S)           delete[] S;
        if(Q)           delete[] Q;
		if(C)			delete[] C;
		if(D)			delete[] D;
		if(fixed)		delete[] fixed;
        if(fixed_X)     delete[] fixed_X;
		if(F_Temp)      delete[] F_Temp;
        if(C_Temp)      delete[] C_Temp;
        if(TK)          delete[] TK;
        if(TU)          delete[] TU;
        if(TV)          delete[] TV;

		if(E)			delete[] E;
		if(TE)			delete[] TE;
		if(VV)			delete[] VV;
		if(VE)			delete[] VE;
		if(vv_num)		delete[] vv_num;
        
        if(T_alpha)     delete[] T_alpha;
        if(T_beta)      delete[] T_beta;
        if(T_gamma)     delete[] T_gamma;
        if(T_A)         delete[] T_A;
        if(T_B01)       delete[] T_B01;
        if(T_B02)       delete[] T_B02;
        if(T_B12)       delete[] T_B12;
        
		if(VTT)			delete[] VTT;
		if(vtt_num)		delete[] vtt_num;
        
        if(cg_r)        delete[] cg_r;
        if(cg_z)        delete[] cg_z;
        if(cg_p)        delete[] cg_p;
        if(cg_Ap)       delete[] cg_Ap;

		if(lbfgs_solver)delete lbfgs_solver;
	}

///////////////////////////////////////////////////////////////////////////////////////////
//  Initialize functions
///////////////////////////////////////////////////////////////////////////////////////////
	
	void Initialize()
	{
		TET_MESH<TYPE>::Initialize();
		Build_VTT();
		Build_Edge_From_Tet();
		Build_VV();
        
        for(int t=0; t<tet_number; t++)
        {
            Q[t*4+0]=0;
            Q[t*4+1]=0;
            Q[t*4+2]=0;
            Q[t*4+3]=1;
        }
        memcpy(fixed_X, X, sizeof(TYPE)*number*3);

		lbfgs_solver=new LBFGS<TYPE>(number*3);

		pcg_solver.Create(number, e_number, vv_num, VV, VE);

	}

	void Build_VV()
	{
		//First set vv_num
		memset(vv_num, 0, sizeof(int)*max_number);
		for(int i=0; i<e_number; i++)
		{
			vv_num[E[i*2+0]]++;
			vv_num[E[i*2+1]]++;
		}
		for(int i=1; i<number; i++)
			vv_num[i]+=vv_num[i-1];
		for(int i=number; i>0; i--)
			vv_num[i]=vv_num[i-1];
		vv_num[0]=0;

		//Then set vv
		int *_vv_num=new int[max_number];
		memcpy(_vv_num, vv_num, sizeof(int)*max_number);
		for(int i=0; i<e_number; i++)
		{
			VV[_vv_num[E[i*2+0]]]=E[i*2+1];
			VV[_vv_num[E[i*2+1]]]=E[i*2+0];
			VE[_vv_num[E[i*2+0]]++]=i;
			VE[_vv_num[E[i*2+1]]++]=i;
		}
		delete []_vv_num;
	}

	void Build_Edge_From_Tet()
	{
		int*	_E=new int[tet_number*18];	
		for(int t=0; t<tet_number; t++)
		{
			int v0=Tet[t*4+0];
			int v1=Tet[t*4+1];
			int v2=Tet[t*4+2];
			int v3=Tet[t*4+3];

			_E[t*18+ 0]=v0; _E[t*18+ 1]=v1; _E[t*18+ 2]=t;
			_E[t*18+ 3]=v0; _E[t*18+ 4]=v2; _E[t*18+ 5]=t;
			_E[t*18+ 6]=v0; _E[t*18+ 7]=v3; _E[t*18+ 8]=t;
			_E[t*18+ 9]=v1; _E[t*18+10]=v2; _E[t*18+11]=t;
			_E[t*18+12]=v1; _E[t*18+13]=v3; _E[t*18+14]=t;
			_E[t*18+15]=v2; _E[t*18+16]=v3; _E[t*18+17]=t;
		}

		for(int i=0; i<tet_number*6; i++)		
			if(_E[i*3+0]>_E[i*3+1])	Swap(_E[i*3+0], _E[i*3+1]);
		

		//Quicksort
		Quick_Sort_RE(_E, 0, tet_number*6-1);		

		e_number=0;
		for(int i=0; i<tet_number*6; i++)
		{
			if(i!=0 && _E[i*3]==_E[(i-1)*3] && _E[i*3+1]== _E[(i-1)*3+1])
			{				
				//Add the edge to  TE
				int v0=Tet[_E[i*3+2]*4+0];
				int v1=Tet[_E[i*3+2]*4+1];
				int v2=Tet[_E[i*3+2]*4+2];
				int v3=Tet[_E[i*3+2]*4+3];
				int p0=_E[i*3+0];
				int p1=_E[i*3+1];

					 if(v0==p0 && v1==p1 || v0==p1 && v1==p0)	TE[_E[i*3+2]*6+0]=e_number-1;
				else if(v0==p0 && v2==p1 || v0==p1 && v2==p0)	TE[_E[i*3+2]*6+1]=e_number-1;
				else if(v0==p0 && v3==p1 || v0==p1 && v3==p0)	TE[_E[i*3+2]*6+2]=e_number-1;
				else if(v1==p0 && v2==p1 || v1==p1 && v2==p0)	TE[_E[i*3+2]*6+3]=e_number-1;
				else if(v1==p0 && v3==p1 || v1==p1 && v3==p0)	TE[_E[i*3+2]*6+4]=e_number-1;
				else if(v2==p0 && v3==p1 || v2==p1 && v3==p0)	TE[_E[i*3+2]*6+5]=e_number-1;
				
				else{	printf("ERROR: unknown TEa.\n");		getchar();}
			}
			else
			{
				//Add the edge to E
				E[e_number*2+0]=_E[i*3+0];
				E[e_number*2+1]=_E[i*3+1];

				//Add the edge to  TE
				int v0=Tet[_E[i*3+2]*4+0];
				int v1=Tet[_E[i*3+2]*4+1];
				int v2=Tet[_E[i*3+2]*4+2];
				int v3=Tet[_E[i*3+2]*4+3];
				int p0=_E[i*3+0];
				int p1=_E[i*3+1];

					 if(v0==p0 && v1==p1 || v0==p1 && v1==p0)	TE[_E[i*3+2]*6+0]=e_number;
				else if(v0==p0 && v2==p1 || v0==p1 && v2==p0)	TE[_E[i*3+2]*6+1]=e_number;
				else if(v0==p0 && v3==p1 || v0==p1 && v3==p0)	TE[_E[i*3+2]*6+2]=e_number;
				else if(v1==p0 && v2==p1 || v1==p1 && v2==p0)	TE[_E[i*3+2]*6+3]=e_number;
				else if(v1==p0 && v3==p1 || v1==p1 && v3==p0)	TE[_E[i*3+2]*6+4]=e_number;
				else if(v2==p0 && v3==p1 || v2==p1 && v3==p0)	TE[_E[i*3+2]*6+5]=e_number;
				
				else{	printf("ERROR: unknown TEb (%d, %d, %d, %d) (%d, %d, %d).\n", v0, v1, v2, v3, p0, p1, _E[i*3+2]);		getchar();}
				e_number++;
			}
		}
		delete []_E;

		//for(int t=0; t<t_number; t++)
		//{
		//	int t3=t*3;
		//	printf("V %d: %d, %d, %d\n", t, T[t3+0], T[t3+1], T[t3+2]);
		//	printf("value %d: %d, %d\n", TE[t3+0], RE[TE[t3+0]*2+0], RE[TE[t3+0]*2+1]);
		//	printf("value %d: %d, %d\n", TE[t3+1], RE[TE[t3+1]*2+0], RE[TE[t3+1]*2+1]);
		//	printf("value %d: %d, %d\n", TE[t3+2], RE[TE[t3+2]*2+0], RE[TE[t3+2]*2+1]);
		//}
		//for(int i=0; i<t_number; i++)
		//{
		//	printf("TE count: %d\n", TE_counter[i]);
		//	printf("TE: %d, %d, %d\n", TE[i*3+0], TE[i*3+1], TE[i*3+2]);
		//	printf("Edge: %d, %d; %d, %d; %d, %d\n", 
		//		RE[TE[i*3+0]*2+0], RE[TE[i*3+0]*2+1], 
		//		RE[TE[i*3+1]*2+0], RE[TE[i*3+1]*2+1],
		//		RE[TE[i*3+2]*2+0], RE[TE[i*3+2]*2+1]); 
		//		printf("edge: %d, %d\n", RE[i*3], RE[i*3+1]);
		//		getchar();
		//}
		//for(int i=0; i<te_number; i++)
		//	printf("edge: %d, %d\n", TE[i*2], TE[i*2+1]);
		//delete []TE_counter;
	}
	
	void Quick_Sort_RE( int a[], int l, int r)
	{
		int j;
		if( l < r ) 
		{
			j = Quick_Sort_Partition_RE( a, l, r);
			Quick_Sort_RE( a, l, j-1);
			Quick_Sort_RE( a, j+1, r);
		}
	}
	
	int Quick_Sort_Partition_RE( int a[], int l, int r) 
	{
		int pivot[3], i, j, t[3];
		pivot[0] = a[l*3+0];
		pivot[1] = a[l*3+1];
		pivot[2] = a[l*3+2];	
		i = l; j = r+1;		
		while( 1)
		{
			do ++i; while( (a[i*3]<pivot[0] || a[i*3]==pivot[0] && a[i*3+1]<=pivot[1]) && i <= r );
			do --j; while(  a[j*3]>pivot[0] || a[j*3]==pivot[0] && a[j*3+1]> pivot[1] );
			if( i >= j ) break;
			//Swap i and j			
			t[0]=a[i*3+0];
			t[1]=a[i*3+1];
			t[2]=a[i*3+2];
			a[i*3+0]=a[j*3+0];
			a[i*3+1]=a[j*3+1];
			a[i*3+2]=a[j*3+2];
			a[j*3+0]=t[0];
			a[j*3+1]=t[1];
			a[j*3+2]=t[2];
		}
		//Swap l and j
		t[0]=a[l*3+0];
		t[1]=a[l*3+1];
		t[2]=a[l*3+2];
		a[l*3+0]=a[j*3+0];
		a[l*3+1]=a[j*3+1];
		a[l*3+2]=a[j*3+2];
		a[j*3+0]=t[0];
		a[j*3+1]=t[1];
		a[j*3+2]=t[2];
		return j;
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
   

///////////////////////////////////////////////////////////////////////////////////////////
//  Constraint functions
///////////////////////////////////////////////////////////////////////////////////////////
    TYPE Energy(TYPE I, TYPE II, TYPE III, TYPE stiffness_0, TYPE stiffness_1, TYPE stiffness_2, TYPE stiffness_3)
	{
		TYPE J		= sqrt(III);
		TYPE rcbrtJ = 1.0 / cbrt(J);
		TYPE bar_I	= I *rcbrtJ*rcbrtJ;
		TYPE bar_II = (I*I-II)*rcbrtJ*rcbrtJ*rcbrtJ*rcbrtJ*0.5;
		//return stiffness_0*(bar_I - 3) + stiffness_1*(bar_II - 3);

		//printf("I: %f, %f, %f (%f)\n", I, II, III, expf(stiffness_3*(bar_I-3)));

		return stiffness_0*(bar_I-3) + stiffness_2*(expf(stiffness_3*(bar_I-3))-1);
	}

	TYPE Energy(TYPE Sigma[3], TYPE stiffness_0, TYPE stiffness_1, TYPE stiffness_2, TYPE stiffness_3)
	{
		TYPE J = Sigma[0] * Sigma[1] * Sigma[2];
		TYPE I = Sigma[0] * Sigma[0] + Sigma[1] * Sigma[1] + Sigma[2] * Sigma[2];
		TYPE II = Sigma[0] * Sigma[0] * Sigma[0] * Sigma[0] + Sigma[1] * Sigma[1] * Sigma[1] * Sigma[1] + Sigma[2] * Sigma[2] * Sigma[2] * Sigma[2];
		TYPE III = J*J;
		return Energy(I, II, III, stiffness_0, stiffness_1, stiffness_2, stiffness_3);
	}

	TYPE Numerical_Estimate_dWdI(TYPE I, TYPE II, TYPE III, int id, TYPE stiffness_0, TYPE stiffness_1, TYPE stiffness_2, TYPE stiffness_3)
	{
		TYPE epsilon = 0.01;

		if(id==0)	I  -= epsilon*0.5;
		if(id==1)	II -= epsilon*0.5;
		if(id==2)	III-= epsilon*0.5;
		TYPE value0 = Energy(I, II, III, stiffness_0, stiffness_1, stiffness_2, stiffness_3);

		if(id==0)	I  += epsilon;
		if(id==1)	II += epsilon;
		if(id==2)	III+= epsilon;
		TYPE value1 = Energy(I, II, III, stiffness_0, stiffness_1, stiffness_2, stiffness_3);

		return (value1-value0)/epsilon;
	}

	TYPE Numerical_Estimate_dW2dI2(TYPE I, TYPE II, TYPE III, int id0, int id1, TYPE stiffness_0, TYPE stiffness_1, TYPE stiffness_2, TYPE stiffness_3)
	{
		TYPE epsilon = 0.01;
		if (id0 == 0)	I   -= epsilon*0.5;
		if (id0 == 1)	II  -= epsilon*0.5;
		if (id0 == 2)	III -= epsilon*0.5;
		TYPE value0 = Numerical_Estimate_dWdI(I, II, III, id1, stiffness_0, stiffness_1, stiffness_2, stiffness_3);

		if (id0 == 0)	I   += epsilon;
		if (id0 == 1)	II  += epsilon;
		if (id0 == 2)	III += epsilon;
		TYPE value1 = Numerical_Estimate_dWdI(I, II, III, id1, stiffness_0, stiffness_1, stiffness_2, stiffness_3);

		//printf("value: %f, %f\n", value0, value1);

		return (value1 - value0) / epsilon;
	}

	TYPE Numerical_Estimate_dWdS(TYPE Sigma[3], int id, TYPE stiffness_0, TYPE stiffness_1, TYPE stiffness_2, TYPE stiffness_3)
	{
		TYPE epsilon = 0.00001;

		Sigma[id] -= epsilon*0.5;
		TYPE value0 = Energy(Sigma, stiffness_0, stiffness_1, stiffness_2, stiffness_3);
		Sigma[id] +=epsilon;
		TYPE value1 = Energy(Sigma, stiffness_0, stiffness_1, stiffness_2, stiffness_3);

		return (value1-value0)/epsilon;
	}

    TYPE Numerical_Estimate_K(int t, int i, int j)
    {
        TYPE epsilon = 0.0001;
        int pi      = Tet[t*4+i/3]*3+i%3;
        TYPE old_Xi = X[pi];
        X[pi]+=epsilon*0.5;
        Compute_FM(t, stiffness_0, stiffness_1, stiffness_2, stiffness_3, false);
        TYPE value0=F_Temp[t*12+j];
        X[pi]-=epsilon;
        Compute_FM(t, stiffness_0, stiffness_1, stiffness_2, stiffness_3, false);
        TYPE value1=F_Temp[t*12+j];
        X[pi]=old_Xi;
        return (value0-value1)/epsilon;
    }
    

    TYPE Compute_FM(int t, TYPE stiffness_0, TYPE stiffness_1, TYPE stiffness_2, TYPE stiffness_3, const bool has_matrix=true)
    {
        stiffness_0=-Vol[t]*stiffness_0;
        stiffness_1=-Vol[t]*stiffness_1;
        stiffness_2=-Vol[t]*stiffness_2;

        //No velocity access in this function
        int p0=Tet[t*4+0]*3;
        int p1=Tet[t*4+1]*3;
        int p2=Tet[t*4+2]*3;
        int p3=Tet[t*4+3]*3;
        
        TYPE Ds[9];
        Ds[0]=X[p1+0]-X[p0+0];
        Ds[3]=X[p1+1]-X[p0+1];
        Ds[6]=X[p1+2]-X[p0+2];
        Ds[1]=X[p2+0]-X[p0+0];
        Ds[4]=X[p2+1]-X[p0+1];
        Ds[7]=X[p2+2]-X[p0+2];
        Ds[2]=X[p3+0]-X[p0+0];
        Ds[5]=X[p3+1]-X[p0+1];
        Ds[8]=X[p3+2]-X[p0+2];
        
        TYPE F[9], U[9], Sigma[3], V[9];
        Matrix_Product_3(Ds, &inv_Dm[t*9], F);
        
        memcpy(U, F, sizeof(TYPE)*9);
        SVD3((TYPE (*)[3])U, Sigma, (TYPE (*)[3])V);
        int small_id;
        if(fabsf(Sigma[0])<fabsf(Sigma[1]) && fabsf(Sigma[0])<fabsf(Sigma[2]))	small_id=0;
        else if(fabsf(Sigma[1])<fabsf(Sigma[2]))                                small_id=1;
        else                                                                    small_id=2;
        if(U[0]*(U[4]*U[8]-U[7]*U[5])+U[3]*(U[7]*U[2]-U[1]*U[8])+U[6]*(U[1]*U[5]-U[4]*U[2])<0)
        {
            U[0+small_id]	=-U[0+small_id];
            U[3+small_id]	=-U[3+small_id];
            U[6+small_id]	=-U[6+small_id];
			Sigma[small_id]	=-Sigma[small_id];
        }
        if(V[0]*(V[4]*V[8]-V[7]*V[5])+V[3]*(V[7]*V[2]-V[1]*V[8])+V[6]*(V[1]*V[5]-V[4]*V[2])<0)
        {
            V[0+small_id]=-V[0+small_id];
            V[3+small_id]=-V[3+small_id];
            V[6+small_id]=-V[6+small_id];
			Sigma[small_id] = -Sigma[small_id];
        }
        
        //SVD3x3((TYPE (*)[3])F, (TYPE (*)[3])U, Sigma, &Q[t*4], (TYPE (*)[3])V, svd_iterations);

		float interpolator0, interpolator1, interpolator2;
		if (Sigma[0] < 1)	interpolator0 = (0.15 - Sigma[0])*10.0;
		if (Sigma[1] < 1)	interpolator1 = (0.15 - Sigma[1])*10.0;
		if (Sigma[2] < 1)	interpolator2 = (0.15 - Sigma[2])*10.0;

		float interpolator = MAX(MAX(interpolator0, interpolator1), interpolator2);
		if(interpolator<0)		interpolator=0;
		if(interpolator>1)		interpolator=1;

		//lambda[t] = MAX(lambda[t], interpolator);
		//lambda[t] = 0;
		//interpolator = lambda[t];

		interpolator = 0;
		stiffness_0 *= (1 - interpolator);
		stiffness_1 *= (1 - interpolator);

		//printf("sigma; %f, %f, %f\n", Sigma[0], Sigma[1], Sigma[2]);


        TYPE I			= Sigma[0]*Sigma[0]+Sigma[1]*Sigma[1]+Sigma[2]*Sigma[2];
        TYPE J			= Sigma[0]*Sigma[1]*Sigma[2];
        TYPE II			= Sigma[0]*Sigma[0]*Sigma[0]*Sigma[0]+Sigma[1]*Sigma[1]*Sigma[1]*Sigma[1]+Sigma[2]*Sigma[2]*Sigma[2]*Sigma[2];
		TYPE III		= J*J;
		TYPE rcbrt_III	= 1.0/cbrt(III);
		TYPE factor_1	= ONE_THIRD*rcbrt_III/III;
		TYPE factor_2	= ONE_THIRD* factor_1/III;
		TYPE dEdI		= 0;
		TYPE dEdII		= 0;
		TYPE dEdIII		= 0;
		TYPE H[3][3];
		memset(&H[0][0], 0, sizeof(TYPE) * 9);
		TYPE energy		= 0;

		
		// StVK
	//	dEdI	+= stiffness_0*(I-3)*0.25-stiffness_1*0.5;
	//	dEdII	+= stiffness_1*0.25;
	//	dEdIII	+= 0;
	//	H[0][0]	+= stiffness_0*0.25;
	//	energy-=stiffness_0*(I-3)*(I-3)*0.125+stiffness_1*(II-2*I+3)*0.25;

        
		//Neo-Hookean
    	dEdI	+= rcbrt_III*stiffness_0;
    	dEdII	+= 0;
    	dEdIII	+= -factor_1*stiffness_0*I;
    	H[0][2]	+= -factor_1*stiffness_0;
    	H[2][2]	+=  factor_2*stiffness_0*I*4;
     
		energy=-(stiffness_0*(I *rcbrt_III-3) + stiffness_1*(J-1)*(J-1));


		// Mooney-Rivlin		
	//	TYPE two_term_a	= stiffness_2*rcbrt_III*I;
	//	TYPE two_term_b	= stiffness_2*rcbrt_III*(I*I-II);
	//	dEdI	+= rcbrt_III*two_term_a;
	//	dEdII	+= -0.5*stiffness_2*rcbrt_III*rcbrt_III;
	//	dEdIII	+= -factor_1*two_term_b;
	//	H[0][0]	+= stiffness_2*rcbrt_III*rcbrt_III;
	//	H[0][2]	+= -factor_1*two_term_a*2;
	//	H[1][2]	+= factor_1*stiffness_2*rcbrt_III;
	//	H[2][2]	+= factor_2*two_term_b*5;
	//	energy-=stiffness_2*(0.5*rcbrt_III*rcbrt_III*(I*I-II)-3);
		
		// Fung
	//	TYPE exp_term	= expf(stiffness_3*(rcbrt_III*I-3))*stiffness_2*stiffness_3;
	//	dEdI	+= exp_term*rcbrt_III;
	//	dEdIII	+= -factor_1*I*exp_term;
	//	H[0][0]	+= rcbrt_III*stiffness_3*rcbrt_III*exp_term;
	//	H[0][2]	+= -factor_1*exp_term*(1+stiffness_3*rcbrt_III*I);
	//	H[2][2]	+= factor_2*I*exp_term* (4 + stiffness_3*I*rcbrt_III);
	//	energy-=stiffness_2*(exp(stiffness_3*(rcbrt_III*I-3))-1);

		// Volume correction
		dEdIII	+=stiffness_1*(J-1)/J;
		H[2][2]	+=stiffness_1/(2*III*J);

		// Make H symmetric
		H[1][0]	= H[0][1];
		H[2][0]	= H[0][2];
		H[2][1]	= H[1][2];

		TYPE P0 = 2 * (dEdI*Sigma[0] + 2 * dEdII*Sigma[0] * Sigma[0] * Sigma[0] + dEdIII*III / Sigma[0]);
		TYPE P1 = 2 * (dEdI*Sigma[1] + 2 * dEdII*Sigma[1] * Sigma[1] * Sigma[1] + dEdIII*III / Sigma[1]);
		TYPE P2 = 2 * (dEdI*Sigma[2] + 2 * dEdII*Sigma[2] * Sigma[2] * Sigma[2] + dEdIII*III / Sigma[2]);


        //printf("update: %d\n", has_matrix);
		//printf("check0: %f, %f\n", P0, Numerical_Estimate_dWdS(Sigma, 0, stiffness_0, stiffness_1, stiffness_2, stiffness_3));
		//printf("check1: %f, %f\n", P1, Numerical_Estimate_dWdS(Sigma, 1, stiffness_0, stiffness_1, stiffness_2, stiffness_3));
		//printf("check2: %f, %f\n", P2, Numerical_Estimate_dWdS(Sigma, 2, stiffness_0, stiffness_1, stiffness_2, stiffness_3));
		//printf("C: %f, %f\n", H[0][0], Numerical_Estimate_dW2dI2(I, II, III, 0, 0, stiffness_0, stiffness_1, stiffness_2, stiffness_3));
		//printf("C: %f, %f\n", H[0][2], Numerical_Estimate_dW2dI2(I, II, III, 0, 2, stiffness_0, stiffness_1, stiffness_2, stiffness_3));
		//printf("C: %f, %f\n", H[1][2], Numerical_Estimate_dW2dI2(I, II, III, 1, 2, stiffness_0, stiffness_1, stiffness_2, stiffness_3));
		//printf("C: %f, %f\n", H[2][2], Numerical_Estimate_dW2dI2(I, II, III, 2, 2, stiffness_0, stiffness_1, stiffness_2, stiffness_3));
		//getchar();


        TYPE PV_transpose[9], P[9], force[9];
        PV_transpose[0]=P0*V[0];
        PV_transpose[1]=P0*V[3];
        PV_transpose[2]=P0*V[6];
        PV_transpose[3]=P1*V[1];
        PV_transpose[4]=P1*V[4];
        PV_transpose[5]=P1*V[7];
        PV_transpose[6]=P2*V[2];
        PV_transpose[7]=P2*V[5];
        PV_transpose[8]=P2*V[8];
        Matrix_Product_3(U, PV_transpose, P);
        Matrix_Product_T_3(P, &inv_Dm[t*9], force);
        
        F_Temp[t*12+ 0]=-(force[0]+force[1]+force[2]);
        F_Temp[t*12+ 1]=-(force[3]+force[4]+force[5]);
        F_Temp[t*12+ 2]=-(force[6]+force[7]+force[8]);
        F_Temp[t*12+ 3]=force[0];
        F_Temp[t*12+ 4]=force[3];
        F_Temp[t*12+ 5]=force[6];
        F_Temp[t*12+ 6]=force[1];
        F_Temp[t*12+ 7]=force[4];
        F_Temp[t*12+ 8]=force[7];
        F_Temp[t*12+ 9]=force[2];
        F_Temp[t*12+10]=force[5];
        F_Temp[t*12+11]=force[8];
        

		//Force[p0+0]+=-(force[0]+force[1]+force[2]);
        //Force[p0+1]+=-(force[3]+force[4]+force[5]);
        //Force[p0+2]+=-(force[6]+force[7]+force[8]);
        //Force[p1+0]+=force[0];
        //Force[p1+1]+=force[3];
        //Force[p1+2]+=force[6];
        //Force[p2+0]+=force[1];
        //Force[p2+1]+=force[4];
        //Force[p2+2]+=force[7];
        //Force[p3+0]+=force[2];
        //Force[p3+1]+=force[5];
        //Force[p3+2]+=force[8];
        
        //Copy UV into TU and TV
        memcpy(&TU[t*9], U, sizeof(TYPE)*9);
        memcpy(&TV[t*9], V, sizeof(TYPE)*9);
        

		//Strain limiting part
		TYPE new_R[9];
		Matrix_Product_T_3(U, V, new_R);
		const TYPE* idm = &inv_Dm[t * 9];
		TYPE half_matrix[3][4], result_matrix[3][4];
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

		Matrix_Substract_3(new_R, F, new_R);
		Matrix_Product(new_R, &half_matrix[0][0], &result_matrix[0][0], 3, 3, 4);

		TYPE pd_stiffness = Vol[t] * stiffness_p * interpolator;
		//Force[p0+0] += result_matrix[0][0] * pd_stiffness;
		//Force[p0+1] += result_matrix[1][0] * pd_stiffness;
		//Force[p0+2] += result_matrix[2][0] * pd_stiffness;
		//Force[p1+0] += result_matrix[0][1] * pd_stiffness;
		//Force[p1+1] += result_matrix[1][1] * pd_stiffness;
		//Force[p1+2] += result_matrix[2][1] * pd_stiffness;
		//Force[p2+0] += result_matrix[0][2] * pd_stiffness;
		//Force[p2+1] += result_matrix[1][2] * pd_stiffness;
		//Force[p2+2] += result_matrix[2][2] * pd_stiffness;
		//Force[p3+0] += result_matrix[0][3] * pd_stiffness;
		//Force[p3+1] += result_matrix[1][3] * pd_stiffness;
		//Force[p3+2] += result_matrix[2][3] * pd_stiffness;
	
		//PD energy
		TYPE pd_energy=0;
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



		TYPE value0 = half_matrix[0][0]*half_matrix[0][0]+half_matrix[1][0]*half_matrix[1][0]+half_matrix[2][0]*half_matrix[2][0];
		TYPE value1 = half_matrix[0][1]*half_matrix[0][1]+half_matrix[1][1]*half_matrix[1][1]+half_matrix[2][1]*half_matrix[2][1];
		TYPE value2 = half_matrix[0][2]*half_matrix[0][2]+half_matrix[1][2]*half_matrix[1][2]+half_matrix[2][2]*half_matrix[2][2];
		TYPE value3 = half_matrix[0][3]*half_matrix[0][3]+half_matrix[1][3]*half_matrix[1][3]+half_matrix[2][3]*half_matrix[2][3];
		ext_C[Tet[t * 4 + 0]]+=value0*pd_stiffness;
		ext_C[Tet[t * 4 + 1]]+=value1*pd_stiffness;
		ext_C[Tet[t * 4 + 2]]+=value2*pd_stiffness;
		ext_C[Tet[t * 4 + 3]]+=value3*pd_stiffness;



        if(has_matrix==false)	return energy;

        
		//Now compute the stiffness matrix
        TYPE alpha[3][3], beta[3][3], gamma[3][3];
        for(int i=0; i<3; i++)
        for(int j=i; j<3; j++)
        {
            alpha[i][j]=2*dEdI+4*(Sigma[i]*Sigma[i]+Sigma[j]*Sigma[j])*dEdII;
            beta[i][j]=4*Sigma[i]*Sigma[j]*dEdII-2*III*dEdIII/(Sigma[i]*Sigma[j]);
                
            TYPE vi[3]={2*Sigma[i], 4*Sigma[i]*Sigma[i]*Sigma[i], 2*III/Sigma[i]};
            TYPE vj[3]={2*Sigma[j], 4*Sigma[j]*Sigma[j]*Sigma[j], 2*III/Sigma[j]};
            TYPE r[3];
            Matrix_Vector_Product_3(&H[0][0], vj, r);
            gamma[i][j]=DOT(vi, r)+4*III*dEdIII/(Sigma[i]*Sigma[j]);
        }
        
        // Save alpha, beta, and gamma
        memcpy(&T_alpha[t*9], &alpha[0][0], sizeof(TYPE)*9);
        memcpy(&T_beta [t*9], &beta [0][0], sizeof(TYPE)*9);
        memcpy(&T_gamma[t*9], &gamma[0][0], sizeof(TYPE)*9);
        
        T_A[t*9+0]=alpha[0][0]+beta[0][0]+gamma[0][0];
        T_A[t*9+1]=gamma[0][1];
        T_A[t*9+2]=gamma[0][2];
        T_A[t*9+3]=gamma[0][1];
        T_A[t*9+4]=alpha[1][1]+beta[1][1]+gamma[1][1];
        T_A[t*9+5]=gamma[1][2];
        T_A[t*9+6]=gamma[0][2];
        T_A[t*9+7]=gamma[1][2];
        T_A[t*9+8]=alpha[2][2]+beta[2][2]+gamma[2][2];
        
        T_B01[t*4+0]=alpha[0][1];
        T_B01[t*4+1]= beta[0][1];
        T_B01[t*4+2]= beta[0][1];
        T_B01[t*4+3]=alpha[0][1];        
        T_B02[t*4+0]=alpha[0][2];
        T_B02[t*4+1]= beta[0][2];
        T_B02[t*4+2]= beta[0][2];
        T_B02[t*4+3]=alpha[0][2];        
        T_B12[t*4+0]=alpha[1][2];
        T_B12[t*4+1]= beta[1][2];
        T_B12[t*4+2]= beta[1][2];
        T_B12[t*4+3]=alpha[1][2];      

		//Fix...
		//printf("t: %d; %f, %f, %f\n", t, Sigma[0], Sigma[1], Sigma[2]);		
		//eigen_project(&T_A[t*9+0], 3);
		//eigen_project(&T_B01[t*4+0], 2);
		//eigen_project(&T_B02[t*4+0], 2);
		//eigen_project(&T_B12[t*4+0], 2);
		//Fix...



        TYPE dGdF[12][9];    //G is related to force, according to [TSIF05], (g0, g3, g6), (g1, g4, g7), (g2, g5, g8), (g9, g10, g11)
        for(int i=0; i<9; i++)
        {
            TYPE dF[9], temp0[9], temp1[9];
            memset(&dF, 0, sizeof(TYPE)*9);
            dF[i]=1;
            
            Matrix_Product_3(dF, V, temp0);
            Matrix_T_Product_3(U, temp0, temp1);
            
            temp0[0]=T_A[t*9+0]*temp1[0]+T_A[t*9+1]*temp1[4]+T_A[t*9+2]*temp1[8];
            temp0[4]=T_A[t*9+3]*temp1[0]+T_A[t*9+4]*temp1[4]+T_A[t*9+5]*temp1[8];
            temp0[8]=T_A[t*9+6]*temp1[0]+T_A[t*9+7]*temp1[4]+T_A[t*9+8]*temp1[8];

            temp0[1]=T_B01[t*4+0]*temp1[1]+T_B01[t*4+1]*temp1[3];
            temp0[3]=T_B01[t*4+2]*temp1[1]+T_B01[t*4+3]*temp1[3];
            
            temp0[2]=T_B02[t*4+0]*temp1[2]+T_B02[t*4+1]*temp1[6];
            temp0[6]=T_B02[t*4+2]*temp1[2]+T_B02[t*4+3]*temp1[6];
            
            temp0[5]=T_B12[t*4+0]*temp1[5]+T_B12[t*4+1]*temp1[7];
            temp0[7]=T_B12[t*4+2]*temp1[5]+T_B12[t*4+3]*temp1[7];
            
            Matrix_Product_T_3(temp0, V, temp1);
            Matrix_Product_3(U, temp1, temp0);
            Matrix_Product_T_3(temp0, &inv_Dm[t*9], &dGdF[i][0]);
        }
        
        //for(int i=0; i<9; i++)
        //    printf("dPdF %d: %f, %f, %f; %f, %f, %f; %f, %f, %f\n", i,
        //           dPdF[t*81+i*9+0], dPdF[t*81+i*9+1], dPdF[t*81+i*9+2], dPdF[t*81+i*9+3], dPdF[t*81+i*9+4], dPdF[t*81+i*9+5], dPdF[t*81+i*9+6], dPdF[t*81+i*9+7], dPdF[t*81+i*9+8]);
        
        //Transpose dGdF
        TYPE temp;
        for(int i=0; i<9; i++) for(int j=i+1; j<9; j++)
            SWAP(dGdF[i][j], dGdF[j][i]);
        
        for(int j=0; j< 9; j++)
        {
            dGdF[ 9][j]=-dGdF[0][j]-dGdF[1][j]-dGdF[2][j];
            dGdF[10][j]=-dGdF[3][j]-dGdF[4][j]-dGdF[5][j];
            dGdF[11][j]=-dGdF[6][j]-dGdF[7][j]-dGdF[8][j];
        }
        
        TYPE new_idm[4][3];
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
        
        C_Temp[t*12+ 0]=-Matrix_Product_T(&new_idm[0][0], &dGdF[ 9][0], 4, 3, 3, 0, 0);
        C_Temp[t*12+ 1]=-Matrix_Product_T(&new_idm[0][0], &dGdF[10][0], 4, 3, 3, 0, 1);
        C_Temp[t*12+ 2]=-Matrix_Product_T(&new_idm[0][0], &dGdF[11][0], 4, 3, 3, 0, 2);
        C_Temp[t*12+ 3]=-Matrix_Product_T(&new_idm[0][0], &dGdF[ 0][0], 4, 3, 3, 1, 0);
        C_Temp[t*12+ 4]=-Matrix_Product_T(&new_idm[0][0], &dGdF[ 3][0], 4, 3, 3, 1, 1);
        C_Temp[t*12+ 5]=-Matrix_Product_T(&new_idm[0][0], &dGdF[ 6][0], 4, 3, 3, 1, 2);
        C_Temp[t*12+ 6]=-Matrix_Product_T(&new_idm[0][0], &dGdF[ 1][0], 4, 3, 3, 2, 0);
        C_Temp[t*12+ 7]=-Matrix_Product_T(&new_idm[0][0], &dGdF[ 4][0], 4, 3, 3, 2, 1);
        C_Temp[t*12+ 8]=-Matrix_Product_T(&new_idm[0][0], &dGdF[ 7][0], 4, 3, 3, 2, 2);
        C_Temp[t*12+ 9]=-Matrix_Product_T(&new_idm[0][0], &dGdF[ 2][0], 4, 3, 3, 3, 0);
        C_Temp[t*12+10]=-Matrix_Product_T(&new_idm[0][0], &dGdF[ 5][0], 4, 3, 3, 3, 1);
        C_Temp[t*12+11]=-Matrix_Product_T(&new_idm[0][0], &dGdF[ 8][0], 4, 3, 3, 3, 2);
                
        //Get K matrix per tetrahedron
        /*Matrix_Product_T(&new_idm[0][0], &dGdF[ 0][0], &TK[t*144+ 3*12], 4, 3, 3);
        Matrix_Product_T(&new_idm[0][0], &dGdF[ 1][0], &TK[t*144+ 6*12], 4, 3, 3);
        Matrix_Product_T(&new_idm[0][0], &dGdF[ 2][0], &TK[t*144+ 9*12], 4, 3, 3);
        Matrix_Product_T(&new_idm[0][0], &dGdF[ 3][0], &TK[t*144+ 4*12], 4, 3, 3);
        Matrix_Product_T(&new_idm[0][0], &dGdF[ 4][0], &TK[t*144+ 7*12], 4, 3, 3);
        Matrix_Product_T(&new_idm[0][0], &dGdF[ 5][0], &TK[t*144+10*12], 4, 3, 3);
        Matrix_Product_T(&new_idm[0][0], &dGdF[ 6][0], &TK[t*144+ 5*12], 4, 3, 3);
        Matrix_Product_T(&new_idm[0][0], &dGdF[ 7][0], &TK[t*144+ 8*12], 4, 3, 3);
        Matrix_Product_T(&new_idm[0][0], &dGdF[ 8][0], &TK[t*144+11*12], 4, 3, 3);
        Matrix_Product_T(&new_idm[0][0], &dGdF[ 9][0], &TK[t*144+ 0*12], 4, 3, 3);
        Matrix_Product_T(&new_idm[0][0], &dGdF[10][0], &TK[t*144+ 1*12], 4, 3, 3);
        Matrix_Product_T(&new_idm[0][0], &dGdF[11][0], &TK[t*144+ 2*12], 4, 3, 3);
        
       
        TYPE G[9];
        TYPE x0[3]={1, 0, 0};
        TYPE x1[3]={0, 0, 0};
        TYPE x2[3]={0, 0, 0};
        TYPE x3[3]={0, 0, 0};
        Apply_K(t, x0, x1, x2, x3, G);
        printf("G: %f; %f\n", G[0], C_Temp[t*12+0]);
        
        printf("Center: %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n",
               C_Temp[t*12+ 0], C_Temp[t*12+ 1], C_Temp[t*12+ 2],
               C_Temp[t*12+ 3], C_Temp[t*12+ 4], C_Temp[t*12+ 5],
               C_Temp[t*12+ 6], C_Temp[t*12+ 7], C_Temp[t*12+ 8],
               C_Temp[t*12+ 9], C_Temp[t*12+10], C_Temp[t*12+11]);
		*/

		// Evaluate TK
		//printf("value: %f, %f\n", Numerical_Estimate_K(t, 5, 6), TK[t*144+5*12+6]);
		return energy;
    }

    
    void A_Times_Vector(TYPE* X, TYPE *Y, TYPE inv_t)
    {
        for(int i=0; i<number; i++)
        {
            TYPE c  = M[i]*inv_t*inv_t+fixed[i];
            Y[i*3+0]=c*X[i*3+0];
            Y[i*3+1]=c*X[i*3+1];
            Y[i*3+2]=c*X[i*3+2];
        }
		//printf("Y0: %f, %f, %f (%f, %f, %f)\n", Y[0], Y[1], Y[2], X[0], X[1], X[2]);

        for(int t=0; t<tet_number; t++)
        {
            //Recover alpha, beta, and gamma
            TYPE alpha[3][3], beta[3][3], gamma[3][3];
            memcpy(&alpha[0][0], &T_alpha[t*9], sizeof(TYPE)*9);
            memcpy(&beta [0][0], &T_beta [t*9], sizeof(TYPE)*9);
            memcpy(&gamma[0][0], &T_gamma[t*9], sizeof(TYPE)*9);
            
            
            int p0=Tet[t*4+0]*3;
            int p1=Tet[t*4+1]*3;
            int p2=Tet[t*4+2]*3;
            int p3=Tet[t*4+3]*3;
            
            TYPE Ds[9];
            Ds[0]=X[p1+0]-X[p0+0];
            Ds[3]=X[p1+1]-X[p0+1];
            Ds[6]=X[p1+2]-X[p0+2];
            Ds[1]=X[p2+0]-X[p0+0];
            Ds[4]=X[p2+1]-X[p0+1];
            Ds[7]=X[p2+2]-X[p0+2];
            Ds[2]=X[p3+0]-X[p0+0];
            Ds[5]=X[p3+1]-X[p0+1];
            Ds[8]=X[p3+2]-X[p0+2];
            TYPE F[9], temp0[9], temp1[9];
            Matrix_Product_3(Ds, &inv_Dm[t*9], F);
            Matrix_Product_3(F, &TV[t*9], temp0);
            Matrix_T_Product_3(&TU[t*9], temp0, temp1);
            
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

			if (p0 == 0 || p1 == 0 || p2 == 0 || p3 == 0)
			{
			//	printf("alpha: %f, %f, %f\n", alpha[0][0], alpha[0][1], alpha[0][2]);
			//	printf("U: %f, %f, %f\n", temp0[0], temp0[1], temp0[2]);
			//	printf("U: %f, %f, %f\n", temp0[3], temp0[4], temp0[5]);
			//	printf("U: %f, %f, %f\n", temp0[6], temp0[7], temp0[8]);
			}

            TYPE G[3][3];
            Matrix_Product_T_3(temp0, &TV[t*9], temp1);
            Matrix_Product_3(&TU[t*9], temp1, temp0);
            Matrix_Product_T_3(temp0, &inv_Dm[t*9], &G[0][0]);

			if(p0==0 || p1==0 || p2==0 || p3==0)
			{				
			//	printf("G: %f, %f, %f\n", G[0][0], G[0][1], G[0][2]);
			//	printf("G: %f, %f, %f\n", G[1][0], G[1][1], G[1][2]);
			//	printf("G: %f, %f, %f\n", G[2][0], G[2][1], G[2][2]);
			}
            
            Y[p0+0]+=G[0][0]+G[0][1]+G[0][2];
            Y[p0+1]+=G[1][0]+G[1][1]+G[1][2];
            Y[p0+2]+=G[2][0]+G[2][1]+G[2][2];
            Y[p1+0]-=G[0][0];
            Y[p1+1]-=G[1][0];
            Y[p1+2]-=G[2][0];
            Y[p2+0]-=G[0][1];
            Y[p2+1]-=G[1][1];
            Y[p2+2]-=G[2][1];
            Y[p3+0]-=G[0][2];
            Y[p3+1]-=G[1][2];
            Y[p3+2]-=G[2][2];
        }

		//printf("Y0: %f, %f, %f\n", Y[0], Y[1], Y[2]);
    }
    

    
    void A_Times(TYPE* X, TYPE *Y, TYPE inv_t)
    {

        for(int i=0; i<number; i++)
        {
            TYPE c  = M[i]*inv_t*inv_t+fixed[i];
            Y[i*3+0]=c*X[i*3+0];
            Y[i*3+1]=c*X[i*3+1];
            Y[i*3+2]=c*X[i*3+2];
        }
        
        for(int t=0; t<tet_number; t++)
        {
            int p0=Tet[t*4+0]*3;
            int p1=Tet[t*4+1]*3;
            int p2=Tet[t*4+2]*3;
            int p3=Tet[t*4+3]*3;
            
            TYPE x[12]={
                X[p0+0], X[p0+1], X[p0+2],
                X[p1+0], X[p1+1], X[p1+2],
                X[p2+0], X[p2+1], X[p2+2],
                X[p3+0], X[p3+1], X[p3+2]};
            TYPE y[12];
            
            Matrix_Product(&TK[t*144], x, y, 12, 12, 1);
            
            Y[p0+0]-=y[ 0];
            Y[p0+1]-=y[ 1];
            Y[p0+2]-=y[ 2];
            Y[p1+0]-=y[ 3];
            Y[p1+1]-=y[ 4];
            Y[p1+2]-=y[ 5];
            Y[p2+0]-=y[ 6];
            Y[p2+1]-=y[ 7];
            Y[p2+2]-=y[ 8];
            Y[p3+0]-=y[ 9];
            Y[p3+1]-=y[10];
            Y[p3+2]-=y[11];
        }
    }
    
    void M_Times(TYPE* X, TYPE *Y)
    {
        memset(Y, 0, sizeof(TYPE)*number*3);
        for(int t=0; t<tet_number; t++)
        {
            int p0=Tet[t*4+0]*3;
            int p1=Tet[t*4+1]*3;
            int p2=Tet[t*4+2]*3;
            int p3=Tet[t*4+3]*3;
            
            TYPE x[12]={
                X[p0+0], X[p0+1], X[p0+2],
                X[p1+0], X[p1+1], X[p1+2],
                X[p2+0], X[p2+1], X[p2+2],
                X[p3+0], X[p3+1], X[p3+2]};
            TYPE y[12];
            
            Matrix_Product(&TK[t*144], x, y, 12, 12, 1);
            
            Y[p0+0]+=y[ 0];
            Y[p0+1]+=y[ 1];
            Y[p0+2]+=y[ 2];
            Y[p1+0]+=y[ 3];
            Y[p1+1]+=y[ 4];
            Y[p1+2]+=y[ 5];
            Y[p2+0]+=y[ 6];
            Y[p2+1]+=y[ 7];
            Y[p2+2]+=y[ 8];
            Y[p3+0]+=y[ 9];
            Y[p3+1]+=y[10];
            Y[p3+2]+=y[11];
        }
    }
    
    void Compute_FM_for_Newton(int t, std::vector<Eigen::Triplet<TYPE>> *coefficients=0)
    {
        TYPE stiffness0=-Vol[t]*stiffness_0;
        TYPE stiffness1=-Vol[t]*stiffness_1;
        
        //No velocity access in this function
        int p0=Tet[t*4+0]*3;
        int p1=Tet[t*4+1]*3;
        int p2=Tet[t*4+2]*3;
        int p3=Tet[t*4+3]*3;        
        
        TYPE Ds[9];
        Ds[0]=X[p1+0]-X[p0+0];
        Ds[3]=X[p1+1]-X[p0+1];
        Ds[6]=X[p1+2]-X[p0+2];
        Ds[1]=X[p2+0]-X[p0+0];
        Ds[4]=X[p2+1]-X[p0+1];
        Ds[7]=X[p2+2]-X[p0+2];
        Ds[2]=X[p3+0]-X[p0+0];
        Ds[5]=X[p3+1]-X[p0+1];
        Ds[8]=X[p3+2]-X[p0+2];
        
        //printf("Ds: %f, %f, %f\n", Ds[0], Ds[1], Ds[2]);
        //printf("Ds: %f, %f, %f\n", Ds[3], Ds[4], Ds[5]);
        //printf("Ds: %f, %f, %f\n", Ds[6], Ds[7], Ds[8]);
        
        TYPE F[9], U[9], Sigma[3], V[9];
        Matrix_Product_3(Ds, &inv_Dm[t*9], F);

		memcpy(U, F, sizeof(TYPE)*9);
        SVD3((TYPE (*)[3])U, Sigma, (TYPE (*)[3])V);
        int small_id;
        if(fabsf(Sigma[0])<fabsf(Sigma[1]) && fabsf(Sigma[0])<fabsf(Sigma[2]))	small_id=0;
        else if(fabsf(Sigma[1])<fabsf(Sigma[2]))                                small_id=1;
        else                                                                    small_id=2;
        if(U[0]*(U[4]*U[8]-U[7]*U[5])+U[3]*(U[7]*U[2]-U[1]*U[8])+U[6]*(U[1]*U[5]-U[4]*U[2])<0)
        {
            U[0+small_id]	=-U[0+small_id];
            U[3+small_id]	=-U[3+small_id];
            U[6+small_id]	=-U[6+small_id];
			Sigma[small_id]	=-Sigma[small_id];
        }
        if(V[0]*(V[4]*V[8]-V[7]*V[5])+V[3]*(V[7]*V[2]-V[1]*V[8])+V[6]*(V[1]*V[5]-V[4]*V[2])<0)
        {
            V[0+small_id]=-V[0+small_id];
            V[3+small_id]=-V[3+small_id];
            V[6+small_id]=-V[6+small_id];
			Sigma[small_id] = -Sigma[small_id];
        }

/*		if(Sigma[0]<0.4)  Sigma[0]=0.4;
        if(Sigma[1]<0.4)  Sigma[1]=0.4;
        if(Sigma[2]<0.4)  Sigma[2]=0.4;        
        if(Sigma[0]>3)    Sigma[0]=3;
        if(Sigma[1]>3)    Sigma[1]=3;
        if(Sigma[2]>3)    Sigma[2]=3;


        {
            TYPE T0[3][3], C[3][3], B[3][3];
            memset(&C[0][0], 0, sizeof(TYPE)*9);
            C[0][0]=Sigma[0];
            C[1][1]=Sigma[1];
            C[2][2]=Sigma[2];
            Matrix_Product_3(U, &C[0][0], &B[0][0]);
            Matrix_Product_T_3(&B[0][0], V, &T0[0][0]);   
            
            printf("F: %f, %f, %f; %f, %f, %f; %f, %f, %f\n",
                   F[0], F[1], F[2],
                   F[3], F[4], F[5],
                   F[6], F[7], F[8]);
            
            printf("T: %f, %f, %f; %f, %f, %f; %f, %f, %f\n",
                   T0[0][0], T0[0][1], T0[0][2],
                   T0[1][0], T0[1][1], T0[1][2],
                   T0[2][0], T0[2][1], T0[2][2]);
            
        }*/        
        

        TYPE I=Sigma[0]*Sigma[0]+Sigma[1]*Sigma[1]+Sigma[2]*Sigma[2];
        TYPE J=Sigma[0]*Sigma[1]*Sigma[2];
        TYPE III=J*J;
        
        TYPE var0=2*stiffness0/cbrt(III);
        TYPE var1=2*stiffness1*(J-1);
        
        TYPE P0=var0*(Sigma[0]-I/(3*Sigma[0]))+var1*Sigma[1]*Sigma[2];
        TYPE P1=var0*(Sigma[1]-I/(3*Sigma[1]))+var1*Sigma[0]*Sigma[2];
        TYPE P2=var0*(Sigma[2]-I/(3*Sigma[2]))+var1*Sigma[0]*Sigma[1];
        
        TYPE PV_transpose[9], P[9], force[9];
        PV_transpose[0]=P0*V[0];
        PV_transpose[1]=P0*V[3];
        PV_transpose[2]=P0*V[6];
        PV_transpose[3]=P1*V[1];
        PV_transpose[4]=P1*V[4];
        PV_transpose[5]=P1*V[7];
        PV_transpose[6]=P2*V[2];
        PV_transpose[7]=P2*V[5];
        PV_transpose[8]=P2*V[8];
        Matrix_Product_3(U, PV_transpose, P);
        Matrix_Product_T_3(P, &inv_Dm[t*9], force);
        Force[p1+0]+=force[0];
        Force[p1+1]+=force[3];
        Force[p1+2]+=force[6];
        Force[p2+0]+=force[1];
        Force[p2+1]+=force[4];
        Force[p2+2]+=force[7];
        Force[p3+0]+=force[2];
        Force[p3+1]+=force[5];
        Force[p3+2]+=force[8];
        Force[p0+0]-=(force[0]+force[1]+force[2]);
        Force[p0+1]-=(force[3]+force[4]+force[5]);
        Force[p0+2]-=(force[6]+force[7]+force[8]);
        
        
        //Now compute the stiffness matrix
        TYPE dEdI   = var0*0.5;
        TYPE dEdII  = 0;
        TYPE dEdIII =-var0*0.5*I/(3*III)+var1*0.5/J;
        TYPE H[3][3];
        memset(&H[0][0], 0, sizeof(TYPE)*9);
        
        H[0][2]=H[2][0]=-var0*0.5/(3*III);
        H[2][2]=var0*I*2.0/(9.0*III*III)+stiffness1*0.5/(J*J*J);
        
        TYPE alpha[3][3], beta[3][3], gamma[3][3];
        for(int i=0; i<3; i++)
        for(int j=i; j<3; j++)
        {
            alpha[i][j]=2*dEdI+4*(Sigma[i]*Sigma[i]+Sigma[j]*Sigma[j])*dEdII;
             beta[i][j]=4*Sigma[i]*Sigma[j]*dEdII-2*III*dEdIII/(Sigma[i]*Sigma[j]);
            
            TYPE vi[3]={2*Sigma[i], 4*Sigma[i]*Sigma[i]*Sigma[i], 2*III/Sigma[i]};
            TYPE vj[3]={2*Sigma[j], 4*Sigma[j]*Sigma[j]*Sigma[j], 2*III/Sigma[j]};
            TYPE r[3];            
            Matrix_Vector_Product_3(&H[0][0], vj, r);
            gamma[i][j]=DOT(vi, r)+4*III*dEdIII/(Sigma[i]*Sigma[j]);
        }

		TYPE dGdF[12][9];    //G is related to force, according to [TSIF05], (g0, g3, g6), (g1, g4, g7), (g2, g5, g8), (g9, g10, g11)
        for(int i=0; i<9; i++)
        {
            TYPE dF[9], temp0[9], temp1[9];
            memset(&dF, 0, sizeof(TYPE)*9);
            dF[i]=1;
            
            Matrix_Product_3(dF, V, temp0);
            Matrix_T_Product_3(U, temp0, temp1);
            
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
            
            Matrix_Product_T_3(temp0, V, temp1);
            Matrix_Product_3(U, temp1, temp0);
            Matrix_Product_T_3(temp0, &inv_Dm[t*9], &dGdF[i][0]);
        }
        
        //for(int i=0; i<9; i++) for(int j=0; j<9; j++)
        //   dGdF[i][j]*=-Vol[t];
        
        //Transpose dGdF
        TYPE temp;
        for(int i=0; i<9; i++) for(int j=i+1; j<9; j++)
            SWAP(dGdF[i][j], dGdF[j][i]);
                
        for(int j=0; j< 9; j++)
        {
            dGdF[ 9][j]=-dGdF[0][j]-dGdF[1][j]-dGdF[2][j];
            dGdF[10][j]=-dGdF[3][j]-dGdF[4][j]-dGdF[5][j];
            dGdF[11][j]=-dGdF[6][j]-dGdF[7][j]-dGdF[8][j];
        }
            
        TYPE new_idm[4][3];
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
        
       
        
        TYPE r[12];
        for(int i=0; i<12; i++)
        {
            Matrix_Product_T(&new_idm[0][0], &dGdF[i][0], r, 4, 3, 3);            
            //printf("r %d: %f, %f, %f; %f, %f, %f; %f, %f, %f; %f, %f, %f\n", i,
            //       r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10], r[11]);
            for(int j=0; j<12; j++)
            {
                int pi, pj;                
                     if(i==0)   pi=p1;
                else if(i==1)   pi=p2;
                else if(i==2)   pi=p3;
                else if(i==3)   pi=p1+1;
                else if(i==4)   pi=p2+1;
                else if(i==5)   pi=p3+1;
                else if(i==6)   pi=p1+2;
                else if(i==7)   pi=p2+2;
                else if(i==8)   pi=p3+2;
                else if(i==9)   pi=p0;
                else if(i==10)  pi=p0+1;
                else if(i==11)  pi=p0+2;
                
                     if(j== 0)  pj=p0+0;
                else if(j== 1)  pj=p0+1;
                else if(j== 2)  pj=p0+2;
                else if(j== 3)  pj=p1+0;
                else if(j== 4)  pj=p1+1;
                else if(j== 5)  pj=p1+2;
                else if(j== 6)  pj=p2+0;
                else if(j== 7)  pj=p2+1;
                else if(j== 8)  pj=p2+2;
                else if(j== 9)  pj=p3+0;
                else if(j==10)  pj=p3+1;
                else if(j==11)  pj=p3+2;                                              
             
                if(coefficients)    coefficients->push_back(Eigen::Triplet<TYPE>(pi, pj, -r[j]));

				//if(pi==pj)
                //{
                //    C[pi]+=-r[j];
                //}
                //M[i][j]=-r[j];
            }
        }

		TYPE K[144];
		Matrix_Product_T(&new_idm[0][0], &dGdF[ 0][0], &K[ 3*12], 4, 3, 3);
        Matrix_Product_T(&new_idm[0][0], &dGdF[ 1][0], &K[ 6*12], 4, 3, 3);
        Matrix_Product_T(&new_idm[0][0], &dGdF[ 2][0], &K[ 9*12], 4, 3, 3);
        Matrix_Product_T(&new_idm[0][0], &dGdF[ 3][0], &K[ 4*12], 4, 3, 3);
        Matrix_Product_T(&new_idm[0][0], &dGdF[ 4][0], &K[ 7*12], 4, 3, 3);
        Matrix_Product_T(&new_idm[0][0], &dGdF[ 5][0], &K[10*12], 4, 3, 3);
        Matrix_Product_T(&new_idm[0][0], &dGdF[ 6][0], &K[ 5*12], 4, 3, 3);
        Matrix_Product_T(&new_idm[0][0], &dGdF[ 7][0], &K[ 8*12], 4, 3, 3);
        Matrix_Product_T(&new_idm[0][0], &dGdF[ 8][0], &K[11*12], 4, 3, 3);
        Matrix_Product_T(&new_idm[0][0], &dGdF[ 9][0], &K[ 0*12], 4, 3, 3);
        Matrix_Product_T(&new_idm[0][0], &dGdF[10][0], &K[ 1*12], 4, 3, 3);
        Matrix_Product_T(&new_idm[0][0], &dGdF[11][0], &K[ 2*12], 4, 3, 3);        

		//printf("K: %f, %f, %f\n", K[0*12+0], K[1*12+1], K[2*12+2]);
		//printf("K: %f, %f, %f\n", K[3*12+3], K[4*12+4], K[5*12+5]);
		//printf("K: %f, %f, %f\n", K[6*12+6], K[7*12+7], K[8*12+8]);
		//printf("K: %f, %f, %f\n", K[9*12+9], K[10*12+10], K[11*12+11]);


		TYPE KK[4][4][3][3];
		for(int ii=0; ii<4; ii++) 
		for(int jj=0; jj<4; jj++)
		for(int  i=0;  i<3;  i++)
		for(int  j=0;  j<3;  j++)
		{
			KK[ii][jj][i][j]=-K[(ii*3+i)*12+jj*3+j];
		}

		// Set vert
		pcg_solver.Set_M_vert(Tet[t*4+0], &KK[0][0][0][0]);
		pcg_solver.Set_M_vert(Tet[t*4+1], &KK[1][1][0][0]);
		pcg_solver.Set_M_vert(Tet[t*4+2], &KK[2][2][0][0]);
		pcg_solver.Set_M_vert(Tet[t*4+3], &KK[3][3][0][0]);

		//Set edge

		if(p0<p1)	pcg_solver.Set_M_edge(TE[t*6+0], &KK[0][1][0][0]);
		else		pcg_solver.Set_M_edge(TE[t*6+0], &KK[1][0][0][0]);

		if(p0<p2)	pcg_solver.Set_M_edge(TE[t*6+1], &KK[0][2][0][0]);
		else		pcg_solver.Set_M_edge(TE[t*6+1], &KK[2][0][0][0]);

		if(p0<p3)	pcg_solver.Set_M_edge(TE[t*6+2], &KK[0][3][0][0]);
		else		pcg_solver.Set_M_edge(TE[t*6+2], &KK[3][0][0][0]);

		if(p1<p2)	pcg_solver.Set_M_edge(TE[t*6+3], &KK[1][2][0][0]);
		else		pcg_solver.Set_M_edge(TE[t*6+3], &KK[2][1][0][0]);

		if(p1<p3)	pcg_solver.Set_M_edge(TE[t*6+4], &KK[1][3][0][0]);
		else		pcg_solver.Set_M_edge(TE[t*6+4], &KK[3][1][0][0]);

		if(p2<p3)	pcg_solver.Set_M_edge(TE[t*6+5], &KK[2][3][0][0]);
		else		pcg_solver.Set_M_edge(TE[t*6+5], &KK[3][2][0][0]);
    }

    void Conjugate_Gradient_Solver(TYPE *X, TYPE* old_X, TYPE *B, TYPE t)
	{
		int number3=number*3;
		pcg_solver.Set_Preconditioner();
		

		pcg_solver.M_Times_Vector(X, pcg_solver.temp);
		for(int i=0; i<number3; i++)
			pcg_solver.pcg_R[i]=B[i]-pcg_solver.temp[i];

		pcg_solver.Solve_Preconditioner_L(pcg_solver.pcg_R, pcg_solver.temp);
		pcg_solver.Solve_Preconditioner_U(pcg_solver.temp, pcg_solver.pcg_Z);

		//p=z
		memcpy(pcg_solver.pcg_P, pcg_solver.pcg_Z, sizeof(TYPE)*number3);

		//Set rz
		TYPE rz, new_rz, alpha, beta;
		rz=Dot(pcg_solver.pcg_R, pcg_solver.pcg_Z, number3);

		for(int tet=0; tet<tet_number; tet++)
			Compute_FM(tet, stiffness_0, stiffness_1, stiffness_2, stiffness_3, true);

		int l;
		for(l=0; l<97; l++)
		{
			//printf("\n");
			//printf("%d: %f\n", l, Norm(pcg_solver.pcg_R, number3));

			TYPE error_mag1=0;
            A_Times(X, pcg_solver.pcg_AP, 1/t);
			for(int i=0; i<number*3; i++)
                error_mag1+=-X[i]*pcg_solver.pcg_AP[i]-2*pcg_solver.pcg_R[i]*X[i];

			for(int i=0; i<number3; i++)
				next_X[i]=old_X[i]+X[i];
			printf("%f\n", Get_Energy(next_X, t), error_mag1, Norm(pcg_solver.pcg_R, number3));

			if(Norm(pcg_solver.pcg_R, number3)<1e-5f)	break;
			pcg_solver.M_Times_Vector(pcg_solver.pcg_P, pcg_solver.pcg_AP);

			//printf("PAP: %f, %f\n", pcg_P[10], pcg_AP[10]);

			alpha=rz/Dot(pcg_solver.pcg_P, pcg_solver.pcg_AP, number3);
			//printf("alphaaaaa %f, %f: %f\n", rz, Dot(pcg_P, pcg_AP, number3), alpha);
			for(int i=0; i<number3; i++)
			{
				X[i]    +=alpha*pcg_solver.pcg_P[i];
				pcg_solver.pcg_R[i]-=alpha*pcg_solver.pcg_AP[i];
			}

			//printf("R: %f\n", pcg_R[10]);

			pcg_solver.Solve_Preconditioner_L(pcg_solver.pcg_R, pcg_solver.temp);
			pcg_solver.Solve_Preconditioner_U(pcg_solver.temp, pcg_solver.pcg_Z);
			new_rz=Dot(pcg_solver.pcg_R, pcg_solver.pcg_Z, number3);

			//printf("rz: %f, %f\n", rz, new_rz);

			beta=new_rz/rz;
			rz=new_rz;

			for(int i=0; i<number3; i++)
				pcg_solver.pcg_P[i]=pcg_solver.pcg_Z[i]+beta*pcg_solver.pcg_P[i];
			//printf("beta: %f\n", beta);
		}
	}


	void Newton_Update(TYPE t, int iterations)
	{
        //Damp the system first
        for(int i=0; i<number; i++)
        {
            V[i*3+0]*=0.999;
            V[i*3+1]*=0.999;
            V[i*3+2]*=0.999;
        }

		//Calculate the expected position: S
		for(int i=0; i<number; i++)
        {
            S[i*3+0]=X[i*3+0]+V[i*3+0]*t;
            S[i*3+1]=X[i*3+1]+V[i*3+1]*t;
            S[i*3+2]=X[i*3+2]+V[i*3+2]*t;

			X[i * 3 + 0] += V[i * 3 + 0]*t;// + (V[i * 3 + 0] - prev_V[i * 3 + 0])*0*t;
			X[i * 3 + 1] += V[i * 3 + 1]*t;// + (V[i * 3 + 1] - prev_V[i * 3 + 1])*0*t;
			X[i * 3 + 2] += V[i * 3 + 2]*t;// + (V[i * 3 + 2] - prev_V[i * 3 + 2])*0*t;
        }

		std::vector<Eigen::Triplet<TYPE>>						coefficients;
	//	Eigen::ConjugateGradient<Eigen::SparseMatrix<TYPE>, Eigen::Lower|Eigen::Upper, Eigen::IncompleteLUT<TYPE> >    solver;
	//	Eigen::ConjugateGradient<Eigen::SparseMatrix<TYPE> >    solver;
	//	Eigen::PardisoLDLT<Eigen::SparseMatrix<TYPE> >	solver;
		
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<TYPE> >	solver;

		Eigen::SparseMatrix<TYPE>								matrix;
		matrix.resize(number*3, number*3);
		Eigen::VectorXd b(number*3), x(number*3);

		for(int l=0; l<iterations; l++)
		{	
			memset(Force, 0, sizeof(TYPE)*number*3);
			coefficients.clear();
			coefficients.reserve(10000000);



			for(int t=0; t<tet_number; t++)			
				Compute_FM_for_Newton(t, &coefficients);
	

			//Add mass diagonal
			for(int i=0; i<number; i++)
			{
				TYPE c=M[i]/(t*t)+fixed[i];
				coefficients.push_back(Eigen::Triplet<TYPE>(i*3+0, i*3+0, c));
				coefficients.push_back(Eigen::Triplet<TYPE>(i*3+1, i*3+1, c));
				coefficients.push_back(Eigen::Triplet<TYPE>(i*3+2, i*3+2, c));

				TYPE m[9];
				m[0]=c; m[1]=0; m[2]=0;
				m[3]=0; m[4]=c; m[5]=0;
				m[6]=0; m[7]=0; m[8]=c;
				pcg_solver.Set_M_vert(i, m);

			}
			matrix.setFromTriplets(coefficients.begin(), coefficients.end());        

			TYPE error_sum=0;
			for(int i=0; i<number; i++)
			{
				TYPE oc=M[i]/(t*t);
				b[i*3+0]=oc*(S[i*3+0]-X[i*3+0])+fixed[i]*(fixed_X[i*3+0]-X[i*3+0])+Force[i*3+0];
                b[i*3+1]=oc*(S[i*3+1]-X[i*3+1])+fixed[i]*(fixed_X[i*3+1]-X[i*3+1])+Force[i*3+1]+M[i]*gravity;
                b[i*3+2]=oc*(S[i*3+2]-X[i*3+2])+fixed[i]*(fixed_X[i*3+2]-X[i*3+2])+Force[i*3+2];
				
				B[i*3+0]=b[i*3+0];
				B[i*3+1]=b[i*3+1];
				B[i*3+2]=b[i*3+2];

				error_sum+=b[i*3+0]*b[i*3+0]+b[i*3+1]*b[i*3+1]+b[i*3+2]*b[i*3+2];
			}
			//printf("%f, %f\n", error_sum, timer.Get_Time());
			
		//	printf("s: %d\n", coefficients.size());
			
		//	TYPE* XX=new TYPE[number*3];
		//	memset(XX, 0, sizeof(TYPE)*number*3);
		//	TIMER timer2;
		//	Conjugate_Gradient_Solver(XX, X, B, t);
		//	printf("timer2: %f\n", timer2.Get_Time());
		//	for(int i=0; i<number*3; i++)
		//		printf("solution %d: %f\n", i, XX[i]);
		//	delete[] XX;
		///	getchar();

		//	solver.setMaxIterations(20);

			solver.analyzePattern(matrix);

			TIMER timer;
			solver.factorize(matrix);
			x = solver.solve(b);
			//printf("time here: %f\n", timer.Get_Time());
			printf("%f, %f, (%f)\n", Get_Energy(X, t), error_sum, timer.Get_Time());
	       // printf("end error: %f; %f\n", Get_Error(t), Get_Energy(X, t));
        
			for(int i=0; i<number*3; i++)			
				X[i]+=x[i]*0.5;


		}

		//printf("it: %d\n", solver.iterations());
        printf("%f\n", Get_Energy(X, t));
		TYPE inv_t=1/t;
		for(int i=0; i<number; i++)
        {
			prev_V[i*3+0]=V[i*3+0];
			prev_V[i*3+1]=V[i*3+1];
			prev_V[i*3+2]=V[i*3+2];

            V[i*3+0]+=(X[i*3+0]-S[i*3+0])*inv_t;
            V[i*3+1]+=(X[i*3+1]-S[i*3+1])*inv_t;
            V[i*3+2]+=(X[i*3+2]-S[i*3+2])*inv_t;
        }


		//printf("time: %f\n", timer.Get_Time());
        //TYPE sum=0;
		//for(int i=0; i<number*3; i++)   sum+=X[i];
		//printf("sum: %f\n", sum);        
        //printf("V: %f, %f, %f\n", V[0]/t, V[1]/t, V[2]/t);
	}
    

    void CG_Update(TYPE t, int iterations)
    {
        //memset(lambda, 0, sizeof(TYPE)*tet_number);

        //Damp the system first
        for(int i=0; i<number; i++) //if(fixed[i]==0)
        {
            V[i*3+0]*=0.999;
            V[i*3+1]*=0.999;
            V[i*3+2]*=0.999;
        }
        
        //Calculate the expected position: S
        for(int i=0; i<number; i++)
        {
            S[i*3+0]=X[i*3+0]+V[i*3+0]*t;
            S[i*3+1]=X[i*3+1]+V[i*3+1]*t;
            S[i*3+2]=X[i*3+2]+V[i*3+2]*t;
            
            //X[i*3+0]+=V[i*3+0]*t;
            //X[i*3+1]+=V[i*3+1]*t;
            //X[i*3+2]+=V[i*3+2]*t;

			X[i * 3 + 0] += (V[i * 3 + 0] + (V[i * 3 + 0] - prev_V[i * 3 + 0])*0)*t;
			X[i * 3 + 1] += (V[i * 3 + 1] + (V[i * 3 + 1] - prev_V[i * 3 + 1])*0)*t;
			X[i * 3 + 2] += (V[i * 3 + 2] + (V[i * 3 + 2] - prev_V[i * 3 + 2])*0)*t;
        }        

        
        TYPE inv_t = 1/t;
        TYPE zr = 0;
        
   /*    for(int tet=0; tet<tet_number; tet++)
				Compute_FM(tet, stiffness_0, stiffness_1, stiffness_2, stiffness_3, true);

        
        for(int i=0; i<number; i++)
        {
            
            TYPE c=M[i]*inv_t*inv_t+fixed[i];
            C[i*3+0]=c;
            C[i*3+1]=c;
            C[i*3+2]=c;
            
            TYPE force[3]={0, 0, 0};
            for(int index=vtt_num[i]; index<vtt_num[i+1]; index++)
            {
                force[0]+=F_Temp[VTT[index]*3+0];
                force[1]+=F_Temp[VTT[index]*3+1];
                force[2]+=F_Temp[VTT[index]*3+2];
                
                C[i*3+0]+=C_Temp[VTT[index]*3+0];
                C[i*3+1]+=C_Temp[VTT[index]*3+1];
                C[i*3+2]+=C_Temp[VTT[index]*3+2];
            }
            
            TYPE oc = M[i]*inv_t*inv_t;
            B[i*3+0]=oc*(S[i*3+0])+fixed[i]*(fixed_X[i*3+0])+force[0];
            B[i*3+1]=oc*(S[i*3+1])+fixed[i]*(fixed_X[i*3+1])+force[1]+M[i]*gravity;
            B[i*3+2]=oc*(S[i*3+2])+fixed[i]*(fixed_X[i*3+2])+force[2];
            
            //printf("%f, %f, %f, ", B[i*3+0], B[i*3+1], B[i*3+2]);
            //printf("f: %f, %f, %f\n", force[0], force[1], force[2]);
            //printf("sx: %f, %f, %f; %f, %f, %f\n", S[i*3+0]+V[i*3+0]*t, S[i*3+1]+V[i*3+1]*t, S[i*3+2]+V[i*3+2]*t,
            //       X[i*3+0], X[i*3+1], X[i*3+2]);
        }*/

        //for(int i=0; i<number*3; i++)
        //    X[i]=RandomFloat()*2-1;
        
       /* printf("\n");
        printf("x is:\n");
        for(int i=0; i<number*3; i++)
            printf("%f, ", X[i]);
        printf("\n");*/
        //            printf("%f\n", Get_Error(t));
		//			getchar();



        for(int l=0; l<iterations; l++)
        {
			memset(ext_C, 0, sizeof(TYPE)*number);
			memset(lambda, 0, sizeof(TYPE)*tet_number);

            //Get r
            for(int tet=0; tet<tet_number; tet++)
				Compute_FM(tet, stiffness_0, stiffness_1, stiffness_2, stiffness_3, true);
            
            for(int i=0; i<number; i++)
            {
                TYPE c=M[i]*inv_t*inv_t+fixed[i];
                C[i*3+0]=c;
                C[i*3+1]=c;
                C[i*3+2]=c;
                
                TYPE force[3]={0, 0, 0};
                for(int index=vtt_num[i]; index<vtt_num[i+1]; index++)
                {
                    force[0]+=F_Temp[VTT[index]*3+0];
                    force[1]+=F_Temp[VTT[index]*3+1];
                    force[2]+=F_Temp[VTT[index]*3+2];
                    
                    C[i*3+0]+=C_Temp[VTT[index]*3+0];
                    C[i*3+1]+=C_Temp[VTT[index]*3+1];
                    C[i*3+2]+=C_Temp[VTT[index]*3+2];
                }
                
                TYPE oc = M[i]*inv_t*inv_t;
                cg_r[i*3+0]=oc*(S[i*3+0]-X[i*3+0])+fixed[i]*(fixed_X[i*3+0]-X[i*3+0])+force[0];
                cg_r[i*3+1]=oc*(S[i*3+1]-X[i*3+1])+fixed[i]*(fixed_X[i*3+1]-X[i*3+1])+force[1]+M[i]*gravity;
                cg_r[i*3+2]=oc*(S[i*3+2]-X[i*3+2])+fixed[i]*(fixed_X[i*3+2]-X[i*3+2])+force[2];
            }
            
            
			//printf("error: %f, %f\n", Get_Error(t), Get_Energy(X, t));
			//printf("cgr: %f\n", cg_r[10]);
            //A_Times(X, cg_r, inv_t);
            //TYPE diff_error=0;
            //for(int i=0; i<number*3; i++)
            //    diff_error+=X[i]*cg_r[i]-2*B[i]*X[i];
            //for(int i=0; i<number*3; i++)   cg_r[i]=B[i]-cg_r[i];
            
            //for(int i=0; i<number*3; i++)
            //    printf("r %d: %f (%f)\n", i, cg_r[i], X[i]);
            
            for(int k=0; k<1; k++)
            {
            TYPE error_mag=0;
            for(int i=0; i<number; i++)
                error_mag+=cg_r[i*3+0]*cg_r[i*3+0]+cg_r[i*3+1]*cg_r[i*3+1]+cg_r[i*3+2]*cg_r[i*3+2];
            
			TYPE error_mag1=0;
            A_Times(X, cg_Ap, inv_t);
			for(int i=0; i<number*3; i++)
                error_mag1+=-X[i]*cg_Ap[i]-2*cg_r[i]*X[i];
            
			printf("%f\n", Get_Energy(X, t), error_mag);
			//printf("%f\n", error_mag1);
            //getchar();

            //Precondition solve
            for(int i=0; i<number*3; i++)   cg_z[i]=cg_r[i];///(C[i]+ext_C[i/3]); //note here here
            
            //Update p
            TYPE old_zr=zr;
            zr=0;
            for(int i=0; i<number*3; i++)   zr+=cg_z[i]*cg_r[i];
            TYPE beta=0;
			if(l!=0)    
			//if(k!=0)	
				beta=zr/old_zr;			
            
		//	printf("beta: %f (%f, %f)\n", beta, old_zr, zr);

            for(int i=0; i<number*3; i++)   cg_p[i]=cg_z[i]+beta*cg_p[i];
            
            //printf("beta: %f (%f, %f)\n", beta, zr, old_zr);
            
            //for(int i=0; i<number*3; i++)
            //    printf("P: %d, %f\n", i, cg_p[i]);
            
            //get alpha
            A_Times(cg_p, cg_Ap, inv_t);

            //for(int i=0; i<number*3; i++)
            //    printf("P: %d, %f, %f\n", i, cg_p[i], cg_Ap[i]);
		//	printf("PAP: %f, %f\n", cg_p[10], cg_Ap[10]);

            
            TYPE alpha=0;
            for(int i=0; i<number*3; i++)   alpha+=cg_p[i]*cg_Ap[i];    //alpha=pAp
            
		//	printf("alpha: %f\n", alpha);

            //for(int i=0; i<number*3; i++)
            //    printf("value: %f, %f, %f, %f\n", cg_z[i], cg_r[i], cg_p[i], cg_Ap[i]);
            
            //printf("value: %f, %f\n", zr, alpha);

            alpha=zr/alpha;                                             //alpha=zr/pAp
            
            //printf("alpha: %ef\n", alpha);
            
            //Update X
            for(int i=0; i<number*3; i++)
            {
                X[i]=X[i]+alpha*cg_p[i];
                cg_r[i]=cg_r[i]-alpha*cg_Ap[i];
            }

            //for(int i=0; i<number*3; i++)   X[i]=X[i]+alpha*cg_p[i];
		//	printf("end: %f\n", cg_r[10]);
            
            //getchar();
            }
        }
        
        for(int i=0; i<number; i++)
        {
            V[i*3+0]+=(X[i*3+0]-S[i*3+0])*inv_t;
            V[i*3+1]+=(X[i*3+1]-S[i*3+1])*inv_t;
            V[i*3+2]+=(X[i*3+2]-S[i*3+2])*inv_t;
        }
    }

	void Single_Newton_Update(TYPE t, int iterations)
    {   
		//t=0.0005;
        //Damp the system first
        for(int i=0; i<number; i++) //if(fixed[i]==0)
        {
            V[i*3+0]*=0.999;
            V[i*3+1]*=0.999;
            V[i*3+2]*=0.999;
        }
        
        //Calculate the expected position: S
        for(int i=0; i<number; i++)
        {
            S[i*3+0]=X[i*3+0]+V[i*3+0]*t;
            S[i*3+1]=X[i*3+1]+V[i*3+1]*t;
            S[i*3+2]=X[i*3+2]+V[i*3+2]*t;
            
			X[i * 3 + 0] += (V[i * 3 + 0] + (V[i * 3 + 0] - prev_V[i * 3 + 0])*0)*t;
			X[i * 3 + 1] += (V[i * 3 + 1] + (V[i * 3 + 1] - prev_V[i * 3 + 1])*0)*t;
			X[i * 3 + 2] += (V[i * 3 + 2] + (V[i * 3 + 2] - prev_V[i * 3 + 2])*0)*t;
        }     

        TYPE inv_t = 1/t;

		//Get r
        for(int tet=0; tet<tet_number; tet++)
			Compute_FM(tet, stiffness_0, stiffness_1, stiffness_2, stiffness_3, true);
            
        for(int i=0; i<number; i++)
        {
			TYPE c=M[i]*inv_t*inv_t+fixed[i];
            C[i*3+0]=c;
			C[i*3+1]=c;
			C[i*3+2]=c;
                
			TYPE force[3]={0, 0, 0};
			for(int index=vtt_num[i]; index<vtt_num[i+1]; index++)
			{
				force[0]+=F_Temp[VTT[index]*3+0];
				force[1]+=F_Temp[VTT[index]*3+1];
				force[2]+=F_Temp[VTT[index]*3+2];                    
				C[i*3+0]+=C_Temp[VTT[index]*3+0];
				C[i*3+1]+=C_Temp[VTT[index]*3+1];
				C[i*3+2]+=C_Temp[VTT[index]*3+2];
			}
                
			TYPE oc = M[i]*inv_t*inv_t;
			cg_r[i*3+0]=oc*(S[i*3+0]-X[i*3+0])+fixed[i]*(fixed_X[i*3+0]-X[i*3+0])+force[0];
			cg_r[i*3+1]=oc*(S[i*3+1]-X[i*3+1])+fixed[i]*(fixed_X[i*3+1]-X[i*3+1])+force[1]+M[i]*gravity;
			cg_r[i*3+2]=oc*(S[i*3+2]-X[i*3+2])+fixed[i]*(fixed_X[i*3+2]-X[i*3+2])+force[2];
		}


		TYPE stepping=0.001;
		TYPE omega, rho;

		for(int l=0; l<iterations; l++)
		{
			TYPE error=0;
			for(int i=0; i<number*3; i++)	error+=cg_r[i]*cg_r[i];
			printf("%d: %f\n", l, Get_Energy(X, t));

			for(int i=0; i<number*3; i++)	cg_p[i]=cg_r[i]/C[i];

			for(int i=0; i<number*3; i++)	next_X[i]=X[i]+cg_p[i]*stepping;
			
	/*		if (l == 0)				omega=1;
			if (l == 1)				{ rho = 0.96;	omega = 2 / (2 - rho*rho);		}
			if (l > 1 && l < 7)		{ rho = 0.96;	omega = 4 / (4 - rho*rho*omega);}
			if(l==6)				omega=1;
			if (l == 7)				{ rho = 0.99;	omega = 2 / (2 - rho*rho);		}
			if (l > 7 && l < 12)	{ rho = 0.99;	omega = 4 / (4 - rho*rho*omega);}
			if (l == 11)			omega=1;
			if (l == 12)			{ rho = 0.999;  omega = 2 / (2 - rho*rho);		}	//0.9992
			if (l > 12 && l<20)		{ rho = 0.999;  omega = 4 / (4 - rho*rho*omega);}	//0.9992
			if (l == 20)			omega=1;
			if (l == 21)			{ rho = 0.99996;  omega = 2 / (2 - rho*rho);	}
			if (l > 21)				{ rho = 0.99996;  omega = 4 / (4 - rho*rho*omega);}	//0.9999
*/

			if (l == 0)				omega=1;
			if (l == 1)				{ rho = 0.96;	omega = 2 / (2 - rho*rho);		}
			if (l > 1 && l < 7)		{ rho = 0.96;	omega = 4 / (4 - rho*rho*omega);}
			if(l==6)				omega=1;
			if (l == 7)				{ rho = 0.99;	omega = 2 / (2 - rho*rho);		}
			if (l > 7 && l < 12)	{ rho = 0.99;	omega = 4 / (4 - rho*rho*omega);}
			if (l == 11)			omega=1;
			if (l == 12)			{ rho = 0.991;  omega = 2 / (2 - rho*rho);		}	
			if (l > 12 && l<20)		{ rho = 0.991;  omega = 4 / (4 - rho*rho*omega);}	
			if (l == 20)			omega=1;
			if (l == 21)			{ rho = 0.99;  omega = 2 / (2 - rho*rho);		}
			if (l > 21)				{ rho = 0.99;  omega = 4 / (4 - rho*rho*omega);}	//0.99996

			for(int i=0; i<number*3; i++)
            {
                next_X[i]=omega*(next_X[i]-prev_X[i])+prev_X[i];
                prev_X[i]=X[i];
                X[i]=next_X[i];
				cg_p[i]=X[i]-prev_X[i];
            }	

			//Update r
			A_Times(cg_p, cg_Ap, inv_t);
			for(int i=0; i<number*3; i++)	cg_r[i]-=cg_Ap[i];
		}

		for(int i=0; i<number; i++)
        {
            V[i*3+0]=(X[i*3+0]-S[i*3+0])*inv_t;
            V[i*3+1]=(X[i*3+1]-S[i*3+1])*inv_t;
            V[i*3+2]=(X[i*3+2]-S[i*3+2])*inv_t;
        }

	}


	TYPE Get_Error(TYPE t)
	{
		TYPE error_sum=0;
		for (int tet = 0; tet < tet_number; tet++)
			Compute_FM(tet, stiffness_0, stiffness_1, stiffness_2, stiffness_3, false);
		for (int i = 0; i < number; i++)
		{
			TYPE oc = M[i] / (t*t);
			TYPE c = oc + fixed[i];

			TYPE force[3] = { 0, 0, 0 };
			for (int index = vtt_num[i]; index < vtt_num[i + 1]; index++)
			{
				force[0] += F_Temp[VTT[index] * 3 + 0];
				force[1] += F_Temp[VTT[index] * 3 + 1];
				force[2] += F_Temp[VTT[index] * 3 + 2];
			}
			B[i * 3 + 0] = oc*(S[i * 3 + 0] - X[i * 3 + 0]) + fixed[i] * (fixed_X[i * 3 + 0] - X[i * 3 + 0]) + force[0];
			B[i * 3 + 1] = oc*(S[i * 3 + 1] - X[i * 3 + 1]) + fixed[i] * (fixed_X[i * 3 + 1] - X[i * 3 + 1]) + force[1] + M[i] * gravity;
			B[i * 3 + 2] = oc*(S[i * 3 + 2] - X[i * 3 + 2]) + fixed[i] * (fixed_X[i * 3 + 2] - X[i * 3 + 2]) + force[2];
			error_sum += B[i * 3 + 0] * B[i * 3 + 0] + B[i * 3 + 1] * B[i * 3 + 1] + B[i * 3 + 2] * B[i * 3 + 2];
		}
		return error_sum;
	}
    
	TYPE Get_Tet_Energy(int t0, int t1, TYPE stiffness_0, TYPE stiffness_1, TYPE stiffness_2, TYPE stiffness_3)
	{
		TYPE energy=0;
		for(int tet=t0; tet<t1; tet++)
			energy+=Compute_FM(tet, stiffness_0, stiffness_1, stiffness_2, stiffness_3, false);
		return energy;
	}

	TYPE Get_Energy(TYPE* _X, TYPE t)
	{
		Swap(X, _X);

        TYPE inv_t	= 1/t;

		TYPE _energy[8];

#pragma omp parallel for
		for(int i=0; i<8; i++)
		{
			int t0=tet_number*((i+0)/8.0);
			int t1=tet_number*((i+1)/8.0);
			_energy[i]=Get_Tet_Energy(t0, t1, stiffness_0, stiffness_1, stiffness_2, stiffness_3);
		}

		TYPE energy	= 0;
		for(int i=0; i<8; i++)
			energy+=_energy[i];

		//for(int tet=0; tet<tet_number; tet++)
		//	energy+=Compute_FM(tet, stiffness_0, stiffness_1, stiffness_2, stiffness_3, false);

		for(int i=0; i<number; i++)
		{
			TYPE oc = M[i]*inv_t*inv_t;
            TYPE c  = oc+fixed[i];
         
			energy+=oc*(S[i*3+0]-X[i*3+0])*(S[i*3+0]-X[i*3+0])*0.5;					// kinetic energy
			energy+=oc*(S[i*3+1]-X[i*3+1])*(S[i*3+1]-X[i*3+1])*0.5;					// kinetic energy
			energy+=oc*(S[i*3+2]-X[i*3+2])*(S[i*3+2]-X[i*3+2])*0.5;					// kinetic energy
			energy+=(c-oc)*(fixed_X[i*3+0]-X[i*3+0])*(fixed_X[i*3+0]-X[i*3+0])*0.5;	// fixed energy
			energy+=(c-oc)*(fixed_X[i*3+1]-X[i*3+1])*(fixed_X[i*3+1]-X[i*3+1])*0.5;	// fixed energy
			energy+=(c-oc)*(fixed_X[i*3+2]-X[i*3+2])*(fixed_X[i*3+2]-X[i*3+2])*0.5;	// fixed energy
	        energy+=-gravity*M[i]*X[i*3+1];
		}

		Swap(X, _X);

		return energy;
	}

	TYPE Wolfe_Stepping(TYPE* D, TYPE t, TYPE value)
	{
		TYPE stepping=2;
		TYPE e0=Get_Energy(X, t);
		TYPE e1;
		int l=0;
		for(l=0; l<1000; l++)
		{
			for(int i=0; i<number*3; i++)
				next_X[i]=X[i]-stepping*D[i];

			e1=Get_Energy(next_X, t);
			
		//	printf("value: %f, %f, %f\n", e1, e0, 0.01*stepping*value);

		//	if(e1<e0+0.001*stepping*value)	break;
			if(e1<e0)	break;


			stepping*=0.7;
		}
		printf("stepping: %ef: %ef (%ef, %ef)\n", stepping, value, e0, e1);
		//getchar();
		return stepping;
	}

	void LBFGS_Update(TYPE t, int iterations)
	{
		//Damp the system first
        for(int i=0; i<number; i++)
        {
			V[i * 3 + 0] *= 0.999;
			V[i * 3 + 1] *= 0.999;
			V[i * 3 + 2] *= 0.999;

			S[i * 3 + 0] = X[i * 3 + 0] + V[i * 3 + 0] * t;
			S[i * 3 + 1] = X[i * 3 + 1] + V[i * 3 + 1] * t;
			S[i * 3 + 2] = X[i * 3 + 2] + V[i * 3 + 2] * t;	
        }
		
		for (int i = 0; i < number; i++)
		{
			X[i*3+0]+=V[i*3+0]*t;
			X[i*3+1]+=V[i*3+1]*t;
			X[i*3+2]+=V[i*3+2]*t;
		//	X[i * 3 + 0] += (V[i * 3 + 0] + (V[i * 3 + 0] - prev_V[i * 3 + 0])*0.2)*t;
		//	X[i * 3 + 1] += (V[i * 3 + 1] + (V[i * 3 + 1] - prev_V[i * 3 + 1])*0.2)*t;
		//	X[i * 3 + 2] += (V[i * 3 + 2] + (V[i * 3 + 2] - prev_V[i * 3 + 2])*0.2)*t;
		}

        TYPE inv_t              = 1/t;
		TYPE last_error_sum     = 1;
        TYPE error_sum          = 1;
		TYPE rho                = 0.9992;
        TYPE omega;
		TYPE theta				= 1;

		for(int l=0; l<iterations; l++)
		{
			bool update_C=l%8==0;

			memset(lambda, 0, sizeof(TYPE)*tet_number);
			memset(ext_C, 0, sizeof(TYPE)*number);
			memset(Force, 0, sizeof(TYPE)*number*3);
			for(int tet=0; tet<tet_number; tet++)	
				Compute_FM(tet, stiffness_0, stiffness_1, stiffness_2, stiffness_3, update_C);
            
			last_error_sum  = error_sum;
            error_sum       = 0;
            for(int i=0; i<number; i++)
            {
                TYPE oc = M[i]*inv_t*inv_t;
                TYPE c  = oc+fixed[i];        
				                
                if(update_C)
                {
                    C[i*3+0]=c;
                    C[i*3+1]=c;
                    C[i*3+2]=c;
                }                
                TYPE force[3]={0, 0, 0};
                for(int index=vtt_num[i]; index<vtt_num[i+1]; index++)
                {
                    if(update_C)
                    {
                        C[i*3+0]+=C_Temp[VTT[index]*3+0];
                        C[i*3+1]+=C_Temp[VTT[index]*3+1];
                        C[i*3+2]+=C_Temp[VTT[index]*3+2];
                    }
                    force[0]+=F_Temp[VTT[index]*3+0];
                    force[1]+=F_Temp[VTT[index]*3+1];
                    force[2]+=F_Temp[VTT[index]*3+2];
                }             
                
				B[i*3+0]=oc*(S[i*3+0]-X[i*3+0])+fixed[i]*(fixed_X[i*3+0]-X[i*3+0])+Force[i*3+0];
                B[i*3+1]=oc*(S[i*3+1]-X[i*3+1])+fixed[i]*(fixed_X[i*3+1]-X[i*3+1])+Force[i*3+1]+M[i]*gravity;
                B[i*3+2]=oc*(S[i*3+2]-X[i*3+2])+fixed[i]*(fixed_X[i*3+2]-X[i*3+2])+Force[i*3+2];
                error_sum+=B[i*3+0]*B[i*3+0]+B[i*3+1]*B[i*3+1]+B[i*3+2]*B[i*3+2];

				B[i*3+0]=-B[i*3+0];
				B[i*3+1]=-B[i*3+1];
				B[i*3+2]=-B[i*3+2];
			}

			printf("%f\n", Get_Energy(X, t));

			//lbfgs_solver->Run(X, B, D, C);

			memcpy(D, B, sizeof(TYPE)*number*3);


			TYPE stepping=0.00000002f;//0.00000002f

			if(l!=0)
			{
				//stepping=1.0;
				//TYPE value=0;
				//for(int i=0; i<number*3; i++)	value-=B[i]*D[i];
				//stepping=Wolfe_Stepping(D, t, value);				
			}

			for(int i=0; i<number*3; i++)
			{
				//if(l!=0)	printf("check %d: %f, %f\n", l, B[i], D[i]);
				next_X[i]=X[i]-stepping*D[i];
			}

			
			if (l == 0)				omega=1;
			if (l == 1)				{ rho = 0.9;	omega = 2 / (2 - rho*rho);		}
			if (l > 1 && l < 7)		{ rho = 0.9;	omega = 4 / (4 - rho*rho*omega);}
			if(l==6)				omega=1;
			if (l == 7)				{ rho = 0.99;	omega = 2 / (2 - rho*rho);		}
			if (l > 7 && l < 32)	{ rho = 0.99;	omega = 4 / (4 - rho*rho*omega);}
			if (l == 31)			omega=1;
			if (l == 32)			{ rho = 0.999;  omega = 2 / (2 - rho*rho);		}	//0.9992
			if (l > 32)				{ rho = 0.999;  omega = 4 / (4 - rho*rho*omega);}	//0.9992

            

		/*	memcpy(prev_Y, Y, sizeof(TYPE)*number*3);
			memcpy(Y, next_X, sizeof(TYPE)*number*3);
			TYPE new_theta=(sqrtf(theta*theta*theta*theta+4*theta*theta)-theta*theta)*0.5;
			TYPE beta=theta*(1-theta)/(theta*theta+new_theta);
			theta=new_theta;
			omega=beta+1;

			omega=1;*/

            for(int i=0; i<number*3; i++)
            {
                if(omega!=1)	next_X[i]=omega*(next_X[i]-prev_X[i])+prev_X[i];

				//if(omega!=1)	next_X[i]=omega*(next_X[i]-prev_Y[i])+prev_Y[i];

                prev_X[i]=X[i];
                X[i]=next_X[i];
            }
		}

		for(int i=0; i<number; i++)
        {
			prev_V[i*3+0]=V[i*3+0];
			prev_V[i*3+1]=V[i*3+1];
			prev_V[i*3+2]=V[i*3+2];

            V[i*3+0]+=(X[i*3+0]-S[i*3+0])*inv_t;
            V[i*3+1]+=(X[i*3+1]-S[i*3+1])*inv_t;
            V[i*3+2]+=(X[i*3+2]-S[i*3+2])*inv_t;
        }
	}

	void Jacobi_Update(TYPE t, int iterations)
	{


        //Damp the system first
        for(int i=0; i<number; i++)
        {
			V[i * 3 + 0] *= 0.999;
			V[i * 3 + 1] *= 0.999;
			V[i * 3 + 2] *= 0.999;

			S[i * 3 + 0] = X[i * 3 + 0] + V[i * 3 + 0] * t;
			S[i * 3 + 1] = X[i * 3 + 1] + V[i * 3 + 1] * t;
			S[i * 3 + 2] = X[i * 3 + 2] + V[i * 3 + 2] * t;	
        }

		
		for (int i = 0; i < number; i++)
		{
		//	X[i*3+0]+=V[i*3+0]*t;
		//	X[i*3+1]+=V[i*3+1]*t;
		//	X[i*3+2]+=V[i*3+2]*t;

			//printf("V: %f, %f, %f; %f, %f, %f\n", V[i*3+0], V[i*3+1], V[i*3+2], prev_V[i*3+0], prev_V[i*3+1], prev_V[i*3+2]);

			X[i * 3 + 0] += (V[i * 3 + 0] + (V[i * 3 + 0] - prev_V[i * 3 + 0])*0)*t;
			X[i * 3 + 1] += (V[i * 3 + 1] + (V[i * 3 + 1] - prev_V[i * 3 + 1])*0)*t;
			X[i * 3 + 2] += (V[i * 3 + 2] + (V[i * 3 + 2] - prev_V[i * 3 + 2])*0)*t;
		}

	/*	
		for (int tet = 0; tet < tet_number; tet++)
				Compute_FM(tet, stiffness_0, stiffness_1, stiffness_2, stiffness_3, true);
		for (int i = 0; i < number; i++)
		{
			TYPE oc = M[i] / (t*t);
			TYPE c = oc + fixed[i];

			TYPE force[3] = { 0, 0, 0 };
			for (int index = vtt_num[i]; index < vtt_num[i + 1]; index++)
			{
				force[0] += F_Temp[VTT[index] * 3 + 0];
				force[1] += F_Temp[VTT[index] * 3 + 1];
				force[2] += F_Temp[VTT[index] * 3 + 2];
			}
			B[i * 3 + 0] = oc*(S[i * 3 + 0] - X[i * 3 + 0]) + fixed[i] * (fixed_X[i * 3 + 0] - X[i * 3 + 0]) + force[0];
			B[i * 3 + 1] = oc*(S[i * 3 + 1] - X[i * 3 + 1]) + fixed[i] * (fixed_X[i * 3 + 1] - X[i * 3 + 1]) + force[1] + M[i] * gravity;
			B[i * 3 + 2] = oc*(S[i * 3 + 2] - X[i * 3 + 2]) + fixed[i] * (fixed_X[i * 3 + 2] - X[i * 3 + 2]) + force[2];
		}
		TYPE* Y = new TYPE[number * 3];
		TYPE* Z = new TYPE[number * 3];
		A_Times_Vector(V, Y, 1/t);

		TYPE value0=0, value1=0;
		for(int i=0; i<number*3; i++)
		{
			value0+=V[i]*Y[i];
			value1+=V[i]*B[i];
			//printf("vY: %f, %f\n",V[i], Y[i]);
		}

		TYPE tt;
		if(value0==0)	tt=0;
		else
		{
			tt=value1/value0;
			if(tt<-t)	tt=-t;
		//	if(tt> t)	tt= t;
		}		

		//printf("t: %f\n", value1/value0);
		//printf("ERROR: %f\n", Get_Error(t));		
		TYPE start_error=Get_Energy(X, t);


		for (int i = 0; i < number; i++)
		{
			X[i * 3 + 0] = X[i * 3 + 0] + V[i * 3 + 0] * t;
			X[i * 3 + 1] = X[i * 3 + 1] + V[i * 3 + 1] * t;
			X[i * 3 + 2] = X[i * 3 + 2] + V[i * 3 + 2] * t;
		}
		TYPE end_error=Get_Energy(X, t);
		for (int i = 0; i < number; i++)
		{
			X[i * 3 + 0] = X[i * 3 + 0] - V[i * 3 + 0] * t;
			X[i * 3 + 1] = X[i * 3 + 1] - V[i * 3 + 1] * t;
			X[i * 3 + 2] = X[i * 3 + 2] - V[i * 3 + 2] * t;
		}

		for (int i = 0; i < number; i++)
		{
			X[i * 3 + 0] = X[i * 3 + 0] + V[i * 3 + 0] * tt;
			X[i * 3 + 1] = X[i * 3 + 1] + V[i * 3 + 1] * tt;
			X[i * 3 + 2] = X[i * 3 + 2] + V[i * 3 + 2] * tt;
		}
		TYPE mid_error = Get_Energy(X, t);	
		*/
		//printf("tt: %f (%f, %f) (%ef, %ef, %ef)\n", tt, value1, value0, start_error, end_error, mid_error);
		
        TYPE sum	= 0;
		TYPE time_a	= 0;
        
        TYPE inv_t              = 1/t;
        TYPE rho                = 0.9998;//9995;//0.5;//0.9995;//0.999
        TYPE under_relaxation   = 0.1;//38; (init) //0.43
        TYPE omega;
        TYPE last_error_sum     = 1;
        TYPE error_sum          = 1;
		TYPE theta				= 1;

		memset(lambda, 0, sizeof(TYPE)*tet_number);

        for(int l=0; l<iterations; l++)
        {
 			TIMER timer_a;
           
            //Compute force (and central matrix), store them into Temps
			bool update_C=l==0;//l%32==0;
			
			memset(lambda, 0, sizeof(TYPE)*tet_number);

			memset(ext_C, 0, sizeof(TYPE)*number);
			memset(Force, 0, sizeof(TYPE)*number*3);
	

#pragma omp parallel for
			for(int tet=0; tet<tet_number; tet++)	
				Compute_FM(tet, stiffness_0, stiffness_1, stiffness_2, stiffness_3, update_C);
            

			/*int index = 12834*3+1;
			TYPE gap = 0.0001;
			printf("f: %f\n", Force[index]);
			X[index]-=gap*0.5;
			TYPE e0=0;
			for(int tet=0; tet<tet_number; tet++)
				e0+=Compute_FM(tet, stiffness_0, stiffness_1, stiffness_2, stiffness_3);
			printf("e0: %f\n", e0);
			X[index]+=gap;
			TYPE e1=0;
			for(int tet=0; tet<tet_number; tet++)
				e1+=Compute_FM(tet, stiffness_0, stiffness_1, stiffness_2, stiffness_3);
			printf("e1: %f\n", e1);
			printf("? %f;\n", -(e1-e0)/gap);
			getchar();*/

			last_error_sum  = error_sum;
            error_sum       = 0;


#pragma omp parallel for
            for(int i=0; i<number; i++)
            {
                TYPE oc = M[i]*inv_t*inv_t;
                TYPE c  = oc+fixed[i];
                
                if(update_C)
                {
                    C[i*3+0]=c;
                    C[i*3+1]=c;
                    C[i*3+2]=c;
                }
                
                TYPE force[3]={0, 0, 0};
                for(int index=vtt_num[i]; index<vtt_num[i+1]; index++)
                {
                    if(update_C)
                    {
                        C[i*3+0]+=C_Temp[VTT[index]*3+0];
                        C[i*3+1]+=C_Temp[VTT[index]*3+1];
                        C[i*3+2]+=C_Temp[VTT[index]*3+2];
                    }
                    force[0]+=F_Temp[VTT[index]*3+0];
                    force[1]+=F_Temp[VTT[index]*3+1];
                    force[2]+=F_Temp[VTT[index]*3+2];
                }
               
                
                B[i*3+0]=oc*(S[i*3+0]-X[i*3+0])+fixed[i]*(fixed_X[i*3+0]-X[i*3+0])+force[0];
                B[i*3+1]=oc*(S[i*3+1]-X[i*3+1])+fixed[i]*(fixed_X[i*3+1]-X[i*3+1])+force[1]+M[i]*gravity;
                B[i*3+2]=oc*(S[i*3+2]-X[i*3+2])+fixed[i]*(fixed_X[i*3+2]-X[i*3+2])+force[2];
                error_sum+=B[i*3+0]*B[i*3+0]+B[i*3+1]*B[i*3+1]+B[i*3+2]*B[i*3+2];

                next_X[i*3+0]=X[i*3+0]+B[i*3+0]*under_relaxation/(C[i*3+0]+ext_C[i]);
                next_X[i*3+1]=X[i*3+1]+B[i*3+1]*under_relaxation/(C[i*3+1]+ext_C[i]);
                next_X[i*3+2]=X[i*3+2]+B[i*3+2]*under_relaxation/(C[i*3+2]+ext_C[i]);
            }

		   if(l!=0)	time_a+=timer_a.Get_Time();


			/*printf("start: %f\n", error_sum);            
            memset(next_X, 0, sizeof(TYPE)*number*3);
            for(int it=0; it<60; it++)
            {                
                Times_A(next_X, Y, inv_t);                
                error_sum=0;
                for(int i=0; i<number; i++)
                {
                    next_X[i*3+0]=next_X[i*3+0]+(B[i*3+0]-Y[i*3+0])*0.0000001;
                    next_X[i*3+1]=next_X[i*3+1]+(B[i*3+1]-Y[i*3+1])*0.0000001;
                    next_X[i*3+2]=next_X[i*3+2]+(B[i*3+2]-Y[i*3+2])*0.0000001;
                    
                    error_sum+=(B[i*3+0]-Y[i*3+0])*(B[i*3+0]-Y[i*3+0]);
                    error_sum+=(B[i*3+1]-Y[i*3+1])*(B[i*3+1]-Y[i*3+1]);
                    error_sum+=(B[i*3+2]-Y[i*3+2])*(B[i*3+2]-Y[i*3+2]);
                }                
                printf("error %d: %f\n", it, error_sum);
            }*/  
            /*for(int i=0; i<number; i++)
            {
                next_X[i*3+0]+=X[i*3+0];
                next_X[i*3+1]+=X[i*3+1];
                next_X[i*3+2]+=X[i*3+2];
            }*/
            
            //if(update_C)
            //printf("error %d %f: %f (%f)\n", l, rho, error_sum, error_sum/last_error_sum);
			//if(l==iterations-1)
		//	printf("%f, %f\n", Get_Energy(X, t), error_sum);

            
	/*		if (l == 0)				omega=1;
			if (l == 1)				{ rho = 0.96;	omega = 2 / (2 - rho*rho);		}
			if (l > 1 && l < 7)		{ rho = 0.96;	omega = 4 / (4 - rho*rho*omega);}
			if(l==6)				omega=1;
			if (l == 7)				{ rho = 0.99;	omega = 2 / (2 - rho*rho);		}
			if (l > 7 && l < 12)	{ rho = 0.99;	omega = 4 / (4 - rho*rho*omega);}
			if (l == 11)			omega=1;
			if (l == 12)			{ rho = 0.991;  omega = 2 / (2 - rho*rho);		}	
			if (l > 12 && l<20)		{ rho = 0.991;  omega = 4 / (4 - rho*rho*omega);}	
			if (l == 20)			omega=1;
			if (l == 21)			{ rho = 0.995;  omega = 2 / (2 - rho*rho);		}
			if (l > 21)				{ rho = 0.995;  omega = 4 / (4 - rho*rho*omega);}	//0.99996
	*/	
			if (l == 0)				omega=1;
			if (l == 1)				{ rho = 0.96;	omega = 2 / (2 - rho*rho);		}
			if (l > 1 && l < 7)		{ rho = 0.96;	omega = 4 / (4 - rho*rho*omega);}
			if(l==6)				omega=1;
			if (l == 7)				{ rho = 0.99;	omega = 2 / (2 - rho*rho);		}
			if (l > 7 && l < 12)	{ rho = 0.99;	omega = 4 / (4 - rho*rho*omega);}
			if (l == 11)			omega=1;
			if (l == 12)			{ rho = 0.999;  omega = 2 / (2 - rho*rho);		}	//0.9992
			if (l > 12 && l<20)		{ rho = 0.999;  omega = 4 / (4 - rho*rho*omega);}	//0.9992
			if (l == 20)			omega=1;
			if (l == 21)			{ rho = 0.99996;  omega = 2 / (2 - rho*rho);		}
			if (l > 21)				{ rho = 0.99996;  omega = 4 / (4 - rho*rho*omega);}	//0.9999
			
			//if(l<=10)		omega=1;
            //else if(l==11)	omega=2/(2-rho*rho);
            //else			omega=4/(4-rho*rho*omega);			
			//omega=1;
            
			/*memcpy(prev_Y, Y, sizeof(TYPE)*number*3);
			memcpy(Y, next_X, sizeof(TYPE)*number*3);
			TYPE new_theta=(sqrtf(theta*theta*theta*theta+4*theta*theta)-theta*theta)*0.5;
			TYPE beta=theta*(1-theta)/(theta*theta+new_theta);
			theta=new_theta;
			omega=beta+1;*/

			//for(int tet=0; tet<tet_number*12; tet++)	
			//	sum+=C_Temp[tet];
			//	continue;


            for(int i=0; i<number*3; i++)
            {
                next_X[i]=omega*(next_X[i]-prev_X[i])+prev_X[i];

				//if(omega!=1)	next_X[i]=omega*(next_X[i]-prev_Y[i])+prev_Y[i];

                prev_X[i]=X[i];
                X[i]=next_X[i];
            }			

        }
        
        
        for(int i=0; i<number; i++)
        {
			prev_V[i*3+0]=V[i*3+0];
			prev_V[i*3+1]=V[i*3+1];
			prev_V[i*3+2]=V[i*3+2];

            V[i*3+0]+=(X[i*3+0]-S[i*3+0])*inv_t;
            V[i*3+1]+=(X[i*3+1]-S[i*3+1])*inv_t;
            V[i*3+2]+=(X[i*3+2]-S[i*3+2])*inv_t;
        }
        
       // delete[] Y;
        //delete[] Z;
  
		TIMER timer;

		printf("%f; %f\n", time_a, Get_Energy(X, t));
		printf("energy time: %f\n", timer.Get_Time());
        //TYPE sum=0;
		//for(int i=0; i<number*3; i++)   sum+=X[i];
		//printf("sum: %f\n", sum);        
        //printf("V: %f, %f, %f\n", V[0]/t, V[1]/t, V[2]/t);
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

    bool Read_Float_File(const char *file_name)
    {
        std::fstream input;
        input.open(file_name,std::ios::in|std::ios::binary);
        if(!input.is_open())	{printf("Error, file not open.\n");	return false;}
        
        float* input_X=new float[number*3];
        float* input_V=new float[number*3];
        
        Read_Binaries(input, input_X, number*3);
        Read_Binaries(input, input_V, number*3);
        
        for(int i=0; i<number*3; i++)
        {
            X[i]=input_X[i];
            V[i]=input_V[i];
        }
        
        delete[] input_X;
        delete[] input_V;
        
        input.close();
        return true;
    }
};


#endif