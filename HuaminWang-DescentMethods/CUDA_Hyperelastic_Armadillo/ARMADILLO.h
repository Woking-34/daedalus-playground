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
//  Class ARMADILLO
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef	__WHMIN_ARMADILLO_H__
#define __WHMIN_ARMADILLO_H__

#include "../lib/CUDA_HYPER_TET_MESH.h"


template <class TYPE>
class ARMADILLO: public CUDA_HYPER_TET_MESH<TYPE> 
{
public:
	ARMADILLO()
	{
		//Read_Original_File("armadillo.1");
		Read_Original_File("armadillo_10k.1");
		//Read_Original_File("sorted_armadillo");

		//Create_A_Tet();
		Scale(0.008);
		Centralize();
		Rotate_X(-0.2);

		printf("N: %d, %d\n", number, tet_number);
		for(int v=0; v<number; v++)
			//if(X[v*3+1]>-0.04 && X[v*3+1]<0)		
			//	if(fabsf(X[v*3+1]+0.01)<1*(X[v*3+2]-0.1))
			if(fabsf(X[v*3+1]+0.01)<2*(X[v*3+2]-0.2))
				fixed[v]=100000;	//10000000
		//10000000
		Rotate_X(1.2);

		//Neo-Hookean
		stiffness_0 = 2000000;	//2000000
		stiffness_1 = 2000000;	//2000000

		//stVK
	//	stiffness_0 = 100000;	//2000000
	//	stiffness_1 = 6000000;	//2000000
	//	stiffness_k = 5000;	//2000

		//Mooney
	//	stiffness_0	= 2000000;	//2000000
     //   stiffness_1	= 2000000;	//2000000
     //   stiffness_2	= 2000000;	//2000000
	//	stiffness_3 = 0.5;



		control_mag	= 100000;  //100000

		//Neo-Hookean
		stiffness_0 = 2000000;	//2000000
		stiffness_1 = 2000000;	//2000000

		//stVK

		model		= NH_MODEL;
		stiffness_0	= 2000000;	//2000000
		stiffness_1	= 20000000;	//2000000
		stiffness_2	= 0;	//2000000
		stiffness_3 = 0.5;
		stiffness_p	= 10000000;

		lower_bound = 0.15;
		upper_bound = 1000.0;



		//Mooney
		// stiffness_0	= 2000000;	//2000000
		// stiffness_1	= 2000000;	//2000000
		// stiffness_2	= 2000000;	//2000000
		// stiffness_3 = 0.5;


		control_mag	= 1000000;  
		profile_v[2] = 0.9997;
	}


};


#endif