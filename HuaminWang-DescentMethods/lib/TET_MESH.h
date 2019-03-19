#ifndef	__WHMIN_TET_H__
#define __WHMIN_TET_H__

#include <iostream>
#include <fstream>
#include <stdio.h>
#include "IO_FUNC.h"
#include "MY_MATH.h"
#include "INTERSECTION.h"
#include "MESH.h"
#include "DISTANCE.h"

template <class TYPE>
class TET_MESH
{
public:
	int		max_number;

	//Vertex
	int		number;
	TYPE*	X;
	TYPE*	M;

	int*	Tet;
	int		tet_number;
	TYPE*	Dm;
	TYPE*	inv_Dm;
	TYPE*	Vol;

	
	// triangles for rendering purpose
	int			t_number;
	int*		T;
	TYPE*		VN;		//Vertex Normal
	TYPE*		TN;		//Triangle Normal

	MESH<TYPE>	mesh;
	int*		v_map;

	TET_MESH(): number(0)
	{

		max_number	= 500000;
		X			= new TYPE	[max_number*3];
		M			= new TYPE	[max_number  ];
		Tet			= new int	[max_number*4];
		Dm			= new TYPE	[max_number*9];
		inv_Dm		= new TYPE	[max_number*9];
		Vol			= new TYPE	[max_number  ];
		T			= new int	[max_number*3];
		VN			= new TYPE	[max_number*3];
		TN			= new TYPE	[max_number*3];
		v_map		= new int	[max_number  ];
	}
	
	~TET_MESH()
	{
		if(X)		delete[] X;
		if(M)		delete[] M;
		if(Tet)		delete[] Tet;
		if(Dm)		delete[] Dm;
		if(inv_Dm)	delete[] inv_Dm;
		if(Vol)		delete[] Vol;
		if(T)		delete[] T;
		if(VN)		delete[] VN;
		if(TN)		delete[] TN;
		if(v_map)	delete[] v_map;
	}

	void Create_A_Tet()
	{
		X[ 0]=0;
		X[ 1]=0;
		X[ 2]=0;
		
		X[ 3]=0.1;
		X[ 4]=0;
		X[ 5]=0;
		
		X[ 6]=0;
		X[ 7]=0;
		X[ 8]=0.1;
		
		X[ 9]=0;
		X[10]=0.1;
		X[11]=0;

		X[12]=0;
		X[13]=-0.1;
		X[14]=0;
		number=4;

		Tet[0]=0;
		Tet[1]=1;
		Tet[2]=2;
		Tet[3]=3;
		Tet[4]=0;
		Tet[5]=1;
		Tet[6]=4;
		Tet[7]=2;
		tet_number=1;

		Build_Boundary_Triangles();
	}

	void Create_Prism(int v0, int v1, int v2, int v3, int v4, int v5)
	{
		Tet[tet_number * 4 +  0] = v0;
		Tet[tet_number * 4 +  1] = v1;
		Tet[tet_number * 4 +  2] = v2;
		Tet[tet_number * 4 +  3] = v3;
		Tet[tet_number * 4 +  4] = v2;
		Tet[tet_number * 4 +  5] = v3;
		Tet[tet_number * 4 +  6] = v1;
		Tet[tet_number * 4 +  7] = v5;
		Tet[tet_number * 4 +  8] = v1;
		Tet[tet_number * 4 +  9] = v4;
		Tet[tet_number * 4 + 10] = v5;
		Tet[tet_number * 4 + 11] = v3;
		tet_number += 3;
	}

	void Create_Block(TYPE min_x, TYPE max_x, int nx, TYPE min_y, TYPE max_y, int ny, TYPE min_z, TYPE max_z, int nz)
	{
		for(int i=0; i<=nx; i++)
		for(int j=0; j<=ny; j++)
		for(int k=0; k<=nz; k++)
		{
			X[number * 3 + 0] = min_x + (max_x - min_x)*((TYPE)(i)) / (nx);
			X[number * 3 + 1] = min_y + (max_y - min_y)*((TYPE)(j)) / (ny);
			X[number * 3 + 2] = min_z + (max_z - min_z)*((TYPE)(k)) / (nz);
			number++;
		}

		for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++)
		for (int k = 0; k < nz; k++)
		{
			int v0 = (i  )*(ny + 1)*(nz + 1) + (j  )*(nz + 1) + (k  );
			int v1 = (i  )*(ny + 1)*(nz + 1) + (j  )*(nz + 1) + (k+1);
			int v2 = (i  )*(ny + 1)*(nz + 1) + (j+1)*(nz + 1) + (k  );
			int v3 = (i  )*(ny + 1)*(nz + 1) + (j+1)*(nz + 1) + (k+1);
			int v4 = (i+1)*(ny + 1)*(nz + 1) + (j  )*(nz + 1) + (k  );
			int v5 = (i+1)*(ny + 1)*(nz + 1) + (j  )*(nz + 1) + (k+1);
			int v6 = (i+1)*(ny + 1)*(nz + 1) + (j+1)*(nz + 1) + (k  );
			int v7 = (i+1)*(ny + 1)*(nz + 1) + (j+1)*(nz + 1) + (k+1);

			Tet[tet_number * 4 +  0] = v0;
			Tet[tet_number * 4 +  1] = v1;
			Tet[tet_number * 4 +  2] = v2;
			Tet[tet_number * 4 +  3] = v4;
			tet_number++;

			Tet[tet_number * 4 + 0] = v2;
			Tet[tet_number * 4 + 1] = v6;
			Tet[tet_number * 4 + 2] = v4;
			Tet[tet_number * 4 + 3] = v5;
			tet_number++;

			Tet[tet_number * 4 + 0] = v1;
			Tet[tet_number * 4 + 1] = v4;
			Tet[tet_number * 4 + 2] = v5;
			Tet[tet_number * 4 + 3] = v2;
			tet_number++;

			Tet[tet_number * 4 + 0] = v5;
			Tet[tet_number * 4 + 1] = v6;
			Tet[tet_number * 4 + 2] = v7;
			Tet[tet_number * 4 + 3] = v3;
			tet_number++;

			Tet[tet_number * 4 + 0] = v2;
			Tet[tet_number * 4 + 1] = v1;
			Tet[tet_number * 4 + 2] = v3;
			Tet[tet_number * 4 + 3] = v5;
			tet_number++;

			Tet[tet_number * 4 + 0] = v2;
			Tet[tet_number * 4 + 1] = v3;
			Tet[tet_number * 4 + 2] = v6;
			Tet[tet_number * 4 + 3] = v5;
			tet_number++;
		}
		Build_Boundary_Triangles();

		
		for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++)
		{
			int k = 0;
			int v0 = (i)*(ny + 1)*(nz + 1) + (j)*(nz + 1) + (k);
			int v1 = (i)*(ny + 1)*(nz + 1) + (j + 1)*(nz + 1) + (k);
			int v2 = (i + 1)*(ny + 1)*(nz + 1) + (j)*(nz + 1) + (k);
			int v3 = (i + 1)*(ny + 1)*(nz + 1) + (j + 1)*(nz + 1) + (k);
			int m0 = mesh.number + (i)*(ny + 1) + (j);
			int m1 = mesh.number + (i)*(ny + 1) + (j + 1);
			int m2 = mesh.number + (i + 1)*(ny + 1) + (j);
			int m3 = mesh.number + (i + 1)*(ny + 1) + (j + 1);

			mesh.T[mesh.t_number * 3 + 0] = m0;
			mesh.T[mesh.t_number * 3 + 1] = m2;
			mesh.T[mesh.t_number * 3 + 2] = m1;
			mesh.T[mesh.t_number * 3 + 3] = m1;
			mesh.T[mesh.t_number * 3 + 4] = m2;
			mesh.T[mesh.t_number * 3 + 5] = m3;
			mesh.t_number += 2;
			v_map[m0] = v0;
			v_map[m1] = v1;
			v_map[m2] = v2;
			v_map[m3] = v3;
		}
		mesh.number += (nx + 1)*(ny + 1);

		
		for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++)
		{
			int k = nz;
			int v0 = (i)*(ny + 1)*(nz + 1) + (j)*(nz + 1) + (k);
			int v1 = (i)*(ny + 1)*(nz + 1) + (j + 1)*(nz + 1) + (k);
			int v2 = (i + 1)*(ny + 1)*(nz + 1) + (j)*(nz + 1) + (k);
			int v3 = (i + 1)*(ny + 1)*(nz + 1) + (j + 1)*(nz + 1) + (k);
			int m0 = mesh.number + (i)*(ny + 1) + (j);
			int m1 = mesh.number + (i)*(ny + 1) + (j + 1);
			int m2 = mesh.number + (i + 1)*(ny + 1) + (j);
			int m3 = mesh.number + (i + 1)*(ny + 1) + (j + 1);

			mesh.T[mesh.t_number * 3 + 0] = m0;
			mesh.T[mesh.t_number * 3 + 1] = m1;
			mesh.T[mesh.t_number * 3 + 2] = m2;
			mesh.T[mesh.t_number * 3 + 3] = m1;
			mesh.T[mesh.t_number * 3 + 4] = m3;
			mesh.T[mesh.t_number * 3 + 5] = m2;
			mesh.t_number += 2;
			v_map[m0] = v0;
			v_map[m1] = v1;
			v_map[m2] = v2;
			v_map[m3] = v3;
		}
		mesh.number += (nx + 1)*(ny + 1);
		

		for (int j = 0; j < ny; j++)
		for (int k = 0; k < nz; k++)
		{
			int i = 0;
			int v0 = (i  )*(ny + 1)*(nz + 1) + (j  )*(nz + 1) + (k  );
			int v1 = (i  )*(ny + 1)*(nz + 1) + (j  )*(nz + 1) + (k+1);
			int v2 = (i  )*(ny + 1)*(nz + 1) + (j+1)*(nz + 1) + (k  );
			int v3 = (i  )*(ny + 1)*(nz + 1) + (j+1)*(nz + 1) + (k+1);

			int m0 = mesh.number + (j)*(nz + 1) + (k);
			int m1 = mesh.number + (j)*(nz + 1) + (k+1);
			int m2 = mesh.number + (j + 1)*(nz + 1) + (k);
			int m3 = mesh.number + (j + 1)*(nz + 1) + (k+1);

			mesh.T[mesh.t_number * 3 + 0] = m0;
			mesh.T[mesh.t_number * 3 + 1] = m2;
			mesh.T[mesh.t_number * 3 + 2] = m1;
			mesh.T[mesh.t_number * 3 + 3] = m1;
			mesh.T[mesh.t_number * 3 + 4] = m2;
			mesh.T[mesh.t_number * 3 + 5] = m3;
			mesh.t_number += 2;
			v_map[m0] = v0;
			v_map[m1] = v1;
			v_map[m2] = v2;
			v_map[m3] = v3;
		}
		mesh.number += (nz + 1)*(ny + 1);


		for (int j = 0; j < ny; j++)
		for (int k = 0; k < nz; k++)
		{
			int i = nx;
			int v0 = (i)*(ny + 1)*(nz + 1) + (j)*(nz + 1) + (k);
			int v1 = (i)*(ny + 1)*(nz + 1) + (j)*(nz + 1) + (k + 1);
			int v2 = (i)*(ny + 1)*(nz + 1) + (j + 1)*(nz + 1) + (k);
			int v3 = (i)*(ny + 1)*(nz + 1) + (j + 1)*(nz + 1) + (k + 1);

			int m0 = mesh.number + (j)*(nz + 1) + (k);
			int m1 = mesh.number + (j)*(nz + 1) + (k + 1);
			int m2 = mesh.number + (j + 1)*(nz + 1) + (k);
			int m3 = mesh.number + (j + 1)*(nz + 1) + (k + 1);

			mesh.T[mesh.t_number * 3 + 0] = m0;
			mesh.T[mesh.t_number * 3 + 1] = m1;
			mesh.T[mesh.t_number * 3 + 2] = m2;
			mesh.T[mesh.t_number * 3 + 3] = m1;
			mesh.T[mesh.t_number * 3 + 4] = m3;
			mesh.T[mesh.t_number * 3 + 5] = m2;
			mesh.t_number += 2;
			v_map[m0] = v0;
			v_map[m1] = v1;
			v_map[m2] = v2;
			v_map[m3] = v3;
		}
		mesh.number += (nz + 1)*(ny + 1);


		for (int i = 0; i < nx; i++)
		for (int k = 0; k < nz; k++)
		{
			int j = 0;
			int v0 = (i  )*(ny + 1)*(nz + 1) + (j)*(nz + 1) + (k);
			int v1 = (i  )*(ny + 1)*(nz + 1) + (j)*(nz + 1) + (k + 1);
			int v2 = (i+1)*(ny + 1)*(nz + 1) + (j)*(nz + 1) + (k);
			int v3 = (i+1)*(ny + 1)*(nz + 1) + (j)*(nz + 1) + (k + 1);

			int m0 = mesh.number + (i    )*(nz + 1) + (k    );
			int m1 = mesh.number + (i    )*(nz + 1) + (k + 1);
			int m2 = mesh.number + (i + 1)*(nz + 1) + (k    );
			int m3 = mesh.number + (i + 1)*(nz + 1) + (k + 1);

			mesh.T[mesh.t_number * 3 + 0] = m0;
			mesh.T[mesh.t_number * 3 + 1] = m1;
			mesh.T[mesh.t_number * 3 + 2] = m2;
			mesh.T[mesh.t_number * 3 + 3] = m1;
			mesh.T[mesh.t_number * 3 + 4] = m3;
			mesh.T[mesh.t_number * 3 + 5] = m2;
			mesh.t_number += 2;
			v_map[m0] = v0;
			v_map[m1] = v1;
			v_map[m2] = v2;
			v_map[m3] = v3;
		}
		mesh.number += (nz + 1)*(nx + 1);


		for (int i = 0; i < nx; i++)
		for (int k = 0; k < nz; k++)
		{
			int j = ny;
			int v0 = (i)*(ny + 1)*(nz + 1) + (j)*(nz + 1) + (k);
			int v1 = (i)*(ny + 1)*(nz + 1) + (j)*(nz + 1) + (k + 1);
			int v2 = (i + 1)*(ny + 1)*(nz + 1) + (j)*(nz + 1) + (k);
			int v3 = (i + 1)*(ny + 1)*(nz + 1) + (j)*(nz + 1) + (k + 1);

			int m0 = mesh.number + (i)*(nz + 1) + (k);
			int m1 = mesh.number + (i)*(nz + 1) + (k + 1);
			int m2 = mesh.number + (i + 1)*(nz + 1) + (k);
			int m3 = mesh.number + (i + 1)*(nz + 1) + (k + 1);

			mesh.T[mesh.t_number * 3 + 0] = m0;
			mesh.T[mesh.t_number * 3 + 1] = m2;
			mesh.T[mesh.t_number * 3 + 2] = m1;
			mesh.T[mesh.t_number * 3 + 3] = m1;
			mesh.T[mesh.t_number * 3 + 4] = m2;
			mesh.T[mesh.t_number * 3 + 5] = m3;
			mesh.t_number += 2;
			v_map[m0] = v0;
			v_map[m1] = v1;
			v_map[m2] = v2;
			v_map[m3] = v3;
		}
		mesh.number += (nz + 1)*(nx + 1);

		Update_Mesh();
	}

	void Update_Mesh()
	{
		for(int v=0; v<mesh.number; v++)
		{
			mesh.X[v * 3 + 0] = X[v_map[v] * 3 + 0];
			mesh.X[v * 3 + 1] = X[v_map[v] * 3 + 1];
			mesh.X[v * 3 + 2] = X[v_map[v] * 3 + 2];
		}
		mesh.Build_VN();
	}

	void Center(TYPE c[])
	{
		c[0]=c[1]=c[2]=0;
		TYPE mass_sum=0;
		for(int i=0; i<number; i++)
		{
			c[0]	 += X[i*3+0];
			c[1]	 += X[i*3+1];
			c[2]	 += X[i*3+2];
			mass_sum += 1;
		}
		c[0]/=mass_sum;
		c[1]/=mass_sum;
		c[2]/=mass_sum;
	}

	void Translate(TYPE tx, TYPE ty, TYPE tz)
	{
		for(int i=0; i<number; i++)
		{
			X[i*3+0]+=tx;
			X[i*3+1]+=ty;
			X[i*3+2]+=tz;
		}
	}

	void Centralize()
	{
		TYPE c[3];
		Center(c);
		Translate(-c[0], -c[1], -c[2]);
	}

	void Read_Original_File(char *name)
	{
		char filename[1024];
		int temp_value;
		int bound;

		sprintf(filename, "%s.node", name);

		std::string rootStr0 = std::string(SAMPLE_NAME) + "/data/";
		std::string filePath0 = filename;

		bool fileFound0 = findFullPath(rootStr0, filePath0);

		FILE *fp;
		fp=fopen(filePath0.c_str(), "r+");
		if(fp==NULL)	{printf("ERROR: file %s not open.\n", filename); return;}
		fscanf(fp, "%d %d %d %d\n", &number, &temp_value, &temp_value, &bound);
		if(bound==0)
			for(int i=0; i<number; i++)
			{
				float temp_x0, temp_x1, temp_x2;
				fscanf(fp, "%d %f %f %f\n", &temp_value, &temp_x0, &temp_x1, &temp_x2);
				X[i*3+0]=temp_x0;
				X[i*3+1]=temp_x1;
				X[i*3+2]=temp_x2;
			}
		else
			for(int i=0; i<number; i++)
			{
				float temp_x0, temp_x1, temp_x2;
				fscanf(fp, "%d %f %f %f %d\n", &temp_value, &temp_x0, &temp_x1, &temp_x2, &temp_value);
				X[i*3+0]=temp_x0;
				X[i*3+1]=temp_x1;
				X[i*3+2]=temp_x2;
				//printf("en %d: %f, %f, %f\n", i, X[i*3], X[i*3+1], X[i*3+2]);
			}

		fclose(fp);
		//for(int i=0; i<number; i++)
		//	printf("v %d: %f, %f, %f\n", i, X[i*3], X[i*3+1], X[i*3+2]);

		sprintf(filename, "%s.ele", name);

		std::string rootStr1 = std::string(SAMPLE_NAME) + "/data/";
		std::string filePath1 = filename;

		bool fileFound1 = findFullPath(rootStr1, filePath1);

		fp=fopen(filePath1.c_str(), "r+");
		if(fp==NULL)	{printf("ERROR: file %s not open.\n", filename); return;}
		fscanf(fp, "%d %d %d\n", &tet_number, &temp_value, &bound);
		
		if(bound==0)
			for(int i=0; i<tet_number; i++)
				fscanf(fp, "%d %d %d %d %d\n", &temp_value, &Tet[i*4+0], &Tet[i*4+1], &Tet[i*4+2], &Tet[i*4+3]);
		else if(bound==1)
			for(int i=0; i<tet_number; i++)
				fscanf(fp, "%d %d %d %d %d %d\n", &temp_value, &Tet[i*4+0], &Tet[i*4+1], &Tet[i*4+2], &Tet[i*4+3], &temp_value);
		fclose(fp);

		for(int i=0; i<tet_number; i++)
		{
			Tet[i*4+0]-=1;
			Tet[i*4+1]-=1;
			Tet[i*4+2]-=1;
			Tet[i*4+3]-=1;
		}

		//printf("Read %s: %d, %d\n", name, number, tet_number);
		Build_Boundary_Triangles();
	}

	void Write_Original_File(char *name)
	{
		char filename[1024];

		sprintf(filename, "%s.node", name);
		FILE *fp=fopen(filename, "w+");
		if(fp==NULL)	{printf("ERROR: file %s not open.\n", filename); return;}
		fprintf(fp, "%d %d %d %d\n", number, 3, 0, 0);
		for(int i=0; i<number; i++)
			fprintf(fp, "%d %f %f %f\n", i+1, X[i*3+0], X[i*3+1], X[i*3+2]);
		fclose(fp);

		sprintf(filename, "%s.ele", name);
		fp=fopen(filename, "w+");
		if(fp==NULL)	{printf("ERROR: file %s not open.\n", filename); return;}
		fprintf(fp, "%d %d %d\n", tet_number, 4, 0);

		for(int i=0; i<tet_number; i++)
			fprintf(fp, "%d %d %d %d %d\n", i+1, Tet[i*4+0]+1, Tet[i*4+1]+1, Tet[i*4+2]+1, Tet[i*4+3]+1);
		fclose(fp);
	}

	void Permutation(char *filename)
	{
		int* p=new int[number];
		int *q=new int[number];
		FILE *fp=fopen(filename, "r+");
		if(fp==NULL)	{printf("ERROR: file %s not open.\n", filename); return;}
		for(int i=0; i<number; i++)
		{
			fscanf(fp, "%d", &p[i]);
			p[i]-=1;
			//printf("read %d: %d\n", i, p[i]);
			q[p[i]]=i;
		}
		fclose(fp);

		TYPE*	new_X=new TYPE[number*3];
		int*	new_Tet=new int[tet_number*4];

		for(int i=0; i<number; i++)
		{
			int old_i=p[i];
			new_X[i*3+0]=X[old_i*3+0];
			new_X[i*3+1]=X[old_i*3+1];
			new_X[i*3+2]=X[old_i*3+2];
		}

		
		for(int t=0; t<tet_number; t++)
		{
			new_Tet[t*4+0]=q[Tet[t*4+0]];
			new_Tet[t*4+1]=q[Tet[t*4+1]];
			new_Tet[t*4+2]=q[Tet[t*4+2]];
			new_Tet[t*4+3]=q[Tet[t*4+3]];
		}

		memcpy(X, new_X, sizeof(TYPE)*number*3);
		memcpy(Tet, new_Tet, sizeof(int)*tet_number*4);
		delete[] new_X;
		delete[] new_Tet;
		delete[] p;
		delete[] q;
		
		Build_Boundary_Triangles();
	}

	void Scale(TYPE s)
	{
		for(int i=0; i<number*3; i++)	X[i]*=s;
	}

	void Scale(TYPE sx, TYPE sy, TYPE sz)
	{
		for(int i=0; i<number; i++)
		{
			X[i*3+0]*=sx;
			X[i*3+1]*=sy;
			X[i*3+2]*=sz;
		}
	}

	void Rotate_X(TYPE angle)
	{
		for(int i=0; i<number; i++)
		{
			TYPE y=X[i*3+1];
			TYPE z=X[i*3+2];
			X[i*3+1]= y*cos(angle)+z*sin(angle);
			X[i*3+2]=-y*sin(angle)+z*cos(angle);
		}
	}

	void Rotate_Y(TYPE angle)
	{
		for(int i=0; i<number; i++)
		{
			TYPE x=X[i*3+0];
			TYPE z=X[i*3+2];
			X[i*3+0]= x*cos(angle)+z*sin(angle);
			X[i*3+2]=-x*sin(angle)+z*cos(angle);
		}
	}

	void Rotate_Z(TYPE angle)
	{
		for(int i=0; i<number; i++)
		{
			TYPE x=X[i*3+0];
			TYPE y=X[i*3+1];
			X[i*3+0]= x*cos(angle)+y*sin(angle);
			X[i*3+1]=-x*sin(angle)+y*cos(angle);
		}
	}

	void Render(int visual_mode=0, int render_mode=0)
	{
		Build_VN();

		if(visual_mode==0)
		{			
			glDisable(GL_LIGHTING);
			glColor3f(1, 0, 0);
			for(int v=0; v<number; v++)
			{
				if(v!=0)	continue;

				glPushMatrix();
				glTranslatef(X[v*3+0], X[v*3+1], X[v*3+2]);
				glutSolidSphere(0.01, 5, 5);
				glPopMatrix();
			}
			
			/*for(int t=0; t<tet_number; t++)
			{
				int v0=Tet[t*4+0];
				int v1=Tet[t*4+1];
				int v2=Tet[t*4+2];
				int v3=Tet[t*4+3];

				glBegin(GL_LINES);
				glVertex3dv(&X[v0*3]);	glVertex3dv(&X[v1*3]);
				glVertex3dv(&X[v0*3]);	glVertex3dv(&X[v2*3]);
				glVertex3dv(&X[v0*3]);	glVertex3dv(&X[v3*3]);
				glVertex3dv(&X[v1*3]);	glVertex3dv(&X[v2*3]);
				glVertex3dv(&X[v1*3]);	glVertex3dv(&X[v3*3]);
				glVertex3dv(&X[v2*3]);	glVertex3dv(&X[v3*3]);
				glEnd();
			}*/


			glEnable(GL_LIGHTING);		
			float diffuse_color[3]={0.8, 0.8, 0.8};
			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse_color);
			//if(0)
			for(int i=0; i<t_number; i++)
			{
				TYPE *p0=&X[T[i*3+0]*3];
				TYPE *p1=&X[T[i*3+1]*3];
				TYPE *p2=&X[T[i*3+2]*3];

				if(render_mode==0)	glNormal3f(TN[i*3+0], TN[i*3+1], TN[i*3+2]);

				glBegin(GL_TRIANGLES);
				if(render_mode)		glNormal3f(VN[T[i*3+0]*3+0], VN[T[i*3+0]*3+1], VN[T[i*3+0]*3+2]);
				glVertex3d(p0[0], p0[1], p0[2]);
				if(render_mode)		glNormal3f(VN[T[i*3+1]*3+0], VN[T[i*3+1]*3+1], VN[T[i*3+1]*3+2]);
				glVertex3d(p1[0], p1[1], p1[2]);
				if(render_mode)		glNormal3f(VN[T[i*3+2]*3+0], VN[T[i*3+2]*3+1], VN[T[i*3+2]*3+2]);
				glVertex3d(p2[0], p2[1], p2[2]);
				glEnd();
			}
		}
	}


	void Build_Boundary_Triangles()
	{
		int *temp_T=new int[tet_number*4*4];

		for(int i=0; i<tet_number; i++)
		{
			temp_T[i*16+0]=Tet[i*4+0];
			temp_T[i*16+1]=Tet[i*4+1];
			temp_T[i*16+2]=Tet[i*4+2];
			temp_T[i*16+3]=1;

			temp_T[i*16+4]=Tet[i*4+0];
			temp_T[i*16+5]=Tet[i*4+2];
			temp_T[i*16+6]=Tet[i*4+3];
			temp_T[i*16+7]=1;

			temp_T[i*16+8]=Tet[i*4+0];
			temp_T[i*16+9]=Tet[i*4+3];
			temp_T[i*16+10]=Tet[i*4+1];
			temp_T[i*16+11]=1;

			temp_T[i*16+12]=Tet[i*4+1];
			temp_T[i*16+13]=Tet[i*4+3];
			temp_T[i*16+14]=Tet[i*4+2];
			temp_T[i*16+15]=1;
		}

		for(int i=0; i<tet_number*4; i++)
		{
			if(temp_T[i*4+1]<temp_T[i*4+0])
			{
				Swap(temp_T[i*4+0], temp_T[i*4+1]);
				temp_T[i*4+3]=(temp_T[i*4+3]+1)%2;
			}
			if(temp_T[i*4+2]<temp_T[i*4+0])
			{
				Swap(temp_T[i*4+0], temp_T[i*4+2]);
				temp_T[i*4+3]=(temp_T[i*4+3]+1)%2;
			}
			if(temp_T[i*4+2]<temp_T[i*4+1])
			{
				Swap(temp_T[i*4+1], temp_T[i*4+2]);
				temp_T[i*4+3]=(temp_T[i*4+3]+1)%2;
			}
		}

		QuickSort(temp_T, 0, tet_number*4-1);

		t_number=0;
		for(int i=0; i<tet_number*4; i++)
		{
			if(i!=tet_number*4-1 && temp_T[i*4+0]==temp_T[i*4+4] && temp_T[i*4+1]==temp_T[i*4+5] && temp_T[i*4+2]==temp_T[i*4+6])
			{
				i++;
				continue;
			}

			if(temp_T[i*4+3]==1)
			{
				T[t_number*3+0]=temp_T[i*4+0];
				T[t_number*3+1]=temp_T[i*4+1];
				T[t_number*3+2]=temp_T[i*4+2];
			}
			else
			{
				T[t_number*3+0]=temp_T[i*4+1];
				T[t_number*3+1]=temp_T[i*4+0];
				T[t_number*3+2]=temp_T[i*4+2];
			}
			t_number++;
		}

		delete []temp_T;
	}

	void QuickSort( int a[], int l, int r)
	{
		if( l < r ) 
		{
			int j=QuickSort_Partition(a, l, r);
			QuickSort(a, l, j-1);
			QuickSort(a, j+1, r);
		}
	}
	
	int QuickSort_Partition( int a[], int l, int r) 
	{
		int pivot[4], i, j, t[4];
		pivot[0] = a[l*4+0];
		pivot[1] = a[l*4+1];
		pivot[2] = a[l*4+2];
		pivot[3] = a[l*4+3];
		i = l; j = r+1;
		
		while( 1)
		{
			do ++i; while( (a[i*4+0]<pivot[0] || a[i*4+0]==pivot[0] && a[i*4+1]<pivot[1] || a[i*4+0]==pivot[0] && a[i*4+1]==pivot[1] && a[i*4+2]<=pivot[2]) && i <= r );
			do --j; while(  a[j*4+0]>pivot[0] || a[j*4+0]==pivot[0] && a[j*4+1]>pivot[1] || a[j*4+0]==pivot[0] && a[j*4+1]==pivot[1] && a[j*4+2]> pivot[2]);
			if( i >= j ) break;
			//Swap i and j
			t[0]=a[i*4+0];
			t[1]=a[i*4+1];
			t[2]=a[i*4+2];
			t[3]=a[i*4+3];
			a[i*4+0]=a[j*4+0];
			a[i*4+1]=a[j*4+1];
			a[i*4+2]=a[j*4+2];
			a[i*4+3]=a[j*4+3];
			a[j*4+0]=t[0];
			a[j*4+1]=t[1];
			a[j*4+2]=t[2];
			a[j*4+3]=t[3];
		}
		//Swap l and j
		t[0]=a[l*4+0];
		t[1]=a[l*4+1];
		t[2]=a[l*4+2];
		t[3]=a[l*4+3];
		a[l*4+0]=a[j*4+0];
		a[l*4+1]=a[j*4+1];
		a[l*4+2]=a[j*4+2];
		a[l*4+3]=a[j*4+3];
		a[j*4+0]=t[0];
		a[j*4+1]=t[1];
		a[j*4+2]=t[2];
		a[j*4+3]=t[3];
		return j;
	}

	void Build_TN()
	{
		memset(TN, 0, sizeof(TYPE)*t_number*3);
		for(int i=0; i<t_number; i++)
		{
			TYPE *p0=&X[T[i*3+0]*3];
			TYPE *p1=&X[T[i*3+1]*3];
			TYPE *p2=&X[T[i*3+2]*3];
			Normal(p0, p1, p2, &TN[i*3]);
		}
	}

	void Build_VN()
	{
		memset(VN, 0, sizeof(TYPE)*number*3);
		Build_TN();

		for(int i=0; i<t_number; i++)
		{
			int v0=T[i*3+0];
			int v1=T[i*3+1];
			int v2=T[i*3+2];
		
			VN[v0*3+0]+=TN[i*3+0];
			VN[v0*3+1]+=TN[i*3+1];
			VN[v0*3+2]+=TN[i*3+2];

			VN[v1*3+0]+=TN[i*3+0];
			VN[v1*3+1]+=TN[i*3+1];
			VN[v1*3+2]+=TN[i*3+2];

			VN[v2*3+0]+=TN[i*3+0];
			VN[v2*3+1]+=TN[i*3+1];
			VN[v2*3+2]+=TN[i*3+2];
		}

		TYPE length2, inv_length;
		for(int i=0; i<number; i++)
		{
			length2=VN[i*3+0]*VN[i*3+0]+VN[i*3+1]*VN[i*3+1]+VN[i*3+2]*VN[i*3+2];
			if(length2<1e-16f)	continue;
			inv_length=1.0f/sqrtf(length2);
			
			VN[i*3+0]*=inv_length;
			VN[i*3+1]*=inv_length;
			VN[i*3+2]*=inv_length;
		}
	}

	void Initialize()	
	{
		TYPE density=50000;

		memset(M, 0, sizeof(TYPE)*number);

		for(int t=0; t<tet_number; t++)
		{
			int p0=Tet[t*4+0]*3;
			int p1=Tet[t*4+1]*3;
			int p2=Tet[t*4+2]*3;
			int p3=Tet[t*4+3]*3;

			Dm[t*9+0]=X[p1+0]-X[p0+0];
			Dm[t*9+3]=X[p1+1]-X[p0+1];
			Dm[t*9+6]=X[p1+2]-X[p0+2];
			Dm[t*9+1]=X[p2+0]-X[p0+0];
			Dm[t*9+4]=X[p2+1]-X[p0+1];
			Dm[t*9+7]=X[p2+2]-X[p0+2];
			Dm[t*9+2]=X[p3+0]-X[p0+0];
			Dm[t*9+5]=X[p3+1]-X[p0+1];
			Dm[t*9+8]=X[p3+2]-X[p0+2];
			Vol[t]=fabs(Matrix_Inverse_3(&Dm[t*9], &inv_Dm[t*9]))/6.0;

			M[p0/3]+=Vol[t]*density;
			M[p1/3]+=Vol[t]*density;
			M[p2/3]+=Vol[t]*density;
			M[p3/3]+=Vol[t]*density;

			//printf("vol: %f\n", Vol[t]);
			//printf("invDm %d: %f, %f, %f; %f, %f, %f; %f, %f, %f\n", t,
			//	inv_Dm[t*9+0], inv_Dm[t*9+1], inv_Dm[t*9+2], 
			//	inv_Dm[t*9+3], inv_Dm[t*9+4], inv_Dm[t*9+5], 
			//	inv_Dm[t*9+6], inv_Dm[t*9+7], inv_Dm[t*9+8]);
		}

		//for(int i=0; i<number; i++)
		//	M[i]=1;
		//	printf("M %f: %f\n", M[i], 1/M[i]);

		
		//printf("12370: %ef\n", Vol[12370]);
	}

	void Select(TYPE p[], TYPE q[], int& select_v)
	{
		TYPE dir[3];
		dir[0]=q[0]-p[0];
		dir[1]=q[1]-p[1];
		dir[2]=q[2]-p[2];
		Normalize(dir);

		TYPE min_t=MY_INFINITE;
		int	 select_t;
		for(int t=0; t<t_number; t++)
		{
			TYPE _min_t=MY_INFINITE;
			if(Ray_Triangle_Intersection(&X[T[t*3+0]*3], &X[T[t*3+1]*3], &X[T[t*3+2]*3], p, dir, _min_t) && _min_t>0 && _min_t<min_t)
			{
				select_t = t;
				min_t = _min_t;
			}
		}

		if(min_t!=MY_INFINITE)	//Selection made
		{
			TYPE r;
			TYPE d0=Squared_VE_Distance(&X[T[select_t*3+0]*3], p, q, r);
			TYPE d1=Squared_VE_Distance(&X[T[select_t*3+1]*3], p, q, r);
			TYPE d2=Squared_VE_Distance(&X[T[select_t*3+2]*3], p, q, r);
			if(d0<d1 && d0<d2)	select_v=T[select_t*3+0];
			else if(d1<d2)		select_v=T[select_t*3+1];
			else				select_v=T[select_t*3+2];
		}
	}

};


#endif