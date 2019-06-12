#include <cutil_math.h>

__constant__ unsigned char c_perm_3d[256];
__shared__ unsigned char s_perm_3d[256]; // shared memory copy of permuation array
//unsigned char* d_perm_parts=NULL; // global memory copy of permutation array
// host version of permutation array
const static unsigned char h_perm[] = { 151, 160, 137, 91, 90, 15,
131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23,
190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33,
88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166,
77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244,
102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196,
135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123,
5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42,
223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9,
129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, 228,
251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107,
49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254,
138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180
};

__device__ inline int perm(int i) { return(c_perm_3d[i & 0xff]); }
__device__ inline float fade(float t) { return t * t * t * (t * (t * 6.f - 15.f) + 10.f); }
__device__ inline float lerpP(float t, float a, float b) { return a + t * (b - a); }
__device__ inline float grad(int hash, float x, float y, float z) {
	int h = hash & 15;                      // CONVERT LO 4 BITS OF HASH CODE
	float u = h<8 ? x : y,                 // INTO 12 GRADIENT DIRECTIONS.
		v = h<4 ? y : h == 12 || h == 14 ? x : z;
	return ((h & 1) == 0 ? u : -u) + ((h & 2) == 0 ? v : -v);
}

__device__ float inoise_parts2(float x, float y, float z) {
	int X = ((int)floorf(x)) & 255, // FIND UNIT CUBE THAT
		Y = ((int)floorf(y)) & 255,   // CONTAINS POINT.
		Z = ((int)floorf(z)) & 255;
	x -= floorf(x);               // FIND RELATIVE X,Y,Z
	y -= floorf(y);               // OF POINT IN CUBE.
	z -= floorf(z);
	float u = fade(x),            // COMPUTE FADE CURVES
		v = fade(y),                // FOR EACH OF X,Y,Z.
		w = fade(z);
	int A = perm(X) + Y, AA = perm(A) + Z, AB = perm(A + 1) + Z, // HASH COORDINATES OF
		B = perm(X + 1) + Y, BA = perm(B) + Z, BB = perm(B + 1) + Z; // THE 8 CUBE CORNERS,

	return lerpP(w, lerpP(v, lerpP(u, grad(perm(AA), x, y, z), // AND ADD
		grad(perm(BA), x - 1.f, y, z)),   // BLENDED
		lerpP(u, grad(perm(AB), x, y - 1.f, z),    // RESULTS
		grad(perm(BB), x - 1.f, y - 1.f, z))),     // FROM  8
		lerpP(v, lerpP(u, grad(perm(AA + 1), x, y, z - 1.f),  // CORNERS
		grad(perm(BA + 1), x - 1.f, y, z - 1.f)),    // OF CUBE
		lerpP(u, grad(perm(AB + 1), x, y - 1.f, z - 1.f),
		grad(perm(BB + 1), x - 1.f, y - 1.f, z - 1.f))));
#ifdef ORIG
	return(perm(X));
#endif

}

__device__ inline float noise1D(float x, float y, float z, int octaves,
	float lacunarity, float gain, float freq, float amp)
{
	float sum = 0.f;
	for (int i = 0; i<octaves; i++) {
		sum += inoise_parts2(x*freq, y*freq, z*freq)*amp;
		freq *= lacunarity;
		amp *= gain;
	}
	return sum;
}

__global__ void dampVelKernel(float3* vel, float damping, float dt, int nParts)
{
    unsigned int x = blockIdx.x*blockDim.x + threadIdx.x;

	int n = x;

	if (n<nParts) {

		vel[n] = (1 - damping*dt)*vel[n];

	}
}

__global__ void addGravityKernel(float3* vel, float3 gravityDir, float gravityStrength, float dt, int nParts)
{
    unsigned int x = blockIdx.x*blockDim.x + threadIdx.x;

	int n = x;

	if (n<nParts) {
		vel[n] = vel[n] + gravityDir*gravityStrength*dt;
	}
}

__global__ void addTurbulenceKernel(float3* vel, float3* pos,
								  float3 noiseAmp, float3 noiseOffset, int noiseOct, float noiseLac, float noiseFreq,
								  float dt, int nParts)
{
    unsigned int x = blockIdx.x*blockDim.x + threadIdx.x;

	int n = x;

	if (n<nParts) {

		float3 currVel = vel[n];
		
		float3 noise = make_float3(0,0,0);
		if (noiseAmp.x != 0)
			noise.x = noiseAmp.x*noise1D(pos[n].x+noiseOffset.x, pos[n].y+noiseOffset.y, pos[n].z+noiseOffset.z,
										noiseOct, noiseLac, 0.5, noiseFreq,1);
		if (noiseAmp.y != 0)
			noise.y  = noiseAmp.y*noise1D(pos[n].x+noiseOffset.x+2000, pos[n].y+noiseOffset.y, pos[n].z+noiseOffset.z,
										noiseOct, noiseLac, 0.5, noiseFreq,1);
		if (noiseAmp.x != 0)
			noise.z  += noiseAmp.z*noise1D(pos[n].x+noiseOffset.x+5000, pos[n].y+noiseOffset.y, pos[n].z+noiseOffset.z,
										noiseOct, noiseLac, 0.5, noiseFreq,1);

		vel[n] = vel[n] + noise*dt;

	}
}

__global__ void integrateParticlesKernel(float3* pos, float3* vel, float* age, float* life,
								  float4* colour, float opacity, float3 col1, float3 col2, float dt, int nParts)
{
    unsigned int x = blockIdx.x*blockDim.x + threadIdx.x;

	int n = x;

	if (n<nParts) {

		pos[n] = pos[n] + vel[n]*dt;

		age[n] = age[n] + dt;
		if (age[n] > life[n])
			age[n] = life[n];

		float ageNorm = age[n]/life[n];

		float3 col = lerp(col1,col2,ageNorm);
		float alpha = opacity*(1-pow(age[n]/life[n],2));

		colour[n] = make_float4(col.x,col.y,col.z,alpha);


	}
}

__global__ void initNewParticlesKernel(float3* pos, float3* vel, float* age, float* life,
										float3 initPos, float3 initVel, float radVelAmp,
										float3 noiseVelAmp, float3 noiseVelOffset, int noiseVelOct, float noiseVelLac, float noiseVelFreq,
										float initLife, float time, int nParts)
{
    unsigned int x = blockIdx.x*blockDim.x + threadIdx.x;

	int n = x;

	//__device__ inline float noise1D(float x, float y, float z, int octaves,
	//	     float lacunarity, float gain, float freq, float amp)

	if (n<nParts) {
						
		float3 radVel = radVelAmp * (pos[n] - initPos);	
		
		vel[n] = initVel + radVel;
							
		if (noiseVelAmp.x != 0)
			vel[n].x += noiseVelAmp.x*noise1D(pos[n].x+noiseVelOffset.x, pos[n].y+noiseVelOffset.y, pos[n].z+noiseVelOffset.z,
										noiseVelOct, noiseVelLac, 0.5, noiseVelFreq,1);
			
		if (noiseVelAmp.y != 0)
			vel[n].y  += noiseVelAmp.y*noise1D(pos[n].x+noiseVelOffset.x+2000, pos[n].y+noiseVelOffset.y, pos[n].z+noiseVelOffset.z,
										noiseVelOct, noiseVelLac, 0.5, noiseVelFreq,1);
	
		if (noiseVelAmp.x != 0)
			vel[n].z  += noiseVelAmp.z*noise1D(pos[n].x+noiseVelOffset.x+5000, pos[n].y+noiseVelOffset.y, pos[n].z+noiseVelOffset.z,
										noiseVelOct, noiseVelLac, 0.5, noiseVelFreq,1);
	
		age[n] = 0.0;
	}

	
}

__global__ void resetParticlesKernel(float3* pos, float3* vel, float* age, float* life, int nParts)
{
    unsigned int x = blockIdx.x*blockDim.x + threadIdx.x;

	int n = x;

	if (n<nParts) {
		pos[n] = make_float3(0.0,0.0,0.0);
		vel[n] =  make_float3(0.0,0.0,0.0);
		age[n] = 1.0;
		life[n] = 1.0;
	}
}


extern "C" void dampVelCu(float3* vel, float damping, float dt, int nParts){

	int nthreads = min(256, nParts);
	int nBlocks = nParts/nthreads + (!(nParts%nthreads)?0:1);
    dim3 blocks(nBlocks, 1,1); dim3 threads(nthreads, 1, 1);

	dampVelKernel<<< blocks, threads>>>(vel, damping, dt, nParts);
}

extern "C" void addGravityCu(float3* vel, float3 gravityDir, float gravityStrength, float dt, int nParts){

	int nthreads = min(256, nParts);
	int nBlocks = nParts/nthreads + (!(nParts%nthreads)?0:1);
    dim3 blocks(nBlocks, 1,1); dim3 threads(nthreads, 1, 1);

	addGravityKernel<<< blocks, threads>>>(vel, gravityDir, gravityStrength, dt, nParts);
}

extern "C" void addTurbulenceCu(float3* vel, float3* pos, float3 noiseAmp, float3 noiseOffset,
								int noiseOct, float noiseLac, float noiseFreq, float dt, int nParts){

	int nthreads = min(256, nParts);
	int nBlocks = nParts/nthreads + (!(nParts%nthreads)?0:1);
    dim3 blocks(nBlocks, 1,1); dim3 threads(nthreads, 1, 1);

	addTurbulenceKernel<<< blocks, threads>>>(vel, pos, noiseAmp, noiseOffset, noiseOct, noiseLac, noiseFreq, dt, nParts);
}

extern "C" void integrateParticlesCu(float3* pos, float3* vel, float* age, float* life,
								  float4* colour, float opacity, float3 col1, float3 col2,
									float dt, int nParts){

	int nthreads = min(256, nParts);
	int nBlocks = nParts/nthreads + (!(nParts%nthreads)?0:1);
    dim3 blocks(nBlocks, 1,1); dim3 threads(nthreads, 1, 1);

	integrateParticlesKernel<<< blocks, threads>>>(pos, vel, age, life, colour, opacity, col1, col2, dt, nParts);

}


extern "C" void resetParticlesCu(float3* pos, float3* vel, float* age, float* life, int nParts) {

	int nthreads = min(256, nParts);

	int nBlocks = nParts/nthreads + (!(nParts%nthreads)?0:1);

    dim3 blocks(nBlocks, 1,1);
    dim3 threads(nthreads, 1, 1);

    resetParticlesKernel<<< blocks, threads>>>(pos, vel, age, life, nParts);
}

extern "C" void initNewParticlesCu(float3* pos, float3* vel, float* age, float* life,
									float3 initPos, float3 initVel, float radVelAmp, 
									float3 noiseVelAmp, float3 noiseVelOffset, int noiseVelOct, float noiseVelLac, float noiseVelFreq,
									float initLife, float time, int nParts) {

	int nthreads = min(256, nParts);

	int nBlocks = nParts/nthreads + (!(nParts%nthreads)?0:1);

    dim3 blocks(nBlocks, 1,1);
    dim3 threads(nthreads, 1, 1);

	cudaMemcpyToSymbol(c_perm_3d, h_perm, sizeof(h_perm),0,cudaMemcpyHostToDevice );

    initNewParticlesKernel<<< blocks, threads>>>(pos, vel, age, life, initPos, initVel, radVelAmp,
												noiseVelAmp, noiseVelOffset, noiseVelOct, noiseVelLac, noiseVelFreq,
												initLife, time, nParts);
}