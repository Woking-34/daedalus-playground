#pragma once

#ifdef __CUDACC__
#include "cutil_math.h"
//#include <cutil_inline.h>    // includes cuda.h and cuda_runtime_api.h

typedef unsigned int       uint32_t;
#else
typedef Vec2 float2;
typedef Vec3 float3;
#endif


struct GrainSystem;

struct GrainParams
{
	float2 mGravity;
	float mDamp;
	
	float mBaumgarte;
	float mFriction;
	float mRestitution;
	float mOverlap;

	float3 mPlanes[8];
	int mNumPlanes;
};

struct GrainTimers
{
	GrainTimers() 
		: mCreateCellIndices(0.0f)
		, mSortCellIndices(0.0f)
		, mCreateGrid(0.0f)
		, mCollide(0.0f)
		, mIntegrate(0.0f)
		, mReorder(0.0f)
	{
	}

	float mCreateCellIndices;
	float mSortCellIndices;
	float mCreateGrid;
	float mCollide;
	float mIntegrate;
	float mReorder;
};

GrainSystem* grainCreateSystem(int numGrains);
void grainDestroySystem(GrainSystem* s);

void grainSetSprings(GrainSystem* s, const uint32_t* indices, const float* restLengths, uint32_t numSprings);

void grainSetPositions(GrainSystem* s, float* p, int n);
void grainSetVelocities(GrainSystem* s, float* v, int n);

void grainGetPositions(GrainSystem* s, float* p);
void grainGetVelocities(GrainSystem* s, float* v);

void grainSetRadii(GrainSystem* s, float* r);
void grainGetRadii(GrainSystem* s, float r);

void grainGetMass(GrainSystem* s, float* r);

void grainSetParams(GrainSystem* s, GrainParams* params);

void grainUpdateSystem(GrainSystem* s, float dt, int iterations, GrainTimers* timers);
