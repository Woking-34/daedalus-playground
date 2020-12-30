#ifndef _MATH_3D_
#define _MATH_3D_

#include <math.h>
#include <memory.h>

typedef float M3DVector3f[3]; // Vector of three floats (x, y, z)
typedef float M3DMatrix33f[9]; // 3x3 Matrix of floats, see below
/*
0 3 6
1 4 7
2 5 8
*/
typedef float M3DMatrix44f[16]; // 4x4 Matrix of floats, see below
/*
0 4  8 12
1 5  9 13
2 6 10 14
3 7 11 15
*/

// Copies vector src to vector dst
void M3DCopyVector3f(M3DVector3f dst, M3DVector3f src);

// Calculates c = a x b
void M3DCrossProduct33f(M3DVector3f c, M3DVector3f a, M3DVector3f b);

// Gets the 3x3 matrix that rotates a point around a line that passes through the origin
// the angle is given in radian
void M3DGetRotationMatrix33f(M3DMatrix33f m, float angle, float x, float y, float z);

// Gets the 4x4 matrix that rotates a point around a line that passes through the origin
// the angle is given in radian
void M3DGetRotationMatrix44f(M3DMatrix44f m, float angle, float x, float y, float z);

// Loads the identify 3x3 matrix of floats
void M3DLoadIdentityMatrix33f(M3DMatrix33f m);

// Loads the identify 4x4 matrix of floats
void M3DLoadIdentityMatrix44f(M3DMatrix44f m);

// Multiply a 4x4 matrix to a 3x3 vector
// v = mu
void M3DMultiplyMatrix44Vector3f(M3DVector3f v, M3DMatrix44f m, M3DVector3f u);

// Multiply a 3x3 matrix to a 3x3 vector
// v = mu
void M3DMultiplyMatrix33Vector3f(M3DVector3f v, M3DMatrix44f m, M3DVector3f u);

void M3DInvertMatrix44f(M3DMatrix44f a, M3DMatrix44f b);

#endif
