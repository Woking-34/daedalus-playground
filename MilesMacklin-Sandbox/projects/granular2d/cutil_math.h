    /*
 * Copyright 1993-2007 NVIDIA Corporation.  All rights reserved.
 *
 * NOTICE TO USER:
 *
 * This source code is subject to NVIDIA ownership rights under U.S. and
 * international Copyright laws.  Users and possessors of this source code
 * are hereby granted a nonexclusive, royalty-free license to use this code
 * in individual and commercial software.
 *
 * NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE
 * CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR
 * IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
 * IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL,
 * OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS,  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION,  ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOURCE CODE.
 *
 * U.S. Government End Users.   This source code is a "commercial item" as
 * that term is defined at  48 C.F.R. 2.101 (OCT 1995), consisting  of
 * "commercial computer  software"  and "commercial computer software
 * documentation" as such terms are  used in 48 C.F.R. 12.212 (SEPT 1995)
 * and is provided to the U.S. Government only as a commercial end item.
 * Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through
 * 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the
 * source code with only those rights set forth herein.
 *
 * Any use of this source code in individual and commercial software must
 * include, in the user documentation and internal comments to the code,
 * the above Disclaimer and U.S. Government End Users Notice.
 */

/*
    This file implements common mathematical operations on vector types
    (float3, float4 etc.) since these are not provided as standard by CUDA.

    The syntax is modelled on the Cg standard library.
*/

#ifndef CUTIL_MATH_H
#define CUTIL_MATH_H

#include "cuda_runtime.h"

////////////////////////////////////////////////////////////////////////////////
typedef unsigned int uint;
typedef unsigned short ushort;

#ifndef __CUDACC__
#include <math.h>

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

inline float fminf(float a, float b)
{
  return a < b ? a : b;
}

inline float fmaxf(float a, float b)
{
  return a < b ? a : b;
}

inline int max(int a, int b)
{
  return a > b ? a : b;
}

inline int min(int a, int b)
{
  return a < b ? a : b;
}
#endif

// float functions
////////////////////////////////////////////////////////////////////////////////

// lerp
inline __device__ __host__ float lerp(float a, float b, float t)
{
    return a + t*(b-a);
}

// clamp
inline __device__ __host__ float clamp(float f, float a, float b)
{
    return fmaxf(a, fminf(f, b));
}

//select
inline __device__ __host__ float select(float a, float b, float c)
{
	return c ? b : a;
}

// int2 functions
////////////////////////////////////////////////////////////////////////////////

// additional constructors
inline __host__ __device__ int2 make_int2(float2 s)
{
	return make_int2(s.x, s.y);
}

inline __host__ __device__ int2 make_int2(int s)
{
	return make_int2(s, s);
}

// negate
inline __host__ __device__ int2 operator-(int2 &a)
{
    return make_int2(-a.x, -a.y);
}

// addition
inline __host__ __device__ int2 operator+(int2 a, int2 b)
{
    return make_int2(a.x + b.x, a.y + b.y);
}
inline __host__ __device__ void operator+=(int2 &a, int2 b)
{
    a.x += b.x; a.y += b.y;
}

// subtract
inline __host__ __device__ int2 operator-(int2 a, int2 b)
{
    return make_int2(a.x - b.x, a.y - b.y);
}
inline __host__ __device__ void operator-=(int2 &a, int2 b)
{
    a.x -= b.x; a.y -= b.y;
}

// multiply
inline __host__ __device__ int2 operator*(int2 a, int2 b)
{
    return make_int2(a.x * b.x, a.y * b.y);
}
inline __host__ __device__ int2 operator*(int2 a, int s)
{
    return make_int2(a.x * s, a.y * s);
}
inline __host__ __device__ int2 operator*(int s, int2 a)
{
    return make_int2(a.x * s, a.y * s);
}
inline __host__ __device__ void operator*=(int2 &a, int s)
{
    a.x *= s; a.y *= s;
}

// div
inline __host__ __device__ int2 operator/(int2 a, int2 b)
{
	return make_int2(a.x / b.x, a.y / b.y);
}

// float2 functions
////////////////////////////////////////////////////////////////////////////////

// additional constructors
inline __host__ __device__ float2 make_float2(float s)
{
    return make_float2(s, s);
}
inline __host__ __device__ float2 make_float2(int2 a)
{
    return make_float2(float(a.x), float(a.y));
}

// min
static __inline__ __host__ __device__ float2 fminf(float2 a, float2 b)
{
	return make_float2(fminf(a.x, b.x), fminf(a.y, b.y));
}

// max
static __inline__ __host__ __device__ float2 fmaxf(float2 a, float2 b)
{
	return make_float2(fmaxf(a.x, b.x), fmaxf(a.y, b.y));
}

// negate
inline __host__ __device__ float2 operator-(float2 &a)
{
    return make_float2(-a.x, -a.y);
}

// addition
inline __host__ __device__ float2 operator+(float2 a, float2 b)
{
    return make_float2(a.x + b.x, a.y + b.y);
}
inline __host__ __device__ void operator+=(float2 &a, float2 b)
{
    a.x += b.x; a.y += b.y;
}

// subtract
inline __host__ __device__ float2 operator-(float2 a, float2 b)
{
    return make_float2(a.x - b.x, a.y - b.y);
}
inline __host__ __device__ void operator-=(float2 &a, float2 b)
{
    a.x -= b.x; a.y -= b.y;
}

// multiply
inline __host__ __device__ float2 operator*(float2 a, float2 b)
{
    return make_float2(a.x * b.x, a.y * b.y);
}
inline __host__ __device__ float2 operator*(float2 a, float s)
{
    return make_float2(a.x * s, a.y * s);
}
inline __host__ __device__ float2 operator*(float s, float2 a)
{
    return make_float2(a.x * s, a.y * s);
}
inline __host__ __device__ void operator*=(float2 &a, float s)
{
    a.x *= s; a.y *= s;
}

// divide
inline __host__ __device__ float2 operator/(float2 a, float2 b)
{
    return make_float2(a.x / b.x, a.y / b.y);
}
inline __host__ __device__ float2 operator/(float2 a, float s)
{
    float inv = 1.0f / s;
    return a * inv;
}
inline __host__ __device__ float2 operator/(float s, float2 a)
{
    float inv = 1.0f / s;
    return a * inv;
}
inline __host__ __device__ void operator/=(float2 &a, float s)
{
    float inv = 1.0f / s;
    a *= inv;
}

// lerp
inline __device__ __host__ float2 lerp(float2 a, float2 b, float t)
{
    return a + t*(b-a);
}

// clamp
inline __device__ __host__ float2 clamp(float2 v, float a, float b)
{
    return make_float2(clamp(v.x, a, b), clamp(v.y, a, b));
}

inline __device__ __host__ float2 clamp(float2 v, float2 a, float2 b)
{
    return make_float2(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y));
}

// dot product
inline __host__ __device__ float dot(float2 a, float2 b)
{ 
    return a.x * b.x + a.y * b.y;
}

// length
inline __host__ __device__ float length(float2 v)
{
    return sqrtf(dot(v, v));
}

// normalize
inline __host__ __device__ float2 normalize(float2 v)
{
    float invLen = 1.0f / sqrtf(dot(v, v));
    return v * invLen;
}

// floor
inline __host__ __device__ float2 floor(const float2 v)
{
    return make_float2(floor(v.x), floor(v.y));
}

// fract
inline __host__ __device__ float2 fract(const float2 v)
{
	return make_float2(v.x - floor(v.x), v.y - floor(v.y));
}

// ceil
inline __host__ __device__ float2 ceil(const float2 v)
{
    return make_float2(ceil(v.x), ceil(v.y));
}

// reflect
inline __host__ __device__ float2 reflect(float2 i, float2 n)
{
	return i - 2.0f * n * dot(n,i);
}

// float3 functions
////////////////////////////////////////////////////////////////////////////////

// additional constructors
inline __host__ __device__ float3 make_float3(float s)
{
    return make_float3(s, s, s);
}
inline __host__ __device__ float3 make_float3(float2 a)
{
    return make_float3(a.x, a.y, 0.0f);
}
inline __host__ __device__ float3 make_float3(float2 a, float s)
{
    return make_float3(a.x, a.y, s);
}
inline __host__ __device__ float3 make_float3(float4 a)
{
    return make_float3(a.x, a.y, a.z);  // discards w
}
inline __host__ __device__ float3 make_float3(int3 a)
{
    return make_float3(float(a.x), float(a.y), float(a.z));
}

// negate
inline __host__ __device__ float3 operator-(float3 &a)
{
    return make_float3(-a.x, -a.y, -a.z);
}

// min
static __inline__ __host__ __device__ float3 fminf(float3 a, float3 b)
{
	return make_float3(fminf(a.x,b.x), fminf(a.y,b.y), fminf(a.z,b.z));
}

// max
static __inline__ __host__ __device__ float3 fmaxf(float3 a, float3 b)
{
	return make_float3(fmaxf(a.x,b.x), fmaxf(a.y,b.y), fmaxf(a.z,b.z));
}

// addition
inline __host__ __device__ float3 operator+(float3 a, float3 b)
{
    return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline __host__ __device__ float3 operator+(float3 a, float b)
{
    return make_float3(a.x + b, a.y + b, a.z + b);
}
inline __host__ __device__ void operator+=(float3 &a, float3 b)
{
    a.x += b.x; a.y += b.y; a.z += b.z;
}

// subtract
inline __host__ __device__ float3 operator-(float3 a, float3 b)
{
    return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}
inline __host__ __device__ float3 operator-(float3 a, float b)
{
    return make_float3(a.x - b, a.y - b, a.z - b);
}
inline __host__ __device__ void operator-=(float3 &a, float3 b)
{
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
}

// multiply
inline __host__ __device__ float3 operator*(float3 a, float3 b)
{
    return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);
}
inline __host__ __device__ float3 operator*(float3 a, float s)
{
    return make_float3(a.x * s, a.y * s, a.z * s);
}
inline __host__ __device__ float3 operator*(float s, float3 a)
{
    return make_float3(a.x * s, a.y * s, a.z * s);
}
inline __host__ __device__ void operator*=(float3 &a, float s)
{
    a.x *= s; a.y *= s; a.z *= s;
}

// divide
inline __host__ __device__ float3 operator/(float3 a, float3 b)
{
    return make_float3(a.x / b.x, a.y / b.y, a.z / b.z);
}
inline __host__ __device__ float3 operator/(float3 a, float s)
{
    float inv = 1.0f / s;
    return a * inv;
}
inline __host__ __device__ float3 operator/(float s, float3 a)
{
    float inv = 1.0f / s;
    return a * inv;
}
inline __host__ __device__ void operator/=(float3 &a, float s)
{
    float inv = 1.0f / s;
    a *= inv;
}

// lerp
inline __device__ __host__ float3 lerp(float3 a, float3 b, float t)
{
    return a + t*(b-a);
}

// clamp
inline __device__ __host__ float3 clamp(float3 v, float a, float b)
{
    return make_float3(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b));
}

inline __device__ __host__ float3 clamp(float3 v, float3 a, float3 b)
{
    return make_float3(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z));
}

// dot product
inline __host__ __device__ float dot(float3 a, float3 b)
{ 
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// cross product
inline __host__ __device__ float3 cross(float3 a, float3 b)
{ 
    return make_float3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x); 
}

// length
inline __host__ __device__ float length(float3 v)
{
    return sqrtf(dot(v, v));
}

// normalize
inline __host__ __device__ float3 normalize(float3 v)
{
    float invLen = 1.0f / sqrtf(dot(v, v));
    return v * invLen;
}

// floor
inline __host__ __device__ float3 floor(const float3 v)
{
    return make_float3(floor(v.x), floor(v.y), floor(v.z));
}

// reflect
inline __host__ __device__ float3 reflect(float3 i, float3 n)
{
	return i - 2.0f * n * dot(n,i);
}

// float4 functions
////////////////////////////////////////////////////////////////////////////////

// additional constructors
inline __host__ __device__ float4 make_float4(float s)
{
    return make_float4(s, s, s, s);
}
inline __host__ __device__ float4 make_float4(float3 a)
{
    return make_float4(a.x, a.y, a.z, 0.0f);
}
inline __host__ __device__ float4 make_float4(float3 a, float w)
{
    return make_float4(a.x, a.y, a.z, w);
}
inline __host__ __device__ float4 make_float4(int4 a)
{
    return make_float4(float(a.x), float(a.y), float(a.z), float(a.w));
}
inline __host__ __device__ float4 make_float4(uchar4 a)
{
	return make_float4(float(a.x), float(a.y), float(a.z), float(a.w));
}

// negate
inline __host__ __device__ float4 operator-(float4 &a)
{
    return make_float4(-a.x, -a.y, -a.z, -a.w);
}

// min
static __inline__ __host__ __device__ float4 fminf(float4 a, float4 b)
{
	return make_float4(fminf(a.x,b.x), fminf(a.y,b.y), fminf(a.z,b.z), fminf(a.w,b.w));
}

// max
static __inline__ __host__ __device__ float4 fmaxf(float4 a, float4 b)
{
	return make_float4(fmaxf(a.x,b.x), fmaxf(a.y,b.y), fmaxf(a.z,b.z), fmaxf(a.w,b.w));
}

// addition
inline __host__ __device__ float4 operator+(float4 a, float4 b)
{
    return make_float4(a.x + b.x, a.y + b.y, a.z + b.z,  a.w + b.w);
}
inline __host__ __device__ void operator+=(float4 &a, float4 b)
{
    a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
}

// subtract
inline __host__ __device__ float4 operator-(float4 a, float4 b)
{
    return make_float4(a.x - b.x, a.y - b.y, a.z - b.z,  a.w - b.w);
}
inline __host__ __device__ void operator-=(float4 &a, float4 b)
{
    a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;
}

// multiply
inline __host__ __device__ float4 operator*(float4 a, float4 b)
{
    return make_float4(a.x * b.x, a.y * b.y, a.z * b.z,  a.w * b.w);
}
inline __host__ __device__ float4 operator*(float4 a, float s)
{
    return make_float4(a.x * s, a.y * s, a.z * s, a.w * s);
}
inline __host__ __device__ float4 operator*(float s, float4 a)
{
    return make_float4(a.x * s, a.y * s, a.z * s, a.w * s);
}
inline __host__ __device__ void operator*=(float4 &a, float s)
{
    a.x *= s; a.y *= s; a.z *= s; a.w *= s;
}

// divide
inline __host__ __device__ float4 operator/(float4 a, float4 b)
{
    return make_float4(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w);
}
inline __host__ __device__ float4 operator/(float4 a, float s)
{
    //float inv = 1.0f / s;
    return a / make_float4(s);
}
inline __host__ __device__ float4 operator/(float s, float4 a)
{
    //float inv = 1.0f / s;
	return make_float4(s) / a;
}

inline __host__ __device__ void operator/=(float4 &a, float s)
{
    //float inv = 1.0f / s;
	a = a / make_float4(s);
}

// lerp
inline __device__ __host__ float4 lerp(float4 a, float4 b, float t)
{
    return a + t*(b-a);
}

// clamp
inline __device__ __host__ float4 clamp(float4 v, float a, float b)
{
    return make_float4(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b), clamp(v.w, a, b));
}

inline __device__ __host__ float4 clamp(float4 v, float4 a, float4 b)
{
    return make_float4(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z), clamp(v.w, a.w, b.w));
}

// dot product
inline __host__ __device__ float dot(float4 a, float4 b)
{ 
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

// cross product
inline __host__ __device__ float4 cross(float4 a, float4 b)
{
	return make_float4(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x, 0.0f);
}

// length
inline __host__ __device__ float length(float4 r)
{
    return sqrtf(dot(r, r));
}

// normalize
inline __host__ __device__ float4 normalize(float4 v)
{
    float invLen = 1.0f / sqrtf(dot(v, v));
    return v * invLen;
}

// floor
inline __host__ __device__ float4 floor(const float4 v)
{
    return make_float4(floor(v.x), floor(v.y), floor(v.z), floor(v.w));
}

// exp
inline __host__ __device__ float4 expf(const float4 v)
{
	return make_float4(expf(v.x), expf(v.y), expf(v.z), expf(v.w));
}

// sqrt
inline __host__ __device__ float4 sqrtf(const float4 v)
{
	return make_float4(sqrtf(v.x), sqrtf(v.y), sqrtf(v.z), sqrtf(v.w));
}

// powf
inline __host__ __device__ float4 powf(const float4 v, float e)
{
	return make_float4(powf(v.x, e), powf(v.y, e), powf(v.z, e), powf(v.w, e));
}

// int3 functions
////////////////////////////////////////////////////////////////////////////////

// additional constructors
inline __host__ __device__ int3 make_int3(int s)
{
    return make_int3(s, s, s);
}
inline __host__ __device__ int3 make_int3(float3 a)
{
    return make_int3(int(a.x), int(a.y), int(a.z));
}

// negate
inline __host__ __device__ int3 operator-(int3 &a)
{
    return make_int3(-a.x, -a.y, -a.z);
}

// min
inline __host__ __device__ int3 min(int3 a, int3 b)
{
    return make_int3(min(a.x,b.x), min(a.y,b.y), min(a.z,b.z));
}

// max
inline __host__ __device__ int3 max(int3 a, int3 b)
{
    return make_int3(max(a.x,b.x), max(a.y,b.y), max(a.z,b.z));
}

// addition
inline __host__ __device__ int3 operator+(int3 a, int3 b)
{
    return make_int3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline __host__ __device__ void operator+=(int3 &a, int3 b)
{
    a.x += b.x; a.y += b.y; a.z += b.z;
}

// subtract
inline __host__ __device__ int3 operator-(int3 a, int3 b)
{
    return make_int3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline __host__ __device__ void operator-=(int3 &a, int3 b)
{
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
}

// multiply
inline __host__ __device__ int3 operator*(int3 a, int3 b)
{
    return make_int3(a.x * b.x, a.y * b.y, a.z * b.z);
}
inline __host__ __device__ int3 operator*(int3 a, int s)
{
    return make_int3(a.x * s, a.y * s, a.z * s);
}
inline __host__ __device__ int3 operator*(int s, int3 a)
{
    return make_int3(a.x * s, a.y * s, a.z * s);
}
inline __host__ __device__ void operator*=(int3 &a, int s)
{
    a.x *= s; a.y *= s; a.z *= s;
}

// divide
inline __host__ __device__ int3 operator/(int3 a, int3 b)
{
    return make_int3(a.x / b.x, a.y / b.y, a.z / b.z);
}
inline __host__ __device__ int3 operator/(int3 a, int s)
{
    return make_int3(a.x / s, a.y / s, a.z / s);
}
inline __host__ __device__ int3 operator/(int s, int3 a)
{
    return make_int3(a.x / s, a.y / s, a.z / s);
}
inline __host__ __device__ void operator/=(int3 &a, int s)
{
    a.x /= s; a.y /= s; a.z /= s;
}

// clamp
inline __device__ __host__ int clamp(int f, int a, int b)
{
    return max(a, min(f, b));
}

inline __device__ __host__ int3 clamp(int3 v, int a, int b)
{
    return make_int3(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b));
}

inline __device__ __host__ int3 clamp(int3 v, int3 a, int3 b)
{
    return make_int3(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z));
}


// int4 functions
////////////////////////////////////////////////////////////////////////////////

// additional constructors
inline __host__ __device__ int4 make_int4(int s)
{
	return make_int4(s, s, s, s);
}
inline __host__ __device__ int4 make_int4(float4 a)
{
	return make_int4(int(a.x), int(a.y), int(a.z), int(a.w));
}
inline __host__ __device__ int4 make_int4(uchar4 a)
{
	return make_int4(int(a.x), int(a.y), int(a.z), int(a.w));
}
inline __host__ __device__ int4 make_int4(uint4 a)
{
	return make_int4(int(a.x), int(a.y), int(a.z), int(a.w));
}

// negate
inline __host__ __device__ int4 operator-(int4 &a)
{
	return make_int4(-a.x, -a.y, -a.z, -a.w);
}

// min
inline __host__ __device__ int4 min(int4 a, int4 b)
{
	return make_int4(min(a.x, b.x), min(a.y, b.y), min(a.z, b.z), min(a.w, b.w));
}

// max
inline __host__ __device__ int4 max(int4 a, int4 b)
{
	return make_int4(max(a.x, b.x), max(a.y, b.y), max(a.z, b.z), max(a.w, b.w));
}

// addition
inline __host__ __device__ int4 operator+(int4 a, int4 b)
{
	return make_int4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}
inline __host__ __device__ void operator+=(int4 &a, int4 b)
{
	a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
}

// subtract
inline __host__ __device__ int4 operator-(int4 a, int4 b)
{
	return make_int4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

inline __host__ __device__ void operator-=(int4 &a, int4 b)
{
	a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;
}

// multiply
inline __host__ __device__ int4 operator*(int4 a, int4 b)
{
	return make_int4(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);
}
inline __host__ __device__ int4 operator*(int4 a, int s)
{
	return make_int4(a.x * s, a.y * s, a.z * s, a.w * s);
}
inline __host__ __device__ int4 operator*(int s, int4 a)
{
	return make_int4(a.x * s, a.y * s, a.z * s, a.w * s);
}
inline __host__ __device__ void operator*=(int4 &a, int s)
{
	a.x *= s; a.y *= s; a.z *= s; a.w *= s;
}

// divide
inline __host__ __device__ int4 operator/(int4 a, int4 b)
{
	return make_int4(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w);
}
inline __host__ __device__ int4 operator/(int4 a, int s)
{
	return make_int4(a.x / s, a.y / s, a.z / s, a.w / s);
}
inline __host__ __device__ int4 operator/(int s, int4 a)
{
	return make_int4(a.x / s, a.y / s, a.z / s, a.w / s);
}
inline __host__ __device__ void operator/=(int4 &a, int s)
{
	a.x /= s; a.y /= s; a.z /= s; a.w /= s;
}

// clamp
inline __device__ __host__ int4 clamp(int4 v, int a, int b)
{
	return make_int4(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b), clamp(v.w, a, b));
}

inline __device__ __host__ int4 clamp(int4 v, int4 a, int4 b)
{
	return make_int4(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z), clamp(v.w, a.w, b.w));
}

// uint2 functions
////////////////////////////////////////////////////////////////////////////////

// additional constructors
inline __host__ __device__ uint2 make_uint2(uint s)
{
	return make_uint2(s, s);
}

inline __host__ __device__ uint2 operator&(uint2 v0, uint2 v1)
{
	return make_uint2(v0.x & v1.x, v0.y & v1.y);
}

inline __host__ __device__ uint2 operator|(uint2 v0, uint2 v1)
{
	return make_uint2(v0.x | v1.x, v0.y | v1.y);
}

inline __host__ __device__ uint2 operator>>(uint2 v0, uint2 shift)
{
	return make_uint2(v0.x >> shift.x, v0.y >> shift.y);
}

inline __host__ __device__ uint2 operator<<(uint2 v0, uint2 shift)
{
	return make_uint2(v0.x << shift.x, v0.y << shift.y);
}

inline __host__ __device__ uint2 operator+(uint2 v0, uint2 v1)
{
	return make_uint2(v0.x + v1.x, v0.y + v1.y);
}

inline __host__ __device__ uint2 operator-(uint2 v0, uint2 v1)
{
	return make_uint2(v0.x - v1.x, v0.y - v1.y);
}

inline __host__ __device__ uint2 operator*(uint2 v0, uint2 v1)
{
	return make_uint2(v0.x * v1.x, v0.y * v1.y);
}

inline __host__ __device__ uint2 operator/(uint2 v0, uint2 v1)
{
	return make_uint2(v0.x / v1.x, v0.y / v1.y);
}

inline __host__ __device__ float __int_as_float(int a)
{
	union { int a; float b; } u;
	u.a = a;
	return u.b;
}

inline __host__ __device__ float2 as_float2(uint2 v0)
{
	return make_float2(__int_as_float((int)v0.x), __int_as_float((int)v0.y));
}


// uint3 functions
////////////////////////////////////////////////////////////////////////////////

// additional constructors
inline __host__ __device__ uint3 make_uint3(uint s)
{
    return make_uint3(s, s, s);
}
inline __host__ __device__ uint3 make_uint3(float3 a)
{
    return make_uint3(uint(a.x), uint(a.y), uint(a.z));
}

// min
inline __host__ __device__ uint3 min(uint3 a, uint3 b)
{
    return make_uint3(min(a.x,b.x), min(a.y,b.y), min(a.z,b.z));
}

// max
inline __host__ __device__ uint3 max(uint3 a, uint3 b)
{
    return make_uint3(max(a.x,b.x), max(a.y,b.y), max(a.z,b.z));
}

// addition
inline __host__ __device__ uint3 operator+(uint3 a, uint3 b)
{
    return make_uint3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline __host__ __device__ void operator+=(uint3 &a, uint3 b)
{
    a.x += b.x; a.y += b.y; a.z += b.z;
}

// subtract
inline __host__ __device__ uint3 operator-(uint3 a, uint3 b)
{
    return make_uint3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline __host__ __device__ void operator-=(uint3 &a, uint3 b)
{
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
}

// multiply
inline __host__ __device__ uint3 operator*(uint3 a, uint3 b)
{
    return make_uint3(a.x * b.x, a.y * b.y, a.z * b.z);
}
inline __host__ __device__ uint3 operator*(uint3 a, uint s)
{
    return make_uint3(a.x * s, a.y * s, a.z * s);
}
inline __host__ __device__ uint3 operator*(uint s, uint3 a)
{
    return make_uint3(a.x * s, a.y * s, a.z * s);
}
inline __host__ __device__ void operator*=(uint3 &a, uint s)
{
    a.x *= s; a.y *= s; a.z *= s;
}

// divide
inline __host__ __device__ uint3 operator/(uint3 a, uint3 b)
{
    return make_uint3(a.x / b.x, a.y / b.y, a.z / b.z);
}
inline __host__ __device__ uint3 operator/(uint3 a, uint s)
{
    return make_uint3(a.x / s, a.y / s, a.z / s);
}
inline __host__ __device__ uint3 operator/(uint s, uint3 a)
{
    return make_uint3(a.x / s, a.y / s, a.z / s);
}
inline __host__ __device__ void operator/=(uint3 &a, uint s)
{
    a.x /= s; a.y /= s; a.z /= s;
}

// clamp
inline __device__ __host__ uint clamp(uint f, uint a, uint b)
{
    return max(a, min(f, b));
}

inline __device__ __host__ uint3 clamp(uint3 v, uint a, uint b)
{
    return make_uint3(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b));
}

inline __device__ __host__ uint3 clamp(uint3 v, uint3 a, uint3 b)
{
    return make_uint3(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z));
}

// uchar4 functions
////////////////////////////////////////////////////////////////////////////////

// additional constructors
inline __host__ __device__ uchar4 make_uchar4(int s)
{
	return make_uchar4((unsigned char)s, (unsigned char)s, (unsigned char)s, (unsigned char)s);
}

inline __host__ __device__ uchar4 make_uchar4(unsigned char s)
{
	return make_uchar4(s, s, s, s);
}

inline __host__ __device__ uchar4 make_uchar4(int4 s)
{
	return make_uchar4((unsigned char)s.x, (unsigned char)s.y, (unsigned char)s.z, (unsigned char)s.w);
}

inline __host__ __device__ uchar4 make_uchar4(float4 s)
{
	return make_uchar4((unsigned char)s.x, (unsigned char)s.y, (unsigned char)s.z, (unsigned char)s.w);
}

inline __host__ __device__ uchar4 operator&(uchar4 v0, uchar4 v1)
{
	return make_uchar4(v0.x & v1.x, v0.y & v1.y, v0.z & v1.z, v0.w & v1.w);
}

inline __host__ __device__ uchar4 operator*(uchar4 v0, unsigned char s)
{
	return make_uchar4(v0.x * s, v0.y * s, v0.z * s, v0.w * s);
}

inline __host__ __device__ uchar4 operator-(uchar4 v0, uchar4 v1)
{
	return make_uchar4(v0.x - v1.x, v0.y - v1.y, v0.z - v1.z, v0.w - v1.w);
}

inline __host__ __device__ uchar4 operator+(uchar4 v0, uchar4 v1)
{
	return make_uchar4(v0.x + v1.x, v0.y + v1.y, v0.z + v1.z, v0.w + v1.w);
}


#endif
