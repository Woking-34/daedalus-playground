#include "math3d.h"

// Copies vector src to vector dst
void M3DCopyVector3f(M3DVector3f dst, M3DVector3f src)
{
  memcpy(dst, src, sizeof(M3DVector3f));
}

// Calculates c = a x b
void M3DCrossProduct33f(M3DVector3f c, M3DVector3f a, M3DVector3f b)
{
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

// Gets the matrix that rotates a point around a line that passes through the origin
// the angle is given in radian
void M3DGetRotationMatrix33f(M3DMatrix33f m, float angle, float x, float y, float z)
{
  float length = sqrt(x * x + y * y + z * z);
  if (length == 0.0f)
  {
    M3DLoadIdentityMatrix33f(m);
    return;
  }

  // Normalize
  x /= length;
  y /= length;
  z /= length;

  float s = (float) sin(angle);
  float c = (float) cos(angle);
  float xx = x * x;
  float yy = y * y;
  float zz = z * z;
  float xy = x * y;
  float yz = y * z;
  float zx = z * x;
  float xs = x * s;
  float ys = y * s;
  float zs = z * s;
  float oc = 1.0f - c;

  m[0] = oc * xx + c;
  m[3] = oc * xy - zs;
  m[6] = oc * zx + ys;

  m[1] = oc * xy + zs;
  m[4] = oc * yy + c;
  m[7] = oc * yz - xs;

  m[2] = oc * zx - ys;
  m[5] = oc * yz + xs;
  m[8] = oc * zz + c;
}

// Gets the 4x4 matrix that rotates a point around a line that passes through the origin
// the angle is given in radian
void M3DGetRotationMatrix44f(M3DMatrix44f m, float angle, float x, float y, float z)
{
  float length = sqrt(x * x + y * y + z * z);
  if (length == 0.0f)
  {
    M3DLoadIdentityMatrix44f(m);
    return;
  }

  // Normalize
  x /= length;
  y /= length;
  z /= length;

  float s = (float) sin(angle);
  float c = (float) cos(angle);
  float xx = x * x;
  float yy = y * y;
  float zz = z * z;
  float xy = x * y;
  float yz = y * z;
  float zx = z * x;
  float xs = x * s;
  float ys = y * s;
  float zs = z * s;
  float oc = 1.0f - c;

  m[0] = oc * xx + c;
  m[4] = oc * xy - zs;
  m[8] = oc * zx + ys;
  m[12] = 0.0f;

  m[1] = oc * xy + zs;
  m[5] = oc * yy + c;
  m[9] = oc * yz - xs;
  m[13] = 0.0f;

  m[2] = oc * zx - ys;
  m[6] = oc * yz + xs;
  m[10] = oc * zz + c;
  m[14] = 0.0f;

  m[3] = 0.0f;
  m[7] = 0.0f;
  m[11] = 0.0f;
  m[15] = 1.0f;
}

// Loads the identify 3x3 matrix of floats
void M3DLoadIdentityMatrix33f(M3DMatrix33f m)
{
  static M3DMatrix33f identity = { 1.0f, 0.0f, 0.0f,
                                   0.0f, 1.0f, 0.0f,
                                   0.0f, 0.0f, 1.0f
                                 };
  memcpy(m, identity, sizeof(M3DMatrix33f));
}

// Loads the identify 4x4 matrix of floats
void M3DLoadIdentityMatrix44f(M3DMatrix44f m)
{
  static M3DMatrix44f identity = { 1.0f, 0.0f, 0.0f, 0.0f,
                                   0.0f, 1.0f, 0.0f, 0.0f,
                                   0.0f, 0.0f, 1.0f, 0.0f,
                                   0.0f, 0.0f, 0.0f, 1.0f
                                 };
  memcpy(m, identity, sizeof(M3DMatrix44f));
}

// Multiply a 4x4 matrix to a 3x3 vector
// v = mu
void M3DMultiplyMatrix44Vector3f(M3DVector3f v, M3DMatrix44f m, M3DVector3f u)
{
  v[0] = m[0] * u[0] + m[4] * u[1] + m[8] * u[2];
  v[1] = m[1] * u[0] + m[5] * u[1] + m[9] * u[2];
  v[2] = m[2] * u[0] + m[6] * u[1] + m[10] * u[2];
}

// Multiply a 3x3 matrix to a 3x3 vector
// v = mu
void M3DMultiplyMatrix33Vector3f(M3DVector3f v, M3DMatrix44f m, M3DVector3f u)
{
  v[0] = m[0] * u[0] + m[3] * u[1] + m[6] * u[2];
  v[1] = m[1] * u[0] + m[4] * u[1] + m[7] * u[2];
  v[2] = m[2] * u[0] + m[5] * u[1] + m[8] * u[2];
}

void M3DMultiplyMatrix44f(M3DMatrix44f a, M3DMatrix44f b, M3DMatrix44f c)
{
  c[0] = a[0] * b[0] + a[4] * b[1] + a[8] * b[2] + a[12] * b[3];
  c[1] = a[1] * b[0] + a[5] * b[1] + a[9] * b[2] + a[13] * b[3];
  c[2] = a[2] * b[0] + a[6] * b[1] + a[10] * b[2] + a[14] * b[3];
  c[3] = a[3] * b[0] + a[7] * b[1] + a[11] * b[2] + a[15] * b[3];
  c[4] = a[0] * b[4] + a[4] * b[5] + a[8] * b[6] + a[12] * b[7];
  c[5] = a[1] * b[4] + a[5] * b[5] + a[9] * b[6] + a[13] * b[7];
  c[6] = a[2] * b[4] + a[6] * b[5] + a[10] * b[6] + a[14] * b[7];
  c[7] = a[3] * b[4] + a[7] * b[5] + a[11] * b[6] + a[15] * b[7];
  c[8] = a[0] * b[8] + a[4] * b[9] + a[8] * b[10] + a[12] * b[11];
  c[9] = a[1] * b[8] + a[5] * b[9] + a[9] * b[10] + a[13] * b[11];
  c[10] = a[2] * b[8] + a[6] * b[9] + a[10] * b[10] + a[14] * b[11];
  c[11] = a[3] * b[8] + a[7] * b[9] + a[11] * b[10] + a[15] * b[11];
  c[12] = a[0] * b[12] + a[4] * b[13] + a[8] * b[14] + a[12] * b[15];
  c[13] = a[1] * b[12] + a[5] * b[13] + a[9] * b[14] + a[13] * b[15];
  c[14] = a[2] * b[12] + a[6] * b[13] + a[10] * b[14] + a[14] * b[15];
  c[15] = a[3] * b[12] + a[7] * b[13] + a[11] * b[14] + a[15] * b[15];
}

void M3DInvertMatrix44f(M3DMatrix44f a, M3DMatrix44f b)
{
  M3DMatrix44f RInverse = {a[0], a[4], a[8], 0.0f,
                           a[1], a[5], a[9], 0.0f,
                           a[2], a[6], a[10], 0.0f,
                           0.0f, 0.0f, 0.0f, 1.0f};
  M3DMatrix44f TInverse = {1.0f, 0.0f, 0.0f, 0.0f,
                           0.0f, 1.0f, 0.0f, 0.0f,
                           0.0f, 0.0f, 1.0f, 0.0f,
                           -a[12], -a[13], -a[14], 1.0f};
  memcpy(b, RInverse, sizeof(M3DMatrix44f));
  b[12] = RInverse[0] * TInverse[12] + RInverse[4] * TInverse[13] + RInverse[8] * TInverse[14];
  b[13] = RInverse[1] * TInverse[12] + RInverse[5] * TInverse[13] + RInverse[9] * TInverse[14];
  b[14] = RInverse[2] * TInverse[12] + RInverse[6] * TInverse[13] + RInverse[10] * TInverse[14];
}