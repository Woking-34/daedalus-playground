#include "frame3d.h"

// Constructors
Frame3D::Frame3D()
{
  Origin[0] = 0.0f;
  Origin[1] = 0.0f;
  Origin[2] = 0.0f;

  Forward[0] = 0.0f;
  Forward[1] = 0.0f;
  Forward[2] = -1.0f;

  Up[0] = 0.0f;
  Up[1] = 1.0f;
  Up[2] = 0.0f;  
}

Frame3D::Frame3D(M3DVector3f origin, M3DVector3f up, M3DVector3f forward)
{
  SetOrigin(origin[0], origin[1], origin[2]);
  SetForward(forward[0], forward[1], forward[2]);
  SetUp(up[0], up[1], up[2]);  
}

// Gets
void Frame3D::GetOrigin(M3DVector3f origin)
{
  M3DCopyVector3f(origin, Origin);
}
void Frame3D::GetForward(M3DVector3f forward)
{
  M3DCopyVector3f(forward, Forward);
}
void Frame3D::GetUp(M3DVector3f up)
{
  M3DCopyVector3f(up, Up);
}

// Sets
void Frame3D::SetOrigin(float x, float y, float z)
{
  Origin[0] = x;
  Origin[1] = y;
  Origin[2] = z;
}

void Frame3D::SetForward(float x, float y, float z)
{
  Forward[0] = x;
  Forward[1] = y;
  Forward[2] = z;
}

void Frame3D::SetUp(float x, float y, float z)
{
  Up[0] = x;
  Up[1] = y;
  Up[2] = z;
}

// Moves the origin forward
void Frame3D::MoveForward(float d)
{
  Origin[0] += Forward[0] * d;
  Origin[1] += Forward[1] * d;
  Origin[2] += Forward[2] * d;
}

// Moves the origin to the right (cross vector's direction)
void Frame3D::MoveRight(float d)
{
  M3DVector3f cross;
  M3DCrossProduct33f(cross, Up, Forward);
  Origin[0] += cross[0] * d;
  Origin[1] += cross[1] * d;
  Origin[2] += cross[2] * d;
}

// Moves the origin up
void Frame3D::MoveUp(float d)
{
  Origin[0] += Up[0] * d;
  Origin[1] += Up[1] * d;
  Origin[2] += Up[2] * d;
}

// Rotates the frame around local axis X (cross vector)
void Frame3D::RotateLocalX(float angle)
{
  M3DVector3f cross;
  M3DCrossProduct33f(cross, Up, Forward);

  M3DMatrix33f m;
  M3DGetRotationMatrix33f(m, angle, cross[0], cross[1], cross[2]);

  // Rotate Up vector
  M3DVector3f newVect;
  M3DMultiplyMatrix33Vector3f(newVect, m, Up);
  M3DCopyVector3f(Up, newVect);

  // Rotate Forward vector
  M3DMultiplyMatrix33Vector3f(newVect, m, Forward);
  M3DCopyVector3f(Forward, newVect);
}

// Rotates the frame around local axis Y (Up vector)
void Frame3D::RotateLocalY(float angle)
{
  M3DMatrix33f m;
  M3DGetRotationMatrix33f(m, angle, Up[0], Up[1], Up[2]);

  // Rotate Forward vector
  M3DVector3f newVect;
  M3DMultiplyMatrix33Vector3f(newVect, m, Forward);
  M3DCopyVector3f(Forward, newVect);
}

// Rotates the frame around local axis Z (Forward vector)
void Frame3D::RotateLocalZ(float angle)
{
  M3DMatrix33f m;
  M3DGetRotationMatrix33f(m, angle, Forward[0], Forward[1], Forward[2]);

  // Rotate Up vector
  M3DVector3f newVect;
  M3DMultiplyMatrix33Vector3f(newVect, m, Up);
  M3DCopyVector3f(Up, newVect);
}

void Frame3D::RotateWorld(float fAngle, float x, float y, float z)
{
  M3DMatrix44f rotMat;

  // Create the Rotation matrix
  M3DGetRotationMatrix44f(rotMat, fAngle, x, y, z);

  M3DVector3f newVect;

  // Transform the up axis (inlined 3x3 rotation)
  newVect[0] = rotMat[0] * Up[0] + rotMat[4] * Up[1] + rotMat[8] *  Up[2];	
  newVect[1] = rotMat[1] * Up[0] + rotMat[5] * Up[1] + rotMat[9] *  Up[2];	
  newVect[2] = rotMat[2] * Up[0] + rotMat[6] * Up[1] + rotMat[10] * Up[2];	
  M3DCopyVector3f(Up, newVect);

  // Transform the forward axis
  newVect[0] = rotMat[0] * Forward[0] + rotMat[4] * Forward[1] + rotMat[8] *  Forward[2];	
  newVect[1] = rotMat[1] * Forward[0] + rotMat[5] * Forward[1] + rotMat[9] *  Forward[2];	
  newVect[2] = rotMat[2] * Forward[0] + rotMat[6] * Forward[1] + rotMat[10] * Forward[2];	
  M3DCopyVector3f(Forward, newVect);
}

// Builds the transformation matrix
void Frame3D::BuildTransformationMatrix(M3DMatrix44f m, bool rotateOnly)
{
  M3DLoadIdentityMatrix44f(m);

  M3DVector3f cross;
  M3DCrossProduct33f(cross, Up, Forward);

  m[0] = cross[0];
  m[1] = cross[1];
  m[2] = cross[2];
  m[3] = 0.0f;

  m[4] = Up[0];
  m[5] = Up[1];
  m[6] = Up[2];
  m[7] = 0.0f;

  m[8] = Forward[0];
  m[9] = Forward[1];
  m[10] = Forward[2];
  m[11] = 0.0f;

  if (rotateOnly)
  {
    m[12] = 0.0f;
    m[13] = 0.0f;
    m[14] = 0.0f;
  }
  else
  {
    m[12] = Origin[0];
    m[13] = Origin[1];
    m[14] = Origin[2];
  }
  m[15] = 1.0f;
}

// Applies camera transformation
void Frame3D::ApplyCameraTransform(bool rotateOnly)
{
  M3DMatrix44f m;
  GetCameraOrientation(m);
  glMultMatrixf(m);

  // Multiply with the translation matrix if not only rotate
  if (!rotateOnly)
    glTranslatef(-Origin[0], -Origin[1], -Origin[2]);
}

// Gets camera orientation
void Frame3D::GetCameraOrientation(M3DMatrix44f m)
{
  M3DVector3f x, z;

  z[0] = -Forward[0];
  z[1] = -Forward[1];
  z[2] = -Forward[2];

  M3DCrossProduct33f(x, Up, z);

  m[0]  = x[0];
  m[4]  = x[1];
  m[8]  = x[2];
  m[12] = 0.0f;
  m[1]  = Up[0];
  m[5]  = Up[1];
  m[9]  = Up[2];
  m[13] = 0.0f;
  m[2]  = z[0];
  m[6]  = z[1];
  m[10] = z[2];
  m[14] = 0.0f;
  m[3]  = 0.0f;
  m[7]  = 0.0f;
  m[11] = 0.0f;
  m[15] = 1.0f;
}

// Applies actor transformation
void Frame3D::ApplyActorTransform(bool rotateOnly)
{
  M3DMatrix44f m;
  BuildTransformationMatrix(m, rotateOnly);
  glMultMatrixf(m);
}
