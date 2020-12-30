#ifndef _FRAME_3D_
#define _FRAME_3D_

#include "math3d.h"
#include "GL/glew.h"
#include "GL/freeglut.h"

class Frame3D
{

protected:
  M3DVector3f Origin;
  M3DVector3f Up;
  M3DVector3f Forward;  

public:
  /* Constructors */
  Frame3D();
  Frame3D(M3DVector3f Origin, M3DVector3f Up, M3DVector3f Forward);

  /* Gets */
  void GetOrigin(M3DVector3f origin);
  void GetUp(M3DVector3f up);
  void GetForward(M3DVector3f forward);  

  /* Sets */
  void SetOrigin(float x, float y, float z);
  void SetUp(float x, float y, float z);
  void SetForward(float x, float y, float z);

  /* Moves the origin forward */
  void MoveForward(float d);

  /* Moves the origin to the right (cross vector's direction) */
  void MoveRight(float d);

  /* Moves the origin up */
  void MoveUp(float d);

  /* Rotates the frame around local axis X (cross vector) */
  void RotateLocalX(float angle);

  /* Rotates the frame around local axis Y (Up vector) */
  void RotateLocalY(float angle);

  /* Rotates the frame around local axis Z (Forward vector) */
  void RotateLocalZ(float angle);
  
  void RotateWorld(float fAngle, float x, float y, float z);

  /* Builds the transformation matrix */
  void BuildTransformationMatrix(M3DMatrix44f m, bool rotateOnly);

  /* Applies camera transformation */
  void ApplyCameraTransform(bool rotateOnly = false);

  /* Gets camera orientation */
  void GetCameraOrientation(M3DMatrix44f m);

  /* Applies camera transformation */
  void ApplyActorTransform(bool rotateOnly = false);  
};

#endif