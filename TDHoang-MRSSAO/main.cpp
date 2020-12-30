/**********************************************************************************
 *                                                                                *
 * Multi-Resolution Screen-Space Ambient Occlusion                                *
 * Author: Thai-Duong Hoang                                                       *
 *                                                                                *
 * You are welcome to use and modify this code for any purpose, but please note   *
 * that it comes without any guarantees.                                          *
 *                                                                                *
 * Please keep in mind that this project should only be considered a prototype to *
 * prove that an idea works. As such, do not expect high code quality =)          *
 *                                                                                *
 **********************************************************************************/

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>

#include "GL/glew.h"
#include "GL/freeglut.h"

#include "frame3d.h"
#include "frameGrab.h"
#include "frameRate.h"
#include "glm.h"
#include "math3d.h"
#include "shader.h"

#define LEVEL_COUNT 5 // number of mip-map levels
#define RESOLUTION 800 // the final output resolution
#define MOVEMENT_SPEED 30.0f // camera's movement speed
#define PI 3.14159265

const float scale = 1.0; // used to scale the scene model
const float dMax = 3.0f; // AO radius of influence
const float rMax = 5.0; // maximum screen-space sampling radius
const float zNear = 0.1f;
const float zFar = 1000.0f;
const float fov = 60.0f;

bool oddFrame = true; // used to ping-pong buffers in the temporal filtering step
bool key_state[256] = { false }; // used to smoothly control camera's movement

// Euler angles to describe camera's orientation
float phi = 0.0f;
float theta = 0.0f;
float psi = 0.0f;

int minResolution = RESOLUTION; // the coarsest resolution

FrameRate frameRate(30); // frame rates are averaged over 30 frames
FrameGrab frameGrab; // used to get screenshots

GLuint modelList; // scene models' display list
GLMmodel *model; // glm object containing scene models

Frame3D camera;

std::string meshFile; // scene models' file name

float gluOrtho[LEVEL_COUNT][16]; // orthogonal projection matrix
float iMVMat[16]; // inversed model-view matrix
float mVMat[16]; // model-view matrix
float projMat[16]; // perspective projection matrix

// timers, used to control camera's movement and animations
int previousTime;
int currentTime;
int elapsedTime;

// buffers and textures ([0] = finest resolution)
GLuint frameBufs[LEVEL_COUNT];
GLuint depthBufs[LEVEL_COUNT];
GLuint posTex[LEVEL_COUNT];
GLuint normTex[LEVEL_COUNT];
GLuint aoTex[LEVEL_COUNT];
GLuint aoTexBlur[LEVEL_COUNT];
float *randRot[LEVEL_COUNT];
GLuint randRotTex[LEVEL_COUNT];
GLuint lastFrameAOTex;
GLuint lastFramePosTex;

// shader programs
GLuint geometryProg;
GLuint downsampleProg[LEVEL_COUNT]; // [0] = finest resolution
GLuint aoProg[LEVEL_COUNT]; // [LEVEL_COUNT - 1] = coarsest resolution
GLuint aoBlurProg[LEVEL_COUNT];

// render targets
GLenum bufs01[2] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1 };
GLenum bufs51[2] = { GL_COLOR_ATTACHMENT5, GL_COLOR_ATTACHMENT1 };

std::string rootData = "";

/* Blurs the AO textures */
void Blur(int size, int index)
{  
  glUseProgram(aoBlurProg[index]);

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_RECTANGLE, aoTex[index]);      
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_RECTANGLE, normTex[index]);
  glActiveTexture(GL_TEXTURE2);
  glBindTexture(GL_TEXTURE_RECTANGLE, posTex[index]);  

  glBindFramebuffer(GL_FRAMEBUFFER, frameBufs[index]);
  glDrawBuffer(GL_COLOR_ATTACHMENT3);   

  glPushAttrib(GL_VIEWPORT_BIT);  
  glViewport(0, 0, size, size);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  
  glBegin(GL_QUADS);
  glVertex2d(0, 0);
  glVertex2d(0, size);
  glVertex2d(size, size);
  glVertex2d(size, 0);
  glEnd();  
  glPopAttrib();
}

/* Creates orthographic projection matrices */
void BuildMatrices()
{
  int size = RESOLUTION;

  for (int i = 0; i < LEVEL_COUNT; ++i)
  {    
    glPushMatrix();

    glLoadIdentity();
    gluOrtho2D(0, size, 0, size);

    glGetFloatv(GL_MODELVIEW_MATRIX, gluOrtho[i]);

    glPopMatrix();

    size /= 2;
  }

  glPushMatrix();
  glLoadIdentity();
  gluPerspective(fov, (GLfloat) RESOLUTION / (GLfloat) RESOLUTION, zNear, zFar);  
  glGetFloatv(GL_MODELVIEW_MATRIX, projMat);
  glPopMatrix();

  glPushMatrix();
  glLoadIdentity();
  camera.ApplyCameraTransform();
  glGetFloatv(GL_MODELVIEW_MATRIX, mVMat);
  M3DInvertMatrix44f(mVMat, iMVMat);
  glPopMatrix();
}

// Frees mem before quiting
void CleanUp()
{
  printf("Exiting...\n");

  glDeleteFramebuffers(LEVEL_COUNT, frameBufs);  
  glDeleteRenderbuffers(LEVEL_COUNT, depthBufs);

  glDeleteTextures(1, &lastFrameAOTex);
  glDeleteTextures(1, &lastFramePosTex);
  glDeleteTextures(LEVEL_COUNT, posTex);
  glDeleteTextures(LEVEL_COUNT, normTex);
  glDeleteTextures(LEVEL_COUNT, aoTex);
  glDeleteTextures(LEVEL_COUNT, aoTexBlur);

  glDeleteLists(modelList, 1);

  for (int i = 0; i < LEVEL_COUNT; ++i)
    delete [] randRot[i];  
}

/* Downsamples the g-buffer */
void Downsample(int size, int index)
{
  glUseProgram(downsampleProg[index]);

  glActiveTexture(GL_TEXTURE0);
  if (index == 1 && !oddFrame)
    glBindTexture(GL_TEXTURE_RECTANGLE, lastFramePosTex);        
  else  
    glBindTexture(GL_TEXTURE_RECTANGLE, posTex[index - 1]);  

  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_RECTANGLE, normTex[index - 1]);

  glBindFramebuffer(GL_FRAMEBUFFER, frameBufs[index]);
  glDrawBuffers(2, bufs01);

  glPushAttrib(GL_VIEWPORT_BIT);
  glViewport(0, 0, size, size);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glBegin(GL_QUADS);
  glVertex2d(0, 0);
  glVertex2d(0, size);
  glVertex2d(size, size);
  glVertex2d(size, 0);
  glEnd();
  glPopAttrib();
}

/* Draws the model */
void DrawModel()
{
  glPushMatrix();
  glScalef(scale, scale, scale);
  glCallList(modelList);
  glPopMatrix();
}

/* Initializes some OpenGL states */
void GLInit()
{
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_TEXTURE_RECTANGLE);

  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
}

/* Initializes the camera's position and orientation */
void InitCamera()
{      
  // sibenik.obj (scale: 1.0)
  camera.SetOrigin(16.610f, 0.759f, 14.594f);        
  camera.RotateWorld(0.850376771f, 1.0f, 0.0f, 0.0f);
  camera.RotateWorld(0.00879645943f, 0.0f, 1.0f, 0.0f); 
  camera.RotateWorld(1.60875724f, 0.0f, 0.0f, 1.0f);

  // conference.obj (scale: 1.0)
  /*camera.SetOrigin(4.624f, 14.788f, 8.425f);
  camera.RotateWorld(1.019780f, 1.0f, 0.0f, 0.0f);
  camera.RotateWorld(-0.016726f, 0.0f, 1.0f, 0.0f);
  camera.RotateWorld(-1.648948f, 0.0f, 0.0f, 1.0f);*/

  // sponza.obj (scale: 1.0)
  /*camera.SetOrigin(-5.505f, 8.170f, 0.086f);
  camera.RotateWorld(2.209693f, 1.0f, 0.0f, 0.0f);
  camera.RotateWorld(-1.549080f, 0.0f, 1.0f, 0.0f);
  camera.RotateWorld(-2.583589f, 0.0f, 0.0f, 1.0f);*/

  // Bedroom.obj (scale: 0.02)
  /*camera.SetOrigin(3.1142983f, 4.5599704f, 5.0486407f);
  camera.SetForward(-0.7761337f, -0.5180085f, -0.3595259f);
  camera.SetUp(-0.4808656f, 0.8550724f, -0.1939232f);*/

  // LocalTrain_Scene.obj (scale: 0.4)
  /*camera.SetOrigin(-1.9109584f, 1.8348960f, -5.8745866f);
  camera.SetForward(0.3675672f, -0.3185951f, -0.8737243f);
  camera.SetUp(0.0899387f, 0.9472651f, -0.3075740f);*/

  // NeonChrome.obj (scale: 0.4)
  /*camera.SetOrigin(-4.6555371f, 2.8841534f, 13.7699118f);
  camera.SetForward(0.6703116f, -0.2015510f, -0.7141659f);
  camera.SetUp(0.1565388f, 0.9791557f, -0.1294100f);*/

  // ChristmasChallenge3.obj (scale: 0.1)
  /*camera.SetOrigin(-7.2839103f, 4.1200843f, 1.0289787f);
  camera.SetForward(0.9957945f, -0.0314440f, -0.0860315f);
  camera.SetUp(0.0296265f, 0.9993098f, -0.0223240f);*/

  // FilmNoirChallenge.obj (scale: 0.5)
  /*camera.SetOrigin(.5294018f, 2.8532739f, 14.9101143f);
  camera.SetForward(-0.7219405f, -0.0246199f, 0.6915044f);
  camera.SetUp(0.0016653f, 0.9992970f, 0.0373177f);*/

  // TheShopGirls.obj (scale: 0.05)
  /*camera.SetOrigin(10.4550095f, 5.0569258f, 6.1843805f);
  camera.SetForward(-0.4004210f, -0.4440890f, -0.8015903f);
  camera.SetUp(-0.2241063f, 0.8956328f, -0.3842444f);*/

  // TheCarnival.obj (scale: 0.005)
  /*camera.SetOrigin(1.9547788f, 2.4798114f, -27.7151527f);
  camera.SetForward(-0.9308686f, -0.2969295f, -0.2128343f);
  camera.SetUp(-0.2859567f, 0.9547774f, -0.0813668f);*/

  // ScienceFictionChallenge.obj (scale: 0.1)
  /*camera.SetOrigin(13.2646532f, 1.7795181f, 0.0235817f);
  camera.SetForward(-0.9973823f, -0.0721106f, -0.0027442f);
  camera.SetUp(-0.0721242f, 0.9973760f, 0.0062320f);*/

  // museumhallRD.obj (scale: 0.5)
  /*camera.SetOrigin(23.4290047f, -2.1287637f, -5.1787529f);
  camera.SetForward(-0.9612807f, -0.0343969f, -0.2733888f);
  camera.SetUp(-0.0351828f, 0.9993753f, -0.0020315f);*/

  // KingsTreasure.obj (scale: 0.1)
  /*camera.SetOrigin(-2.3643525f, 2.0722952f, 2.0264993f);
  camera.SetForward(0.6282180f, -0.3486778f, -0.6955153f);
  camera.SetUp(0.2197803f, 0.9370723f, -0.2712604f);*/
}

/* Keys down */
void KeyDown(unsigned char key, int x, int y)
{
  key_state[key - 'a'] = true;
}

/* Keys up */
void KeyUp(unsigned char key, int x, int y)
{
  key_state[key - 'a'] = false;
}

/* Loads the model from file */
void LoadModel()
{  
  if (!model)
  { 
    char *cstr = new char [meshFile.size() + 1];
    strcpy (cstr, meshFile.c_str());

    model = glmReadOBJ(cstr);
    if (!model)
      exit(0);

    glmFacetNormals(model);
    //glmVertexNormals(model, 90.0f);
    modelList = glmList(model, GLM_SMOOTH);
    glmDelete(model);
  }
}

/* Converts from matrix representation of camera's orientation to Euler angles */
void Matrix2Euler(float *m, float *psi, float *theta, float *phi)
{
  if (m[8] != -1.0f && m[8] != 1.0f)
  {
    *theta = -asin(m[8]);
    *psi = atan2(m[9] / cos(*theta), m[10] / cos(*theta));
    *phi = atan2(m[4] / cos(*theta), m[0] / cos(*theta));
  }
  else
  {
    *phi = 0.0f;
    if (m[8] == -1.0f)
    {
      *theta = PI / 2.0f;
      *psi = *phi + atan2(m[1], m[2]);
    }
    else
    {
      *theta = -PI / 2.0f;
      *psi = -(*phi) + atan2(-m[1], -m[2]);
    }
  }
}

/* Moves or rotates the camera */
void MoveCamera()
{
  if (key_state['w' - 'a'])
    camera.MoveForward(0.1f * MOVEMENT_SPEED * (elapsedTime / 1000.0));   
  if (key_state['s' - 'a'])
    camera.MoveForward(-0.1f * MOVEMENT_SPEED * (elapsedTime / 1000.0));     
  if (key_state[GLUT_KEY_LEFT])
    camera.RotateLocalY(0.02f * MOVEMENT_SPEED * (elapsedTime / 1000.0));
  if (key_state[GLUT_KEY_RIGHT])
    camera.RotateLocalY(-0.02f * MOVEMENT_SPEED * (elapsedTime / 1000.0));
  if (key_state[GLUT_KEY_UP])
    camera.RotateLocalX(0.02f * MOVEMENT_SPEED * (elapsedTime / 1000.0));
  if (key_state[GLUT_KEY_DOWN])
    camera.RotateLocalX(-0.02f * MOVEMENT_SPEED * (elapsedTime / 1000.0));
  if (key_state['a' - 'a'])
    camera.MoveRight(0.1f * MOVEMENT_SPEED * (elapsedTime / 1000.0) );
  if (key_state['d' - 'a'])
    camera.MoveRight(-0.1f * MOVEMENT_SPEED * (elapsedTime / 1000.0));
}

/* glut's idle function */
void Idle()
{
  glutPostRedisplay();

  currentTime = glutGet(GLUT_ELAPSED_TIME);
  elapsedTime = currentTime - previousTime;

  MoveCamera();
}

/* Prints the camera position and orientation */
void PrintCamera()
{
  FILE *file = fopen("./txt/camera.txt","w");

  fprintf(file, "Camera:\n");
  
  M3DVector3f origin, up, forward;
  camera.GetOrigin(origin);
  camera.GetForward(forward);
  camera.GetUp(up);
    
  fprintf(file, "At: %.7f %.7f %.7f\n", origin[0], origin[1], origin[2]);
  fprintf(file, "Forward: %.7f %.7f %.7f\n", forward[0], forward[1], forward[2]);
  fprintf(file, "Up: %.7f %.7f %.7f\n", up[0], up[1], up[2]);

  float m[16];
  camera.GetCameraOrientation(m);  
  Matrix2Euler(m, &psi, &theta, &phi);

  fprintf(file, "Psi: %.7f\n", psi);
  fprintf(file, "Theta: %.7f\n", theta);
  fprintf(file, "Phi: %.7f\n", phi);

  fclose(file);
}

/* Renders the scene for the first time */
void Render2GBuffer()
{  
  glUseProgram(geometryProg);

  glBindFramebuffer(GL_FRAMEBUFFER, frameBufs[0]);
  if (oddFrame)
    glDrawBuffers(2, bufs01);
  else
    glDrawBuffers(2, bufs51);  

  glPushAttrib(GL_VIEWPORT_BIT);  
  glViewport(0, 0, RESOLUTION, RESOLUTION);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  DrawModel();
  glPopAttrib();  
}

/* Renders AO at different resolutions */
void RenderAO(int size, int index)
{  
  glUseProgram(aoProg[index]);  

  glActiveTexture(GL_TEXTURE0);
  if (index == 0 && !oddFrame)
    glBindTexture(GL_TEXTURE_RECTANGLE, lastFramePosTex);
  else
    glBindTexture(GL_TEXTURE_RECTANGLE, posTex[index]);
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_RECTANGLE, normTex[index]);

  if (index < LEVEL_COUNT - 1 || LEVEL_COUNT == 1)
  {    
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_RECTANGLE, aoTexBlur[index + 1]);
    //glBindTexture(GL_TEXTURE_RECTANGLE, aoTex[index + 1]);
    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_RECTANGLE, posTex[index + 1]);
    glActiveTexture(GL_TEXTURE4);
    glBindTexture(GL_TEXTURE_RECTANGLE, normTex[index + 1]);
  }

  if (size == RESOLUTION)
  {        
    if (oddFrame)
    {
      glActiveTexture(GL_TEXTURE5);
      glBindTexture(GL_TEXTURE_RECTANGLE, lastFrameAOTex);
      glActiveTexture(GL_TEXTURE6);
      glBindTexture(GL_TEXTURE_RECTANGLE, lastFramePosTex);
    }
    else
    {
      glActiveTexture(GL_TEXTURE5);
      glBindTexture(GL_TEXTURE_RECTANGLE, aoTex[0]);
      glActiveTexture(GL_TEXTURE6);
      glBindTexture(GL_TEXTURE_RECTANGLE, posTex[0]);
    }        

    glActiveTexture(GL_TEXTURE7);
    glBindTexture(GL_TEXTURE_RECTANGLE, randRotTex[index]);

    glUniformMatrix4fv(glGetUniformLocation(aoProg[0], "mVMat"), 1, false, mVMat);
    glGetFloatv(GL_MODELVIEW_MATRIX, mVMat);    
    M3DInvertMatrix44f(mVMat, iMVMat);
    glUniformMatrix4fv(glGetUniformLocation(aoProg[0], "iMVMat"), 1, false, iMVMat);
  }

  if (size < RESOLUTION)
  {
    glBindFramebuffer(GL_FRAMEBUFFER, frameBufs[index]);
    glDrawBuffer(GL_COLOR_ATTACHMENT2);
  }
  else
  {
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
  }

  glPushAttrib(GL_VIEWPORT_BIT);  
  glViewport(0, 0, size, size);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glBegin(GL_QUADS);
  glVertex2d(0, 0);
  glVertex2d(0, size);
  glVertex2d(size, size);
  glVertex2d(size, 0);
  glEnd();

  glPopAttrib();

  if (size < RESOLUTION)
    Blur(size, index);
}

/* Resizes the display window */
void Reshape(int w, int h)
{
  // avoid dividing by 0
  if (h == 0)
    h = 1;

  glViewport(0, 0, w, h);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(fov, (GLfloat) w / (GLfloat) h, zNear, zFar);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

/* Setups the shader programs used in the AO blur steps */
void SetupAOBlurPrograms()
{
  GLuint vs, fs;  

  for (int i = LEVEL_COUNT - 1; i >= 0; --i)
  {
    vs = LoadShader(EVertexShader, "shaders/ortho.vert", rootData);
    fs = LoadShader(EFragmentShader, "shaders/blur.frag", rootData);

    aoBlurProg[i] = CreateProgram();
    AttachShader(aoBlurProg[i], vs);
    AttachShader(aoBlurProg[i], fs);
    glBindFragDataLocation(aoBlurProg[i], 0, "AO");
    LinkProgram(aoBlurProg[i]);         
    glUseProgram(aoBlurProg[i]);

    glUniform1i(glGetUniformLocation(aoBlurProg[i], "aoTex"), 0); 
    glUniform1i(glGetUniformLocation(aoBlurProg[i], "normTex"), 1); 
    glUniform1i(glGetUniformLocation(aoBlurProg[i], "posTex"), 2); 
    glUniformMatrix4fv(glGetUniformLocation(aoBlurProg[i], "gluOrtho"), 1, false, gluOrtho[i]);
  }
}

/* glut's display function */
void Display()
{  
  frameRate.StartFrame();

  static int counter = 0; // this is used to control the printing of fps

  previousTime = currentTime;

  glPushMatrix();  
  camera.ApplyCameraTransform();
  Render2GBuffer();  
  for (int index = 1, size = RESOLUTION / 2; index < LEVEL_COUNT; ++index)
  {
    Downsample(size, index);
    size /= 2;
  }  
  for (int index = LEVEL_COUNT - 1, size = minResolution; index >= 0; --index)
  {
    RenderAO(size, index);
    size *= 2;
  }
  glPopMatrix();

  if (oddFrame)
    glBindTexture(GL_TEXTURE_RECTANGLE, aoTex[0]);
  else
    glBindTexture(GL_TEXTURE_RECTANGLE, lastFrameAOTex);
  oddFrame = !oddFrame;

  glCopyTexSubImage2D(GL_TEXTURE_RECTANGLE, 0, 0, 0, 0, 0, RESOLUTION, RESOLUTION); 

  glutSwapBuffers();

  frameRate.EndFrame();

  // display the fps once every 30 frames
  counter++;
  if (counter > 30)
  {
    printf("%.2f fps\n", frameRate.GetLastFrameRate());
    counter = 0;
  }
}

/* Setups the shader programs used in the AO computation steps */
void SetupAOPrograms()
{
  GLfloat poissonDisk[] = { -0.6116678f,  0.04548655f, -0.26605980f, -0.6445347f,
    -0.4798763f,  0.78557830f, -0.19723210f, -0.1348270f,
    -0.7351842f, -0.58396650f, -0.35353550f,  0.3798947f,
    0.1423388f,  0.39469180f, -0.01819171f,  0.8008046f,
    0.3313283f, -0.04656135f,  0.58593510f,  0.4467109f,
    0.8577477f,  0.11188750f,  0.03690137f, -0.9906120f,
    0.4768903f, -0.84335800f,  0.13749180f, -0.4746810f,
    0.7814927f, -0.48938420f,  0.38269190f,  0.8695006f };

  GLuint fs, vs;

  float size = (float) RESOLUTION;

  vs = LoadShader(EVertexShader, "shaders/ortho.vert", rootData);
  for (int i = 0; i < LEVEL_COUNT; ++i)
  {    
    if (i == 0)    
      fs = LoadShader(EFragmentShader, "shaders/aoLast.frag", rootData);    
    else if (i == LEVEL_COUNT - 1)       
      fs = LoadShader(EFragmentShader, "shaders/aoFirst.frag", rootData);
    else          
      fs = LoadShader(EFragmentShader, "shaders/ao.frag", rootData);

    aoProg[i] = CreateProgram();
    AttachShader(aoProg[i], vs);
    AttachShader(aoProg[i], fs);
    glBindFragDataLocation(aoProg[i], 0, "AO");
    LinkProgram(aoProg[i]);
    glUseProgram(aoProg[i]);

    glUniform1i(glGetUniformLocation(aoProg[i], "posTex"), 0);
    glUniform1i(glGetUniformLocation(aoProg[i], "normTex"), 1);
    if (i < LEVEL_COUNT - 1 || LEVEL_COUNT == 1)
    {
      glUniform1i(glGetUniformLocation(aoProg[i], "loResAOTex"), 2);
      glUniform1i(glGetUniformLocation(aoProg[i], "loResPosTex"), 3);
      glUniform1i(glGetUniformLocation(aoProg[i], "loResNormTex"), 4);
    }
    if (i == 0)
    {
      glUniform1i(glGetUniformLocation(aoProg[i], "lastFrameAOTex"), 5);
      glUniform1i(glGetUniformLocation(aoProg[i], "lastFramePosTex"), 6);
      glUniformMatrix4fv(glGetUniformLocation(aoProg[i], "projMat"), 1, false, projMat);
      glUniformMatrix4fv(glGetUniformLocation(aoProg[i], "iMVMat"), 1, false, iMVMat);
      glUniformMatrix4fv(glGetUniformLocation(aoProg[i], "mVMat"), 1, false, mVMat);

    }
    glUniformMatrix4fv(glGetUniformLocation(aoProg[i], "gluOrtho"), 1, false, gluOrtho[i]);
    glUniform1f(glGetUniformLocation(aoProg[i], "dMax"), dMax);
    glUniform1f(glGetUniformLocation(aoProg[i], "rMax"), rMax);
    float r = size * dMax / (2.0f * abs(tan(fov * PI / 180.0f / 2.0f)));
    glUniform1f(glGetUniformLocation(aoProg[i], "r"), r);
    glUniform1f(glGetUniformLocation(aoProg[i], "resolution"), size);
    glUniform1i(glGetUniformLocation(aoProg[i], "randRotTex"), 7);
    glUniform1fv(glGetUniformLocation(aoProg[i], "poissonDisk"), 32, poissonDisk);

    size /= 2.0f;
  }
}

/* Setups the shader programs used in the downsampling step */
void SetupDownsamplePrograms()
{ 
  GLuint vs = LoadShader(EVertexShader, "shaders/ortho.vert", rootData);
  GLuint fs = LoadShader(EFragmentShader, "shaders/downsample.frag", rootData);

  for (int i = 1; i < LEVEL_COUNT; ++i)
  {
    downsampleProg[i] = CreateProgram();
    AttachShader(downsampleProg[i], vs);
    AttachShader(downsampleProg[i], fs);

    glBindFragDataLocation(downsampleProg[i], 0, "Pos");
    glBindFragDataLocation(downsampleProg[i], 1, "Norm");

    LinkProgram(downsampleProg[i]);
    glUseProgram(downsampleProg[i]);

    glUniform1i(glGetUniformLocation(downsampleProg[i], "hiResPosTex"), 0);
    glUniform1i(glGetUniformLocation(downsampleProg[i], "hiResNormTex"), 1);
    glUniformMatrix4fv(glGetUniformLocation(downsampleProg[i], "gluOrtho"), 1, false, gluOrtho[i]);    
  }
}

/* Creates the buffers and textures used throughout the program */
void SetupFBOs()
{  
  glGenFramebuffers(LEVEL_COUNT, frameBufs);
  glGenRenderbuffers(LEVEL_COUNT, depthBufs);

  glGenTextures(LEVEL_COUNT, posTex);
  glGenTextures(LEVEL_COUNT, normTex);
  glGenTextures(LEVEL_COUNT, aoTex);
  glGenTextures(LEVEL_COUNT, aoTexBlur);
  glGenTextures(1, &lastFrameAOTex);
  glGenTextures(1, &lastFramePosTex);

  glBindTexture(GL_TEXTURE_RECTANGLE, lastFrameAOTex);
  glTexImage2D(GL_TEXTURE_RECTANGLE, 0, GL_RGBA32F, RESOLUTION, RESOLUTION, 0, GL_RGBA, GL_FLOAT, NULL);
  glTexParameterf(GL_TEXTURE_RECTANGLE, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameterf(GL_TEXTURE_RECTANGLE, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

  glBindTexture(GL_TEXTURE_RECTANGLE, lastFramePosTex);
  glTexImage2D(GL_TEXTURE_RECTANGLE, 0, GL_RGBA32F, RESOLUTION, RESOLUTION, 0, GL_RGBA, GL_FLOAT, NULL);
  glTexParameterf(GL_TEXTURE_RECTANGLE, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameterf(GL_TEXTURE_RECTANGLE, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

  int size = RESOLUTION;

  for (int i = 0; i < LEVEL_COUNT; ++i)
  {     
    glBindTexture(GL_TEXTURE_RECTANGLE, posTex[i]);
    glTexImage2D(GL_TEXTURE_RECTANGLE, 0, GL_RGBA32F, size, size, 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameterf(GL_TEXTURE_RECTANGLE, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameterf(GL_TEXTURE_RECTANGLE, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    glBindTexture(GL_TEXTURE_RECTANGLE, normTex[i]);
    glTexImage2D(GL_TEXTURE_RECTANGLE, 0, GL_RGBA32F, size, size, 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameterf(GL_TEXTURE_RECTANGLE, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameterf(GL_TEXTURE_RECTANGLE, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    glBindTexture(GL_TEXTURE_RECTANGLE, aoTex[i]);
    glTexImage2D(GL_TEXTURE_RECTANGLE, 0, GL_RGBA32F, size, size, 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameterf(GL_TEXTURE_RECTANGLE, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameterf(GL_TEXTURE_RECTANGLE, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MIN_FILTER, GL_NEAREST);    

    glBindTexture(GL_TEXTURE_RECTANGLE, aoTexBlur[i]);
    glTexImage2D(GL_TEXTURE_RECTANGLE, 0, GL_RGBA32F, size, size, 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameterf(GL_TEXTURE_RECTANGLE, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameterf(GL_TEXTURE_RECTANGLE, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    glBindRenderbuffer(GL_RENDERBUFFER, depthBufs[i]);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, size, size);

    glBindFramebuffer(GL_FRAMEBUFFER, frameBufs[i]);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthBufs[i]);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_RECTANGLE, posTex[i], 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_RECTANGLE, normTex[i], 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_RECTANGLE, aoTex[i], 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT3, GL_TEXTURE_RECTANGLE, aoTexBlur[i], 0);
    if (i == 0)
    {
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT4, GL_TEXTURE_RECTANGLE, lastFrameAOTex, 0);
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT5, GL_TEXTURE_RECTANGLE, lastFramePosTex, 0);
    }

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE)
      printf("Framebuffer created successfully.\n");
    else
      printf("Framebuffer created unsuccessfully. Error code: %d.\n", glCheckFramebufferStatus(frameBufs[i]));

    size /= 2;
  }
}

/* Setups the shader program used in the geometry step */
void SetupGeometryProgram()
{  
  GLuint vs = LoadShader(EVertexShader, "shaders/geometry.vert", rootData);
  GLuint fs = LoadShader(EFragmentShader, "shaders/geometry.frag", rootData);

  geometryProg = CreateProgram();
  AttachShader(geometryProg, vs);
  AttachShader(geometryProg, fs);

  glBindFragDataLocation(geometryProg, 0, "Pos");
  glBindFragDataLocation(geometryProg, 1, "Norm");

  LinkProgram(geometryProg);
}

/* Creates a random rotation texture in a 3x3 interleaved pattern */
void SetupRandRotTex()
{  
  float pattern[LEVEL_COUNT][3][3][2];

  srand(927);

  glGenTextures(LEVEL_COUNT, randRotTex);

  for (int index = 0; index < LEVEL_COUNT; ++index)
  {
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        float alpha = (float) rand() / RAND_MAX * 2.0 * PI;
        pattern[index][i][j][0] = sin(alpha);
        pattern[index][i][j][1] = cos(alpha);
      }
    }
  }

  int size = RESOLUTION;

  for (int i = 0; i < LEVEL_COUNT; ++i)
  {
    randRot[i] = new float[size * size * 2];
    size /= 2;
  }

  size = RESOLUTION;

  for (int index = 0; index < LEVEL_COUNT; ++index)
  {
    for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < size; ++j)
      {      
        randRot[index][i * size * 2 + j * 2 + 0] = pattern[index][i % 3][j % 3][0];
        randRot[index][i * size * 2 + j * 2 + 1] = pattern[index][i % 3][j % 3][1];
      }
    }

    glBindTexture(GL_TEXTURE_RECTANGLE_ARB, randRotTex[index]);
    glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_RG32F, size, size, 0, GL_RG, GL_FLOAT, randRot[index]);
    glTexParameterf(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameterf(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    size /= 2;
  }
}

/* Special keys up */
void SpecialUp(int key, int x, int y)
{
  key_state[key] = false;
}

/* Special keys down */
void Special(int key, int x, int y)
{
  key_state[key] = true;

  switch (key)
  {
    case GLUT_KEY_INSERT:      
      PrintCamera();
      break;
    case GLUT_KEY_F11:
      frameGrab.CaptureFrame();
      break;
  };     
}

/** file path helper */
static bool findFullPath(const std::string& root, std::string& filePath)
{
    bool fileFound = false;
    const std::string resourcePath = root;

    filePath = resourcePath + filePath;
    for (unsigned int i = 0; i < 16; ++i)
    {
        std::ifstream file;
        file.open(filePath.c_str());
        if (file.is_open())
        {
            fileFound = true;
            break;
        }

        filePath = "../" + filePath;
    }

    return fileFound;
}

/* Main */
int main(int argc, char* argv[])
{  
  {
    std::string locStr = "filesystem.loc";
    size_t len = locStr.size();

    bool fileFound = findFullPath(std::string(SAMPLE_NAME) + "/data/", locStr);
    rootData = locStr.substr(0, locStr.size() - len);
  }

  meshFile = std::string(rootData + "/scn/");
  if (argc == 2)
  {
    meshFile.append(argv[1]);
  }
  else
  {
    meshFile.append("sibenik.obj");
  }
 
  for (int i = 1; i < LEVEL_COUNT; ++i)
    minResolution /= 2;

  glutInit(&argc, argv);
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInitWindowSize(RESOLUTION, RESOLUTION);
  glutCreateWindow("SSAO");

  GLInit();
  glewInit();

  glutReshapeFunc(Reshape);
  glutDisplayFunc(Display);
  glutKeyboardFunc(KeyDown);
  glutKeyboardUpFunc(KeyUp);
  glutSpecialFunc(Special);
  glutSpecialUpFunc(SpecialUp);
  glutIdleFunc(Idle);

  LoadModel();
  
  InitCamera();

  BuildMatrices();  
  
  SetupFBOs();
  SetupRandRotTex();
  SetupGeometryProgram();  
  SetupDownsamplePrograms();  
  SetupAOPrograms();
  SetupAOBlurPrograms();  
  currentTime = glutGet(GLUT_ELAPSED_TIME);
  
  glutMainLoop();

  CleanUp();
}