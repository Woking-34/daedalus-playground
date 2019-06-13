//---------------------------------------------------------------------
// Program's Name:-bullet_PrimitiveShapes
// API Used:-------OpenGL in Win32 API
// Programmer:-----Krishna Kumar
// Email:----------krishx007@yahoo.com
// WebSite:--------www.gfxguru.net
// Compiler Used:- MS VC++ 2008 Express
//---------------------------------------------------------------------

//-------------Controls------------------------------------------------
// MouseLeftButtonDrag  = Rotate modelview(scene)
// MouseRightButtonDrag = Translate modelview(scene) in Z direction
// Up Key, Down Key     = Translate modelview(scene) in Y direction
// Left Key, Right Key  = Translate modelview(scene) in X direction
// R                    = Reset modelview(scene)
//---------------------------------------------------------------------

// To compile the program you need to download the Bullet Physics API and 
// include the dependencies in the project. You can download latest Bullet
// Physics library from http://code.google.com/p/bullet/downloads/list
// Also visit http://www.bulletphysics.org
//----------------------------------------------------------------------

#include <GL/freeglut.h>

int currWidth = 800;
int currHeight = 600;

int prev_time = 0;
int curr_time = 0;
int delta_time = 0;

bool mouseState[3] = { false, false, false };

int mouseX = 0;
int mouseY = 0;
int mouseX_old = 0;
int mouseY_old = 0;

//----------Bullet Physics--------------
#include <btBulletDynamicsCommon.h>

#include "geometry.h"

//variables for modelview (i.e scene) transformation
float  g_fSpinX = 0.0f;
float  g_fSpinY = 0.0f;

float  g_fPosX  = 0.0f;
float  g_fPosY  = -5.0f;
float  g_fPosZ  = 0.0f;

//variables for opengl materials and lighting 
GLfloat mat_specular[]   = { 1.0f, 1.0f, 1.0f, 1.0f };
GLfloat mat_shininess[]  = { 50.0f };
GLfloat mat_diffuse[]    = {1.0f, 1.0f, 1.0f, 1.0f};

GLfloat light_ambient[]  = { 0.1f, 0.1f, 0.1f, 1.0f };
GLfloat light_position[] = { 1.0f, 10.0f, -10.0f, 0.0f };

//variables for to bullet physics API
btAlignedObjectArray<btCollisionShape*>	m_collisionShapes; //keep the collision shapes, for deletion/cleanup
btBroadphaseInterface*					m_broadphase;
btCollisionDispatcher*					m_dispatcher;
btConstraintSolver*						m_solver;
btDefaultCollisionConfiguration*		m_collisionConfiguration;
btDynamicsWorld*						m_dynamicsWorld; //this is the most important class

void InitGL(); //Initialize and creates OpenGL rendering context.
void SetCamera(); //Creates initial OpenGL camera setup
void ProcessInput(); //Process the keyboard input(s)
void RenderScene(); // This is our main rendering command sequence 
void DrawGrid(int); //Draw grids at origin.
void DrawDimensionLines(float,float); //Draws dimension lines at origin.

void StepBulletPhysics()
{
	if(m_dynamicsWorld)//step the simulation
		m_dynamicsWorld->stepSimulation(delta_time / 1000.0f);
}

void InitBulletPhysics()
{

	m_collisionConfiguration = new btDefaultCollisionConfiguration(); //collision configuration contains default setup for memory, collision setup
	m_dispatcher			 = new btCollisionDispatcher(m_collisionConfiguration); //use the default collision dispatcher 
	m_broadphase		  	 = new btDbvtBroadphase();
	m_solver			 	 = new btSequentialImpulseConstraintSolver; //the default constraint solver
	m_dynamicsWorld			 = new btDiscreteDynamicsWorld(m_dispatcher,m_broadphase,m_solver,m_collisionConfiguration);
	m_dynamicsWorld->setGravity(btVector3(0,-10,0));

	//Creating a static shape which will act as ground
	{
		btCollisionShape* groundShape = new btBoxShape(btVector3(50,50,50));
		// btCollisionShape* groundShape = new btStaticPlaneShape(btVector3(0,1,0),50);
		m_collisionShapes.push_back(groundShape);

		btScalar mass = 0; //rigidbody is static if mass is zero, otherwise dynamic
		btVector3 localInertia(0,0,0);

		groundShape->calculateLocalInertia(mass,localInertia);

		btTransform groundTransform;
		groundTransform.setIdentity();
		groundTransform.setOrigin(btVector3(0,-50,0));

		btDefaultMotionState* myMotionState = new btDefaultMotionState(groundTransform); //motionstate provides interpolation capabilities, and only synchronizes 'active' objects
		btRigidBody::btRigidBodyConstructionInfo rbInfo(mass,myMotionState,groundShape,localInertia);
		btRigidBody* body = new btRigidBody(rbInfo);

		m_dynamicsWorld->addRigidBody(body); //add the body to the dynamics world
	}

	//Creating some dynamic boxShape
	{
		btCollisionShape* boxShape = new btBoxShape(btVector3(1,1,1));
		m_collisionShapes.push_back(boxShape);

		btScalar mass = 1.0f;
		btVector3 localInertia(0,0,0);

		boxShape->calculateLocalInertia(mass,localInertia);

		btTransform startTransform;
		startTransform.setIdentity();
		startTransform.setOrigin(btVector3(5,20,0));

		for(int j=0; j<20; j++)
		{
			startTransform.setOrigin(btVector3(0,(j*5)+10 ,0));

			btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
			btRigidBody::btRigidBodyConstructionInfo rbInfo(mass,myMotionState,boxShape,localInertia);
			btRigidBody* body = new btRigidBody(rbInfo);

			m_dynamicsWorld->addRigidBody(body);
		}
	}

	//Creating some dynamic sphereShape
	{
		btCollisionShape* sphereShape = new btSphereShape(1);
		m_collisionShapes.push_back(sphereShape);

		btScalar mass = 1.0f;
		btVector3 localInertia(0,0,0);

		sphereShape->calculateLocalInertia(mass,localInertia);

		btTransform startTransform;
		startTransform.setIdentity();
		startTransform.setOrigin(btVector3(5,20,0));

		for(int j=0; j<20; j++)
		{
			startTransform.setOrigin(btVector3(-5,(j*5)+10 ,0));

			btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
			btRigidBody::btRigidBodyConstructionInfo rbInfo(mass,myMotionState,sphereShape,localInertia);
			btRigidBody* body = new btRigidBody(rbInfo);

			m_dynamicsWorld->addRigidBody(body);
		}
	}

	//Creating some dynamic CylinderShape
	{
		btCollisionShape* cylinderShapeZ = new btCylinderShapeZ(btVector3(1,1,2));
		m_collisionShapes.push_back(cylinderShapeZ);

		btScalar mass = 1.0f;
		btVector3 localInertia(0,0,0);

		cylinderShapeZ->calculateLocalInertia(mass,localInertia);

		btTransform startTransform;
		startTransform.setIdentity();
		startTransform.setOrigin(btVector3(5,20,0));

		for(int j=0; j<20; j++)
		{
			startTransform.setOrigin(btVector3(5,(j*5)+10 ,0));

			btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
			btRigidBody::btRigidBodyConstructionInfo rbInfo(mass,myMotionState,cylinderShapeZ,localInertia);
			btRigidBody* body = new btRigidBody(rbInfo);

			m_dynamicsWorld->addRigidBody(body);
		}
	}

	//Creating some dynamic coneShape
	{
		btCollisionShape* coneShapeZ = new btConeShapeZ(1.5,3);
		m_collisionShapes.push_back(coneShapeZ);

		btScalar mass = 1.0f;
		btVector3 localInertia(0,0,0);

		coneShapeZ->calculateLocalInertia(mass,localInertia);

		btTransform startTransform;
		startTransform.setIdentity();
		startTransform.setOrigin(btVector3(5,20,0));

		for(int j=0; j<20; j++)
		{
			startTransform.setOrigin(btVector3(10,(j*5)+10 ,0));

			btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
			btRigidBody::btRigidBodyConstructionInfo rbInfo(mass,myMotionState,coneShapeZ,localInertia);
			btRigidBody* body = new btRigidBody(rbInfo);

			m_dynamicsWorld->addRigidBody(body);
		}
	}
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0.35f, 0.53f, 0.6f, 1.0f);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glTranslatef(g_fPosX, g_fPosY, g_fPosZ);
	glRotatef(g_fSpinY, 1.0f, 0.0f, 0.0f);
	glRotatef(-g_fSpinX, 0.0f, 1.0f, 0.0f);

	DrawDimensionLines(5, 5);
	DrawGrid(50);

	//---------------------------------------------------
	StepBulletPhysics();

	btScalar	m[16];
	btMatrix3x3	rot; rot.setIdentity();
	const int	numObjects = m_dynamicsWorld->getNumCollisionObjects();

	glEnable(GL_COLOR_MATERIAL);

	for (int i = 0; i<numObjects; i++)
	{
		btCollisionObject*	colObj = m_dynamicsWorld->getCollisionObjectArray()[i];
		btRigidBody*		body = btRigidBody::upcast(colObj);

		if (body && body->getMotionState())
		{
			btDefaultMotionState* myMotionState = (btDefaultMotionState*)body->getMotionState();
			myMotionState->m_graphicsWorldTrans.getOpenGLMatrix(m);
			rot = myMotionState->m_graphicsWorldTrans.getBasis();

			char *shapeName = (char*)body->getCollisionShape()->getName();

			if (!strcmp(shapeName, "Box"))
			{
				glPushMatrix();
				glColor3f(1, 0, 0);
				glMultMatrixf(m);
				renderSolidCube(2);
				glPopMatrix();
			}
			else if (!strcmp(shapeName, "SPHERE"))
			{
				glPushMatrix();
				glColor3f(0, 1, 0);
				glMultMatrixf(m);
				renderSolidSphere(1, 25, 10);
				glPopMatrix();
			}
			else if (!strcmp(shapeName, "CylinderZ"))
			{
				glPushMatrix();
				glColor3f(0, 0, 1);
				glMultMatrixf(m);
				renderSolidCylinder(1, 2, 20, 10);
				glPopMatrix();
			}
		}
	}

	glFinish();
	glutSwapBuffers();
}

void keyfunc(unsigned char k, int x, int y)
{
	switch (k)
	{
	case 27:
		// escape key
		glutLeaveMainLoop();
		break;
	}
}

void mouseclickfunc(int button, int state, int x, int y)
{
	mouseX = x;
	mouseY = y;

	mouseState[button] = state == GLUT_DOWN;
}

void mousemovefunc(int x, int y)
{
	mouseX_old = mouseX;
	mouseY_old = mouseY;

	mouseX = x;
	mouseY = y;

	if (mouseState[GLUT_LEFT_BUTTON])
	{
		if (mouseX_old != -1 && mouseY_old != -1)
		{
			g_fSpinX -= (mouseX - mouseX_old);
			g_fSpinY -= (mouseY - mouseY_old);
		}

		mouseX_old = mouseX;
		mouseY_old = mouseY;
	}

	if (mouseState[GLUT_RIGHT_BUTTON])
	{
		if (mouseX_old != -1 && mouseY_old != -1)
		{
			g_fPosZ -= (mouseY - mouseY_old);
		}

		mouseX_old = mouseX;
		mouseY_old = mouseY;
	}

	{
		mouseX_old = -1;
		mouseY_old = -1;
	}
}

void reshape(int w, int h)
{
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(45.0,	/*FOV in degree*/ (GLdouble)currWidth / (GLdouble)currHeight /*Aspect ratio*/, 0.1,/*Z near*/   -1000.0 /*Z far*/);
	gluLookAt(0, 5, -25, /*Eye position*/		0, 0, 0, /*Position of object to focus*/     0, 1, 0 /*Up is in positive Y direction*/);

	currWidth = w;
	currHeight = h;
}

void idle(void)
{
	prev_time = curr_time;
	curr_time = glutGet(GLUT_ELAPSED_TIME);
	delta_time = curr_time - prev_time;

	glutPostRedisplay();
}

int main(int argc, char **argv)
{
	glutInit(&argc, argv);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(currWidth, currHeight);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_STENCIL);

	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);

	glutCreateWindow("Main Window        krishx007@yahoo.com    www.gfxguru.net");

	InitBulletPhysics();

	SetCamera();

	glutIdleFunc(idle);
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);

	glutKeyboardFunc(keyfunc);
	glutMouseFunc(mouseclickfunc);
	glutMotionFunc(mousemovefunc);

	glutMainLoop();

	return 0;
}

//-----------------------------------------------------------------------------
// Name:  SetCamera()
// Desc: Sets initial OpenGL camera
//-----------------------------------------------------------------------------
void SetCamera()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(45.0,	/*FOV in degree*/ (GLdouble)currWidth / (GLdouble)currHeight /*Aspect ratio*/, 0.1,/*Z near*/   -1000.0 /*Z far*/);
	gluLookAt(0, 5, -25, /*Eye position*/		0, 0, 0, /*Position of object to focus*/     0, 1, 0 /*Up is in positive Y direction*/);

	glShadeModel(GL_SMOOTH);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);
}

//-----------------------------------------------------------------------------
// Name: DrawGrid()
// Desc: Draw grids at origin. 
//-----------------------------------------------------------------------------
void DrawGrid(int GRID_SIZE)
{
	glDisable(GL_LIGHTING);
	glBegin(GL_LINES);
	glColor3f(0.75f, 0.75f, 0.75f);
	for (int i = -GRID_SIZE; i <= GRID_SIZE; i++)
	{
		glVertex3f((float)i, 0, (float)-GRID_SIZE);
		glVertex3f((float)i, 0, (float)GRID_SIZE);

		glVertex3f((float)-GRID_SIZE, 0, (float)i);
		glVertex3f((float)GRID_SIZE, 0, (float)i);
	}
	glEnd();
	glEnable(GL_LIGHTING);
}

//-----------------------------------------------------------------------------
// Name: DrawDimensionLines()
// Desc: Draws dimension lines at origin. 
//-----------------------------------------------------------------------------
void DrawDimensionLines(float size, float lineWidth)
{
	//Draw dimension Lines at origin 
	glDisable(GL_LIGHTING);
	glPushMatrix();
	glBegin(GL_LINES);
	glLineWidth(lineWidth);

	glColor3f(1, 0, 0);//X-axis
	glVertex3f(-size, 0, 0);
	glVertex3f(size, 0, 0);

	glColor3f(0.0f, 1, 0.0f);//Y-axis
	glVertex3f(0, -size, 0);
	glVertex3f(0, size, 0);

	glColor3f(0, 0, 1);//Z-axis
	glVertex3f(0, 0, -size);
	glVertex3f(0, 0, size);
	glEnd();
	glPopMatrix();
	glColor3f(1, 1, 1);//Reseting color
	glEnable(GL_LIGHTING);
}