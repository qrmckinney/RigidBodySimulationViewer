
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <freeglut.h> 
#include <windows.h>

using namespace std;

// Global variables 
char title[] = "3D Shapes";
float angle_x = 1;
float angle_y = 0;
float angle_z = 0;
ifstream ReadingDude;
int NumComponents;
GLfloat Vertices[50][3];
int kappa = 0;
double CurrentCenterOfMass[3];
double CurrentRotation[3][3];

double testNode[3] = { 0,0,1.25 };
double transformedNode[3]; 

void MatrixOperation(double M[3][3], double V[3], double OutV[3]); 
void VSum(double V1[3], double V2[3]); 
double CommandedPlace[3];

double p1_x[100]; 
double p1_y[100];
double p1_z[100];

double p2_x[100];
double p2_y[100];
double p2_z[100];

double p3_x[100];
double p3_y[100];
double p3_z[100];

double p4_x[100];
double p4_y[100];
double p4_z[100];

double p1_V[3] = {  1, 1,0.505 };
double p2_V[3] = { -1, 1,0.505 };
double p3_V[3] = {  1,-1,0.505 };
double p4_V[3] = { -1,-1,0.505 };


void createPropeller(double V[3], double x[100], double y[100], double z[100]) {

	double t;

	for (int i = 0; i < 100; i++) {
		t = ((double) i) / 100; 

		x[i] = 0.5*cos(79 * 3.14159*t) + V[0]; 
		y[i] = 0.5*sin(79 * 3.14159*t) + V[1];
		z[i] = V[2]; 
	}

}

void transformPropeller(double Vin[3], double Vout[3]) {

	MatrixOperation(CurrentRotation, Vin, Vout);
	VSum(Vout, CurrentCenterOfMass);

}


void drawPropeller(double x[100], double y[100], double z[100],double r, double g, double b) {


	double tmp[3]; 
	double tmpOut[3]; 

	glColor3f(r,g,b);

	glBegin(GL_TRIANGLES);
	for (int i = 0; i < 100; i++)	{
		tmp[0] = x[i]; 
		tmp[1] = y[i];
		tmp[2] = z[i];

		transformPropeller(tmp, tmpOut); 		
		glVertex3f(tmpOut[0],tmpOut[2],tmpOut[1]);

	}
	glEnd();

}



/* Initialize OpenGL Graphics */
void initGL() {

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f); 				// Set background color to black and opaque
	glClearDepth(1.0f);                  				// Set background depth to farthest
	glEnable(GL_DEPTH_TEST);   							// Enable depth testing for z-culling
	glDepthFunc(GL_LEQUAL);    							// Set the type of depth-test
	glShadeModel(GL_SMOOTH);   							// Enable smooth shading
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  	// Nice perspective corrections
}

void drawLine(float x, float y, float z, float x1, float y1, float z1) {

	glBegin(GL_LINES);


	glVertex3d(x, y, z);
	glVertex3d(x1, y1, z1);

	glEnd();

}

void display() {

	MatrixOperation(CurrentRotation, testNode, transformedNode);
	VSum(transformedNode, CurrentCenterOfMass);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 	// Clear color and depth buffers

	glMatrixMode(GL_MODELVIEW);     						// To operate on model-view matrix

	glLoadIdentity();                 						// Reset the model-view matrix
	glTranslatef(0,0 , -8.0f);  						// Move into the screen
	//glRotatef(angle_x, 1, 0, 0);
	//glRotatef(angle_y, 0, 1, 0);
	//glRotatef(angle_z, 0, 0, 1);

	glLineWidth(1.25); 

	drawPropeller(p1_x, p1_y, p1_z, 1, 0, 0); 
	drawPropeller(p2_x, p2_y, p2_z, 1, 0, 0);
	drawPropeller(p3_x, p3_y, p3_z, 1, 0, 0);
	drawPropeller(p4_x, p4_y, p4_z, 1, 0, 0);

	glColor3f(1.0f, 1.0f, 0.0f);

		glBegin(GL_POINTS);
		for (int i = 0; i < NumComponents; i++)
		{
			glVertex3f(Vertices[i][0], Vertices[i][2], Vertices[i][1]);
		}
		glEnd();
		
	int i = 1; 
		int j = 3; 

		drawLine(Vertices[i][0], Vertices[i][2], Vertices[i][1], Vertices[j][0], Vertices[j][2], Vertices[j][1]);
		i = 3; 
		j = 4; 

		drawLine(Vertices[i][0], Vertices[i][2], Vertices[i][1], Vertices[j][0], Vertices[j][2], Vertices[j][1]);
		i = 4; 
		j = 6; 

		drawLine(Vertices[i][0], Vertices[i][2], Vertices[i][1], Vertices[j][0], Vertices[j][2], Vertices[j][1]);

		i = 6; 
		j = 1; 

		drawLine(Vertices[i][0], Vertices[i][2], Vertices[i][1], Vertices[j][0], Vertices[j][2], Vertices[j][1]);

		i = 0; 
		j = 2; 

		drawLine(Vertices[i][0], Vertices[i][2], Vertices[i][1], Vertices[j][0], Vertices[j][2], Vertices[j][1]);

		i = 2;
		j = 5;

		drawLine(Vertices[i][0], Vertices[i][2], Vertices[i][1], Vertices[j][0], Vertices[j][2], Vertices[j][1]);

		i = 5;
		j = 7;

		drawLine(Vertices[i][0], Vertices[i][2], Vertices[i][1], Vertices[j][0], Vertices[j][2], Vertices[j][1]);

		i = 7;
		j = 0;

		drawLine(Vertices[i][0], Vertices[i][2], Vertices[i][1], Vertices[j][0], Vertices[j][2], Vertices[j][1]);

		i = 0;
		j = 1;

		drawLine(Vertices[i][0], Vertices[i][2], Vertices[i][1], Vertices[j][0], Vertices[j][2], Vertices[j][1]);

		i = 2;
		j = 3;

		drawLine(Vertices[i][0], Vertices[i][2], Vertices[i][1], Vertices[j][0], Vertices[j][2], Vertices[j][1]);

		i = 5;
		j = 4;

		drawLine(Vertices[i][0], Vertices[i][2], Vertices[i][1], Vertices[j][0], Vertices[j][2], Vertices[j][1]);

		i = 7;
		j = 6;

		drawLine(Vertices[i][0], Vertices[i][2], Vertices[i][1], Vertices[j][0], Vertices[j][2], Vertices[j][1]);

		

		glBegin(GL_LINES);
			glVertex3f(CurrentCenterOfMass[0], CurrentCenterOfMass[2], CurrentCenterOfMass[1]);
			glVertex3f(transformedNode[0], transformedNode[2], transformedNode[1]);
		glEnd();

		glColor3f(0.0f, 0.0f, 1.0f);

		glBegin(GL_POINTS);
			glVertex3f(CommandedPlace[0], CommandedPlace[2], CommandedPlace[1]);
			glEnd();

			//Insert Windows Sleep Function?
			Sleep(.5);
			

			//glReadPixels(0, 0, 640, 480, GL_RGB, GL_FLOAT,  );

		


		glutSwapBuffers();  // Swap the front and back frame buffers (double buffering)

}

/* Handler for window re-size event. Called back when the window first appears and
   whenever the window is re-sized with its new width and height */
void reshape(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
   // Compute aspect ratio of the new window
	if (height == 0) height = 1;                // To prevent divide by 0
	GLfloat aspect = (GLfloat)width / (GLfloat)height;

	// Set the viewport to cover the new window
	glViewport(0, 0, width, height);

	// Set the aspect ratio of the clipping volume to match the viewport
	glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
	glLoadIdentity();             // Reset
	// Enable perspective projection with fovy, aspect, zNear and zFar
	gluPerspective(45.0f, aspect, 0.1f, 10000.0f);
}

float t = 0; 

void idle() {
	
	for (int i = 0; i < NumComponents; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			ReadingDude >> Vertices[i][j];
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			ReadingDude >> CurrentRotation[i][j];
		}
	}

	for (int i = 0; i < 3; i++)
	{
		ReadingDude >> CurrentCenterOfMass[i];
	}

	for (int i = 0; i < 3; i++)
	{
		ReadingDude >> CommandedPlace[i];
	}

	MatrixOperation(CurrentRotation, testNode, transformedNode); 
	VSum(transformedNode, CurrentCenterOfMass); 
	

	glutPostRedisplay(); 
}


void MatrixOperation(double M[3][3], double V[3], double OutV[3])
{
	for (int i = 0; i < 3; i++)
	{
		OutV[i] = 0;
		for (int j = 0; j < 3; j++)
		{
			OutV[i] += M[i][j] * V[j];
		}
	}
}

void VSum(double V1[3], double V2[3])
{
	for (int i = 0; i < 3; i++)
	{
		V1[i] += V2[i];
	}
}

/* Main function: GLUT runs as a console application starting at main() */
int main(int argc, char** argv) {

	createPropeller(p1_V, p1_x, p1_y, p1_z); 
	createPropeller(p2_V, p2_x, p2_y, p2_z);
	createPropeller(p3_V, p3_x, p3_y, p3_z);
	createPropeller(p4_V, p4_x, p4_y, p4_z);

	ReadingDude.open("ParticlePositions.txt");
	NumComponents = 8;

	glutInit(&argc, argv);            // Initialize GLUT
	glutInitDisplayMode(GLUT_DOUBLE| GLUT_RGB); // Enable double buffered mode
	glutInitWindowSize(640, 480);   // Set the window's initial width & height
	glutInitWindowPosition(0, 0); // Position the window's initial top-left corner
	glutCreateWindow(title);          // Create window with the given title
	glutDisplayFunc(display);       // Register callback handler for window re-paint event
	glutReshapeFunc(reshape); // Register callback handler for window re-size event
	system("pause");
	glutIdleFunc(idle);
	initGL();                       // Our own OpenGL initialization
	glutMainLoop();                 // Enter the infinite event-processing loop
	return 0;
}


