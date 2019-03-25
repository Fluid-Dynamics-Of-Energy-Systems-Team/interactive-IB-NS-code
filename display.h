


#ifndef DISPLAY_H_
#define DISPLAY_H_




#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "math.h"
#include <ctime>

#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include <algorithm>
using namespace std;

#ifdef LINUX
  #include "GL/gl.h"
  #include "GL/glu.h"
  #include "GL/glut.h"
#else
  #include "/opt/X11/include/GL/gl.h"
  #include "/opt/X11/include/GL/glu.h"
  #include "/opt/X11/include/GL/glut.h"
#endif


#include "point.h"
#include "ib.h"
#include "globals.h"
#include "main.h"



//     x---------- 3/1
//     |            |
//     |            |
//   -3/-1 -------- x
//
//            0      3   *2
//            1      3   *2
//



bool codeCrash = false;

float scaleDisplay = 2.0;
float translateX = -LoH/scaleDisplay;
float translateY = -1.0/scaleDisplay;
//float scaleDisplay = 1.0;
//float translateX = 0.0;
//float translateY = 0.0;


int width = 600*3;
int height = 600;
float aspectRatio = float(height)/float(width);
int mouseDownX, mouseDownY, modifierKey;
bool rightButton, leftButton, middleButton, ctrKey;
//bool isDragging;
bool drawing = false;


deque <POINT> ptsInteractive;

bool postpro = false;



struct MyColor {
	float r,g,b;
};
class ColorMap {
public:
	string name;
	vector<MyColor> c;
public:
	MyColor getColor(float v) {

		if (v >= 1.0) return c[c.size()-1];
		if (v <= 0.0) return c[0];

		float  ci = v*(c.size()-1);
		int   ici = int(ci);
		float dci = (ci - ici);

		MyColor cret;
		cret.r = c[ici].r + dci*(c[ici+1].r-c[ici].r);
		cret.g = c[ici].g + dci*(c[ici+1].g-c[ici].g);
		cret.b = c[ici].b + dci*(c[ici+1].b-c[ici].b);
		//        cout << ci << "\t" << ci - int(ci) << "\t" << endl;
		return cret;
	}
};
deque<ColorMap> cmap;




//#define MAXOBJS 10000
//#define MAXSELECT 100
//#define MAXFEED 300
//GLuint selectBuf[MAXSELECT];
//GLint viewport[4];




void initDisplay(int argc, char** argv);
void display();
void drawContour();
void drawIBObjects(int mode);
void drawArrow(POINT &p1, POINT &t);
void keyboard(unsigned char key, int x, int y);
void SpecialInput(int key, int x, int y);
void reshape(GLsizei w, GLsizei h);
void mouseButtonEvent(int button, int state, int x, int y);
void mouseButtonEventMove(int x, int y);
void print_bitmap_string(void* font, char* s);
void displayInfo();
//GLint DoSelect(GLint x, GLint y);
void loadMyColorMap(string fname);
void writeTableEntry(string s1, string s2, int entry, float offsetX, float offset);




void initDisplay(int argc, char** argv)  {

	glutInit(&argc, argv);                      // Initialize GLUT
	glutInitWindowSize(width, height);               // Set the window's initial width & height - non-square
	glutInitWindowPosition(50, 50);             // Position the window's initial top-left corner
	glutCreateWindow("NS solver");              // Create window with the given title
	glutDisplayFunc(display);                     // Register display callback handler for window re-paint
	glutReshapeFunc(reshape);                   // Register callback handler for window re-size event
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(SpecialInput);
	glutIdleFunc(run);
	glutMouseFunc(mouseButtonEvent);
	glutMotionFunc(mouseButtonEventMove);

	loadMyColorMap("colormaps/pink.txt");
	loadMyColorMap("colormaps/bone.txt");
	//    loadMyColorMap("colormaps/bone2.txt");
	//    loadMyColorMap("colormaps/bone3.txt");
	//    loadMyColorMap("colormaps/bone4.txt");
	//    loadMyColorMap("colormaps/bone5.txt");
	loadMyColorMap("colormaps/copper.txt");
	loadMyColorMap("colormaps/jet.txt");
	loadMyColorMap("colormaps/cool.txt");

}


void loadMyColorMap(string fname) {

	FILE *fp;
	if ((fp = fopen(fname.c_str(), "rt")) != NULL) {
		cout << "reading file: " << fname << endl;
		POINT pt;
		ColorMap map;
		map.name = fname;
		MyColor c;
		while (fscanf(fp, "%f,%f,%f\n", &c.r, &c.g, &c.b) != EOF) map.c.push_back(c);
		cmap.push_back(map);
	}
	else { cout << "file not found: " << fname << endl; }
}



void reshape(GLsizei w, GLsizei h)
{
	if (h == 0) height = 1;
	aspectRatio = (GLfloat)w / (GLfloat)h;

	// Set the viewport to cover the new window
	glViewport(0, 0, w, h);

	// Set the aspect ratio of the clipping area to match the viewport
	glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
	glLoadIdentity();             // Reset the projection matrix

	gluOrtho2D(-1.0*aspectRatio, 1.0*aspectRatio, -1.0, 1.0);

	//    glGetIntegerv(GL_VIEWPORT, viewport);

	width = w;
	height = h;
}


void print_bitmap_string(void* font, char* s)
{
	if (s && strlen(s)) {
		while (*s) {
			glutBitmapCharacter(font, *s);
			s++;
		}
	}
}



void display() {

	glClearColor(178.0f/255.0f, 178.0f/255.0f, 178.0f/255.0f, 0.0f);

	glClear(GL_COLOR_BUFFER_BIT);   // Clear the color buffer with current clearing color

	glPushMatrix();
	//    glTranslatef(-2.8333, 0.0, 0.0); // for loh=3
	//    glTranslatef(-1.5+translateX, -0.5+translateY, 0.0); // for loh=3
	glScalef(scaleDisplay, scaleDisplay, 1.0);
	glTranslatef(translateX, translateY, 0.0); // for loh=3


	//    glTranslatef(-1.1, -0.5, 0.0);

	drawContour();
	drawIBObjects(GL_RENDER);

	displayInfo();


	glPopMatrix();
	glFlush();  // Render now
}


void writeTableEntry(char *s1, char *s2, float offsetX, float offsetY) {
	glRasterPos2f(-aspectRatio + 30.0/height,             0.9 - offsetY);
	print_bitmap_string(GLUT_BITMAP_HELVETICA_12, s1);

	glRasterPos2f(-aspectRatio + 30.0/height + offsetX, 0.9 - offsetY);
	print_bitmap_string(GLUT_BITMAP_HELVETICA_12, s2);
}

void displayInfo() {

	char name[200];
	char s1[200];
	char s2[200];

	glTranslatef(-translateX, -translateY, 0.0); // for loh=3
	glScalef(1./scaleDisplay, 1./scaleDisplay, 1.0);

	//    glColor4f(0.0, 0.5, 1.0, 0.0);
	glColor4f(1.0, 1.0, 1.0, 0.0);
	float offsetX = 65.0/height;
	float offsetY = 28.0/height;
	string str;


	sprintf(s1, "Use these keys to change settings:");
	writeTableEntry(s1, "\0", offsetX, 0.0*offsetY-0.02);
	sprintf(s1, "................................................................");
	writeTableEntry(s1, "\0", offsetX, 0.0*offsetY-0.003);


	sprintf(s1, "1-9");
	sprintf(s2, "|  load predefined cases");
	writeTableEntry(s1, s2, offsetX, 1.0*offsetY);


	sprintf(s1, "0");
	sprintf(s2, "|  delete immersed boundaries");
	writeTableEntry(s1, s2, offsetX, 2.0*offsetY);


	sprintf(s1, "R, r");
	sprintf(s2, "|  Reynolds number = %.0f", Re);
	writeTableEntry(s1, s2, offsetX, 3.0*offsetY);


	sprintf(s1, "A, a");
	switch (advecScheme) {
	case 1: str = "central"; break;
	case 2: str = "central (every 40th step QUICK)"; break;
	case 3: str = "central (every 20th step QUICK)"; break;
	case 4: str = "central (every 5th step QUICK)"; break;
	case 5: str = "QUICK"; break;     }
	sprintf(s2, "|  advection: %s", str.c_str());
	writeTableEntry(s1, s2, offsetX, 4.0*offsetY);


	sprintf(s1, "p");
	str = (periodicX == false) ? "inlet/outlet" : "periodic";
	sprintf(s2, "|  boundary condition in X: %s", str.c_str());
	writeTableEntry(s1, s2, offsetX, 5.0*offsetY);


	sprintf(s1, "w");
	switch (wall) {
	case 1: str = "slip"; break;
	case 2: str = "no-slip"; break;
	case 3: str = "periodic"; break;
	}
	sprintf(s2, "|  boundary condition in Y: %s", str.c_str());
	writeTableEntry(s1, s2, offsetX, 6.0*offsetY);


	sprintf(s1, "i");
	sprintf(s2, "|  initialize flow");
	writeTableEntry(s1, s2, offsetX, 7.0*offsetY);


	sprintf(s1, "v");
	sprintf(s2, "|  toggle vectors");
	writeTableEntry(s1, s2, offsetX, 8.0*offsetY);


	sprintf(s1, "d");
	sprintf(s2, "|  begin/end drawing immersed boundary");
	writeTableEntry(s1, s2, offsetX, 9.0*offsetY);


	sprintf(s1, "f");
	sprintf(s2, "|  begin/end drawing averaging line");
	writeTableEntry(s1, s2, offsetX, 10.0*offsetY);


	sprintf(s1, "+,-");
	sprintf(s2, "|  scale immersed objects");
	writeTableEntry(s1, s2, offsetX, 11.0*offsetY);


	//////




	//    string str = (periodicX == false) ? "false" : "true";
	sprintf(name, "Velocity (max): %.2f", umax);
	glRasterPos2f(-aspectRatio+30.0/height, -0.9+3*offsetY);
	print_bitmap_string(GLUT_BITMAP_HELVETICA_12, name);

	sprintf(name, "Delta t: %.2e", dt);
	glRasterPos2f(-aspectRatio+30.0/height, -0.9+2*offsetY);
	print_bitmap_string(GLUT_BITMAP_HELVETICA_12, name);

	sprintf(name, "Time: %.2f", simtime);
	glRasterPos2f(-aspectRatio+30.0/height, -0.9+1*offsetY);
	print_bitmap_string(GLUT_BITMAP_HELVETICA_12, name);

	sprintf(name, "Tot. divergence: %.2e", divVolTot);
	glRasterPos2f(-aspectRatio+30.0/height, -0.9+0*offsetY);
	print_bitmap_string(GLUT_BITMAP_HELVETICA_12, name);



	if (codeCrash == true) {
		sprintf(name, "You broke the code! Press \"i\" to reinitialize the flow");
		glRasterPos2f(-0.2, 0.0 );
		print_bitmap_string(GLUT_BITMAP_HELVETICA_18, name);
	}

}


void drawContour() {


	codeCrash = false;
	if (u[ind(1,1)] != u[ind(1,1)]) {
		codeCrash = true;
		return;
	}

	float uNode[imax-1][jmax-1], vNode[imax-1][jmax-1], vel[imax-1][jmax-1], vort[imax-1][jmax-1];

	float velMaxMin[2]  = {-1.0e10, 1.0e10};
	float vortMaxMin[2] = {-1.0e10, 1.0e10};

	for (int i=0; i<imax-1; i++)
		for (int j=0; j<jmax-1; j++)
		{
			uNode[i][j] = 0.5*(u[ind(i,j+1)] + u[ind(i,j)]);
			vNode[i][j] = 0.5*(v[ind(i+1,j)] + v[ind(i,j)]);

			vel[i][j] = sqrt(pow(uNode[i][j], 2.0) + pow(vNode[i][j], 2.0));
			velMaxMin[0] = MAX(velMaxMin[0], fabs(vel[i][j]));
			velMaxMin[1] = MIN(velMaxMin[1], fabs(vel[i][j]));
		}
	velMaxMin[0] = velMaxMin[0]+1.0e-6;

	for (int i=0; i<imax-1; i++)
		for (int j=0; j<jmax-1; j++)
			vort[i][j] = (vel[i][j]-velMaxMin[1])/(velMaxMin[0]-velMaxMin[1]);



	//    for (int i=1; i<imax-1; i++)
		//    for (int j=1; j<jmax-1; j++)
			//    {
	//        vort[i-1][j-1] = (-(u[ind(i,j)] - u[ind(i,j-1)])/dx + (v[ind(i,j)] - v[ind(i-1,j)])/dx);
	//        vortMaxMin[0] = max(vortMaxMin[0], vort[i-1][j-1]);
	//        vortMaxMin[1] = min(vortMaxMin[1], vort[i-1][j-1]);
	//    }
	//
	//    vortMaxMin[0] = contourMax;
	//    vortMaxMin[1] = contourMin;
	//
	//    for (int i=0; i<imax-2; i++)
	//    for (int j=0; j<jmax-2; j++)
	//        vort[i][j] = (vort[i][j]-vortMaxMin[1])/(vortMaxMin[0]-vortMaxMin[1]);




	//    cout << vortMaxMin[0] << "/" << vortMaxMin[1] << endl;



	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glBegin(GL_QUADS);
	for (int i=0; i<imax-2; i++)
		for (int j=0; j<jmax-2; j++)
		{
			MyColor c = cmap[colorMap].getColor(vort[i][j]);
			glColor3f(c.r, c.g, c.b);
			glVertex3f(dx*i,     dx*j,     0.0);

			c = cmap[colorMap].getColor(vort[i+1][j]);
			glColor3f(c.r, c.g, c.b);
			glVertex3f(dx*(i+1), dx*j,     0.0);

			c = cmap[colorMap].getColor(vort[i+1][j+1]);
			glColor3f(c.r, c.g, c.b);
			glVertex3f(dx*(i+1), dx*(j+1), 0.0);

			c = cmap[colorMap].getColor(vort[i][j+1]);
			glColor3f(c.r, c.g, c.b);
			glVertex3f(dx*i,     dx*(j+1), 0.0);
		}
	glEnd();


	if (vectorToggle == true) {
		float scaleVec = 0.04;
		glColor3f(0.0, 0.0, 0.0);
		//        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		//        glEnable (GL_BLEND);
		//        glEnable (GL_LINE_SMOOTH);
		//        glHint (GL_LINE_SMOOTH_HINT, GL_NICEST);
		glLineWidth(0.5f);
		int step = 8;
		for (int i=0; i<imax-1; i=i+step)
			for (int j=0; j<jmax-1; j=j+step)
			{
				float x = scaleVec*(vel[i][j]-velMaxMin[1])/(velMaxMin[0]-velMaxMin[1]);
				POINT p1(dx*i, dx*j, 0.01);
				POINT t(x*uNode[i][j]/vel[i][j], x*vNode[i][j]/vel[i][j], 0.0);

				drawArrow(p1, t);
			}
	}
}


void drawArrow(POINT &p1, POINT &t) {

	POINT p2 = p1 + t;
	t.rotateDegZ(20);
	POINT p3 = p2 - 0.4*t;
	t.rotateDegZ(-40);
	POINT p4 = p2 - 0.4*t;

	glBegin(GL_LINE_STRIP);
	glVertex3f(p1.x, p1.y, 0.05);
	glVertex3f(p2.x, p2.y, 0.05);
	glEnd();

	glBegin(GL_LINE_STRIP);
	glVertex3f(p3.x, p3.y, 0.05);
	glVertex3f(p2.x, p2.y, 0.05);
	glVertex3f(p4.x, p4.y, 0.05);
	glEnd();
}



void drawIBObjects(int mode) {


	//    for (int obj=0; obj<ibObj.size(); obj++) {
	//
	////        GLint i;
	////        if (mode == GL_SELECT)  glLoadName(i);
	//
	//        glBegin(GL_QUADS);
	//        for (int ib=0; ib<ibObj[obj].iIB.size(); ib++) {
	//            int i = ibObj[obj].iIB[ib].i-1;
	//            int j = ibObj[obj].iIB[ib].j-1;
	//
	////            if (colorMap == 7)
	//            glColor3f(0.0, 1.0, 0.0);
	////            else                     glColor3f(1.0, 1.0, 1.0);
	//
	//            glVertex3f(dx*(i-0.0), dx*(j-0.0), 0.0);
	//            glVertex3f(dx*(i+1.0), dx*(j-0.0), 0.0);
	//            glVertex3f(dx*(i+1.0), dx*(j+1.0), 0.0);
	//            glVertex3f(dx*(i-0.0), dx*(j+1.0), 0.0);
	//        }
	//        glEnd();
	//    }

	//    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);


	for (int obj=0; obj<ibObj.size(); obj++) {
		for (int el=0; el<ibObj[obj].elem.size(); el++) {
			Triangle tri = ibObj[obj].elem[el];
			glLineWidth(1.0);
			glColor3f(0.0, 0.7, 0.7);
			glBegin(GL_TRIANGLES);
			glVertex3f(ibObj[obj].node[tri.n1].x, ibObj[obj].node[tri.n1].y, 0.01);
			glVertex3f(ibObj[obj].node[tri.n2].x, ibObj[obj].node[tri.n2].y, 0.01);
			glVertex3f(ibObj[obj].node[tri.n3].x, ibObj[obj].node[tri.n3].y, 0.01);
			glEnd();
		}
	}

	for (int obj=0; obj<ibObj.size(); obj++) {
		for (int el=0; el<ibObj[obj].elem.size(); el++) {
			Triangle tri = ibObj[obj].elem[el];
			glLineWidth(1.0);
			glColor3f(1.0, 0.0, 0.0);
			glBegin(GL_LINE_STRIP);
			glVertex3f(ibObj[obj].node[tri.n1].x, ibObj[obj].node[tri.n1].y, 0.01);
			glVertex3f(ibObj[obj].node[tri.n2].x, ibObj[obj].node[tri.n2].y, 0.01);
			glVertex3f(ibObj[obj].node[tri.n3].x, ibObj[obj].node[tri.n3].y, 0.01);
			glVertex3f(ibObj[obj].node[tri.n1].x, ibObj[obj].node[tri.n1].y, 0.01);
			glEnd();
		}
	}


	for (int obj=0; obj<ibObj.size(); obj++) {
		for (int ib=0; ib<ibObj[obj].pts.size()-1; ib++) {
			glLineWidth(3.0);
			glColor3f(1.0, 0.0, 0.0);
			glBegin(GL_LINES);
			glVertex3f(ibObj[obj].pts[ib  ].x, ibObj[obj].pts[ib  ].y, 0.01);
			glVertex3f(ibObj[obj].pts[ib+1].x, ibObj[obj].pts[ib+1].y, 0.01);
			glEnd();
		}
	}

	//            glLineWidth(3.0);
	//            glColor3f(1.0, 0.0, 0.0);
	//    for (int obj=0; obj<ibObj.size(); obj++) {
		//        for (int ib=0; ib<ibObj[obj].iIB.size(); ib++) {
			//            POINT p1(ibObj[obj].iIB[ib].i*dx, ibObj[obj].iIB[ib].j*dx);
	//            POINT p2(ibObj[obj].iIB[ib].nvec);
	////            POINT p2(0.2, 0.2);
	////            p2 = p2 * 100.0;
	//            drawArrow(p1, p2);
	//        }
	//    }


	glLineWidth(2.0);
	glColor3f(1.0, 1.0, 0.0);
	glBegin(GL_LINE_STRIP);
	for (int i=0; i<ptsInteractive.size(); i++)
		glVertex3f(ptsInteractive[i].x, ptsInteractive[i].y, 0.01);
	glEnd();


	glPointSize(5.0);
	glColor3f(1.0, 0.0, 1.0);
	glBegin(GL_POINTS);
	for (int i=0; i<ptsInteractive.size(); i++)
		glVertex3f(ptsInteractive[i].x, ptsInteractive[i].y, 0.01);
	glEnd();


	glColor3f(0.0, 1.0, 0.0);
	glLineWidth(0.5f);
	for (int i=0; i<objPostPro.size(); i++) {
		for (int p=0; p<objPostPro[i].iPP.size(); p++) {

			Index ix = objPostPro[i].iPP[p];
			POINT p1(dx*(ix.i+0.5), dx*(ix.j+0.5), 0.01);
			POINT p2 = objPostPro[i].velAverage[p]*0.08/objPostPro[i].umax;
			drawArrow(p1, p2);
		}
	}
}





void keyboard(unsigned char key, int x, int y)
{
	cout << "key pressed: " << key << endl;

	switch (key)
	{
	//----------------------------------------------------------------------
	case '1':
	{
		ibObj.clear();
		IBObject obj = IBObject("objects/nozzle1.txt", dx, imax, jmax);

		obj.translate(POINT(-0.5, 0.0));
		ibObj.push_back(obj);
		//            obj.mirrorY();
		//            ibObj.push_back(obj);
		//            vInlet = 0.0;
	}

	break;

	case '2':
	{
		ibObj.clear();
		IBObject obj = IBObject("objects/nozzle2.txt", dx, imax, jmax);

		ibObj.push_back(obj);
		obj.mirrorY();
		ibObj.push_back(obj);
		vInlet = 0.0;
	}
	break;

	case '3':
	{
		ibObj.clear();
		IBObject obj = IBObject("objects/circ.txt", dx, imax, jmax);
		obj.scale(0.5);
		ibObj.push_back(obj);
		vInlet = 0.0;
	}
	break;

	case '4':
	{
		ibObj.clear();
		IBObject obj = IBObject("objects/circ.txt", dx, imax, jmax);
		obj.scale(0.3);
		obj.translate(POINT(-0.1, 0.0));
		ibObj.push_back(obj);

		obj.translate(POINT(0.5, 0.2));
		ibObj.push_back(obj);

		obj.translate(POINT(0.0, -2.0*0.2));
		ibObj.push_back(obj);
		vInlet = 0.0;
	}
	break;

	case '5':
	{
		ibObj.clear();
		IBObject obj = IBObject("objects/triangle.txt", dx, imax, jmax);
		obj.scale(2);
		ibObj.push_back(obj);
		vInlet = 0.0;
	}
	break;

	case '6':
	{
		ibObj.clear();
		ibObj.push_back(IBObject("objects/blade.txt", dx, imax, jmax));
		ibObj.push_back(ibObj[0]);
		ibObj.push_back(ibObj[0]);
		ibObj[1].translate(POINT(0.0,  0.5));
		ibObj[2].translate(POINT(0.0, -0.5));

		vInlet = 1.2;
		wall = 3;
	}
	break;

	case '7':
	{
		ibObj.clear();
		IBObject obj = IBObject("objects/naca.txt", dx, imax, jmax);
		obj.translate(POINT(0.3, 0.3));
		obj.scale(0.7);
		ibObj.push_back(obj);
		vInlet = 0.0;

		vInlet = 0.1;
		wall = 3;
	}
	break;

	case '8':
	{
		ibObj.clear();
		IBObject obj = IBObject("objects/porsche.txt", dx, imax, jmax);
		obj.scale(0.2, false);
		obj.translate(POINT(0.4, 0.015));
		ibObj.push_back(obj);

		obj = IBObject("objects/circ.txt", dx, imax, jmax);
		obj.scale(0.4);
		obj.translate(POINT(0.148, -0.42));
		ibObj.push_back(obj);

		obj.translate(POINT(0.562, 0.0));
		ibObj.push_back(obj);

		vInlet = 0.0;
		wall = 2;
	}
	break;

	case '9':
	{
		ibObj.clear();
		IBObject obj = IBObject("objects/ferrari_fxxk.txt", dx, imax, jmax);
		obj.scale(1.40, false);
		obj.translate(POINT(0.32, 0.03));
		ibObj.push_back(obj);

		obj = IBObject("objects/circ.txt", dx, imax, jmax);
		obj.scale(0.4);
		obj.translate(POINT(0.136, -0.42));
		ibObj.push_back(obj);

		obj.translate(POINT(0.672, 0.0));
		ibObj.push_back(obj);

		vInlet = 0.0;
		wall = 2;
	}
	break;


	case '0':
		cout << "deleting IB" << endl;
		ibObj.clear();
		objPostPro.clear();
		break;

		//----------------------------------------------------------------------
	case 'R':
		Re += pow(10.0, int(log10(Re)));
		cout << "Re = " << Re << endl;
		break;

	case 'r':
		Re = MAX(Re - pow(10.0, int(log10(Re-0.001))), 1.0);
		cout << "Re = " << Re << endl;
		break;
		//----------------------------------------------------------------------
	case 'C':
		CFL += 0.05;
		cout << "CFL = " << CFL << endl;
		break;

	case 'c':
		CFL -= 0.05;
		cout << "CFL = " << CFL << endl;
		break;
		//----------------------------------------------------------------------
	case 'M':
		colorMap = MIN(colorMap+1, cmap.size()-1);
		cout << "colorMap = " << colorMap << endl;
		break;

	case 'm':
		colorMap = MAX(colorMap-1, 0);
		cout << "colorMap = " << colorMap << endl;
		break;
		//----------------------------------------------------------------------
	case 'T':
		timeInt = min(timeInt+1, 4);
		cout << "timeInt = " << timeInt << endl;
		break;

	case 't':
		timeInt = max(timeInt-1, 1);
		cout << "timeInt = " << timeInt << endl;
		break;
		//----------------------------------------------------------------------
	case 'A':
		advecScheme = min(advecScheme+1, 5);
		cout << "advecScheme = " << advecScheme << endl;
		break;

	case 'a':
		advecScheme = max(advecScheme-1, 1);
		cout << "advecScheme = " << advecScheme << endl;
		break;
		//----------------------------------------------------------------------
	case '+':
		for (int i=0; i<ibObj.size(); i++)
			ibObj[i].scale(1.1, true);
		cout << "scale = 1.1" << endl;
		//        contourMax += 1;
		break;

	case '-':
		for (int i=0; i<ibObj.size(); i++)
			ibObj[i].scale(0.9, true);
		cout << "scale = 0.9" << endl;
		break;
		////----------------------------------------------------------------------
		//    case '=':
		//        contourMin += 1;
		//        cout << "contourMin = " << contourMin << endl;
		//        break;
		//
		//    case '-':
		//        contourMin -= 1;
		//        cout << "contourMin = " << contourMin << endl;
		//        break;
		//----------------------------------------------------------------------
	case 'w':
		wall = (wall+1)%3+1;
		cout << "wall = " << wall << endl;
		break;
		break;


		//----------------------------------------------------------------------
	case ']':
		value += 0.025; //min(wall+1, 3);
		cout << "value = " << value << endl;

		ibObj.clear();
		ibObj.push_back(IBObject("objects/blade.txt", dx, imax, jmax));
		ibObj[0].scale(value);
		ibObj.push_back(ibObj[0]);
		ibObj.push_back(ibObj[0]);
		ibObj[1].translate(POINT(0.0,  0.5));
		ibObj[2].translate(POINT(0.0, -0.5));

		vInlet = 1.2;
		wall = 3;
		break;

	case '\\':
		value -= 0.025;
		cout << "value = " << value << endl;

		ibObj.clear();
		ibObj.push_back(IBObject("objects/blade.txt", dx, imax, jmax));
		ibObj[0].scale(value);
		ibObj.push_back(ibObj[0]);
		ibObj.push_back(ibObj[0]);
		ibObj[1].translate(POINT(0.0,  0.5));
		ibObj[2].translate(POINT(0.0, -0.5));

		vInlet = 1.2;
		wall = 3;
		break;

		//----------------------------------------------------------------------
	case '[':
		vInlet += 0.025; //min(wall+1, 3);
		cout << "vInlet = " << vInlet << endl;
		break;

	case '\'':
		vInlet -= 0.025;
		cout << "vInlet = " << vInlet << endl;
		break;


		//----------------------------------------------------------------------
	case 'p':
		periodicX = !periodicX;
		cout << "periodicX = " << periodicX << endl;
		break;

	case 'z':
		change = !change;
		cout << "change = " << change << endl;
		break;

		//----------------------------------------------------------------------
	case 27:
		runSim = !runSim;
		cout << " runSim = " << runSim << endl;
		break;
		//----------------------------------------------------------------------
	case 'v':
		vectorToggle = !vectorToggle;
		cout << "vectorToggle = " << vectorToggle << endl;
		break;
		//----------------------------------------------------------------------
	case 'd':

		if ((drawing == true) && (ptsInteractive.size()>0)) {
			ibObj.push_back(IBObject(ptsInteractive, dx, imax, jmax));
			ptsInteractive.clear();

			glutSetCursor(GLUT_CURSOR_INHERIT);
		}
		else
			glutSetCursor(GLUT_CURSOR_CROSSHAIR);


		drawing = !drawing;
		cout << "drawing = " << drawing << endl;
		break;

		//----------------------------------------------------------------------
	case 'f':

		if ((postpro == true) && (ptsInteractive.size()>0)) {
			objPostPro.push_back(PostProcess(ptsInteractive, dx, imax, jmax, umax));
			ptsInteractive.clear();

			glutSetCursor(GLUT_CURSOR_INHERIT);
		}
		else
			glutSetCursor(GLUT_CURSOR_INFO);


		postpro = !postpro;
		cout << "postpro = " << postpro << endl;
		break;

		//----------------------------------------------------------------------
	case 'i':
		initFlowField();
		dt = calcDT();
		cout << "init flow field" << endl;
		break;

		//----------------------------------------------------------------------
		//    case 32:
		//        run();
		//        display();
		//        cout << "run" << endl;
		//        break;
		//----------------------------------------------------------------------
	}

	glutPostRedisplay();

}





void mouseButtonEvent(int button, int state, int x, int y)
{
	GLint hit;

	//    if (state == GLUT_DOWN) {
	//        hit = DoSelect((GLint) x, (GLint) y);
	//
	//        if (hit != -1) {
	//            if (button == GLUT_LEFT_BUTTON)
	//                cout << hit << endl;
	//                //{
	////                RecolorTri(hit);
	////            } else if (button == GLUT_MIDDLE_BUTTON) {
	////                GrowTri(hit);
	////            } else if (button == GLUT_RIGHT_BUTTON) {
	////                DeleteTri(hit);
	////            }
	////            glutPostRedisplay();
	//        }
	//    }


	modifierKey = glutGetModifiers();

	mouseDownX = x;
	mouseDownY = y;

	float coordX = (2.0*float(x)/float(width)*aspectRatio - aspectRatio);
	float coordY = (2.0*float(height-y)/height) - 1.0;
	coordX = (coordX - translateX*scaleDisplay)/scaleDisplay;
	coordY = (coordY - translateY*scaleDisplay)/scaleDisplay;

	//    cout << "mouse: "<<width<<"/"<<height<<  "\t"  <<winx<<"/"<<winy<<  "\t"  <<x<<"/"<<y<<  "\t";
	//    cout << translateX<<"/"<<translateY << "\t"<< scaleDisplay <<"\t"<<coordX<<"/"<<coordY<<endl;

	//    if (modifierKey == GLUT_ACTIVE_CTRL)
	//        ptsInteractive.push_back(POINT(coordX, coordY));

	if (state == GLUT_DOWN) {
		if ((drawing == true) || (postpro == true)) {
			ptsInteractive.push_back(POINT(coordX, coordY));
		}
	}

	leftButton   = (button == GLUT_LEFT_BUTTON); //   && (state == GLUT_DOWN));
	rightButton  = (button == GLUT_RIGHT_BUTTON); //  && (state == GLUT_DOWN));
	middleButton = (button == GLUT_MIDDLE_BUTTON); // && (state == GLUT_DOWN));


	//    if (state == )

	//    isDragging = true;                      // Prepare For Dragging

	//    if (modifierKey == GLUT_ACTIVE_CTRL)
	//        glutSetCursor(GLUT_CURSOR_CROSSHAIR);

	glutPostRedisplay();

}




void mouseButtonEventMove(int x, int y)
{
	if (leftButton) // && (modifierKey != GLUT_ACTIVE_CTRL))
	{
		translateX += (float)(x-mouseDownX)/width*aspectRatio*2.0/scaleDisplay;
		translateY += (float)(mouseDownY-y)/width*aspectRatio*2.0/scaleDisplay;
	} // translate

	if (rightButton) // && (modifierKey == GLUT_ACTIVE_CTRL))
	{
		scaleDisplay += (float)(y-mouseDownY)/width*aspectRatio*scaleDisplay;
		scaleDisplay = MAX(scaleDisplay, 0.1);
	} // scale

	mouseDownX = x;
	mouseDownY = y;

	glutPostRedisplay();
}



void SpecialInput(int key, int x, int y) {

	int iIBstart = 50;

	switch (key) {
	//    case GLUT_ACTIVE_CTRL:
	//        if ((drawing == true) && (ptsInteractive.size()>0)) {
	////            ibObj.push_back(IBObject(ptsInteractive, dx, imax, jmax));
	////            ptsInteractive.clear();
	//
	//            glutSetCursor(GLUT_CURSOR_INHERIT);
	//        }
	//        else
	//            glutSetCursor(GLUT_CURSOR_CROSSHAIR);
	//
	//        break;
	//    case GLUT_KEY_UP:
	//            for (int i=iIBstart; i<iIB.size(); i++)   iIB[i].ind[1] += 1;
	//            printf("%d\t%d\n",x,y);
	//            break;
	//    case GLUT_KEY_DOWN:
	//            for (int i=iIBstart; i<iIB.size(); i++)   iIB[i].ind[1] -= 1;
	//            break;
	//    case GLUT_KEY_LEFT:
	//            for (int i=iIBstart; i<iIB.size(); i++)   iIB[i].ind[0] -= 1;
	//            break;
	//    case GLUT_KEY_RIGHT:
	//            for (int i=iIBstart; i<iIB.size(); i++)   iIB[i].ind[0] += 1;
	//            break;
	}
}



#endif







//GLint DoSelect(GLint x, GLint y)
//{
//    GLint hits;
//return hits;
//    glSelectBuffer(MAXSELECT, selectBuf);
//    glRenderMode(GL_SELECT);
//    glInitNames();
//    glPushName(~0);
//
//    glPushMatrix();
//
//    glMatrixMode(GL_PROJECTION);
//    glLoadIdentity();
//    gluPickMatrix(x, height - y, 4, 4, viewport);
//
//    GLfloat aspect = (GLfloat)width / (GLfloat)height;
//
//    if (width >= height) {
//            // aspect >= 1, set the height from -1 to 1, with larger width
//            gluOrtho2D(-1.0 * aspect, 1.0 * aspect, -1.0, 1.0);
//        } else {
//            //     aspect < 1, set the width to -1 to 1, with larger height
//            gluOrtho2D(-1.0, 1.0, -1.0 / aspect, 1.0 / aspect);
//        }
//    glMatrixMode(GL_MODELVIEW);
//
//    glClearColor(1.0, 1.0, 1.0, 1.0);
//    glClear(GL_COLOR_BUFFER_BIT);
//
//    float scale = 1.8*scaleDisplay;
////    glScalef(scale, scale, scale);
////    glTranslatef(-1.5+translateX, -0.5+translateY, 0.0); // for loh=3
//
//
//    drawContour();
//    drawIBObjects(GL_SELECT);
//
////    glScalef(1./scale, 1./scale, 1./scale);
////    glTranslatef(1.5-translateX, 0.5-translateY, 0.0); // for loh=3
//
//    glPopMatrix();
//
////    glPopMatrix();
////    glFlush();  // Render now
//
//
//    hits = glRenderMode(GL_RENDER);
//    if (hits <= 0) {
//        return -1;
//    }
//    return selectBuf[(hits - 1) * 4 + 3];
//}



