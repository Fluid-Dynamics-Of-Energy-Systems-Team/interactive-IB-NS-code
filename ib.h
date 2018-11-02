
#ifndef IB_H_
#define IB_H_






#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include <iostream>
#include <string>
#include <deque>
#include <algorithm>
using namespace std;



#include "point.h"
#include "defines.h"

#define REAL float
#define VOID int
#define ANSI_DECLARATORS
#include "triangle.h"






// class index contains i and j values, used to locate immersed boundary
class Index {
public:
	int i;
	int j;

	POINT nvec;

public:
	Index(): i(0), j(0) {}
	Index(int i_, int j_): i(i_), j(j_) {}
	Index(int i_, int j_, POINT n): i(i_), j(j_) {nvec = n;}
	Index(const Index &ind2) {i = ind2.i; j = ind2.j;}
	Index(const Index &ind2, POINT n) {i = ind2.i; j = ind2.j; nvec = n;}

	bool operator ==(Index &ind2) { return ((i==ind2.i) && (j==ind2.j)); }
	bool operator !=(Index &ind2) { return ((i!=ind2.i) && (j!=ind2.j)); }
	void print(string text="") {cout << text << ",\t" << i << ",\t" << j << endl;}
};

bool equalIndex(Index &ind1, Index &ind2) { return ((ind1.i==ind2.i) && (ind1.j==ind2.j)); }

bool lessThanIndex(const Index &ind1, const Index &ind2) {
	if (ind1.i < ind2.i) return true;
	if (ind2.i < ind1.i) return false;
	return ind2.j > ind1.j;
}





// base Object class
class Object {
public:
	int imax, jmax;
	float dx;
	deque<POINT> pts;    // boundary of the IB

public:
	Object(string fname, float dx_, int imax_, int jmax_) {
		imax = imax_;
		jmax = jmax_;
		dx = dx_;
		readIB(fname);
	}

	Object(deque <POINT> &pts_, float dx_, int imax_, int jmax_) {
		imax = imax_;
		jmax = jmax_;
		dx = dx_;
		pts = pts_;
	}

	Object(const Object& obj2) {
		imax = obj2.imax;
		jmax = obj2.jmax;
		dx = obj2.dx;
		pts = obj2.pts;
	}


	virtual ~Object() {
		pts.clear();
	}


	virtual void translate(POINT pt) {
		for (int i=0; i<pts.size(); i++)   	pts[i] += pt;
	}


	virtual POINT scale(float fact, bool center = true) {

		POINT pt(0.0, 0.0);

		if (center == true) {
			for (int i=0; i<pts.size(); i++) pt += pts[i];
			pt = pt/pts.size();
		}
		else pt = pts[0];

		for (int i=0; i<pts.size(); i++)   	pts[i] = fact*(pts[i] - pt) + pt;
		return pt;
	}


	virtual void scale(float fact, POINT pt) {
		for (int i=0; i<pts.size(); i++)   	pts[i] = fact*(pts[i] - pt) + pt;
	}


	virtual void mirrorY() {
		for (int i=0; i<pts.size(); i++)   	pts[i].y = 1.0 - pts[i].y;
	}


	void readIB(string fname) {

		if (pts.size() != 0) {
			pts.clear();
		}

		FILE *fp;
		if ((fp = fopen(fname.c_str(), "rt")) != NULL) {
			cout << "reading file: " << fname << endl;
			POINT pt;
			while (fscanf(fp, "%f,%f\n", &pt.x, &pt.y) != EOF) pts.push_back(pt);
			printBounds("after reading: ");
		}
		else { cout << "file not found: " << fname << endl; }
	}


	void printBounds(string str) {

		float bounds[4] = {-1.0e6, 1.0e6, -1.0e6, 1.0e6};

		for (int i=0; i<pts.size(); i++) {
			bounds[0] = max(pts[i].x, bounds[0]);
			bounds[1] = min(pts[i].x, bounds[1]);
			bounds[2] = max(pts[i].y, bounds[2]);
			bounds[3] = min(pts[i].y, bounds[3]);
		}

		cout << "IB bounds [xmax, xmin, ymax, ymin], " << str << ": ";
		cout << bounds[0] << ", " << bounds[1] << ", ";
		cout << bounds[2] << ", " << bounds[3] << endl;
	}

};


class PostProcess : public Object {

public:
	deque<Index> iPP;
	deque<POINT> velAverage;
	int count;
	float umax;

public:

	PostProcess(deque<POINT> &pts, float dx, int imax, int jmax, float umax_) : Object(pts, dx, imax, jmax) {
		intersectIBwMesh();
		count = 0;
		umax = umax_;
	}

	~PostProcess() {
		iPP.clear();
	}


	void intersectIBwMesh() {

		if (pts.size() == 0) return;
		if (pts.size() == 1) pts.push_back(pts[0]);

		cout << pts.size() << endl;

		for (int p=1; p<pts.size(); p++) {
			POINT p1(pts[p-1].x, pts[p-1].y);
			POINT p2(pts[p].x, pts[p].y);

			POINT vec = p2-p1;

			for (int i=0; i<=10; i++) {
				float fact = float(i)/10.0;
				POINT p = p1 + fact*vec;
				iPP.push_back(Index(int(p.x/dx)+1, int(p.y/dx)+1));
			}
		}

//		POINT p1(pts[0].x, pts[0].y);
//		POINT p2(pts[1].x, pts[1].y);
//
//		POINT vec = p2-p1;
//
//		for (int i=0; i<=10; i++) {
//			float fact = float(i)/10.0;
//			POINT p = p1 + fact*vec;
//			iPP.push_back(Index(int(p.x/dx)+1, int(p.y/dx)+1));
//		}

		deque<Index> iPPtmp;
		for (int i=0; i<iPP.size(); i++)   {
			if ((iPP[i].i < 1) || (iPP[i].i > imax-1) || (iPP[i].j < 1) || (iPP[i].j > jmax-1)) {}
			else iPPtmp.push_back(iPP[i]);
		}

		iPP.clear();
		iPP = iPPtmp;

		velAverage.resize(iPP.size(), 0.0);


		for (int i=0; i<iPP.size(); i++)   {
			cout << iPP[i].i << ", " << iPP[i].j << endl;
		}
//		sort(iIB.begin(), iIB.end(), lessThanIndex);
//		iIB.erase( unique( iIB.begin(), iIB.end(), equalIndex), iIB.end() );


	}

};



// used for triangulation of the immersed boundary
struct Triangle { int n1, n2, n3; };


// main immersed boundary class
class IBObject : public Object {

public:
	deque<Index> iIB;    		// index of the IB in the mesh, therefore dx, imax and jmax from the solver need to be known
	deque<POINT> node;			// node of triangulated IB
	deque<Triangle> elem;		// element of triangulated IB


	IBObject(string fname, float dx, int imax, int jmax) : Object(fname, dx, imax, jmax) {
		init();
	}


	IBObject(deque<POINT> &pts, float dx, int imax, int jmax) : Object(pts, dx, imax, jmax) {
		init();
	}

	void init() {
		if (node.size() > 0) {
			node.clear();
			elem.clear();
		}

		intersectIBwMesh();
		if ((pts.size() > 2) && (fabs(mag(pts[0] - pts[pts.size()-1])) < 1.0e-2)) tesselate();
	}

	IBObject(const IBObject& obj2) : Object(obj2) {
		iIB = obj2.iIB;
		node = obj2.node;
		elem = obj2.elem;
	}


	~IBObject() {
		iIB.clear();
		node.clear();
		elem.clear();
	}


	void translate(POINT pt) {
		Object::translate(pt);

		iIB.clear();
		intersectIBwMesh();

		for (int i=0; i<node.size(); i++) node[i] += pt;
	}

	POINT scale(float fact, bool center = true) {
		POINT pt = Object::scale(fact, center);

		iIB.clear();
		intersectIBwMesh();

		for (int i=0; i<node.size(); i++)   node[i] = fact*(node[i] - pt) + pt;
		return pt;
	}

	void scale(float fact, POINT pt) {
		Object::mirrorY();

		iIB.clear();
		intersectIBwMesh();

		for (int i=0; i<node.size(); i++)   node[i] = fact*(node[i] - pt) + pt;
	}

	void mirrorY() {
		Object::mirrorY();

		iIB.clear();
		intersectIBwMesh();

		for (int i=0; i<node.size(); i++)   node[i].y = 1.0 - node[i].y;
	}


	void intersectIBwMesh() {

		POINT p1(pts[0].x, pts[0].y);
		Index ind1(int(p1.x/dx), int(p1.y/dx));

		iIB.push_back(Index(ind1.i+1, ind1.j+1));

		for (int i=1; i<pts.size(); i++) {

			POINT p2(pts[i].x, pts[i].y);
			Index ind2(int(p2.x/dx), int(p2.y/dx));

			if (ind2 == ind1) { p1 = p2; }
			else {

				POINT dir = (p2-p1);

				float smin = 0.0;
				int maxiter = 8000, iter = 0;
				do {
					POINT face;
					face.x = (dir.x > 0) ? dx*(ind1.i+1) : dx*(ind1.i);
					face.y = (dir.y > 0) ? dx*(ind1.j+1) : dx*(ind1.j);

					POINT delta = face - p1;   // distance from point1 to faces x and y

					float s1, s2;
					s1 = (dir.x != 0.0) ? delta.x/dir.x : 1000.0;
					s2 = (dir.y != 0.0) ? delta.y/dir.y : 1000.0;

					float smin = min(s1, s2);

					if (s1 < s2) {
						ind1.i += sgn(dir.x);
						p1 += s1*dir;
					}
					else {
						ind1.j += sgn(dir.y);
						p1 += s2*dir;
					}
					dir = (p2-p1);

					POINT nvec = POINT(dir.y, dir.x)/mag(dir);
					iIB.push_back(Index(ind1.i+1, ind1.j+1, nvec));
//					nvec.print();

					iter ++;
				} while (!(ind1==ind2) && (iter < maxiter));
			}
		}


		deque<Index> iIBtmp;

		for (int i=0; i<iIB.size(); i++)
		{
			if ((iIB[i].i < 1) || (iIB[i].i > imax-1) || (iIB[i].j < 1) || (iIB[i].j > jmax-1)) {}
			else iIBtmp.push_back(iIB[i]);
		}
		iIB.clear();
		iIB = iIBtmp;

		sort(iIB.begin(), iIB.end(), lessThanIndex);
		iIB.erase( unique( iIB.begin(), iIB.end(), equalIndex), iIB.end() );
	}

	void tesselate() {

		float avgPointDist = 0.0;
		for (int i=1; i<pts.size(); i++)   avgPointDist += mag(pts[i] - pts[i-1]);
		avgPointDist /= pts.size();


		char triangParameters[200];
		sprintf(triangParameters, "pq30a%fFDY", avgPointDist);

		struct triangulateio in, out, vorout;

//		int Next = iIB.size();
		int Next = pts.size();
		in.numberofpoints = Next;

		// ****************************
		// Assign points [x y x y ....]
		in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));

		for (int i=0; i<2*Next; i=i+2)
		{
			// assign all the points
			in.pointlist[i]   = pts[i/2].x;
			in.pointlist[i+1] = pts[i/2].y;
		}

		// *****************************************************
		// Assign markers to points (1 means they are a boundary
		in.pointmarkerlist = (int *) malloc(in.numberofpoints * sizeof(int));
		for (int i=0; i<in.numberofpoints; i++)
			in.pointmarkerlist[i] = 1;

		// ******************************************************************
		// Assign points attributes (can specify area for example). note used
		// point attributes
		in.numberofpointattributes = 1;
		in.pointattributelist = (REAL *) malloc(in.numberofpoints * in.numberofpointattributes * sizeof(REAL));
		for (int i=0; i<in.numberofpoints*in.numberofpointattributes; i++)
			in.pointattributelist[i] = 0.0;

		// *****************************************
		// Define segments [start end start end ...]
		in.numberofsegments = Next;

		printf("number of segments %d\n", in.numberofsegments);
		in.segmentlist = (int *) malloc(in.numberofsegments*2*sizeof(int));

		int vert = 1;
		for (int i=0; i<Next*2; i=i+2)
		{
			in.segmentlist[i]   = vert;
			in.segmentlist[i+1] = vert+1;
			vert++;
		}
		in.segmentlist[2*Next-1] = 1; // connect last point to first

		in.numberofholes = 0;


		// ******************************************
		// Segment markers (1 means it is a boundary)
		in.segmentmarkerlist = (int *) malloc(in.numberofsegments * sizeof(int));
		for (int i=0; i<in.numberofsegments; i++)
			in.segmentmarkerlist[i] = 1;
		printf("Segments markers assigned\n");

		in.numberofregions = 1;
		in.regionlist = (REAL *) malloc(in.numberofregions * 4 * sizeof(REAL));


		// **********************************
		// Initialize out and vorout as NULL
		out.pointlist                 = (float *) NULL;
		out.pointattributelist        = (float *) NULL;
		out.pointmarkerlist           = (int *)   NULL;
		out.trianglelist              = (int *)   NULL;
		out.triangleattributelist     = (float *) NULL;
		out.neighborlist              = (int *)   NULL;
		out.segmentlist               = (int *)   NULL;
		out.segmentmarkerlist         = (int *)   NULL;
		out.edgelist                  = (int *)   NULL;
		out.edgemarkerlist            = (int *)   NULL;

		vorout.pointlist              = (float *) NULL;
		vorout.pointattributelist     = (float *) NULL;
		vorout.edgelist               = (int *)   NULL;
		vorout.normlist               = (float *) NULL;


		// *******************************************************
		// Call to external function to generate unstructured mesh
		triangulate(triangParameters,&in, &out, &vorout);

		for (int i=0; i<2*out.numberofpoints; i=i+2)
			node.push_back(POINT(out.pointlist[i], out.pointlist[i+1], 0.01));

		for (int i=0; i<3*out.numberoftriangles; i=i+3) {
			Triangle tri = {out.trianglelist[i]-1, out.trianglelist[i+1]-1, out.trianglelist[i+2]-1};
			elem.push_back(tri);
		}

//		for (int i=0; i<node.size(); i++)
//			cout << i << "\t" << node[i].x << "\t" << node[i].y << endl;
	}
};




#endif



