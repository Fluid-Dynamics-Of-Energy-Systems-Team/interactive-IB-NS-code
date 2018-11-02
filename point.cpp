/*
 * point.cpp
 *
 *  Created on: Mar 19, 2012
 *      Author: renep
 */

#include "point.h"


float POINT::radX() const {return sqrt(y*y + z*z);}
float POINT::radZ() const {return sqrt(x*x + y*y);}
float POINT::phiRadX() const {return atan2(y, z);}
float POINT::phiRadZ() const {return atan2(x, y);}

float POINT::phiDegX() const {return 180.0/PI*atan2(y, z);}
float POINT::phiDegZ() const {return 180.0/PI*atan2(x, y);}

void POINT::rotateDegX(float phi)
{
	phi *= 4.*atan(1.0)/180.0;
	float tmp = y*cos(phi)-z*sin(phi);
	z = y*sin(phi)+z*cos(phi);
	y = tmp;
}

void POINT::rotateDegZ(float phi)
{
	phi *= 4.*atan(1.0)/180.0;
	float tmp = x*cos(phi)-y*sin(phi);
	y = x*sin(phi)+y*cos(phi);
	x = tmp;
}

float POINT::dot(const POINT &p) {return (x*p.x + y*p.y + z*p.z);}


float mag(const POINT &pt)
{
	return (sqrt(pt.x*pt.x + pt.y*pt.y + pt.z*pt.z));
}

POINT crossProd(const POINT &p1, const POINT &p2)
{
	POINT res;
	res.x = p1.y*p2.z - p1.z-p2.y;
	res.y = p1.z*p2.x - p1.x-p2.z;
	res.z = p1.x*p2.y - p1.y-p2.x;
	return res;
}

POINT operator+(const float val, const POINT &pt)
{
	POINT res;
	res.x = val + pt.x;
	res.y = val + pt.y;
	res.z = val + pt.z;
	return res;
}

POINT operator+(const POINT &pt, const float val)
{
	return (val+pt);
}


POINT operator+(const POINT &p1, const POINT &p2)
{
	POINT res;
	res.x = p1.x + p2.x;
	res.y = p1.y + p2.y;
	res.z = p1.z + p2.z;
	return res;
}

POINT operator-(const POINT &p1, const POINT &p2)
{
	POINT res;
	res.x = p1.x - p2.x;
	res.y = p1.y - p2.y;
	res.z = p1.z - p2.z;
	return res;
}

POINT operator-(const float val, const POINT &pt)
{
	POINT res;
	res.x = val - pt.x;
	res.y = val - pt.y;
	res.z = val - pt.z;
	return res;
}

POINT operator-(const POINT &pt, const float val)
{
	return (val-pt);
}



POINT operator+=(POINT &p1, const POINT &p2)
		{
	p1.x += p2.x;
	p1.y += p2.y;
	p1.z += p2.z;
	return p1;
		}

POINT operator-=(POINT &p1, const POINT &p2)
		{
	p1.x -= p2.x;
	p1.y -= p2.y;
	p1.z -= p2.z;
	return p1;
		}


bool operator==(POINT &p1, const POINT &p2)
		{
	float eps = 1.0e-8;
	return (fabs(p1.x - p2.x) < eps && fabs(p1.y - p2.y) < eps && fabs(p1.z - p2.z) < eps);
		}

bool operator!=(POINT &p1, const POINT &p2)
		{
	return !(p1 == p2);
		}

//bool operator>(const POINT &p1, const POINT &p2)
//{
//  return (p1.x > p2.x);
//}

POINT operator*(const float fact, const POINT &pt)
{
	POINT res;
	res.x = fact*pt.x;
	res.y = fact*pt.y;
	res.z = fact*pt.z;
	return res;
}

POINT operator*(const POINT &pt, const float fact)
{
	return (fact*pt);
}

POINT operator*(const POINT &pt1, const POINT &pt2)
{
	POINT res;
	res.x = pt1.x*pt2.x;
	res.y = pt1.y*pt2.y;
	res.z = pt1.z*pt2.z;
	return res;
}

POINT operator/(const POINT &pt1, const POINT &pt2)
{
	POINT res;

	if (pt2.x != 0.0) res.x = pt1.x/pt2.x;
	if (pt2.y != 0.0) res.y = pt1.y/pt2.y;
	if (pt2.z != 0.0) res.z = pt1.z/pt2.z;
	return res;
}

POINT operator/(const POINT &pt1, const float val)
{
	POINT res;

	if (val != 0.0)
	{
		res.x = pt1.x/val;
		res.y = pt1.y/val;
		res.z = pt1.z/val;
	}
	return res;
}

POINT operator/(const float val, const POINT &pt1)
{
	POINT res;

	if (pt1.x != 0.0) res.x = val/pt1.x;
	if (pt1.y != 0.0) res.y = val/pt1.y;
	if (pt1.z != 0.0) res.z = val/pt1.z;

	return res;
}

POINT rotateDegZ(const float phi, const POINT &p1, const POINT &p2)
{
	POINT p3 = p2-p1;
	p3.rotateDegZ(phi);
	return (p3+p1);
}


// read points from file
void readPts(deque<POINT> &pts, const string &name)
{
	FILE *fp = fopen(name.c_str(), "rt");
	if (fp == NULL)
	{
		printf("couldn't open %s\n", name.c_str());
		throw(-1);
	}

	POINT pt;
	while(fscanf(fp, "%lf%lf%lf", &pt.x, &pt.y, &pt.z) != EOF)  { pts.push_back(pt); }
	fclose(fp);
}

