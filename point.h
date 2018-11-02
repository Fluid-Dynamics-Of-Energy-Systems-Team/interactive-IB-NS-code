/*
 * point.h
 *
 *  Created on: Jan 27, 2012
 *      Author: renep
 */



#ifndef POINT_H_
#define POINT_H_

#include "stdio.h"
#include "math.h"
#include <string>
#include <iostream>
#include <list>
#include <deque>
#include <algorithm>
using namespace std;


#define PI (4.*atan(1.0))

class POINT
{
public:
  float x, y, z, s, dummy;

public:
  POINT(): x(0.0), y(0.0), z(0.0), s(0.0), dummy(0) {}
  POINT(float value): x(value), y(value), z(value), s(0.0), dummy(0) {}
  POINT(float x_, float y_): x(x_), y(y_), z(0), s(0.0), dummy(0) {}
  POINT(float x_, float y_, float z_): x(x_), y(y_), z(z_), s(0.0), dummy(0) {}

  float radX() const;
  float radZ() const;
  float phiRadX() const;
  float phiRadZ() const;

  float phiDegX() const;
  float phiDegZ() const;

  void rotateDegX(float phi);

  void rotateDegZ(float phi);

  float dot(const POINT &p);

  void print(string text = "") { cout << text << ",\txyz: " << x << ",\t" << y << ",\t" << z << endl; }

  bool operator<(const POINT& p1) const{
    return x < p1.x;     // GUSTAVO 20.06.16 comparison
  }

};


#define VEC POINT

float mag(const POINT &pt);

POINT crossProd(const POINT &p1, const POINT &p2);


POINT operator+(const float val, const POINT &pt);

POINT operator+(const POINT &pt, const float val);

POINT operator+(const POINT &p1, const POINT &p2);

POINT operator-(const POINT &p1, const POINT &p2);

POINT operator-(const float val, const POINT &pt);

POINT operator-(const POINT &pt, const float val);

POINT operator+=(POINT &p1, const POINT &p2);

POINT operator-=(POINT &p1, const POINT &p2);


bool operator==(POINT &p1, const POINT &p2);

bool operator!=(POINT &p1, const POINT &p2);

POINT operator*(const float fact, const POINT &pt);

POINT operator*(const POINT &pt, const float fact);

POINT operator*(const POINT &pt1, const POINT &pt2);

POINT operator/(const POINT &pt1, const POINT &pt2);

POINT operator/(const POINT &pt1, const float val);

POINT operator/(const float val, const POINT &pt1);

POINT rotateDegZ(const float phi, const POINT &p1, const POINT &p2);

// read points from file
void readPts(deque<POINT> &pts, const string &name);



#endif /* POINT_H_ */



