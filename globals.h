
#ifndef GLOBALS_H_
#define GLOBALS_H_


#include "ib.h"
#include <fftw3.h>



//const float H = 1.0;
const float LoH = 3.0;//6.25; //3.0;
const int jmax = 194;   // 38, 66, 98, 130, 146, 194, 218, [258], 290, 322, 386, 514, 770
const int imax = LoH*(jmax-2)+2;


const float dx = 1.0/(jmax-2);

int colorMap = 0;
int timeInt = 3;
int advecScheme = 2;
float contourMax =  30.0;
float contourMin = -30.0;
int wall = 3;
bool vectorToggle = false;
int ibcase = 1;
bool periodicX = false;
bool runSim = true;
float value = 1.0;
float uInlet = 1.0;
float vInlet = 0.0;

float umax = 0.0, vmax = 0.0;

int istep = 0;
double timing = 0.0;

float CFL = 0.9;
float dt;
float simtime = 0.0;
float Re = 10000;

float divVolMax, divVolTot;

bool change = false;

// allocate memory

deque<IBObject> ibObj;
deque<PostProcess> objPostPro;


float u[imax*jmax];
float v[imax*jmax];
float p[imax*jmax];

float us[imax*jmax];
float vs[imax*jmax];
float ps[imax*jmax];

float rhsX1[imax*jmax];
float rhsY1[imax*jmax];
float rhsX2[imax*jmax];
float rhsY2[imax*jmax];
float rhsX3[imax*jmax];
float rhsY3[imax*jmax];
float rhsX4[imax*jmax];
float rhsY4[imax*jmax];


float pp[(imax-2)*(jmax-2)];

float evNy[jmax-2];
fftwf_plan dctNy, idctNy;

float evPy[jmax-2];
fftwf_plan dctPy, idctPy;

float evPxNy[(imax-2)*(jmax-2)];
fftwf_plan dctPxNy, idctPxNy;

float evPxPy[(imax-2)*(jmax-2)];
fftwf_plan dctPxPy, idctPxPy;


float profile[jmax-2];
fftwf_plan dst, idst;





#endif
