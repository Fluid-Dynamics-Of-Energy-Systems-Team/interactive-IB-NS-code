

#ifndef MAIN_H_
#define MAIN_H_




void initFlowField();
float calcDT();
void calcDiv(float &divVolMax, float &divVolTot);
void setBCpressure(float *p);
void setBCvelocity(float *u, float *v);
void calcRhs(float *rhsX, float *rhsY, float *u, float *v);
void calcRhsCentral(float *rhsX, float *rhsY, float *u, float *v);
void calcRhsQUICK(float *rhsX, float *rhsY, float *u, float *v);

void initPressureSolver();
void calcPoissonRHS();
void solvePressureDCTPxNy();
void solvePressureDCTPxPy();
void solvePressureTDMAxDCTy();

void run();








#endif
