// Incompressible unsteady Navier-Stokes solver on a staggered mesh
//
//       0         1         2         3         4         5         6   (imax = 7)
//           (0)       (1)       (2)       (3)       (4)       (5)       (6)
//
// 0     p    u    p    u    p    u    p    u    p    u    p    u    p    u
//
//  (0)  v    -----v---------v---------v---------v---------v-----    v
//            |         |         |         |         |         |
// 1     p    u    p    u ---p--- u    p    u    p    u    p    u    p    u
//            |         |\ \ \ \ \|         |         |         |
//  (1)  v    -----v----| \  V  \ |----v---------v---------v-----    v
//            |         |\ \ \ \ \|         |         |         |
// 2     p    u    p    u ---p--- u    p    u    p    u    p    u    p    u
//            |         |         |         |         |         |
//  (2)  v    -----v---------v---------v---------v---------v-----    v
//            |    |\ \ \ \ \|    | CUT-OUT / / / / / |         |
// 3     p    u    p \  U  \ p    u  / p /  u  / p /  u    p    u    p    u
//            |    |\ \ \ \ \|    | / / / / / / / / / |         |
//  (3)  v    -----v---------v-----/-/-v-/-/-/-/-v-/-/-----v33---    v
//            |         |         | / / / / / / / / / |         |
// 4     p    u    p    u    p    u  / p /  u  / p /  u    p    u    p    u
//            |         |         | / / / / / / / / / |         |
//  (4)  v    -----v---------v---------v---------v---------v-----    v
//
// 5     p    u    p    u    p    u    p    u    p    u    p    u    p    u
//
//  (5)  v         v         v         v         v         v         v
// (jmax = 6)




#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include <ctime>
#include <fftw3.h>
#include "omp.h"

#include <iostream>
#include <string>
#include <deque>
#include <algorithm>
using namespace std;



#include "point.h"
#include "ib.h"
#include "defines.h"
#include "display.h"
#include "globals.h"
#include "main.h"






int main(int argc, char** argv)
{

    initDisplay(argc, argv);

	initPressureSolver();

	initFlowField();

	IBObject obj = IBObject("objects/nozzle1.txt", dx, imax, jmax);
	obj.translate(POINT(-0.5, 0.0));
	ibObj.push_back(obj);

//	float umax, vmax;
	dt = calcDT();

	setBCvelocity(u, v);

//	float start = omp_get_wtime();
//	for (int i=0; i<2000; i++)
		timing = omp_get_wtime();

		run();

//		run2();
//	cout << "time = " << timing << "[s]" << endl;


	glutMainLoop();


	fftwf_cleanup_threads();
	fftwf_cleanup();

	return 0;
}



void initPressureSolver() {

	int m = imax-2;
	int n = jmax-2;

	int rank = 1;
	fftw_iodim dims[rank];
	int howmany_rank = 1;
	fftw_iodim howmany_dims[howmany_rank];
	fftw_r2r_kind kind[2];

	int ret = fftwf_init_threads();
	fftwf_plan_with_nthreads(NUM_THREADS_FFT);

	dims[0].n  = n;
	dims[0].is = 1;
	dims[0].os = 1;

	howmany_dims[0].n  = m;
	howmany_dims[0].is = n;
	howmany_dims[0].os = n;

	// init 1d dct in y direction Neumann
	kind[0] = FFTW_REDFT10;
	kind[1] = FFTW_REDFT01;
	 dctNy = fftwf_plan_guru_r2r(rank, dims, howmany_rank, howmany_dims, pp, pp, &kind[0], FFTW_ESTIMATE);
	idctNy = fftwf_plan_guru_r2r(rank, dims, howmany_rank, howmany_dims, pp, pp, &kind[1], FFTW_ESTIMATE);
	for(int j=0; j<n; j++)  evNy[j] = -4.0*pow(sin(M_PI/2.0*j/n)/dx, 2.0);


	// init 1d dct in y direction Periodic
	kind[0] = FFTW_R2HC;
	kind[1] = FFTW_HC2R;
	 dctPy = fftwf_plan_guru_r2r(rank, dims, howmany_rank, howmany_dims, pp, pp, &kind[0], FFTW_ESTIMATE);
	idctPy = fftwf_plan_guru_r2r(rank, dims, howmany_rank, howmany_dims, pp, pp, &kind[1], FFTW_ESTIMATE);
    for(int j=0; j<n; j++)  evPy[j] = -4.0*pow(sin(M_PI*j/n)/dx, 2.0);


	// init 2d dct Periodic in x, Neumann in y
	kind[0] = FFTW_R2HC;
	kind[1] = FFTW_REDFT10;
	dctPxNy   = fftwf_plan_r2r_2d(m, n, pp, pp, kind[0], kind[1], FFTW_MEASURE);

	kind[0] = FFTW_HC2R;
	kind[1] = FFTW_REDFT01;
	idctPxNy  = fftwf_plan_r2r_2d(m, n, pp, pp, kind[0], kind[1], FFTW_MEASURE);

	float evPx[m];
	for(int i=0; i<m; i++)  evPx[i] = -4.0*pow(sin(M_PI*i/m)*LoH/dx, 2.0);

	for(int i=0; i<m; i++)
		for(int j=0; j<n; j++)
			evPxNy[i*n+j] = evPx[i] + evNy[j];
	evPxNy[0] = 1.0;


	// init 2d dct Periodic in x, Periodic in y
	kind[0] = FFTW_R2HC;
	kind[1] = FFTW_R2HC;
	 dctPxPy = fftwf_plan_r2r_2d(m, n, pp, pp, kind[0], kind[1], FFTW_MEASURE);

	kind[0] = FFTW_HC2R;
	kind[1] = FFTW_HC2R;
	idctPxPy = fftwf_plan_r2r_2d(m, n, pp, pp, kind[0], kind[1], FFTW_MEASURE);

	for(int i=0; i<m; i++)
		for(int j=0; j<n; j++)
			evPxPy[i*n+j] = evPx[i] + evPy[j];
	evPxPy[0] = 1.0;



	 dst = fftwf_plan_r2r_1d(jmax-2, profile, profile, FFTW_RODFT10, FFTW_MEASURE);
	idst = fftwf_plan_r2r_1d(jmax-2, profile, profile, FFTW_RODFT01, FFTW_MEASURE);
}



void initFlowField() {
	for (int i=0; i<imax*jmax; i++) {
//		float y = dx*(j-0.5);
//		 u[i] = 4*y*(1-y); //uInlet;
	 	 u[i] = 1.0;
 		 v[i] = 0.0;
		 p[i] = 0.0;
		us[i] = 0.0;
		vs[i] = 0.0;
		ps[i] = 0.0;
	}
}



float calcDT()
{
	umax = 0.0, vmax = 0;

//	for (int i=0; i<imax*jmax; i++) u_max = max(fabs(u[i]), u_max);
//	for (int i=0; i<imax*jmax; i++) v_max = max(fabs(v[i]), v_max);

	for (int i=0; i<imax-1; i++)
	for (int j=1; j<jmax-1; j++)
		umax = MAX(fabs(u[ind(i,j)]), umax);

	for (int i=1; i<imax-1; i++)
	for (int j=0; j<jmax-1; j++)
		vmax = MAX(fabs(v[ind(i,j)]), vmax);

    return CFL/( 4.0/(dx*dx*Re) + (umax+vmax)/dx);
}



void calcDiv(float &divVolMax, float &divVolTot)
{
    divVolMax = -1.0e10;
    divVolTot =  0.0;

	for (int i=1; i<imax-1; i++)
	for (int j=1; j<jmax-1; j++) {
		float divVol = (u[ind(i,j)]-u[ind(i-1,j)] + v[ind(i,j)]-v[ind(i,j-1)])*dx;
		divVolTot = divVolTot + divVol;
		divVolMax = MAX(divVolMax, divVol);
	}
}



void setBCpressure(float *p)
{
	// BC TOP and BOTTOM
	switch (wall) {
	case 1:
	case 2:
		for (int i=1; i<imax-1; i++) {
			p[ind(i,0)]      = p[ind(i,1)];
			p[ind(i,jmax-1)] = p[ind(i,jmax-2)];
		}
		break;
	case 3:
		for (int i=1; i<imax-1; i++) {
			p[ind(i,0)]      = p[ind(i,jmax-2)];
			p[ind(i,jmax-1)] = p[ind(i,1)];
		}
		break;
	}

	// BC IN and OUTFLOW
	if (periodicX == false) {
		for (int j=0; j<jmax; j++) {
			p[ind(0,j)]      = p[ind(1,j)];
			p[ind(imax-1,j)] = p[ind(imax-2,j)];
		}
	}
	else {
		for (int j=0; j<jmax; j++) {
			p[ind(0,j)]      = p[ind(imax-2,j)];
			p[ind(imax-1,j)] = p[ind(1,j)];
		}
	}
}



void setBCvelocity(float *u, float *v)
{
    // BC TOP and BOTTOM
	switch (wall) {
	case 1:
		for (int i=0; i<imax; i++) {
			u[ind(i,0)]      = u[ind(i,1)];
			u[ind(i,jmax-1)] = u[ind(i,jmax-2)];

			v[ind(i,0)]      = 0.0;
			v[ind(i,jmax-2)] = 0.0;
			v[ind(i,jmax-1)] = 0.0;
		}
		break;
	case 2:
		for (int i=0; i<imax; i++) {
			u[ind(i,0)]      = -u[ind(i,1)];
			u[ind(i,jmax-1)] = -u[ind(i,jmax-2)];

			v[ind(i,0)]      = 0.0;
			v[ind(i,jmax-2)] = 0.0;
			v[ind(i,jmax-1)] = 0.0;
		}
		break;
	case 3:
		for (int i=0; i<imax; i++) {
			u[ind(i,0)]      = u[ind(i,jmax-2)];
			u[ind(i,jmax-1)] = u[ind(i,1)];

			v[ind(i,0)]      = v[ind(i,jmax-2)];
			v[ind(i,jmax-1)] = v[ind(i,1)];
		}
		break;
	}


	// BC INFLOW AND OUTFLOW
	if (periodicX == false) {

//        bool flagBackFlow = false;

		float umean = 0.0;
		for (int j=1; j<jmax-1; j++)  umean += 0.5*(u[ind(imax-2,j)]+u[ind(imax-3,j)])*dx;


		for (int j=1; j<jmax-1; j++) {
			float y = dx*(j-0.5);
			u[ind(0,j)] = 1.0; //4*y*(1-y); //uInlet;
//			u[ind(imax-2,j)] = u[ind(imax-3,j)] - (v[ind(imax-2,j)]-v[ind(imax-2,j-1)]);

			u[ind(imax-2,j)] = u[ind(imax-2,j)] - dt/dx*umean*(u[ind(imax-2,j)] - u[ind(imax-3,j)]); // - (v[ind(imax-2,j)]-v[ind(imax-2,j-1)]);
//			u[ind(imax-2,j)] = u[ind(imax-2,j)] - dt/dx*u[ind(imax-3,j)]*(u[ind(imax-2,j)] - u[ind(imax-3,j)]); // - (v[ind(imax-2,j)]-v[ind(imax-2,j-1)]);

//			if (u[ind(imax-2,j)] < 0.0) flagBackFlow = true;

//			profile[j-1] = u[ind(imax-2,j)];
//			u[ind(imax-1,j)] = u[ind(imax-2,j)];

			v[ind(0,j)] = vInlet; //-v[ind(1,j)];
			v[ind(imax-1,j)] = v[ind(imax-2,j)];
		}

//		if (flagBackFlow == true) {
//			fftwf_execute(dst);
//			for (int j=jmax-10; j<jmax-2; j++)  profile[j] = 0.0;
//			fftwf_execute(idst);
//			for (int j=1; j<jmax-1; j++)   u[ind(imax-2,j)] = profile[j-1]/(2*(jmax-2));
//		}

//		if (flagBackFlow == true) {
//
//		}
		
//		if (flagBackFlow == true) {
//  		  double weight = 0.7;
//		  for (int j=1; j<jmax-1; j++)
//			  u[ind(imax-2,j)] = weight*u[ind(imax-2,j)] + (1-weight)/2.0*(u[ind(imax-2,j+1)] + u[ind(imax-2,j-1)]);
//		}

		float MFin = 0.0;
		float MFout = 0.0;
		float dMF;
		for (int j=1; j<jmax-1; j++) {
		   MFin  = MFin  + u[j]*dx;
		   MFout = MFout + u[(imax-2)*jmax+j]*dx;
		}

		 dMF = MFin - MFout;
		 //dMF = MFin/MFout;

//	//	if (fabs(MFout) < 1.0e-15) 		MFout = 1e20;

		for (int j=1; j<jmax-1; j++) {
			u[(imax-2)*jmax+j] += dMF;   //u[(imax-2)*jmax+j] + u[(imax-2)*jmax+j] * dx * dMF/MFout;
			u[(imax-1)*jmax+j] = u[(imax-2)*jmax+j];
		}

//		MFin = 0.0;
//		MFout = 0.0;
//		for (int j=1; j<jmax-1; j++) {
//		   MFin  = MFin  + u[j]*dx;
//		   MFout = MFout + u[ind(imax-2,j)]*dx;
////		   printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\n", j, u[ind(0,j)], u[ind(imax-3,j)], u[ind(imax-2,j)], v[ind(imax-2,j)], v[ind(imax-1,j)]);
//		}

//		dMF = MFin - MFout;
//		printf("%lf\t%lf\t%lf\n", MFin, MFout, dMF);

	}
	else {
		for (int j=1; j<jmax-1; j++) {
			u[ind(0,j)] = u[ind(imax-2,j)];
			u[ind(imax-1,j)] = u[ind(1,j)];

			v[ind(0,j)] = v[ind(imax-2,j)];
			v[ind(imax-1,j)] = v[ind(1,j)];
		}
	}


    // set immersed boundaries
	for (int obj=0; obj<ibObj.size(); obj++) {
		for (int ib=0; ib<ibObj[obj].iIB.size(); ib++) {
			int i = ibObj[obj].iIB[ib].i;
			int j = ibObj[obj].iIB[ib].j;

			u[ind(i-1,j)] = 0.0;
			u[ind(i,  j)] = 0.0;
			v[ind(i,j-1)] = 0.0;
			v[ind(i,j  )] = 0.0;
		}
	}
}



void run()
{
	if (runSim == false) return;

	// advance velocity
	switch (timeInt) {

	// EULER ----------------------------------------------------------------------------------------
	case 1:
		calcRhs(rhsX1, rhsY1, u, v);
#pragma omp parallel for num_threads(NUM_THREADS)
		for (int i=jmax+1; i<(imax-1)*jmax-1; i++) us[i] = u[i] + dt*rhsX1[i];
#pragma omp parallel for num_threads(NUM_THREADS)
		for (int i=jmax+1; i<(imax-1)*jmax-1; i++) vs[i] = v[i] + dt*rhsY1[i];
		setBCvelocity(us, vs);
		break;

	// AB2 ----------------------------------------------------------------------------------------
	case 2:
		calcRhs(rhsX1, rhsY1, u, v);

		if (istep == 0) {
#pragma omp parallel for num_threads(NUM_THREADS)
			for (int i=jmax+1; i<(imax-1)*jmax-1; i++) us[i] = u[i] + dt*rhsX1[i];
#pragma omp parallel for num_threads(NUM_THREADS)
			for (int i=jmax+1; i<(imax-1)*jmax-1; i++) vs[i] = v[i] + dt*rhsY1[i];
			setBCvelocity(us, vs);
		}
		else {
			if (change == false) {
#pragma omp parallel for num_threads(NUM_THREADS)
			for (int i=jmax+1; i<(imax-1)*jmax-1; i++) us[i] = u[i] + dt*(1.5*rhsX1[i] - 0.5*rhsX2[i]);
#pragma omp parallel for num_threads(NUM_THREADS)
			for (int i=jmax+1; i<(imax-1)*jmax-1; i++) vs[i] = v[i] + dt*(1.5*rhsY1[i] - 0.5*rhsY2[i]);
			setBCvelocity(us, vs);
			}
			else {
#pragma omp parallel for num_threads(NUM_THREADS)
			for (int i=jmax+1; i<(imax-1)*jmax-1; i++) us[i] = u[i] + dt*(1.51*rhsX1[i] - 0.51*rhsX2[i]);
#pragma omp parallel for num_threads(NUM_THREADS)
			for (int i=jmax+1; i<(imax-1)*jmax-1; i++) vs[i] = v[i] + dt*(1.51*rhsY1[i] - 0.51*rhsY2[i]);
			setBCvelocity(us, vs);
			}

		}

#pragma omp parallel for num_threads(NUM_THREADS)
		for (int i=jmax+1; i<(imax-1)*jmax-1; i++) rhsX2[i] = rhsX1[i];
#pragma omp parallel for num_threads(NUM_THREADS)
		for (int i=jmax+1; i<(imax-1)*jmax-1; i++) rhsY2[i] = rhsY1[i];
		break;

	// RK3 ----------------------------------------------------------------------------------------
	case 3:
		calcRhs(rhsX1,  rhsY1,  u,  v);
#pragma omp parallel for num_threads(NUM_THREADS)
		for (int i=jmax+1; i<(imax-1)*jmax-1; i++) us[i] = u[i] + 0.5*dt*rhsX1[i];
#pragma omp parallel for num_threads(NUM_THREADS)
		for (int i=jmax+1; i<(imax-1)*jmax-1; i++) vs[i] = v[i] + 0.5*dt*rhsY1[i];
		setBCvelocity(us, vs);

		calcRhs(rhsX2, rhsY2, us, vs);
#pragma omp parallel for num_threads(NUM_THREADS)
		for (int i=jmax+1; i<(imax-1)*jmax-1; i++) us[i] = u[i] + dt*(-rhsX1[i] + 2.0*rhsX2[i]);
#pragma omp parallel for num_threads(NUM_THREADS)
		for (int i=jmax+1; i<(imax-1)*jmax-1; i++) vs[i] = v[i] + dt*(-rhsY1[i] + 2.0*rhsY2[i]);
		setBCvelocity(us, vs);

		calcRhs(rhsX3, rhsY3, us, vs);
#pragma omp parallel for num_threads(NUM_THREADS)
		for (int i=jmax+1; i<(imax-1)*jmax-1; i++) us[i] = u[i] + dt*(rhsX1[i] + 4.0*rhsX2[i] + rhsX3[i])/6.0;
#pragma omp parallel for num_threads(NUM_THREADS)
		for (int i=jmax+1; i<(imax-1)*jmax-1; i++) vs[i] = v[i] + dt*(rhsY1[i] + 4.0*rhsY2[i] + rhsY3[i])/6.0;
		setBCvelocity(us, vs);
		break;

	// RK4 ----------------------------------------------------------------------------------------
	case 4:
		calcRhs(rhsX1,  rhsY1,  u,  v);
#pragma omp parallel for num_threads(NUM_THREADS)
		for (int i=jmax+1; i<(imax-1)*jmax-1; i++) us[i] = u[i] + 0.5*dt*rhsX1[i];
#pragma omp parallel for num_threads(NUM_THREADS)
		for (int i=jmax+1; i<(imax-1)*jmax-1; i++) vs[i] = v[i] + 0.5*dt*rhsY1[i];
		setBCvelocity(us, vs);

		calcRhs(rhsX2, rhsY2, us, vs);
#pragma omp parallel for num_threads(NUM_THREADS)
		for (int i=jmax+1; i<(imax-1)*jmax-1; i++) us[i] = u[i] + 0.5*dt*rhsX2[i];
#pragma omp parallel for num_threads(NUM_THREADS)
		for (int i=jmax+1; i<(imax-1)*jmax-1; i++) vs[i] = v[i] + 0.5*dt*rhsY2[i];
		setBCvelocity(us, vs);

		calcRhs(rhsX3, rhsY3, us, vs);
#pragma omp parallel for num_threads(NUM_THREADS)
		for (int i=jmax+1; i<(imax-1)*jmax-1; i++) us[i] = u[i] + dt*rhsX3[i];
#pragma omp parallel for num_threads(NUM_THREADS)
		for (int i=jmax+1; i<(imax-1)*jmax-1; i++) vs[i] = v[i] + dt*rhsY3[i];
		setBCvelocity(us, vs);

		calcRhs(rhsX4, rhsY4, us, vs);
#pragma omp parallel for num_threads(NUM_THREADS)
		for (int i=jmax+1; i<(imax-1)*jmax-1; i++) us[i] = u[i] + dt*(rhsX1[i]+2*rhsX2[i]+2*rhsX3[i]+rhsX4[i])/6;
#pragma omp parallel for num_threads(NUM_THREADS)
		for (int i=jmax+1; i<(imax-1)*jmax-1; i++) vs[i] = v[i] + dt*(rhsY1[i]+2*rhsY2[i]+2*rhsY3[i]+rhsY4[i])/6;
		setBCvelocity(us, vs);
		break;
	}


	calcPoissonRHS();

	// Poisson equation for pressure using direct solver
	if (periodicX == false) {
		solvePressureTDMAxDCTy();
	}
	else {
		switch (wall) {
		case 1:
		case 2:
			solvePressureDCTPxNy();
			break;
		case 3:
			solvePressureDCTPxPy();
			break;
		}
	}

	// correct velocity
#pragma omp parallel for num_threads(NUM_THREADS)
	for (int i=1; i<imax-1; i++)
		for (int j=1; j<jmax-1; j++)
			u[ind(i,j)] = us[ind(i,j)] - dt/dx*(ps[ind(i+1,j)] - ps[ind(i,j)]);

#pragma omp parallel for num_threads(NUM_THREADS)
	for (int i=1; i<imax-1; i++)
		for (int j=1; j<jmax-1; j++)
			v[ind(i,j)] = vs[ind(i,j)] - dt/dx*(ps[ind(i,j+1)] - ps[ind(i,j)]);

	// update pressure and set boundary condition
#pragma omp parallel for num_threads(NUM_THREADS)
	for (int i=0; i<imax*jmax; i++)
		p[i] += ps[i];

	setBCpressure(p);

#pragma omp for simd
	for (int i=0; i<imax*jmax; i++)
		p[i] -= p[imax*jmax-1];

	setBCvelocity(u, v);

	if ((istep%15) == 0) display();

//	float umax = 0.0, vmax = 0.0;
	if ((istep%5) == 0)  dt = calcDT();

	if ((istep%1000) == 0) {
		cout << "time = " << (omp_get_wtime() - timing) << "[s]" << endl;
		timing = omp_get_wtime();

		calcDiv(divVolMax, divVolTot);
		printf("%d, %le, %le, %le, %le, %le\n", istep, dt, divVolMax, divVolTot, umax, vmax);
	}

	for (int i=0; i<objPostPro.size(); i++) {
		objPostPro[i].count++;
		for (int p=0; p<objPostPro[i].iPP.size(); p++) {
			Index ix = objPostPro[i].iPP[p];
			objPostPro[i].velAverage[p].x += (u[ix.i*jmax+ix.j] - objPostPro[i].velAverage[p].x)/float(objPostPro[i].count);
			objPostPro[i].velAverage[p].y += (v[ix.i*jmax+ix.j] - objPostPro[i].velAverage[p].y)/float(objPostPro[i].count);
		}
	}


//	if ((istep%50000) == 0) {
//		FILE *fp = fopen("out.txt", "wt");
//		for (int i=1; i<imax-2; i++) {
//
//			float mass = 0.0, mom = 0.0, pre = 0.0;
//			for (int j=1; j<jmax-1; j++) {
//				mass += 0.5*(u[ind(i-1,j)] + u[ind(i,j)])*dx;
//				mom += pow(0.5*(u[ind(i-1,j)] + u[ind(i,j)]), 2.0)*dx;
//				pre += p[ind(i,j)]*dx;
//			}
//
//			fprintf(fp, "%f\t%f\t", dx*(i-1), u[ind(i-1,jmax/2-1)]);
//			fprintf(fp, "%f\t%f\t%f\t", dx*(i-0.5), (u[ind(i-1,jmax/2-1)] + u[ind(i,jmax/2-1)])/2.0, p[ind(i,jmax/2-1)]);
//			fprintf(fp, "%f\t%f\t%f\n", mass, mom, pre);
//		}
//		fclose(fp);
////
////		fp = fopen("integral.txt", "wt");
////		fclose(fp);
//	}

	istep = istep + 1;
	simtime += dt;


//	float freq = (float)istep*0.005*2.0*M_PI;
//
//	for (int ib = 0; ib<nIB; ib++)
//	{
//		float xBody = 40.0*pow(cos(freq),333.0);
//		float yBody = 20.0*pow(sin(freq),333.0);
//		indexIB[ib][0] = indexIBOrig[ib][0] + xBody+60;
//		indexIB[ib][1] = indexIBOrig[ib][1] + yBody;
//	}
//
//	cout << "time = " << (clock() - start)/(float)(CLOCKS_PER_SEC) << "[s]" << endl;
}



void calcRhs(float *rhsX, float *rhsY, float *u, float *v) {

	switch (advecScheme) {
	case 1:
		calcRhsCentral(rhsX, rhsY, u, v);
		break;
	case 2:
		if (istep%40 == 0)
			calcRhsQUICK(rhsX, rhsY, u, v);
		else
			calcRhsCentral(rhsX, rhsY, u, v);
		break;
	case 3:
		if (istep%20 == 0)
			calcRhsQUICK(rhsX, rhsY, u, v);
		else
			calcRhsCentral(rhsX, rhsY, u, v);
		break;
	case 4:
		if (istep%5 == 0)
			calcRhsQUICK(rhsX, rhsY, u, v);
		else
			calcRhsCentral(rhsX, rhsY, u, v);
		break;
	case 5:
		calcRhsQUICK(rhsX, rhsY, u, v);
		break;

	}

	if (periodicX == true)
		for (int i=0; i<imax*jmax; i++)
			rhsX[i] += 0.7;
}



void calcRhsCentral(float *rhsX, float *rhsY, float *u, float *v) {

    float dxi = 1.0/dx;
    float dx2Rei = dxi*dxi/Re;
    float uuip05[imax*jmax], vvjp05[imax*jmax], uvp05[imax*jmax];

    // calculate nonlinear terms of velocity
#pragma omp parallel for num_threads(NUM_THREADS)
	for (int i=0; i<(imax-1)*jmax; i++)
		uuip05[i] = 0.25*powf(u[i+jmax] + u[i], 2.0);

#pragma omp parallel for num_threads(NUM_THREADS)
	for (int i=0; i<(imax-1)*jmax; i++)
		vvjp05[i] = 0.25*powf(v[i+1] + v[i], 2.0);

#pragma omp parallel for num_threads(NUM_THREADS)
	for (int i=0; i<(imax-1)*jmax; i++)
		uvp05[i] = 0.25*(u[i+1]+u[i])*(v[i+jmax]+v[i]);


    // x - momentum
#pragma omp parallel for num_threads(NUM_THREADS)
    for (int i=jmax+1; i<(imax-1)*jmax-1; i++)
		rhsX[i] = -dxi*(uuip05[i]-uuip05[i-jmax] + uvp05[i]-uvp05[i-1] + p[i+jmax]-p[i]) + dx2Rei*(u[i+jmax]+u[i-jmax]+u[i+1]+u[i-1]-4.0*u[i]);

    // y - momentum
#pragma omp parallel for num_threads(NUM_THREADS)
    for (int i=jmax+1; i<(imax-1)*jmax-1; i++)
        rhsY[i] = -dxi*(uvp05[i]-uvp05[i-jmax] + vvjp05[i]-vvjp05[i-1] + p[i+1]-p[i]) + dx2Rei*(v[i+jmax]+v[i-jmax]+v[i+1]+v[i-1]-4.0*v[i]);

    // IB forcing
    for (int obj=0; obj<ibObj.size(); obj++) {
		for (int ib=0; ib<ibObj[obj].iIB.size(); ib++) {
			int i = ibObj[obj].iIB[ib].i-1;
			int j = ibObj[obj].iIB[ib].j-1;     rhsX[ind(i,j)] += (dxi*uvp05[ind(i,j)]   - dx2Rei*(u[ind(i,j+1)]+u[ind(i,j)]));
				i = ibObj[obj].iIB[ib].i;       rhsX[ind(i,j)] += (dxi*uvp05[ind(i,j)]   - dx2Rei*(u[ind(i,j+1)]+u[ind(i,j)]));
				j = ibObj[obj].iIB[ib].j+1;     rhsX[ind(i,j)] -= (dxi*uvp05[ind(i,j-1)] + dx2Rei*(u[ind(i,j)]+u[ind(i,j-1)]));
				i = ibObj[obj].iIB[ib].i-1;     rhsX[ind(i,j)] -= (dxi*uvp05[ind(i,j-1)] + dx2Rei*(u[ind(i,j)]+u[ind(i,j-1)]));

				j = ibObj[obj].iIB[ib].j-1;     rhsY[ind(i,j)] += (dxi*uvp05[ind(i,j)]   - dx2Rei*(v[ind(i+1,j)]+v[ind(i,j)]));
				j = ibObj[obj].iIB[ib].j;       rhsY[ind(i,j)] += (dxi*uvp05[ind(i,j)]   - dx2Rei*(v[ind(i+1,j)]+v[ind(i,j)]));
				i = ibObj[obj].iIB[ib].i+1;     rhsY[ind(i,j)] -= (dxi*uvp05[ind(i-1,j)] + dx2Rei*(v[ind(i,j)]+v[ind(i-1,j)]));
				j = ibObj[obj].iIB[ib].j-1;     rhsY[ind(i,j)] -= (dxi*uvp05[ind(i-1,j)] + dx2Rei*(v[ind(i,j)]+v[ind(i-1,j)]));
		}
    }
}



void calcRhsQUICK(float *rhsX, float *rhsY, float *u, float *v) {

    float dxi = 1.0/dx;
    float dx2Rei = dxi*dxi/Re;
    float vAtU[imax*jmax], uAtV[imax*jmax];
    float phi = 1.0/8.0;
//    float phi = 0.0;

#pragma omp parallel for num_threads(NUM_THREADS)
    for (int i=jmax+1; i<(imax-1)*jmax-1; i++)
		vAtU[i] = 0.25*(v[i] + v[i-1] + v[i+jmax] + v[i+jmax-1]);

#pragma omp parallel for num_threads(NUM_THREADS)
    for (int i=jmax+1; i<(imax-1)*jmax-1; i++)
		uAtV[i] = 0.25*(u[i] + u[i-jmax] + u[i+1] + u[i-jmax+1]);



    // x - momentum
#pragma omp parallel for num_threads(NUM_THREADS)
	for (int i=1; i<imax-1; i++)
	for (int j=1; j<jmax-1; j++)
	{
		rhsX[ind(i,j)] = -0.5*phi*( u[ind(i,j)]*(u[ind(i+1,j)]-u[ind(i-1,j)]) + vAtU[ind(i,j)]*(u[ind(i,j+1)]-u[ind(i,j-1)]) ) - (p[ind(i+1,j)] - p[ind(i,j)]);

		if ((u[ind(i,j)] > 0.0) && (i>1))          rhsX[ind(i,j)] -= u[ind(i,j)]*(1.0-phi)*(1.5*u[ind(i,j)]-2.0*u[ind(i-1,j)]+0.5*u[ind(i-2,j)]);
		if ((u[ind(i,j)] < 0.0) && (i<imax-2))     rhsX[ind(i,j)] += u[ind(i,j)]*(1.0-phi)*(1.5*u[ind(i,j)]-2.0*u[ind(i+1,j)]+0.5*u[ind(i+2,j)]);

		if ((vAtU[ind(i,j)] > 0.0) && (j>1))       rhsX[ind(i,j)] -= vAtU[ind(i,j)]*(1.0-phi)*(1.5*u[ind(i,j)]-2.0*u[ind(i,j-1)]+0.5*u[ind(i,j-2)]);
		if ((vAtU[ind(i,j)] < 0.0) && (j<jmax-2))  rhsX[ind(i,j)] += vAtU[ind(i,j)]*(1.0-phi)*(1.5*u[ind(i,j)]-2.0*u[ind(i,j+1)]+0.5*u[ind(i,j+2)]);

		rhsX[ind(i,j)] = rhsX[ind(i,j)]*dxi + dx2Rei*(u[ind(i+1,j)] + u[ind(i-1,j)] + u[ind(i,j+1)] + u[ind(i,j-1)] - 4*u[ind(i,j)]);
	}


    // y - momentum
#pragma omp parallel for num_threads(NUM_THREADS)
	for (int i=1; i<imax-1; i++)
	for (int j=1; j<jmax-1; j++)
	{
		rhsY[ind(i,j)] = -0.5*phi*( uAtV[ind(i,j)]*(v[ind(i+1,j)]-v[ind(i-1,j)]) + v[ind(i,j)]*(v[ind(i,j+1)]-v[ind(i,j-1)]) ) - (p[ind(i,j+1)] - p[ind(i,j)]);

		if ((uAtV[ind(i,j)] > 0.0) && (i>1))       rhsY[ind(i,j)] -= uAtV[ind(i,j)]*(1.0-phi)*(1.5*v[ind(i,j)]-2.0*v[ind(i-1,j)]+0.5*v[ind(i-1,j)]);
		if ((uAtV[ind(i,j)] < 0.0) && (i<imax-2))  rhsY[ind(i,j)] += uAtV[ind(i,j)]*(1.0-phi)*(1.5*v[ind(i,j)]-2.0*v[ind(i+1,j)]+0.5*v[ind(i+2,j)]);

		if ((v[ind(i,j)] > 0.0) && (j>1))          rhsY[ind(i,j)] -= v[ind(i,j)]*(1.0-phi)*(1.5*v[ind(i,j)]-2.0*v[ind(i,j-1)]+0.5*v[ind(i,j-2)]);
		if ((v[ind(i,j)] < 0.0) && (j<jmax-2))     rhsY[ind(i,j)] += v[ind(i,j)]*(1.0-phi)*(1.5*v[ind(i,j)]-2.0*v[ind(i,j+1)]+0.5*v[ind(i,j+2)]);

		rhsY[ind(i,j)] = rhsY[ind(i,j)]*dxi + dx2Rei*(v[ind(i+1,j)] + v[ind(i-1,j)] + v[ind(i,j+1)] + v[ind(i,j-1)] - 4*v[ind(i,j)]);
	}




//    // x - momentum
////	for (int k=1; k<imax-1; k++)
////	for (int l=1; l<jmax-1; l++)
//#pragma omp parallel for num_threads(NUM_THREADS)
//    for (int i=jmax+1; i<(imax-1)*jmax-1; i++)
//	{
//		if (u[i] > 0.0) 		rhsX[i] = -u[i]*(0.5*phi*(u[ip1[i]]-u[im1[i]]) + (1.0-phi)*(1.5*u[i]-2.0*u[im1[i]]+0.5*u[im2[i]]))*dxi;
//		else					rhsX[i] = +u[i]*(0.5*phi*(u[ip1[i]]-u[im1[i]]) + (1.0-phi)*(1.5*u[i]-2.0*u[ip1[i]]+0.5*u[ip2[i]]))*dxi;
//
//		if (vAtU[i] > 0.0) 		rhsX[i] -= vAtU[i]*(0.5*phi*(u[jp1[i]]-u[jm1[i]]) + (1.0-phi)*(1.5*u[i]-2.0*u[jm1[i]]+0.5*u[jm2[i]]))*dxi;
//		else					rhsX[i] += vAtU[i]*(0.5*phi*(u[jp1[i]]-u[jm1[i]]) + (1.0-phi)*(1.5*u[i]-2.0*u[jp1[i]]+0.5*u[jp2[i]]))*dxi;
//
//		rhsX[i] += dx2Rei*(u[ip1[i]] + u[im1[i]] + u[jp1[i]] + u[jm1[i]] - 4*u[i]);
//		rhsX[i] -= (p[ip1[i]] - p[i])*dxi;
//	}
//    // y - momentum
////	for (int k=1; k<imax-1; k++)
////	for (int l=1; l<jmax-1; l++)
//#pragma omp parallel for num_threads(NUM_THREADS)
//    for (int i=jmax+1; i<(imax-1)*jmax-1; i++)
//	{
//		if (uAtV[i] > 0.0) 		rhsY[i] = -uAtV[i]*(0.5*phi*(v[ip1[i]]-v[im1[i]])+(1.0-phi)*(1.5*v[i]-2.0*v[im1[i]]+0.5*v[im2[i]]))*dxi;
//		else					rhsY[i] = +uAtV[i]*(0.5*phi*(v[ip1[i]]-v[im1[i]])+(1.0-phi)*(1.5*v[i]-2.0*v[ip1[i]]+0.5*v[ip2[i]]))*dxi;
//
//		if (v[i] > 0.0) 		rhsY[i] -= v[i]*(0.5*phi*(v[jp1[i]]-v[jm1[i]])+(1.0-phi)*(1.5*v[i]-2.0*v[jm1[i]]+0.5*v[jm2[i]]))*dxi;
//		else					rhsY[i] += v[i]*(0.5*phi*(v[jp1[i]]-v[jm1[i]])+(1.0-phi)*(1.5*v[i]-2.0*v[jp1[i]]+0.5*v[jp2[i]]))*dxi;
//
//		rhsY[i] += dx2Rei*(v[ip1[i]] + v[im1[i]] + v[jp1[i]] + v[jm1[i]] - 4*v[i]);
//		rhsY[i] -= (p[jp1[i]] - p[i])*dxi;
//	}

    // IB forcing
    for (int obj=0; obj<ibObj.size(); obj++) {
		for (int ib=0; ib<ibObj[obj].iIB.size(); ib++) {
			int i = ibObj[obj].iIB[ib].i-1;
			int j = ibObj[obj].iIB[ib].j-1;  	  rhsX[ind(i,j)] += (- dx2Rei*(u[ind(i,j+1)]+u[ind(i,j)]));
				i = ibObj[obj].iIB[ib].i;         rhsX[ind(i,j)] += (- dx2Rei*(u[ind(i,j+1)]+u[ind(i,j)]));
				j = ibObj[obj].iIB[ib].j+1;       rhsX[ind(i,j)] -= (+ dx2Rei*(u[ind(i,j)]+u[ind(i,j-1)]));
				i = ibObj[obj].iIB[ib].i-1;       rhsX[ind(i,j)] -= (+ dx2Rei*(u[ind(i,j)]+u[ind(i,j-1)]));
				j = ibObj[obj].iIB[ib].j-1;       rhsY[ind(i,j)] += (- dx2Rei*(v[ind(i+1,j)]+v[ind(i,j)]));
				j = ibObj[obj].iIB[ib].j;         rhsY[ind(i,j)] += (- dx2Rei*(v[ind(i+1,j)]+v[ind(i,j)]));
				i = ibObj[obj].iIB[ib].i+1;       rhsY[ind(i,j)] -= (+ dx2Rei*(v[ind(i,j)]+v[ind(i-1,j)]));
				j = ibObj[obj].iIB[ib].j-1;       rhsY[ind(i,j)] -= (+ dx2Rei*(v[ind(i,j)]+v[ind(i-1,j)]));
		}
    }
}



void calcPoissonRHS() {
    float dtdxi = 1.0/dt/dx;
#pragma omp parallel for num_threads(NUM_THREADS)
	for(int i=1; i<imax-1; i++)
	for(int j=1; j<jmax-1; j++)
	    pp[indp(i-1,j-1)] = dtdxi*(us[ind(i,j)] - us[ind(i-1,j)] + vs[ind(i,j)] - vs[ind(i,j-1)]);
}



void solvePressureDCTPxNy()
{
    int m = imax-2;
    int n = jmax-2;

	fftwf_execute(dctPxNy);
	for (int i=0; i<m*n; i++)   pp[i] /= evPxNy[i];
	fftwf_execute(idctPxNy);

	// normalize result
	float norm = 1.0/(1.0*m*2.0*n);
#pragma omp parallel for num_threads(NUM_THREADS)
	for(int i=1; i<imax-1; i++)
		for(int j=1; j<jmax-1; j++)
			ps[ind(i,j)] = pp[indp(i-1,j-1)]*norm;

	setBCpressure(ps);
}



void solvePressureDCTPxPy()
{
    int m = imax-2;
    int n = jmax-2;

	fftwf_execute(dctPxPy);
	for (int i=0; i<m*n; i++)   pp[i] /= evPxPy[i];
	fftwf_execute(idctPxPy);

	// normalize result
	float norm = 1.0/(1.0*m*1.0*n);
#pragma omp parallel for num_threads(NUM_THREADS)
	for(int i=1; i<imax-1; i++)
		for(int j=1; j<jmax-1; j++)
			ps[ind(i,j)] = pp[indp(i-1,j-1)]*norm;

	setBCpressure(ps);
}



void solvePressureTDMAxDCTy()
{
    int m = imax-2;
    int n = jmax-2;

    float evyLoc[n];

//	for(int i=0; i<m; i++)
//	for(int j=0; j<n; j++)
//	    pp[i*n+j] = cos(2.0*M_PI*(j+0.5)/n);
//
//FILE *fp = fopen("dummy","wt");
////	for(int i=0; i<m; i++) {
//	for(int j=0; j<n; j++) {
//		int i = 10;
//		fprintf(fp, "%lf\t%lf\n", 1.0*(j+0.5)/n, pp[i*n+j]);
//	}
//fclose(fp);

	// perform dct in y-direction
	switch (wall) {
	case 1:
	case 2:
		fftwf_execute(dctNy);
	    for (int j=0; j<n; j++)  evyLoc[j] = evNy[j];
		break;
	case 3:
		fftwf_execute(dctPy);
	    for (int j=0; j<n; j++)  evyLoc[j] = evPy[j];
	    break;
	}

	// use Thomas algorithm to solve in x-direction
	float a[m], b[m], c[m], d[m*n];

	for (int i=0; i<m; i++) {
		a[i] =  1.0/dx/dx;
		c[i] =  a[i];
		b[i] = -a[i]-c[i];
	}
	b[0] = b[m-1] = -1.0/dx/dx;

	for (int j=0; j<n; j++) {
		float z = 1.0/(b[0]+evyLoc[j]);
		d[j]  = c[0]*z;
		pp[j] = pp[j]*z;
	}

//#pragma omp parallel for num_threads(NUM_THREADS)
	for (int i=1; i<m-1; i++) {
		for (int j=0; j<n; j++) {
			float z = 1./(b[i]+evyLoc[j]-a[i]*d[(i-1)*n+j]);
			d[i*n+j]  = c[i]*z;
			pp[i*n+j] = (pp[i*n+j]-a[i]*pp[(i-1)*n+j])*z;
		}
	}

	for (int j=0; j<n; j++) {
		float z = b[m-1]+evyLoc[j]-a[m-1]*d[(m-2)*n+j];
		if (z != 0.0)
			pp[(m-1)*n+j] = (pp[(m-1)*n+j]-a[m-1]*pp[(m-2)*n+j])/z;
		else
			pp[(m-1)*n+j] = 0.0;
	}

//#pragma omp parallel for num_threads(NUM_THREADS)
	for (int i=m-2; i>=0; i--)
		for (int j=0; j<n; j++)
			pp[i*n+j] = pp[i*n+j]-d[i*n+j]*pp[(i+1)*n+j];

	// perform idct in y-direction
	float norm = 1.0;
	switch (wall) {
	case 1:
	case 2:
		fftwf_execute(idctNy);
		norm = 1.0/(2.0*n);
		break;
	case 3:
		fftwf_execute(idctPy);
		norm = 1.0/(n);
		break;
	}

	// normalize result
#pragma omp parallel for num_threads(NUM_THREADS)
	for(int i=1; i<imax-1; i++)
		for(int j=1; j<jmax-1; j++)
			ps[ind(i,j)] = pp[indp(i-1,j-1)]*norm;

	setBCpressure(ps);



//	fp = fopen("dummy2","wt");
////		for(int i=1; i<imax-1; i++) {
//		for(int j=1; j<jmax-1; j++) {
//			int i = 10;
//			fprintf(fp, "%lf\t%lf\n", 1.0*(j-0.5)/n, 1./pow(2.0*M_PI, -2.0)*ps[ind(i,j)]-0);
//		}
//	fclose(fp);
//	exit(0);
}






