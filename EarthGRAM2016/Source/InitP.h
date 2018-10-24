//header file for Initial Perturbation

#ifndef InitPert_h_
#define InitPert_h_

#include "Init.h"
#include "RRA.h"
#include "MET.h"
#include "AuxProf.h"
#include "NCEP.h"

#include <string>
#include <fstream>

using namespace std;
class InitPert
{
public:
	InitPert();
	~InitPert();
	void initializeMemberVariables();

	void initpert(double h, double phi, double thet, double elt, int nr1);
	void rterp(double h, double phi, double *x, double *y, double *z);
	void intr25(float **xarray, float **yarray, double h, double phi,
		double *suh, double *svh);
	void intruv(float **uarr, float **varr, double h, double phi, double *su,
		double *sv);
	void interz(double p1, double d1, double t1, double z1, double p2,
		double d2, double t2, double z2, double *p, double *d, double *t,
		double z);
	void intrw(double sarr[], double h, double *sarrw);
	void interw(double u1, double v1, double z1, double u2, double v2,
		double z2, double *u, double *v, double z);
	void getsigw(double h, double phi, double thet, double alpha,
		double alphav, double vlls, double dphiu, double dphiv, double elt);
	void rcarry(double *rvec, int lenv);
	double ppnd(double p, int *ifault);
	double correl(double x);
	double ptail(double x1);
	void zinterp(double clat, double clon);
	int mi24, mj24, nr1;
	double dx, dz, dtime, xbar1, zbar1, xl1, zl1, sxl1, szl1;
	double prh1, drh1, trh1, urh1, vrh1, sph1, sp1l, sp1s, sd1l, sd1s, st1l,
		st1s, su1l, su1s, sv1l, sv1s, uds1, vds1;
	double rp1s, rd1s, rt1s, ru1s, rv1s, waverand, phidens, ampfact;
	double rp1l, rd1l, rt1l, ru1l, rv1l;
	double prh, drh, trh, urh, vrh, sd1, sp1, st1, mseeds[24], mcarry;
	double su1, sv1, hsrf1, usrf1, vsrf1, susrf1, svsrf1, tsrf1, uvtsrf1,
		stsrf1, spdavsrf1, spdsdsrf1;
	double twopi, pi, pi180;
	double swh, sw1, rw1, udlsrf, ustar, uvth;


	Init inits1;
	RRA rras;
	Met mets;
	AuxProf auxs;
	NCEPmods ncps;


private:


	void rcarin(int ijkl);

};

#endif
