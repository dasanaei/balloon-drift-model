#ifndef Pert_h_
#define Pert_h_

#include <string>
#include <fstream>

#include "InitP.h"


using namespace std;

class Pert
{
public:
	Pert();
	~Pert();
	void initializeMemberVariables();
	void pert1(double h, double phi, double thet, int iupdate, int nm);
	void rig(double ch, double phir, double *g0, double *g, double *re, double *r0, double *ri);
	void initsigs();
	double dz, dx, dtime, sd1s, sp1s, st1s, su1s, sv1s, sp1l, sd1l, st1l, sv1l, su1l,
		uds1, vds1, rp1s, rd1s, rt1s, ru1s, rv1s, prh, drh, trh, urh, vrh, pi, pi180;
	double sphs, sdhs, sths, suhs, svhs, sphl, sdhl, sthl, suhl, svhl, prhs, drhs,
		trhs, urhs, vrhs, prhl, drhl, trhl, urhl, vrhl;
	double rp1l, rd1l, rt1l, ru1l, rv1l;
	double phi, thet, thet1, phi1, dphi, dthet, sph, sdh, sth, suh, svh,
		wrh, elt, xlh, zlh, uvt2;
	int isev;
	InitPert iperts;

private:

	void coeff(double sd1, double sd2, double sp1, double sp2, double st1,
		double st2, double su1, double su2, double sv1, double sv2, double ud1,
		double ud2, double vd1, double vd2, double rpd, double rp, double rd,
		double rv, double *a, double *b, double *c, double *d, double *e,
		double *f, double *g, double *h, double *ai, double *aj, double *ak);
	void fair(int i, double hg, double pg, double dg, double tg, double ug, double vg, double wg,
		double hj, double pj, double dj, double tj, double uj, double vj, double wj, double h,
		double *p, double *d, double *t, double *u, double *v, double *w, double *czi);
	double radll(double phi1, double thet1, double phi2, double thet2);

};

#endif