//Header file for Middle Atm Program class

#ifndef Map_h_
#define Map_h_

#include <string>
#include <fstream>

#include "Pert.h"

using namespace std;

class Map {
public:
	Map();
	void mapmod(double h, double phi, double thet, double ri, double g, double *pmm, double *dmm, double *tmm, double *umm, double *vmm, double *wmm, double *dtz);
	void concvals(double z, double phi, int iyr, double pgh, double waterg,
		double *oxygen);
	void concmod(double h, double phi, double pgh, double dgh, double tgh,
		double tdg, double stdg, int iyr);
	double dedt(double t, double p);
	double wexler(double t, double p);
	double d2edt2(double t, double p);
	double valint(int nlat, int nz, int ilat, int iz, double alpha,
		double beta, float **value);
	double valz(int nz, int iz, double beta, double *value);
	Pert perts;

private:
	void gterp(int ih, double phi, double *p, double *d, double *t, double *u, double *dpy, double *dty);
	void pdtuv(int ih, double clat, double clon, double *ps, double *ds, double *ts, double *us, double *vs, double *dpx, double *dpy, double *dtx, double *dty);
	void inter2(double p1, double d1, double t1, double z1, double p2, double d2, double t2, double z2, double *p, double *d, double *t, double z);
	void afglconc(double z, double phi, double *watera);
	void mapconc(double z, double pgh, double phi, double *o3m, double *h2om,
		double *sh2om, double *n2om, double *ch4m, double *oxym);
	void larcwat(double z, double phi, double *waterl);
};
#endif
