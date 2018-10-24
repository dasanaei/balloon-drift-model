//header file for RRA class
#ifndef RRA_h_
#define RRA_h_

#include <string>
#include <fstream>
using namespace std;
class RRA
{
public:
	RRA();
	~RRA();
	void init(string namef);
	void initializeMemberVariables();
	void rrasigs(double h, double phi, double thet, int mn, double hsrf1, double usrf1,
		double vsrf1, double tsrf1, double uvtsrf1, double susrf1, double svsrf1,
		double stsrf1, double spdavsrf1, double spdsdsrf1, double r0, double uvth,
		double *sph, double *sdh, double *sth, double *suh, double *svh,
		string &rra_id, double *rrawt);
	void rramean(double h, double phi, double thet, int mn, double r0,
		double *pgh, double *dgh, double *tgh, double *ugh, double *vgh,
		string &rra_id, double *rrawt);
	void namelist();
	void geocenttogeodet(double r, double zin, double *fi, double *h, double a,
		double b);
	string rrasite[99], rrapath, rralist, rrayr1, namef1;
	double rralat[100], rralon[100], rrasfc[100], zmax[100];
	double z1[300];
	double z2[300];
	double z3[300];
	double pi, pi180, rrascale, sitelim, sitenear, rrawt3, avspdrra, sdspdrra, rrawt1,
		rrawt2;
	int num1, num2, num3, nrra, isitered, lenpath, iyrrra;
	const char *endsep;
	int nsrf1[99], nsrf2[99], nsrf3[99];
	double rrsrf, usrf, vsrf, susrf, svsrf, uvtsrf, spdavsrf, spdsdsrf,
		shsrf, psrf, dsrf, spsrf, sdsrf, tsrf, strrasrf, stsrf, uvrra,
		vprra, tdrra, svprra, stdrra, tdsrf, stdsrf, strra, trra, drra, prra,
		sitewgt;
	float *u, *su, *v, *sv, *ruvt, *avspd, *sdspd, *p, *sp, *t, *st, *d,
		*sd, *vp, *svp, *td, *std;

private:

	void readrra1(string &rra_id, double xlat, double xlon, int m, int isite);
	void readrra2(string &rra_id, double xlat, double xlon, int m, int isite);
	void readrra3(string &rra_id, double xlat, double xlon, int m, int isite);
	double heightwt(double h, double *z, int num);
	int numfind(double h, double *z, int num);
	double weightfact(double delta, double sitelim, double sitenear);
	double radll(double phi1, double thet1, double phi2, double thet2);
	void rig(double ch, double phir, double *g0, double *g, double *re,
		double *r0, double *ri, double *a, double *b);
	void geodettogeocent(double fidet, double h, double *ficent, double *rtot,
		double *xy, double *z, double a, double b);
	void interw(double u1, double v1, double z1, double u2, double v2, double z2,
		double *u, double *v, double z);
	void inter2(double p1, double d1, double t1, double z1, double p2, double d2,
		double t2, double z2, double *p, double *d, double *t, double z);
	void interz(double p1, double d1, double t1, double z1, double p2, double d2,
		double t2, double z2, double *p, double *d, double *t, double z);
};

#endif