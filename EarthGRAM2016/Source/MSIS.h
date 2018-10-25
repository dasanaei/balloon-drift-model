//Header for MSIS model
//PWhite

#include <string>
#include <fstream>

#include "JB2008.h"
#include "Map.h"

using namespace std;


class MSIS
{

public:
	MSIS();
	~MSIS();
	void intializeMemberVariables();

	void msismod(double h, double phi, double thet, double elt, double ri, double g, double *pmj,
		double *dmj, double *tmj, double *umj, double *vmj, double *wmj, double *n2nd,
		double *o2nd, double *ond, double *arnd, double *hend, double *hnd,
		double *wtmol, double *dmdz, double *nnd);

	void namelist(string namef);

	double seco, f10, f10b, ap, sec, delt, dtz;
	int iyr, mn, ida, xmin, ihro, mino, imin, ihr;
	int imr;
	string name[2], isdate[3], istime[2];
	void gtd7bk();

	Map maps;
	JB2008 jb2008;

private:
	void ghp7(int iyd, double sec, double *alt, double glat, double glong, double stl, double f107a, double f107, double *ap, double *d, double *t, double press);
	void gtd7(int iyd, double sec, double alt, double glat, double glong, double stl, double f107a, double f107, double ap, double mass, double *d, double *t);
	void gts7(int iyd, double sec, double alt, double glat, double glong, double stl, double f107a, double f107, double *ap, int mass, double *dm28, double *d, double *t);
	void msishwm(int iyr, int mn, int ida, int ihr, double xmin, double z, double xlat,
		double xlon, double f10b, double f10, double *ap, double *ph, double *dh,
		double *th, double *n2nd, double *o2nd, double *ond, double *arnd, double *hend,
		double *hnd, double *nnd, double *wtmol, double *tex, double *uh, double *vh);
	double globe7(double yrd, double sec, double lat, double lon, double tloc, double f107a, double f107, double *ap, double *p);
	double densm(double alt, double d0, double xm, double *tz, int mn3, double *zn3, double *tn3, double tgn3[2], int mn2, double *zn2, double *tn2, double tgn2[2]);
	double densu(double alt, double dlb, double tinf, double tlb, double xm, double alpha, double *tz, double zlb, double s2, int mn1, double *zn1, double *tn1, double *tgn1);
	double glob7s(double *p);
	double dnet(double dd, double dm, double zhm, double xmm, double xm);
	double ccor(double alt, double r, double h1, double zh);
	double ccor2(double alt, double r, double h1, double zh, double h2);
	double zeta(double zz, double zl, double re) { return (zz - zl)*(re + zl) / (re + zz); }
	double scalh(double alt, double xm, double temp);
	double g0(double a, double *p)
	{
		return (a - 4.0 + (p[25] - 1.0)*(a - 4.0 + (exp(-fabs(p[24])*(a - 4.0)) - 1.0) / fabs(p[24])));
	}
	double sg0(double ex, double *ap, double *p)
	{
		return (g0(ap[1], p) + (g0(ap[2], p)*ex + g0(ap[3], p)*ex*ex
			+ g0(ap[4], p)*pow(ex, 3) + (g0(ap[5], p)*pow(ex, 4)
			+ g0(ap[6], p)*pow(ex, 12))*(1.0 - pow(ex, 8)) / (1.0 - ex))) / sumex(ex);
	}
	double sumex(double ex) { return (1.0 + (1.0 - pow(ex, 19)) / (1.0 - ex)*pow(ex, 0.5)); }
	double vtst7(int iyd, double sec, double glat, double glong, double stl, double f107a, double f107, double *ap, int ic);
	void meters(bool meter);
	void splini(double *xa, double *ya, double *y2a, int n, double x, double *yi);
	void glatf(double lat, double *gv, double *refff);
	double pi, pi180, plg[9][4], dd, ctloc, stloc, c2tloc, s2tloc, c3tloc;
	double s3tloc, day, dfa, apdf, apt[4], xlong, tn1[5], tn2[4], tn3[5];
	double tgn1[2], tgn2[2], tgn3[2], gsurf, re, pt[150], pd[9][150];
	double ps[150], ptl[4][100], pma[10][100], sam[100];
	double pavgm[10], ptm[10], pdm[10][8];

	float **pdl;

};