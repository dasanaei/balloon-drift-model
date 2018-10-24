//Header for Marshall Engineering Thermosphere Model (MET)
//P White

#include <string>
#include <fstream>

using namespace std;

class Met{

public:
	Met();
	void initializeMemberVariables();
	double molwt(double a);
	void caltojul(int iy, int im, int id, int ihour, int imin, double sec, double *xjd);
	void tme(int mn, int ida, int iyr, int ihr, double xmin, double xlng, double *xlat, double *sda, double *sha, double *dd, double *dy, double *sra, double *raloc, double *xmjd);
	void jacmod(double h, double phi, double thet, double elt, double *pmj, double *dmj, double *tmj, double *umj, double *vmj, double *wmj, double *arnd, double *henda, double *hnda, double *o2nda, double *n2nda, double *onda, double *wtmola);
	double h, phi, thet, pi, pi180, sec, wtmol1, seco, f10, f10b, ap, delt, dtz, dmdz;
	int iyr, mn, ida, xmin, ihro, mino, imin, ihr;
	void namelist(string namef);
private:
	void gauss(int nmin, double z2, double tx, double t1, double t3, double t4, double a2, double gphi, double rphi, double *r);
	double temp(double alt, double tx, double t1, double t3, double t4, double a2);
	void met_07_jac(double z, double t, double *tz, double *an, double *ao2, double *ao, double *aa, double *ahe, double *ah, double *em, double *dens, double *dl, double gphi, double rphi);
	void slv(double alt, double xlat, double day, double *den);
	void slvh(double xlat, double sda, double *den, double *denhe);
	void met07(double *indata, double *outdata);
	void tinf(int i1, double f10, double f10b, double gi, double xlat, double sda, double sha, double dy, double *te);
	void jacch(int iyr, int m, int ida, int ihr, double xmin, double z, double phi, double thet, double f10, double f10b, double ap, double *ph, double *dh, double *th, double *n2nd, double *o2nd, double *ond, double *arnd, double *hend, double *hnd, double *wtmol, double *tex);
	void wind(double ph, double dh, double th, double h, double g, double phid, double ri, double dpx, double dpy, double dtx, double dty, double dtz, double *ugh, double *vgh, double *wgh);
	void rig(double ch, double phir, double *g0, double *g, double *re, double *r0, double *ri, double *a, double *b);
	void fair5(double dhel1, double dhel2, double dlg1, double dlg2, double h, double *fdhel, double *fdlg);

};