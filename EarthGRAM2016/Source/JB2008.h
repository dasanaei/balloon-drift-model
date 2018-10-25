//Header for JB2008
//P. White

#include <fstream>
#include <string>
#include "HWM.h"

using namespace std;

class JB2008
{
public:
	JB2008();
	void intializeMemberVariables();

	void JB08mod(double h, double phi, double thet, double elt, double ri, double g, double *pmj, double *dmj, double *tmj, double *umj, double *vmj, double *wmj, double *n2nd,
		double *o2nd, double *ond, double *arnd, double *hend, double *hnd, double *wtmol,
		double *dmdz, double *nnd);
	void namelist(string namef);

	int mn, ida, iyr, ihro, mino, imin, ihr;
	double seco, f10, f10b, ap, delt, sec, s10, s10b, xm10, xm10b, y10, y10b, dstdtc,
		dtz;

	HWM hwm;


private:
	double xambar(double z);
	void semian08(int iyr, double day, double ht, double f10b, double s10b, double xm10b, double *fzz, double *gtz, double *drlog);
	void dtsub(double f10, double xlst, double xlat, double zht, double *dtc);
	void tmoutd(double d1950, int *iyr, double *day);
	double xlocal(double z, double tc[]);
	double xgrav(double z);
	void JB08(double amjd, double sun[], double sat[], double f10, double f10b, double s10, double s10b, double xm10, double xm10b, double y10, double y10b, double dstdtc, double temp[], double *rho, double *pres, double *avgmw, double and1[], double *sumn);
	void tme(int mn, int ida, int iyr, int ihr, double xmin, double xlng, double *xlat, double *sda, double *sha, double *dd, double *dy, double *sra, double *raloc, double *xmjd);
	void caltojul(int iy, int im, int id, int ihour, int imin, double sec, double *xjd);
	void JB08HWM(double xmin, double z, double xlat, double xlon, double *ph, double *dh, double *th,
		double *n2nd, double *o2nd, double *ond, double *arnd, double *hend, double *hnd, double *nnd, double *wtmol,
		double *tex, double *uh, double *vh);

};