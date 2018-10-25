//header file for Atmod


#ifndef Atmod_h_
#define Atmod_h_



#include "MSIS.h"


#include <string>
#include <fstream>

using namespace std;

class Atm
{
public:

	Atm();
	~Atm();

	ofstream out1;
	void initdata();
	void initpert();
	void pert();
	void map();
	void ncep();
	void atmos();
	void initializeMemberVariables();
	void met();
	void jb2008();
	void msis();
	void rra();
	void prof();
	void speconc();
	void fair(int i, double hg, double pg, double dg, double tg, double ug, double vg, double wg,
		double hj, double pj, double dj, double tj, double uj, double vj, double wj, double h,
		double *p, double *d, double *t, double *u, double *v, double *w, double *czi);
	void stdatm(double z, double *t, double *p, double *d);
	double wexler(double t, double p);
	double mixrat(double tx, double p);
	double tdbuck(double t, double rh, double p);
	double h, thet, phi, dphi, dthet, dz, delt, pmean, dmean, tmean, umean, vmean, wmean;
	double ppert, dpert, tpert, upert, vpert, pmpert, dmpert, tmpert, umpert, vmpert;
	double pstd, dstd, tstd, ustd, vstd, pi, piby2, pmean1, dmean1, tmean1, umean1, vmean1, wmean1,
		pmean2, dmean2, tmean2, umean2, vmean2, wmean2, hj1, hj2;
	double prramean, drramean, trramean, urramean, vrramean, prrastd,
		drrastd, trrastd, urrastd, vrrastd, pnmean, dnmean, tnmean, unmean,
		vnmean, wnmean, wpert, wmpert, wstd, vpn, rhn, tdmean, stdd, stg, svpn, srhn,
		wtmol, waterg, eoft, seoft, airmw;
	double wth2o, wto3, wtn2o, wtco, wtch4, wtco2, wtn2, wto2, wto,
		wtar, wthe, wth, wtn, wtair;
	double pstd1, dstd1, tstd1, ustd1, vstd1, wstd1, rhp, srhp, eps, pi180;
	double ppert1, dpert1, tpert1, upert1, vpert1, wpert1, elt;
	double ppmo, ppmh, ppmhe, ppmar, ppmo2, ppmn2, ppmn, ppmo3, ppmn2o,
		ppmco, ppmch4, ppmco2, ppmh2o;
	double arnd, hend, hnd, o2nd, n2nd, nnd, ond, o3nd, h2ond, n2ond, cond,
		ch4nd, co2nd, ppmtond, avn, phi0, thet0, spdavg, spdsd, spg, sdg,
		sug, svg, Umean, Vmean, rhov, srhov, dtz1, dtzm, dmdz1, dmdzm;
	int nr1, yr, nm;
	int mn;



	MSIS msisa;

private:



};

#endif