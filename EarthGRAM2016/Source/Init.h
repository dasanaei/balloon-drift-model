//header file for Initial data

#ifndef Init_h_
#define Init_h_

#include <string>
#include <fstream>


using namespace std;
class Init
{
public:
	Init();
	~Init();
	ifstream traj1;
	ofstream ibl;
	ofstream out;
	void initializeMemberVariables();
	void init1(int month);
	void traj();
	void trajopen();
	void namelist();
	void outopen();
	float ***psp, ***dsp, ***tsp, ***usp, ***vsp, **pg, **dg, **tg, **ug;
	float **pr, **dr, **tr, **ur, **vr, **plp, **dlp, **tlp, **ulp, **vlp;
	float **uds, **vds, **udl, **uvt;
	float **h2ol, **h2oa, **o3, **n2o, **co, **ch4, **o3map, **h2omap, **sh2omap,
		**n2omap, **ch4map, **oxymap;

	double xlbar[29], zlbar[29], xsigl[29], zsigl[29], xlmin[29], zlmin[29], xscale[29], zscale[29], wr[29];
	double hgtl[35], hgta[50];
	double po3map[24], ph2omap[19];
	double pmap31[17], mapzox[19];
	double co2[50], o2[50], n2[50];
	double water, ozone, nitrous, carbmon, methane, carbdiox, moloxy,
		nitrogen, sigwater, wtmol, ond;
	double time, hgt, lat, lon, lat1, lon1, delt1, dtime1;
	double dz1, hgt1, dphi1, dthet1;
	int f;
	int ** ztopo;
	int ** landcd;

	string NCEPpath1, atmpath, trapath, prtpath, nprpath, conpath,
		rrapath, rralist, profile, home, namelst, therm, NCEPpath, patch, Ryear;

	double h1, phi1, thet1, f10, f10b, ap, s10, s10b, xm10, xm10b, y10,
		y10b, dstdtc, seco, dphi, dthet, dhgt, delt, rpscale, ruscale,
		rwscale, sitelim, sitenear, rdinit, rtinit, ruinit, rvinit, rwinit,
		z0in, patchy, time1;
	const char * endsep;
	int mn, ida, iyr, ihro, mino, nmax, iopt, iopp, iu0, iup, ius, iuc, iug,
		NCEPyr, NCEPhr, iopr, nr1, iun, iurra, iyrra, initpert,
		itherm, ibltest, mc, iaux;


private:

	void rtran1(ifstream &atmosdat, string *i1, int *i2, int *i3, int *ndata, int *i5);
	void rtran2(ifstream &atmosdat, string *i1, int *i2, int *i3, int *i4, int *ndata1);
	void rtran3(ifstream &atmosdat, string *i1, int *i2, int *i3, int *ndata2);
	void rtran4(ifstream &atmosdat, string *i1, int *i2, int *i3, int *ndata3);
	void concinit(ifstream &atmosdat, int imon);
	void mapinit(ifstream &atmosdat, int imon);
	void topo();

};

#endif
