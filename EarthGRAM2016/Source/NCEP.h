//Header for NCEP data (lower atmosphere data)
//P. White

#ifndef NCEP_h_
#define NCEP_h_

#include <string>
#include <fstream>


using namespace std;

class NCEPmods
{

public:

	NCEPmods();
	~NCEPmods();

	void initializeMemberVariables();
	void NCEPread();

	void ncepmd(double z, double phi, double thet, double *pz, double *rhoz,
		double *tz, double *uz, double *vz, double *wz, double *tdz,
		double *spz, double *srhoz, double *stz, double *suz, double *svz,
		double *stdz, double *rhn, double *srhn, double *vpn, double *svpn,
		double *spdavz, double *spdsdz, double *uvcorrz);

	void namelist(string namef);

	void gethgs(double phi, double thet, double *hg1, double *hg2);

	float ***temp, ***dens, ***dewp, ***uwnd, ***vwnd, ***geop, **slpn,
		**sfcp, ***sden, ***sdewp, ***suwd, ***svwd, ***sprs, **sslp,
		**ssfcp, ***rhum, ***srhum, ***vprs, ***svprs, ***spdav, ***spdsd,
		***uvcor, ***stmp;

	double pi, pi180, dtz;
	string NCEPpath1;
	int NCEPhr;


private:

	void ncepterp(double alpha, double beta, double value[], double *val);

	double gascon(double td, double p);

	double wexler(double t, double p);

	double ztoh(double g0, double r, double z);

	void rig(double ch, double phir, double *g0, double *g, double *re,
		double *r0, double *ri, double *a, double *b);

	void sortlevel(int lsort[][19], double zsort[][19]);

	void ncepzterp(double g0, double re, double z, double *pzo,
		double *tzo, double *tdzo, double *rhozo, double *uzo,
		double *vzo, double *wzo, double *spzo, double *stzo,
		double *stdzo, double *srhozo, double *suzo, double *svzo,
		double *dtdzo, double *rhno, double *srhno, double *vpno,
		double *svpno, double *spdavzo, double *spdsdzo,
		double *uvcorrzo, double alpha, double beta);

	float **t, **p, **h, **td, **u, **v, **w, **st, **std, **su, **sv,
		**sp, **rho, **srho, **vp, **svp, **spdavl, **spdsdl, **uvcorrl,
		**rh, **srh;


};

#endif