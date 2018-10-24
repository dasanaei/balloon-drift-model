//header file for Auxilliary Profile class

#include <string>
#include <fstream>
#include <vector>

using namespace std;

class AuxProf {

public:

	AuxProf();
	void rdprof(double sitenear, double sitelim, string home1, string profile1);


	void profsigs(double chgt, double clat, double clon, double tin,
		double pin, double din, double uin, double vin, double *ptemp,
		double *ppres, double *pdens, double *puwin, double *pvwin,
		double *profwgt);

	void profterp(double chgt, double clat, double clon, double tin,
		double pin, double din, double uin, double vin, double *ptemp,
		double *ppres, double *pdens, double *puwin, double *pvwin,
		double *profwgt);



private:

	void rig(double ch, double phir, double *g0, double *g, double *re,
		double *r0, double *ri, double *a, double *b);
	void namelist();
	void initializeMemberVariables();
	double radll(double phi1, double thet1, double phi2, double thet2);
	vector<double> phgt, plat, plon, ptmp, pprs, pden, puwn, pvwn,
		pstmp, psprs, psden, psuwn, psvwn;
	double profnear, proffar, pi, pi180;
	int nprof;
	string profile;
};