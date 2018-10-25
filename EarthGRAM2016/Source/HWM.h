//Header for Harmonic Wind Model
//P. White

#include <fstream>
#include <string>

using namespace std;

class HWM
{
public:
	HWM();

	void gws5(int iyd, double sec, double alt, double glat, double glon, double stl,
		double f107a, double f107, double *ap, double *w);
	void spline(double *x, double *y, int n, double yp1, double ypn, double *y2);
	void splint(double *xa, double *ya, double *y2, int n, double x, double *y);
	void tselec(double *sv);
	//Getters
	double getCswSw(int index);
	int getCswIsw();
	double getCswSwc(int index);

private:

	void intializeMemberVariables();

	void gwsbk5();
	void initw5();
	void legpl1(double c, double s, int l, int m, double plg[][20], int lmax);
	void vsphr1(double c, double s, int l, int m, double bt[][20], double bp[][20],
		int lmax);
	void glbw5s(int iyd, double lat, double lon, double stl, double *pb, double *pc,
		double *ww);
	void glbw5e(double yrd, double sec, double lat, double lon, double stl,
		double f107a, double f107, double *ap, double *pb, double *pc, double *ww);
	void glbw5m(double yrd, double sec, double lat, double lon, double stl,
		double f107a, double f107, double *ap, double *pb, double *pc, double *ww);
	double wprof(double z, double zl, double s, double uinf, double ulb, double ulbd,
		int mn1, double *zn1, double *un1, double *ugn1, int mn2, double *zn2,
		double *un2, double *ugn2);
	double pwb[200], pwc[200], pwbl[150], pwcl[150], pwbld[150], pwcld[150],
		pb12[150];
	double pc12[150], pb13[150], pc13[150], pb14[150], pc14[150], pb15[150],
		pc15[150];
	double pb15d[150], pc15d[150], pwp[26][100], clat, slat, clong, slong, c2long,
		s2long;
	double cstl, sstl, c2stl, s2stl, c3stl, s3stl, sw[25], swc[25], bt[20][20],
		bp[20][20];
	double xvl, lvl, mvl, tll, nsvl, xll, ngvl, pi, pi180;
	int isw;
	string isdate[3], istime[2], name[2];

};


