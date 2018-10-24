//P. White
//Harmonic Wind Model Class

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include "HWM.h"

using namespace std;

HWM::HWM()
{

	intializeMemberVariables();

}

void HWM::intializeMemberVariables()
{

	//Initialize member variables for HWM class

	clat = 0.0;         // -- Cosine of the latitude
	slat = 0.0;         // -- Sine of the latitude
	pi = 3.1415926535897931;
	pi180 = pi / 180.0;

	for (int i = 0; i < 20; ++i) {
		for (int j = 0; j < 20; ++j) {
			bt[i][j] = 0.0;
			bp[i][j] = 0.0;
		}
	}

	clong = 0.0;       // -- Cosine of longitude
	slong = 0.0;       // -- Sine of longitude
	c2long = 0.0;       // -- Cosine of 2 * longitude
	s2long = 0.0;       // -- Sine of 2 * longitude

	cstl = 0.0;        // -- Cosine of apparent local solar time
	sstl = 0.0;        // -- Sine of apparent local solar time
	c2stl = 0.0;        // -- Cosine of 2 * apparent local solar time
	s2stl = 0.0;        // -- Sine of 2 * apparent local solar time
	c3stl = 0.0;        // -- Cosine of 3 * apparent local solar time
	s3stl = 0.0;        // -- Sine of 3 * apparent local solar time


	for (int i = 0; i < 25; ++i) {
		sw[i] = 1.0;   // -- Main terms. Hard coded to 1.0.
		swc[i] = 1.0;   // -- Cross terms. Hard coded to 1.0.
	}

	// BLOCK DATA initializers
	gwsbk5();
	initw5();

	return;
}

double HWM::wprof(double z, double zl, double s, double uinf, double ulb, double ulbd, int mn1, double *zn1, double *un1, double *ugn1, int mn2, double *zn2, double *un2, double *ugn2)
{

	//wprof member function from HWM class
	//Compute winds at altitude

	double x, f, z1, z2, zdif, xs[15], ys[15], yd1, yd2, y2out[15], y;
	int mn, k;

	double wprof_out = 0.0;

	if (z >= zl) {

		x = s*(z - zl);
		f = exp(-x);

		wprof_out = uinf + (ulb - uinf)*f + (ulb - uinf + ulbd)*x*f;

		return wprof_out;
	}

	if ((z >= *(zn1 + mn1 - 1)) && (z < *(zn1 + 0))) {

		mn = mn1;
		z1 = *(zn1 + 0);
		z2 = *(zn1 + mn - 1);
		zdif = z2 - z1;

		for (k = 0; k < mn; ++k) {
			xs[k] = (*(zn1 + k) - z1) / zdif;
			ys[k] = *(un1 + k);
		}

		yd1 = *(ugn1 + 0)*zdif;
		yd2 = *(ugn1 + 1)*zdif;

		spline(xs, ys, mn, yd1, yd2, y2out);


		x = (z - z1) / zdif;


		splint(xs, ys, y2out, mn, x, &y);

		wprof_out = y;

		return wprof_out;
	}

	if (z < *(zn2 + 0)) {

		mn = mn2;
		z1 = *(zn2 + 0);
		z2 = *(zn2 + mn - 1);
		zdif = z2 - z1;

		for (k = 0; k < mn; ++k) {
			xs[k] = (*(zn2 + k) - z1) / zdif;
			ys[k] = *(un2 + k);
		}

		yd1 = *(ugn2 + 0);
		if (1.0e30 > *(ugn2 + 0)) {
			yd1 = *(ugn2 + 0)*zdif;
		}

		yd2 = *(ugn2 + 1);
		if (1.0e30 > *(ugn2 + 1)) {
			yd2 = *(ugn2 + 1)*zdif;
		}

		spline(xs, ys, mn, yd1, yd2, y2out);


		x = (z - z1) / zdif;


		splint(xs, ys, y2out, mn, x, &y);

		wprof_out = y;

		return wprof_out;
	}

	return wprof_out;

}



void HWM::legpl1(double c, double s, int l, int m, double plg[][20], int lmax)
{

	//legpl1 member function from HWM class
	//Calculate legendre polynomials plg[l+1][m+1] through order l, m
	//for cosine c and sine s of colatitude.

	int mm, mmx, ll;

	if ((m > l) || (l > lmax - 1)){

		cout << "Illegal indices to legpol" << '\n';

		return;

	}
	plg[0][0] = 1.0;

	if ((l == 0) && (m == 0)) {
		return;
	}

	//Calculate L=M case and L=M+1

	for (mm = 0; mm < m + 1; mm++){
		if (mm > 0){
			plg[mm][mm] = plg[mm - 1][mm - 1] * (2.0*mm - 1.0)*s;
		}

		if (l > mm){
			plg[mm + 1][mm] = plg[mm][mm] * (2.0*mm + 1)*c;
		}
	}

	if (l == 1){
		return;
	}

	mmx = min(m, l - 2);

	for (mm = 0; mm < mmx + 1; mm++){
		for (ll = mm + 2; ll < l + 1; ll++){
			plg[ll][mm] = ((2.0*ll - 1.0)*c*plg[ll - 1][mm] -
				(ll + mm - 1.0)*plg[ll - 2][mm]) / (ll - mm);
		}
	}

	return;
}

void HWM::vsphr1(double c, double s, int l, int m, double bt[][20], double bp[][20], int lmax)
{

	//vsphr1 member function from HWM class 
	//Calculates vector spherical harmonic b field theta and phi functions bt, 
	//bp through order l, m for colatitude (theta) with cosine c and sine s
	//of colatitude

	double plg[20][20] = { { 0.0 } };
	double ic, sqt;
	int lmx, ll, mm;

	if ((m > l) || (l > lmax - 1)){
		return;
	}

	bt[0][0] = 0.0;
	bp[0][0] = 0.0;

	if ((l == 0) && (m == 0)){
		return;
	}

	legpl1(c, s, l + 1, m, plg, 20);

	if (fabs(s) < 1.0e-5){
		ic = copysign(1.0, s);
		s = 0;
	}

	for (ll = 1; ll < l + 1; ll++){
		sqt = sqrt(double(ll)*(double(ll) + 1.0));
		lmx = min(ll, m);

		for (mm = 0; mm < lmx + 1; mm++){
			if (s == 0.0){
				if (mm != 1){
					bt[ll][mm] = 0.0;
					bp[ll][mm] = 0.0;
				}
				else {
					bt[ll][mm] = (ll*(ll + 1.0)*(ll + 2.0)*0.5*pow(ic, ll + 2.0) - (ll + 1.0)*c*ll*
						(ll + 1.0)*0.5*pow(ic, ll + 1.0)) / sqt;
					bp[ll][mm] = mm*ll*(ll + 1.0)*0.5*pow(ic, ll + 1.0) / sqt;

				}
			}
			else {
				bt[ll][mm] = ((ll - mm + 1.0)*plg[ll + 1][mm] - (ll + 1.0)*c*plg[ll][mm]) / (s*sqt);
				bp[ll][mm] = mm*plg[ll][mm] / (s*sqt);
			}
		}
	}

	return;
}


void HWM::glbw5s(int iyd, double lat, double lon, double stl, double *pb, double *pc, double *ww)
{
	//glbw5s member function from HWM class

	const double hr = 0.2618, dr = 1.72142e-02, pset = 5.0;
	const int nsw = 14;
	double pb14 = -1.0, pb18 = -1.0, pc14 = -1.0, pc18 = -1.0, wb[2][15] = { { 0.0 } };
	double wc[2][15] = { { 0.0 } }, day, dayl = 0.0, f1b, f1c, f5b, f5c, f7b, f75b, f7c, f75c;
	double f8b, f8c, pb19 = 0.0, cd19b, cd14b, cd14c, cd18b, cd18c, cd25b, cd26b, cd32c, cd39c;
	double cd64c, cd87c, pb25 = 0.0, pb26 = 0.0, pc32 = 0.0, pc39 = 0.0, pc64 = 0.0, pc87 = 0.0, wbt[2], wct[2];
	int j, k, iyr, lv, mv, nsv, ngv;

	//Confirm parameter set
	if (*(pb + 99) == 0.0){
		*(pb + 99) = pset;
	}

	if (*(pb + 99) != pset){
		return;
	}

	for (j = 0; j < nsw; j++){
		wb[0][j] = 0.0;
		wb[1][j] = 0.0;
		wc[0][j] = 0.0;
		wc[1][j] = 0.0;
	}

	iyr = iyd / 1000;
	day = iyd - iyr*1000.0;

	lv = 7;
	mv = 2;

	if ((xvl != lat) || (lv > lvl) || (mv > mvl)){
		slat = sin(pi180*lat);
		clat = cos(pi180*lat);

		vsphr1(slat, clat, lv, mv, bt, bp, 20);

		xvl = lat;
		lvl = lv;
		mvl = mv;
	}

	nsv = 2;

	if ((tll != stl) || (nsv > nsvl)){
		sstl = sin(hr*stl);
		cstl = cos(hr*stl);
		s2stl = sin(2.0*hr*stl);
		c2stl = cos(2.0*hr*stl);
		tll = stl;
		nsvl = nsv;
	}

	if ((day != dayl) || (*(pb + 13) != pb14)){
		cd14b = cos(dr*(day - *(pb + 13)));
	}

	if ((day != dayl) || (*(pc + 13) != pc14)) {
		cd14c = cos(dr*(day - *(pc + 13)));
	}

	if ((day != dayl) || (*(pb + 17) != pb18)) {
		cd18b = cos(2.0*dr*(day - *(pb + 17)));
	}

	if ((day != dayl) || (*(pc + 17) != pc18)) {
		cd18c = cos(2.0*dr*(day - *(pc + 17)));
	}

	if ((day != dayl) || (*(pb + 18) != pb19)) {
		cd19b = cos(2.0*dr*(day - *(pb + 18)));
	}

	if ((day != dayl) || (*(pb + 24) != pb25)) {
		cd25b = cos(dr*(day - *(pb + 24)));
	}

	if ((day != dayl) || (*(pb + 25) != pb26)) {
		cd26b = cos(dr*(day - *(pb + 25)));
	}

	if ((day != dayl) || (*(pc + 31) != pc32)) {
		cd32c = cos(dr*(day - *(pc + 31)));
	}

	if ((day != dayl) || (*(pc + 38) != pc39)) {
		cd39c = cos(2.0*dr*(day - *(pc + 38)));
	}

	if ((day != dayl) || (*(pc + 63) != pc64)) {
		cd64c = cos(dr*(day - *(pc + 63)));
	}

	if ((day != dayl) || (*(pc + 86) != pc87)) {
		cd87c = cos(2.0*dr*(day - *(pc + 86)));
	}

	dayl = day;
	pb14 = *(pb + 13);
	pc14 = *(pc + 13);
	pb18 = *(pb + 17);
	pc18 = *(pc + 17);
	pb19 = *(pb + 18);
	pb25 = *(pb + 24);
	pb26 = *(pb + 25);
	pc32 = *(pc + 31);
	pc39 = *(pc + 38);
	pc64 = *(pc + 63);
	pc87 = *(pc + 86);

	ngv = 1;

	if (xll != lon || ngv > ngvl) {
		slong = sin(pi180*lon);
		clong = cos(pi180*lon);
		xll = lon;
		ngvl = ngv;
	}

	// Time independent
	f1b = 1.0;

	if (ww[0] != 9898.0) {
		wb[0][1] = (*(pb + 1)*bt[2][0] + *(pb + 2)*bt[4][0] + *(pb + 22)*bt[6][0])*f1b;
	}

	wb[1][1] = 0.0;
	wc[0][1] = 0.0;

	f1c = 1.0;

	if (ww[1] != 9898.0) {
		wc[1][1] = -(*(pc + 1)*bt[1][0] + *(pc + 2)*bt[3][0] + *(pc + 22)*bt[5][0])*f1c
			- (*(pc + 26)*bt[2][0] + *(pc + 14)*bt[4][0] + *(pc + 59)*bt[6][0])*f1c;
	}

	// symmetrical annual
	if (ww[1] != 9898.0) {
		wc[1][2] = -(*(pc + 47)*bt[1][0] + *(pc + 29)*bt[3][0])*cd32c;
	}

	// symmetrical semiannual
	if (ww[0] != 9898.0) {
		wb[0][3] = (*(pb + 16)*bt[2][0] + *(pb + 30)*bt[4][0])*cd18b;
	}

	wb[1][3] = 0.0;
	wc[0][3] = 0.0;

	if (ww[1] != 9898.0) {
		wc[1][3] = -(*(pc + 16)*bt[1][0] + *(pc + 30)*bt[3][0] + *(pc + 49)*bt[5][0])*cd18c;
	}


	// Asymmetrical annual
	f5b = 1.0;

	if (ww[0] != 9898.0) {
		wb[0][4] = (*(pb + 9)*bt[1][0] + *(pb + 10)*bt[3][0])*cd14b*f5b;
	}

	wb[1][4] = 0.0;
	wc[0][4] = 0.0;

	f5c = 1.0;

	if (ww[1] != 9898.0) {
		wc[1][4] = -(*(pc + 9)*bt[2][0] + *(pc + 10)*bt[4][0] + *(pc + 20)*bt[6][0])*cd14c*f5c;
	}

	// Asymmetrical semiannual
	if (ww[1] != 9898.0) {
		wc[1][5] = -(*(pc + 37)*bt[2][0] + *(pc + 98)*bt[4][0])*cd39c;
	}

	// Diurnal
	if (sw[6] != 0.0) {
		f7b = 1.0;
		f75b = 1.0;

		if (ww[0] != 9898.0) {
			wb[0][6] = (*(pb + 6)*bt[1][1] + *(pb + 7)*bt[3][1] + *(pb + 28)*bt[5][1]
				+ *(pb + 88)*bt[2][1])*sstl*f7b + (*(pb + 12)*bt[2][1])*cd25b*sstl*f75b*swc[4]
				+ (*(pb + 3)*bt[1][1] + *(pb + 4)*bt[3][1] + *(pb + 27)*bt[5][1]
				+ *(pb + 87)*bt[2][1])*cstl*f7b + (*(pb + 11)*bt[2][1])*cd25b*cstl*f75b*swc[4];
		}

		if (ww[1] != 9898.0) {
			wb[1][6] = -(*(pb + 3)*bp[1][1] + *(pb + 4)*bp[3][1] + *(pb + 27)*bp[5][1]
				+ *(pb + 87)*bp[2][1])*sstl*f7b - (*(pb + 11)*bp[2][1]) *cd25b*sstl*f75b*swc[4]
				+ (*(pb + 6)*bp[1][1] + *(pb + 7)*bp[3][1] + *(pb + 28)*bp[5][1]
				+ *(pb + 88)*bp[2][1])*cstl*f7b + (*(pb + 12)*bp[2][1])*cd25b*cstl*f75b*swc[4];
		}

		f7c = 1.0;
		f75c = 1.0;

		if (ww[0] != 9898.0) {
			wc[0][6] = -(*(pc + 3)*bp[2][1] + *(pc + 4)*bp[4][1] + *(pc + 27)*bp[6][1]
				+ *(pc + 87)*bp[1][1])*sstl*f7c - (*(pc + 11)*bp[1][1])*cd25b*sstl*f75c*swc[4]
				+ (*(pc + 6)*bp[2][1] + *(pc + 7)*bp[4][1] + *(pc + 28)*bp[6][1]
				+ *(pc + 88)*bp[1][1])*cstl*f7c + (*(pc + 12)*bp[1][1])*cd25b*cstl*f75c*swc[4];
		}

		if (ww[1] != 9898.0) {
			wc[1][6] = -(*(pc + 6)*bt[2][1] + *(pc + 7)*bt[4][1] + *(pc + 28)*bt[6][1]
				+ *(pc + 88)*bt[1][1])*sstl*f7c - (*(pc + 12)*bt[1][1])*cd25b*sstl*f75c*swc[4]
				- (*(pc + 3)*bt[2][1] + *(pc + 4)*bt[4][1] + *(pc + 27)*bt[6][1]
				+ *(pc + 87)*bt[1][1])*cstl*f7c - (*(pc + 11)*bt[1][1])*cd25b*cstl*f75c*swc[4];
		}
	}

	// semidiurnal
	if (sw[7] != 0.0) {
		f8b = 1.0;

		if (ww[0] != 9898.0) {
			wb[0][7] = (*(pb + 8)*bt[2][2] + *(pb + 42)*bt[4][2] + *(pb + 34)*bt[6][2]
				+ *(pb + 97)*bt[3][2] + (*(pb + 33)*bt[3][2])*cd26b*swc[4]
				+ (*(pb + 36)*bt[3][2])*cd19b*swc[5])*s2stl*f8b
				+ (*(pb + 5)*bt[2][2] + *(pb + 41)*bt[4][2] + *(pb + 32)*bt[6][2]
				+ *(pb + 95)*bt[3][2] + (*(pb + 23)*bt[3][2])*cd26b*swc[4]
				+ (*(pb + 35)*bt[3][2])*cd19b*swc[5])*c2stl*f8b;
		}

		if (ww[1] != 9898.0) {
			wb[1][7] = -(*(pb + 5)*bp[2][2] + *(pb + 41)*bp[4][2] + *(pb + 32)*bp[6][2]
				+ *(pb + 95)*bp[3][2] + (*(pb + 23)*bp[3][2])*cd26b*swc[4]
				+ (*(pb + 35)*bp[3][2])*cd19b*swc[5])*s2stl*f8b
				+ (*(pb + 8)*bp[2][2] + *(pb + 42)*bp[4][2] + *(pb + 34)*bp[6][2]
				+ *(pb + 97)*bp[3][2] + (*(pb + 33)*bp[3][2])*cd26b*swc[4]
				+ (*(pb + 36)*bp[3][2])*cd19b*swc[5])*c2stl*f8b;
		}

		f8c = 1.0;

		if (ww[0] != 9898.0) {
			wc[0][7] = -(*(pc + 5)*bp[3][2] + *(pc + 41)*bp[5][2] + *(pc + 32)*bp[7][2]
				+ *(pc + 95)*bp[2][2] + (*(pc + 23)*bp[2][2])*cd26b*swc[4]
				+ (*(pc + 35)*bp[2][2])*cd19b*swc[5])*s2stl*f8c
				+ (*(pc + 8)*bp[3][2] + *(pc + 42)*bp[5][2] + *(pc + 34)*bp[7][2]
				+ *(pc + 97)*bp[2][2] + (*(pc + 33)*bp[2][2])*cd26b*swc[4]
				+ (*(pc + 36)*bp[2][2])*cd19b*swc[5])*c2stl*f8c;
		}

		if (ww[1] != 9898.0) {
			wc[1][7] = -(*(pc + 8)*bt[3][2] + *(pc + 42)*bt[5][2] + *(pc + 34)*bt[7][2]
				+ *(pc + 97)*bt[2][2] + (*(pc + 33)*bt[2][2])*cd26b*swc[4]
				+ (*(pc + 36)*bt[2][2])*cd19b*swc[5])*s2stl*f8c
				- (*(pc + 5)*bt[3][2] + *(pc + 41)*bt[5][2] + *(pc + 32)*bt[7][2]
				+ *(pc + 95)*bt[2][2] + (*(pc + 23)*bt[2][2])*cd26b*swc[4]
				+ (*(pc + 35)*bt[2][2])*cd19b*swc[5])*c2stl*f8c;
		}
	}

	// Longitudinal
	if ((sw[9] != 0.0) || (sw[10] != 0.0)) {

		if (ww[0] != 9898.0) {
			wc[0][10] = -(*(pc + 64)*bp[1][1] + *(pc + 65)*bp[3][1] + *(pc + 66)*bp[5][1]
				+ *(pc + 74)*bp[2][1] + *(pc + 75)*bp[4][1] + *(pc + 76)*bp[6][1]
				+ (*(pc + 56)*bp[1][1] + *(pc + 58)*bp[3][1] + *(pc + 61)*bp[5][1]
				+ *(pc + 50)*bp[2][1] + *(pc + 52)*bp[4][1] + *(pc + 54)*bp[6][1])*cd64c*swc[2]
				+ (*(pc + 73)*bp[1][1] + *(pc + 81)*bp[3][1] + *(pc + 84)*bp[5][1]
				+ *(pc + 67)*bp[2][1] + *(pc + 69)*bp[4][1]
				+ *(pc + 71)*bp[6][1])*cd87c*swc[3])*slong
				+ (*(pc + 90)*bp[1][1] + *(pc + 91)*bp[3][1] + *(pc + 92)*bp[5][1]
				+ *(pc + 77)*bp[2][1] + *(pc + 78)*bp[4][1] + *(pc + 79)*bp[6][1]
				+ (*(pc + 57)*bp[1][1] + *(pc + 60)*bp[3][1] + *(pc + 62)*bp[5][1]
				+ *(pc + 51)*bp[2][1] + *(pc + 53)*bp[4][1] + *(pc + 55)*bp[6][1])*cd64c*swc[2]
				+ (*(pc + 80)*bp[1][1] + *(pc + 83)*bp[3][1] + *(pc + 85)*bp[5][1]
				+ *(pc + 68)*bp[2][1] + *(pc + 70)*bp[4][1]
				+ *(pc + 72)*bp[6][1])*cd87c*swc[3])*clong;
		}

		if (ww[1] != 9898.0) {
			wc[1][10] = -(*(pc + 90)*bt[1][1] + *(pc + 91)*bt[3][1] + *(pc + 92)*bt[5][1]
				+ *(pc + 77)*bt[2][1] + *(pc + 78)*bt[4][1] + *(pc + 79)*bt[6][1]
				+ (*(pc + 57)*bt[1][1] + *(pc + 60)*bt[3][1] + *(pc + 62)*bt[5][1]
				+ *(pc + 51)*bt[2][1] + *(pc + 53)*bt[4][1] + *(pc + 55)*bt[6][1])*cd64c*swc[2]
				+ (*(pc + 80)*bt[1][1] + *(pc + 83)*bt[3][1] + *(pc + 85)*bt[5][1]
				+ *(pc + 68)*bt[2][1] + *(pc + 70)*bt[4][1]
				+ *(pc + 72)*bt[6][1])*cd87c*swc[3])*slong
				- (*(pc + 64)*bt[1][1] + *(pc + 65)*bt[3][1] + *(pc + 66)*bt[5][1]
				+ *(pc + 74)*bt[2][1] + *(pc + 75)*bt[4][1] + *(pc + 76)*bt[6][1]
				+ (*(pc + 56)*bt[1][1] + *(pc + 58)*bt[3][1] + *(pc + 61)*bt[5][1]
				+ *(pc + 50)*bt[2][1] + *(pc + 52)*bt[4][1] + *(pc + 54)*bt[6][1])*cd64c*swc[2]
				+ (*(pc + 73)*bt[1][1] + *(pc + 81)*bt[3][1] + *(pc + 84)*bt[5][1]
				+ *(pc + 67)*bt[2][1] + *(pc + 69)*bt[4][1]
				+ *(pc + 71)*bt[6][1])*cd87c*swc[3])*clong;
		}
	}

	wbt[0] = 0.0;
	wbt[1] = 0.0;
	wct[0] = 0.0;
	wct[1] = 0.0;

	//sum winds and change meridional sign to + North
	for (k = 0; k < nsw; ++k) {
		wbt[0] = wbt[0] - fabs(sw[k])*wb[0][k];
		wct[0] = wct[0] - fabs(sw[k])*wc[0][k];
		wbt[1] = wbt[1] + fabs(sw[k])*wb[1][k];
		wct[1] = wct[1] + fabs(sw[k])*wc[1][k];
	}

	if (ww[0] != 9898.0) {
		ww[0] = wbt[0] * sw[23] + wct[0] * sw[24];
	}

	if (ww[1] != 9898.0) {
		ww[1] = wbt[1] * sw[23] + wct[1] * sw[24];
	}

	return;

}


void HWM::glbw5e(double yrd, double sec, double lat, double lon, double stl, double f107a, double f107, double *ap, double *pb, double *pc, double *ww)
{

	//glbw5e member function from HWM class

	const double sr = 7.2722e-05, hr = 0.2618, dr = 1.72142e-02, pset = 3.0;
	const int nsw = 14, lv = 12, mv = 3, nsv = 3;
	double pb14 = -1.0, pb18 = -1.0, sw9 = 1.0;
	double wb[2][15] = { { 0.0 } }, wc[2][15] = { { 0.0 } };
	int j, k, iyr;
	double day, dayl = 0.0, cd14, cd18, df, dfa, dfc, f1b, f1c, f5b, f5c;
	double f7b, f75b, f7c, f75c, f8b, f8c, f14b, f14c, apd, apdf, apdfc;
	double apt, dbasy1, dbasy2, f11b, dcasy1, dcasy2, f11c, utbasy, f12b, utcasy, f12c;
	double wbt[2], wct[2];


	//Confirm parameter set
	if (pb[99] == 0.0){
		pb[99] = pset;
	}

	if (pb[99] != pset){
		return;
	}

	for (j = 0; j < nsw; j++){
		wb[0][j] = 0.0;
		wb[1][j] = 0.0;
		wc[0][j] = 0.0;
		wc[1][j] = 0.0;
	}

	if (sw[8] > 0.0){
		sw9 = 1.0;
	}

	if (sw[8] < 0.0){
		sw9 = -1.0;
	}

	iyr = int(yrd) / 1000;
	day = yrd - iyr*1000.0;

	if ((xvl != lat) || (lv > lvl) || (mv > mvl)){
		slat = sin(pi180*lat);
		clat = cos(pi180*lat);

		vsphr1(slat, clat, lv, mv, bt, bp, 20);

		xvl = lat;
		lvl = lv;
		mvl = mv;

	}

	if ((tll != stl) || (nsv > nsvl)) {
		sstl = sin(hr*stl);
		cstl = cos(hr*stl);
		s2stl = sin(2.0*hr*stl);
		c2stl = cos(2.0*hr*stl);
		s3stl = sin(3.0*hr*stl);
		c3stl = cos(3.0*hr*stl);
		tll = stl;
		nsvl = nsv;
	}

	if ((day != dayl) || (pb[13] != pb14)) {
		cd14 = cos(dr*(day - pb[13]));
		//sd14 = sin(dr*(day - pb[13]));
	}

	if ((day != dayl) || (pb[17] != pb18)) {
		cd18 = cos(2.0*dr*(day - pb[17]));
	}

	dayl = day;
	pb14 = pb[13];
	pb18 = pb[17];

	if (xll != lon) {
		slong = sin(pi180*lon);
		clong = cos(pi180*lon);
		s2long = sin(2.0*pi180*lon);
		c2long = cos(2.0*pi180*lon);
		xll = lon;
		ngvl = 2;
	}

	// f10.7 effect
	df = f107 - f107a;
	dfa = f107a - 150.0;
	dfc = dfa + pb[19] * df;

	// Time independent
	f1b = 1.0 + pb[21] * dfc*swc[0];
	if (ww[0] != 9898.0) {
		wb[0][1] = (pb[1] * bt[2][0] + pb[2] * bt[4][0] + pb[22] * bt[6][0])*f1b;
	}

	wb[1][1] = 0.0;
	wc[0][1] = 0.0;

	f1c = 1.0 + pc[21] * dfc*swc[0];

	if (ww[1] != 9898.0) {
		wc[1][1] = -(pc[1] * bt[1][0] + pc[2] * bt[3][0] + pc[22] * bt[5][0])*f1c
			- (pc[26] * bt[2][0] + pc[14] * bt[4][0] + pc[59] * bt[6][0]
			+ pc[160] * bt[8][0] + pc[161] * bt[10][0] + pc[162] * bt[12][0])*f1c;
	}

	// symmetrical annual
	// symmetrical semiannual
	if (ww[0] != 9898.0) {
		wb[0][3] = (pb[16] * bt[2][0] + pb[30] * bt[4][0])*cd18;
	}

	wb[1][3] = 0.0;
	wc[0][3] = 0.0;

	if (ww[1] != 9898.0) {
		wc[1][3] = -(pc[16] * bt[1][0] + pc[30] * bt[3][0])*cd18;
	}

	// Asymmetrical annual
	f5b = 1.0 + pb[47] * dfc*swc[0];

	if (ww[0] != 9898.0) {
		wb[0][4] = (pb[9] * bt[1][0] + pb[10] * bt[3][0])*cd14*f5b;
	}

	wb[1][4] = 0.0;
	wc[0][4] = 0.0;

	f5c = 1.0 + pc[47] * dfc*swc[0];

	if (ww[1] != 9898.0) {
		wc[1][4] = -(pc[9] * bt[2][0] + pc[10] * bt[4][0])*cd14*f5c;
	}

	// Asymmetrical semiannual
	//     none
	// Diurnal
	if (sw[6] != 0.0) {
		f7b = 1.0 + pb[49] * dfc*swc[0];
		f75b = 1.0 + pb[82] * dfc*swc[0];

		if (ww[0] != 9898.0) {
			wb[0][6] = (pb[6] * bt[1][1] + pb[7] * bt[3][1] + pb[28] * bt[5][1] + pb[141] * bt[7][1]
				+ pb[143] * bt[9][1] + pb[181] * bt[2][1] + pb[183] * bt[4][1])*sstl*f7b
				+ (pb[12] * bt[2][1] + pb[145] * bt[4][1])*cd14*sstl*f75b*swc[4]
				+ (pb[170] * bt[1][1] + pb[172] * bt[3][1])*cd18*sstl*f75b*swc[3]
				+ (pb[3] * bt[1][1] + pb[4] * bt[3][1] + pb[27] * bt[5][1] + pb[140] * bt[7][1]
				+ pb[142] * bt[9][1] + pb[180] * bt[2][1] + pb[182] * bt[4][1])*cstl*f7b
				+ (pb[11] * bt[2][1] + pb[144] * bt[4][1])*cd14*cstl*f75b*swc[4]
				+ (pb[169] * bt[1][1] + pb[171] * bt[3][1])*cd18*cstl*f75b*swc[3];
		}

		if (ww[1] != 9898.0) {
			wb[1][6] = -(pb[3] * bp[1][1] + pb[4] * bp[3][1] + pb[27] * bp[5][1] + pb[140] * bp[7][1]
				+ pb[142] * bp[9][1] + pb[180] * bp[2][1] + pb[182] * bp[4][1])*sstl*f7b
				- (pb[11] * bp[2][1] + pb[144] * bp[4][1])*cd14*sstl*f75b*swc[4]
				- (pb[169] * bp[1][1] + pb[171] * bp[3][1])*cd18*sstl*f75b*swc[3]
				+ (pb[6] * bp[1][1] + pb[7] * bp[3][1] + pb[28] * bp[5][1] + pb[141] * bp[7][1]
				+ pb[143] * bp[9][1] + pb[181] * bp[2][1] + pb[183] * bp[4][1])*cstl*f7b
				+ (pb[12] * bp[2][1] + pb[145] * bp[4][1])*cd14*cstl*f75b*swc[4]
				+ (pb[170] * bp[1][1] + pb[172] * bp[3][1])*cd18*cstl*f75b*swc[3];
		}

		f7c = 1.0 + pc[49] * dfc*swc[0];
		f75c = 1.0 + pc[82] * dfc*swc[0];

		if (ww[0] != 9898.0) {
			wc[0][6] = -(pc[3] * bp[2][1] + pc[4] * bp[4][1] + pc[27] * bp[6][1] + pc[140] * bp[8][1]
				+ pc[142] * bp[10][1] + pc[180] * bp[1][1] + pc[182] * bp[3][1]
				+ pc[184] * bp[5][1] + pc[186] * bp[7][1] + pc[188] * bp[9][1])*sstl*f7c
				- (pc[11] * bp[1][1] + pc[144] * bp[3][1])*cd14*sstl*f75c*swc[4]
				- (pc[169] * bp[2][1] + pc[171] * bp[4][1])*cd18*sstl*f75c*swc[3]
				+ (pc[6] * bp[2][1] + pc[7] * bp[4][1] + pc[28] * bp[6][1] + pc[141] * bp[8][1]
				+ pc[143] * bp[10][1] + pc[181] * bp[1][1] + pc[183] * bp[3][1] + pc[185] * bp[5][1]
				+ pc[187] * bp[7][1] + pc[189] * bp[9][1])*cstl*f7c
				+ (pc[12] * bp[1][1] + pc[145] * bp[3][1])*cd14*cstl*f75c*swc[4]
				+ (pc[170] * bp[2][1] + pc[172] * bp[4][1])*cd18*cstl*f75c*swc[3];
		}

		if (ww[1] != 9898.0) {
			wc[1][6] = -(pc[6] * bt[2][1] + pc[7] * bt[4][1] + pc[28] * bt[6][1] + pc[141] * bt[8][1]
				+ pc[143] * bt[10][1] + pc[181] * bt[1][1] + pc[183] * bt[3][1]
				+ pc[185] * bt[5][1] + pc[187] * bt[7][1] + pc[189] * bt[9][1])*sstl*f7c
				- (pc[12] * bt[1][1] + pc[145] * bt[3][1])*cd14*sstl*f75c*swc[4]
				- (pc[170] * bt[2][1] + pc[172] * bt[4][1])*cd18*sstl*f75c*swc[3]
				- (pc[3] * bt[2][1] + pc[4] * bt[4][1] + pc[27] * bt[6][1] + pc[140] * bt[8][1]
				+ pc[142] * bt[10][1] + pc[180] * bt[1][1] + pc[182] * bt[3][1] + pc[184] * bt[5][1]
				+ pc[186] * bt[7][1] + pc[188] * bt[9][1])*cstl*f7c
				- (pc[11] * bt[1][1] + pc[144] * bt[3][1])*cd14*cstl*f75c*swc[4]
				- (pc[169] * bt[2][1] + pc[171] * bt[4][1])*cd18*cstl*f75c*swc[3];
		}
	}

	// semidiurnal
	if (sw[7] != 0.0) {
		f8b = 1.0 + pb[89] * dfc*swc[0];

		if (ww[0] != 9898.0) {
			wb[0][7] = (pb[8] * bt[2][2] + pb[42] * bt[4][2] + pb[110] * bt[6][2]
				+ (pb[33] * bt[3][2] + pb[147] * bt[5][2])*cd14*swc[4]
				+ (pb[133] * bt[2][2])*cd18*swc[3] + pb[151] * bt[3][2]
				+ pb[153] * bt[5][2] + pb[155] * bt[7][2] + pb[157] * bt[9][2])*s2stl*f8b
				+ (pb[5] * bt[2][2] + pb[41] * bt[4][2] + pb[109] * bt[6][2]
				+ (pb[23] * bt[3][2] + pb[146] * bt[5][2])*cd14*swc[4]
				+ (pb[134] * bt[2][2])*cd18*swc[3] + pb[150] * bt[3][2]
				+ pb[152] * bt[5][2] + pb[154] * bt[7][2] + pb[156] * bt[9][2])*c2stl*f8b;
		}

		if (ww[1] != 9898.0) {
			wb[1][7] = -(pb[5] * bp[2][2] + pb[41] * bp[4][2] + pb[109] * bp[6][2]
				+ (pb[23] * bp[3][2] + pb[146] * bp[5][2])*cd14*swc[4]
				+ (pb[134] * bp[2][2])*cd18*swc[3] + pb[150] * bp[3][2]
				+ pb[152] * bp[5][2] + pb[154] * bp[7][2] + pb[156] * bp[9][2])*s2stl*f8b
				+ (pb[8] * bp[2][2] + pb[42] * bp[4][2] + pb[110] * bp[6][2]
				+ (pb[33] * bp[3][2] + pb[147] * bp[5][2])*cd14*swc[4]
				+ (pb[133] * bp[2][2])*cd18*swc[3] + pb[151] * bp[3][2]
				+ pb[153] * bp[5][2] + pb[155] * bp[7][2] + pb[157] * bp[9][2])*c2stl*f8b;
		}

		f8c = 1.0 + pc[89] * dfc*swc[0];

		if (ww[0] != 9898.0) {
			wc[0][7] = -(pc[5] * bp[3][2] + pc[41] * bp[5][2] + pc[109] * bp[7][2]
				+ (pc[23] * bp[2][2] + pc[146] * bp[4][2])*cd14*swc[4]
				+ (pc[134] * bp[3][2])*cd18*swc[3] + pc[150] * bp[2][2]
				+ pc[152] * bp[4][2] + pc[154] * bp[6][2] + pc[156] * bp[8][2])*s2stl*f8c
				+ (pc[8] * bp[3][2] + pc[42] * bp[5][2] + pc[110] * bp[7][2]
				+ (pc[33] * bp[2][2] + pc[147] * bp[4][2])*cd14*swc[4]
				+ (pc[133] * bp[3][2])*cd18*swc[3] + pc[151] * bp[2][2]
				+ pc[153] * bp[4][2] + pc[155] * bp[6][2] + pc[157] * bp[8][2])*c2stl*f8c;
		}

		if (ww[1] != 9898.0) {
			wc[1][7] = -(pc[8] * bt[3][2] + pc[42] * bt[5][2] + pc[110] * bt[7][2]
				+ (pc[33] * bt[2][2] + pc[147] * bt[4][2])*cd14*swc[4]
				+ (pc[133] * bt[3][2])*cd18*swc[3] + pc[151] * bt[2][2]
				+ pc[153] * bt[4][2] + pc[155] * bt[6][2] + pc[157] * bt[8][2])*s2stl*f8c
				- (pc[5] * bt[3][2] + pc[41] * bt[5][2] + pc[109] * bt[7][2]
				+ (pc[23] * bt[2][2] + pc[146] * bt[4][2])*cd14*swc[4]
				+ (pc[134] * bt[3][2])*cd18*swc[3] + pc[150] * bt[2][2]
				+ pc[152] * bt[4][2] + pc[154] * bt[6][2] + pc[156] * bt[8][2])*c2stl*f8c;
		}
	}

	// Terdiurnal
	if (sw[13] != 0.0) {
		f14b = 1.0;

		if (ww[0] != 9898.0) {
			wb[0][13] = (pb[39] * bt[3][3] + pb[148] * bt[5][3] + pb[113] * bt[7][3]
				+ (pb[93] * bt[4][3] + pb[46] * bt[6][3])*cd14*swc[4])*s3stl*f14b
				+ (pb[40] * bt[3][3] + pb[149] * bt[5][3] + pb[114] * bt[7][3]
				+ (pb[94] * bt[4][3] + pb[48] * bt[6][3])*cd14*swc[4])*c3stl*f14b;
		}

		if (ww[1] != 9898.0) {
			wb[1][13] = -(pb[40] * bp[3][3] + pb[149] * bp[5][3] + pb[114] * bp[7][3]
				+ (pb[94] * bp[4][3] + pb[48] * bp[6][3])*cd14*swc[4])*s3stl*f14b
				+ (pb[39] * bp[3][3] + pb[148] * bp[5][3] + pb[113] * bp[7][3]
				+ (pb[93] * bp[4][3] + pb[46] * bp[6][3])*cd14*swc[4])*c3stl*f14b;
		}

		f14c = 1.0;

		if (ww[0] != 9898.0) {
			wc[0][13] = -(pc[40] * bp[4][3] + pc[149] * bp[6][3] + pc[114] * bp[8][3]
				+ (pc[94] * bp[3][3] + pc[48] * bp[5][3])*cd14*swc[4])*s3stl*f14c
				+ (pc[39] * bp[4][3] + pc[148] * bp[6][3] + pc[113] * bp[8][3]
				+ (pc[93] * bp[3][3] + pc[46] * bp[5][3])*cd14*swc[4])*c3stl*f14c;
		}

		if (ww[1] != 9898.0) {
			wc[1][13] = -(pc[39] * bt[4][3] + pc[148] * bt[6][3] + pc[113] * bt[8][3]
				+ (pc[93] * bt[3][3] + pc[46] * bt[5][3])*cd14*swc[4])*s3stl*f14c
				- (pc[40] * bt[4][3] + pc[149] * bt[6][3] + pc[114] * bt[8][3]
				+ (pc[94] * bt[3][3] + pc[48] * bt[5][3])*cd14*swc[4])*c3stl*f14c;
		}
	}

	// Magnetic activity
	if (sw[8] != 0.0) {
		if (sw9 != -1.0) {

			// Daily ap
			apd = ap[0] - 4.0;
			apdf = (apd + (pb[44] - 1.0)*(apd + (exp(-pb[43] * apd) - 1.0) / pb[43]));
			// apdfc=(apd + (pc[44] - 1.0)*(apd + (exp( - pc[43]*apd) - 1.0)/pc[43]));
			apdfc = apdf;

			if (apd != 0) {

				if (ww[0] != 9898.0) {
					wb[0][8] = (pb[45] * bt[2][0] + pb[34] * bt[4][0] + pb[32] * bt[6][0])*apdf
						+ (pb[174] * bt[2][2] + pb[176] * bt[4][2])*s2stl*apdf
						+ (pb[173] * bt[2][2] + pb[175] * bt[4][2])*c2stl*apdf;
				}

				if (ww[1] != 9898.0) {
					wb[1][8] = -(pb[173] * bp[2][2] + pb[175] * bp[4][2])*s2stl*apdf
						+ (pb[174] * bp[2][2] + pb[176] * bp[4][2])*c2stl*apdf;
				}

				if (ww[0] != 9898.0) {
					wc[0][8] = swc[6] * wc[0][6] * pc[121] * apdfc
						- (pc[173] * bp[3][2] + pc[175] * bp[5][2])*s2stl*apdfc
						+ (pc[174] * bp[3][2] + pc[176] * bp[5][2])*c2stl*apdfc;
				}

				wc[1][8] = -(pc[45] * bt[1][0] + pc[34] * bt[3][0]
					+ pc[32] * bt[5][0])*apdfc + swc[6] * wc[1][6] * pc[121] * apdfc
					- (pc[174] * bt[3][2] + pc[176] * bt[5][2])*s2stl*apdfc
					- (pc[173] * bt[3][2] + pc[175] * bt[5][2])*c2stl*apdfc;
			}
		}				if (ww[1] != 9898.0) {

		}
		else {

			if (pb[24] < 1.0e-04) {
				pb[24] = 1.0e-04;
			}

			// apt = g0(ap[1], pb);
			apt = 0.0;
			if (apt != 0.0) {

				if (ww[0] != 9898.0) {
					wb[0][8] = (pb[96] * bt[2][0] + pb[54] * bt[4][0] + pb[50] * bt[6][0])*apt
						+ (pb[159] * bt[2][2] + pb[178] * bt[4][2])*s2stl*apt
						+ (pb[158] * bt[2][2] + pb[177] * bt[4][2])*c2stl*apt;
				}

				if (ww[1] != 9898.0) {
					wb[1][8] = -(pb[158] * bp[2][2] + pb[177] * bp[4][2])*s2stl*apt
						+ (pb[159] * bp[2][2] + pb[178] * bp[4][2])*c2stl*apt;
				}

				if (ww[0] != 9898.0) {
					wc[0][8] = swc[6] * wc[0][6] * pc[128] * apt
						- (pc[158] * bp[3][2] + pc[177] * bp[5][2])*s2stl*apt
						+ (pc[159] * bp[3][2] + pc[178] * bp[5][2])*c2stl*apt;
				}

				if (ww[1] != 9898.0) {
					wc[1][8] = -(pc[96] * bt[1][0] + pc[54] * bt[3][0] + pc[50] * bt[5][0])*apt
						+ swc[6] * wc[1][6] * pc[128] * apt
						- (pc[159] * bt[3][2] + pc[178] * bt[5][2])*s2stl*apt
						- (pc[158] * bt[3][2] + pc[177] * bt[5][2])*c2stl*apt;
				}
			}
		}
	}

	if (sw[9] != 0.0) {

		// Longitudinal
		dbasy1 = 1.0 + pb[198] * slat;
		dbasy2 = 1.0 + pb[199] * slat;
		f11b = 1.0 + pb[80] * dfc*swc[0];

		if (sw[10] != 0.0) {

			if (ww[0] != 9898.0) {
				wb[0][10] = (pb[90] * bt[2][1] + pb[91] * bt[4][1] + pb[92] * bt[6][1])*slong*dbasy1*f11b
					+ (pb[64] * bt[2][1] + pb[65] * bt[4][1] + pb[66] * bt[6][1])*clong*dbasy1*f11b
					+ (pb[190] * bt[2][2] + pb[192] * bt[4][2]
					+ pb[194] * bt[6][2] + pb[196] * bt[8][2])*s2long*dbasy2*f11b
					+ (pb[191] * bt[2][2] + pb[193] * bt[4][2]
					+ pb[195] * bt[6][2] + pb[197] * bt[8][2])*c2long*dbasy2*f11b;
			}

			if (ww[1] != 9898.0) {
				wb[1][10] = -(pb[64] * bp[2][1] + pb[65] * bp[4][1] + pb[66] * bp[6][1])*slong*dbasy1*f11b
					+ (pb[90] * bp[2][1] + pb[91] * bp[4][1] + pb[92] * bp[6][1])*clong*dbasy1*f11b
					- (pb[191] * bp[2][2] + pb[193] * bp[4][2]
					+ pb[195] * bp[6][2] + pb[197] * bp[8][2])*s2long*dbasy2*f11b
					+ (pb[190] * bp[2][2] + pb[192] * bp[4][2]
					+ pb[194] * bp[6][2] + pb[196] * bp[8][2])*c2long*dbasy2*f11b;
			}

			dcasy1 = 1.0 + pc[198] * slat;
			dcasy2 = 1.0 + pc[199] * slat;
			f11c = 1.0 + pc[80] * dfc*swc[0];

			if (ww[0] != 9898.0) {
				wc[0][10] = -(pc[64] * bp[1][1] + pc[65] * bp[3][1] + pc[66] * bp[5][1] + pc[72] * bp[7][1]
					+ pc[73] * bp[9][1])*slong*dcasy1*f11c
					+ (pc[90] * bp[1][1] + pc[91] * bp[3][1] + pc[92] * bp[5][1] + pc[86] * bp[7][1]
					+ pc[87] * bp[9][1])*clong*dcasy1*f11c
					- (pc[191] * bp[3][2] + pc[193] * bp[5][2] + pc[195] * bp[7][2]
					+ pc[197] * bp[9][2])*s2long*dcasy2*f11c
					+ (pc[190] * bp[3][2] + pc[192] * bp[5][2] + pc[194] * bp[7][2]
					+ pc[196] * bp[9][2])*c2long*dcasy2*f11c;
			}

			if (ww[1] != 9898.0) {
				wc[1][10] = -(pc[90] * bt[1][1] + pc[91] * bt[3][1] + pc[92] * bt[5][1] + pc[86] * bt[7][1]
					+ pc[87] * bt[9][1])*slong*dcasy1*f11c
					- (pc[64] * bt[1][1] + pc[65] * bt[3][1] + pc[66] * bt[5][1] + pc[72] * bt[7][1]
					+ pc[73] * bt[9][1])*clong*dcasy1*f11c
					- (pc[190] * bt[3][2] + pc[192] * bt[5][2] + pc[194] * bt[7][2]
					+ pc[196] * bt[9][2])*s2long*dcasy2*f11c
					- (pc[191] * bt[3][2] + pc[193] * bt[5][2] + pc[195] * bt[7][2]
					+ pc[197] * bt[9][2])*c2long*dcasy2*f11c;
			}
		}

		// Ut & mixed ut/lon
		utbasy = 1.0;
		f12b = 1.0 + pb[81] * dfc*swc[0];

		if (sw[11] != 0.0) {

			if (ww[0] != 9898.0) {
				wb[0][11] = (pb[68] * bt[1][0] + pb[69] * bt[3][0] + pb[70] * bt[5][0]
					+ pb[115] * bt[7][0] + pb[116] * bt[9][0]
					+ pb[117] * bt[11][0])*cos(sr*(sec - pb[71]))*utbasy*f12b
					+ (pb[76] * bt[3][2] + pb[77] * bt[5][2]
					+ pb[78] * bt[7][2])*cos(sr*(sec - pb[79])
					+ 2.0*pi180*lon)*utbasy*f12b*swc[10];
			}

			if (ww[1] != 9898.0) {
				wb[1][11] = (pb[76] * bp[3][2] + pb[77] * bp[5][2]
					+ pb[78] * bp[7][2])*cos(sr*(sec - pb[79] + 21600.0)
					+ 2.0*pi180*lon)*utbasy*f12b*swc[10];
			}

			utcasy = 1.0;
			f12c = 1.0 + pc[81] * dfc*swc[0];

			if (ww[0] != 9898.0) {
				wc[0][11] = (pc[76] * bp[2][2] + pc[77] * bp[4][2] + pc[78] * bp[6][2] + pc[164] * bp[8][2]
					+ pc[165] * bp[10][2] + pc[166] * bp[12][2])*cos(sr*(sec - pc[79])
					+ 2.0*pi180*lon)*utcasy*f12c*swc[10];
			}

			if (ww[1] != 9898.0) {
				wc[1][11] = -(pc[68] * bt[2][0] + pc[69] * bt[4][0] + pc[70] * bt[6][0]
					+ pc[115] * bt[8][0] + pc[116] * bt[10][0]
					+ pc[117] * bt[12][0])*cos(sr*(sec - pc[71]))*utcasy*f12c
					+ (pc[76] * bt[2][2] + pc[77] * bt[4][2] + pc[78] * bt[6][2]
					+ pc[164] * bt[8][2] + pc[165] * bt[10][2]
					+ pc[166] * bt[12][2])*cos(sr*(sec - pc[79] + 21600.0)
					+ 2.0*pi180*lon)*utcasy*f12c*swc[10];
			}
		}

		// Mixed lon, ut, ap
		if ((sw[12] != 0.0) || (apd != 0.0)) {

			if (sw9 != -1.0) {

				if (ww[0] != 9898.0) {
					wb[0][12] = (pb[60] * bt[2][1] + pb[61] * bt[4][1]
						+ pb[62] * bt[6][1])*cos(pi180*(lon - pb[63]))*apdf*swc[10]
						+ (pb[83] * bt[1][0] + pb[84] * bt[3][0]
						+ pb[85] * bt[5][0])*cos(sr*(sec - pb[75]))*apdf*swc[11];
				}

				if (ww[1] != 9898.0) {
					wb[1][12] = (pb[60] * bp[2][1] + pb[61] * bp[4][1]
						+ pb[62] * bp[6][1])*cos(pi180*(lon - pb[63] + 90.))*apdf*swc[10];
				}

				if (ww[0] != 9898.0) {
					wc[0][12] = swc[10] * wc[0][10] * pc[60] * apdfc + swc[11] * wc[0][11] * pc[83] * apdfc;
				}

				if (ww[1] != 9898.0) {
					wc[1][12] = swc[10] * wc[1][10] * pc[60] * apdfc + swc[11] * wc[1][11] * pc[83] * apdfc;
				}
			}
			else {

				if (apt != 0.0) {

					if (ww[0] != 9898.0) {
						wb[0][12] = (pb[52] * bt[2][1] + pb[98] * bt[4][1]
							+ pb[67] * bt[6][1])*cos(pi180*(lon - pb[97]))*apt*swc[10]
							+ (pb[55] * bt[1][0] + pb[56] * bt[3][0]
							+ pb[57] * bt[5][0])*cos(sr*(sec - pb[58]))*apt*swc[11];
					}

					if (ww[1] != 9898.0) {
						wb[1][12] = (pb[52] * bp[2][1] + pb[98] * bp[4][1]
							+ pb[67] * bp[6][1])*cos(pi180*(lon - pb[97] + 90.0))*apt*swc[10];
					}

					if (ww[0] != 9898.0) {
						wc[0][12] = swc[10] * wc[0][10] * pc[52] * apt + swc[11] * wc[0][11] * pc[55] * apt;
					}

					if (ww[1] != 9898.0) {
						wc[1][12] = swc[10] * wc[1][10] * pc[52] * apt + swc[11] * wc[1][11] * pc[55] * apt;
					}
				}
			}
		}
	}

	wbt[0] = 0.0;
	wbt[1] = 0.0;
	wct[0] = 0.0;
	wct[1] = 0.0;

	//sum winds and change meridional sign to + North
	for (k = 0; k < nsw; ++k) {
		wbt[0] = wbt[0] - fabs(sw[k])*wb[0][k];
		wct[0] = wct[0] - fabs(sw[k])*wc[0][k];
		wbt[1] = wbt[1] + fabs(sw[k])*wb[1][k];
		wct[1] = wct[1] + fabs(sw[k])*wc[1][k];
	}

	if (ww[0] != 9898.0) {
		ww[0] = wbt[0] * sw[23] + wct[0] * sw[24];
	}

	if (ww[1] != 9898.0) {
		ww[1] = wbt[1] * sw[23] + wct[1] * sw[24];
	}

	return;
}

void HWM::glbw5m(double yrd, double sec, double lat, double lon, double stl, double f107a, double f107, double *ap, double *pb, double *pc, double *ww)
{

	//glbw5m member function from HWM class

	const double hr = 0.2618, dr = 1.72142e-02, pset = 4.0;
	double pb14 = -1.0, pb18 = -1.0, sw9 = 1.0, day, dayl = 0.0, cd14, cd18, df, dfa, dfc, f1b, f1c, f5b, f5c, f7b, f75b;
	double f7c, f75c, f8b, f8c, apd, apdf, apdfc, apt, pb19 = 0.0, cd19b, wbt[2], wct[2], wb[2][15] = { { 0.0 } }, wc[2][15] = { { 0.0 } };
	const int nsw = 14, lv = 10, mv = 2, nsv = 2;
	int j, k, iyr;

	// Confirm parameter set
	if (pb[99] == 0.0) {
		pb[99] = pset;
	}

	if (pb[99] != pset){
		return;
	}

	for (j = 0; j < nsw; ++j) {
		wb[0][j] = 0.0;
		wb[1][j] = 0.0;
		wc[0][j] = 0.0;
		wc[1][j] = 0.0;
	}

	if (sw[8] > 0.0) {
		sw9 = 1.0;
	}

	if (sw[8] < 0.0) {
		sw9 = -1.0;
	}

	iyr = int(yrd) / 1000;
	day = yrd - iyr*1000.0;

	if ((xvl != lat) || (lv > lvl) || (mv > mvl)) {
		slat = sin(pi180*lat);
		clat = cos(pi180*lat);

		vsphr1(slat, clat, lv, mv, bt, bp, 20);

		xvl = lat;
		lvl = lv;
		mvl = mv;
	}

	if ((tll != stl) || (nsv > nsvl)) {
		sstl = sin(hr*stl);
		cstl = cos(hr*stl);
		s2stl = sin(2.0*hr*stl);
		c2stl = cos(2.0*hr*stl);
		tll = stl;
		nsvl = nsv;
	}

	if ((day != dayl) || (pb[13] != pb14)) {
		cd14 = cos(dr*(day - pb[13]));
	}

	if ((day != dayl) || (pb[17] != pb18)) {
		cd18 = cos(2.0*dr*(day - pb[17]));
	}

	if ((day != dayl) || (pb[18] != pb19)) {
		cd19b = cos(2.0*dr*(day - pb[18]));
	}

	dayl = day;
	pb14 = pb[13];
	pb18 = pb[17];
	pb19 = pb[18];

	// f10.7 effect
	df = f107 - f107a;
	dfa = f107a - 150.0;
	dfc = dfa + pb[19] * df;

	// Time independent
	f1b = 1.0;

	if (ww[0] != 9898.0) {
		wb[0][1] = (pb[1] * bt[2][0] + pb[2] * bt[4][0] + pb[22] * bt[6][0])*f1b;
	}

	wb[1][1] = 0.0;
	wc[0][1] = 0.0;

	f1c = 1.0;

	if (ww[1] != 9898.0) {
		wc[1][1] = -(pc[1] * bt[1][0] + pc[2] * bt[3][0] + pc[22] * bt[5][0])*f1c
			- (pc[26] * bt[2][0] + pc[14] * bt[4][0] + pc[59] * bt[6][0])*f1c;
	}

	// symmetrical annual
	// symmetrical semiannual
	if (ww[0] != 9898.0) {
		wb[0][3] = (pb[16] * bt[2][0] + pb[30] * bt[4][0])*cd18;
	}

	wb[1][3] = 0.0;
	wc[0][3] = 0.0;

	if (ww[1] != 9898.0) {
		wc[1][3] = -(pc[16] * bt[1][0] + pc[30] * bt[3][0])*cd18;
	}

	// Asymmetrical annual
	f5b = 1.0;

	if (ww[0] != 9898.0) {
		wb[0][4] = (pb[9] * bt[1][0] + pb[10] * bt[3][0])*cd14*f5b;
	}

	wb[1][4] = 0.0;
	wc[0][4] = 0.0;

	f5c = 1.0;

	if (ww[1] != 9898.0) {
		wc[1][4] = -(pc[9] * bt[2][0] + pc[10] * bt[4][0])*cd14*f5c;
	}

	// Asymmetrical semiannual
	// Diurnal
	if (sw[6] != 0.0) {
		f7b = 1.0;
		f75b = 1.0;

		if (ww[0] != 9898.0) {
			wb[0][6] = (pb[6] * bt[1][1] + pb[7] * bt[3][1]
				+ pb[28] * bt[5][1] + pb[88] * bt[2][1])*sstl*f7b
				+ (pb[12] * bt[2][1] + pb[145] * bt[4][1])*cd14*sstl*f75b*swc[4]
				+ (pb[3] * bt[1][1] + pb[4] * bt[3][1]
				+ pb[27] * bt[5][1] + pb[87] * bt[2][1])*cstl*f7b
				+ (pb[11] * bt[2][1] + pb[144] * bt[4][1])*cd14*cstl*f75b*swc[4];
		}

		if (ww[1] != 9898.0) {
			wb[1][6] = -(pb[3] * bp[1][1] + pb[4] * bp[3][1]
				+ pb[27] * bp[5][1] + pb[87] * bp[2][1])*sstl*f7b
				- (pb[11] * bp[2][1] + pb[144] * bp[4][1])*cd14*sstl*f75b*swc[4]
				+ (pb[6] * bp[1][1] + pb[7] * bp[3][1]
				+ pb[28] * bp[5][1] + pb[88] * bp[2][1])*cstl*f7b
				+ (pb[12] * bp[2][1] + pb[145] * bp[4][1])*cd14*cstl*f75b*swc[4];
		}

		f7c = 1.0;
		f75c = 1.0;

		if (ww[0] != 9898.0) {
			wc[0][6] = -(pc[3] * bp[2][1] + pc[4] * bp[4][1] + pc[27] * bp[6][1]
				+ pc[87] * bp[1][1] + pc[140] * bp[8][1] + pc[142] * bp[10][1])*sstl*f7c
				- (pc[11] * bp[1][1] + pc[144] * bp[3][1])*cd14*sstl*f75c*swc[4]
				+ (pc[6] * bp[2][1] + pc[7] * bp[4][1] + pc[28] * bp[6][1]
				+ pc[88] * bp[1][1] + pc[141] * bp[8][1] + pc[143] * bp[10][1])*cstl*f7c
				+ (pc[12] * bp[1][1] + pc[145] * bp[3][1])*cd14*cstl*f75c*swc[4];
		}

		if (ww[1] != 9898.0) {
			wc[1][6] = -(pc[6] * bt[2][1] + pc[7] * bt[4][1] + pc[28] * bt[6][1]
				+ pc[88] * bt[1][1] + pc[141] * bt[8][1] + pc[143] * bt[10][1])*sstl*f7c
				- (pc[12] * bt[1][1] + pc[145] * bt[3][1])*cd14*sstl*f75c*swc[4]
				- (pc[3] * bt[2][1] + pc[4] * bt[4][1] + pc[27] * bt[6][1]
				+ pc[87] * bt[1][1] + pc[140] * bt[8][1] + pc[142] * bt[10][1])*cstl*f7c
				- (pc[11] * bt[1][1] + pc[144] * bt[3][1])*cd14*cstl*f75c*swc[4];
		}
	}

	// semidiurnal
	if (sw[7] != 0.0) {
		f8b = 1.0 + pb[89] * dfc*swc[0];

		if (ww[0] != 9898.0) {
			wb[0][7] = (pb[8] * bt[2][2] + pb[42] * bt[4][2] + pb[110] * bt[6][2] + pb[97] * bt[3][2]
				+ (pb[33] * bt[3][2] + pb[147] * bt[5][2])*cd14*swc[4]
				+ (pb[36] * bt[3][2])*cd19b*swc[5])*s2stl*f8b
				+ (pb[5] * bt[2][2] + pb[41] * bt[4][2] + pb[109] * bt[6][2] + pb[95] * bt[3][2]
				+ (pb[23] * bt[3][2] + pb[146] * bt[5][2])*cd14*swc[4]
				+ (pb[35] * bt[3][2])*cd19b*swc[5])*c2stl*f8b;
		}

		if (ww[1] != 9898.0) {
			wb[1][7] = -(pb[5] * bp[2][2] + pb[41] * bp[4][2] + pb[109] * bp[6][2] + pb[95] * bp[3][2]
				+ (pb[23] * bp[3][2] + pb[146] * bp[5][2])*cd14*swc[4]
				+ (pb[35] * bp[3][2])*cd19b*swc[5])*s2stl*f8b
				+ (pb[8] * bp[2][2] + pb[42] * bp[4][2] + pb[110] * bp[6][2] + pb[97] * bp[3][2]
				+ (pb[33] * bp[3][2] + pb[147] * bp[5][2])*cd14*swc[4]
				+ (pb[36] * bp[3][2])*cd19b*swc[5])*c2stl*f8b;
		}

		f8c = 1.0 + pc[89] * dfc*swc[0];

		if (ww[0] != 9898.0) {
			wc[0][7] = -(pc[5] * bp[3][2] + pc[41] * bp[5][2] + pc[109] * bp[7][2] + pc[95] * bp[2][2]
				+ (pc[23] * bp[2][2] + pc[146] * bp[4][2])*cd14*swc[4]
				+ (pc[35] * bp[2][2])*cd19b*swc[5])*s2stl*f8c
				+ (pc[8] * bp[3][2] + pc[42] * bp[5][2] + pc[110] * bp[7][2] + pc[97] * bp[2][2]
				+ (pc[33] * bp[2][2] + pc[147] * bp[4][2])*cd14*swc[4]
				+ (pc[36] * bp[2][2])*cd19b*swc[5])*c2stl*f8c;
		}

		if (ww[1] != 9898.0) {
			wc[1][7] = -(pc[8] * bt[3][2] + pc[42] * bt[5][2] + pc[110] * bt[7][2] + pc[97] * bt[2][2]
				+ (pc[33] * bt[2][2] + pc[147] * bt[4][2])*cd14*swc[4]
				+ (pc[36] * bt[2][2])*cd19b*swc[5])*s2stl*f8c
				- (pc[5] * bt[3][2] + pc[41] * bt[5][2] + pc[95] * bt[2][2] + pc[109] * bt[7][2]
				+ (pc[23] * bt[2][2] + pc[146] * bt[4][2])*cd14*swc[4]
				+ (pc[35] * bt[2][2])*cd19b*swc[5])*c2stl*f8c;
		}
	}

	// Terdiurnal
	// Magnetic activity
	if (sw[8] != 0.0) {
		if (sw9 != -1.0) {

			// Daily ap
			apd = ap[0] - 4.0;
			apdf = (apd + (pb[44] - 1.0)*(apd + (exp(-pb[43] * apd) - 1.0) / pb[43]));
			// apdfc=(apd  +  (pc[44]  -  1.0)*(apd  +  (exp(  -  pc[43]*apd)  -  1.0)/pc[43]));
			apdfc = apdf;

			if (apd != 0) {

				if (ww[0] != 9898.0) {
					wb[0][8] = (pb[45] * bt[2][0] + pb[34] * bt[4][0])*apdf
						+ (pb[121] * bt[1][1] + pb[122] * bt[3][1]
						+ pb[123] * bt[5][1])*cos(hr*(stl - pb[124]))*apdf*swc[6];
				}

				if (ww[1] != 9898.0) {
					wb[1][8] = (pb[121] * bp[1][1] + pb[122] * bp[3][1]
						+ pb[123] * bp[5][1])*cos(hr*(stl - pb[124] + 6.0))*apdf*swc[6];
				}

				if (ww[0] != 9898.0) {
					wc[0][8] = (pc[121] * bp[2][1] + pc[122] * bp[4][1]
						+ pc[123] * bp[6][1])*cos(hr*(stl - pc[124]))*apdfc*swc[6];
				}

				if (ww[1] != 9898.0) {
					wc[1][8] = -(pc[45] * bt[1][0] + pc[34] * bt[3][0])*apdfc
						+ (pc[121] * bt[2][1] + pc[122] * bt[4][1]
						+ pc[123] * bt[6][1])*cos(hr*(stl - pc[124] + 6.0))*apdfc*swc[6];
				}
			}
		}
		else {

			if (pb[24] < 1.0e-04) {
				pb[24] = 1.0e-04;
			}

			//apt = g0(ap[1], pb);
			apt = 0.0;
			if (apt != 0.0) {

				if (ww[0] != 9898.0) {
					wb[0][8] = (pb[96] * bt[2][0] + pb[54] * bt[4][0])*apt
						+ (pb[128] * bt[1][1] + pb[129] * bt[3][1]
						+ pb[130] * bt[5][1])*cos(hr*(stl - pb[131]))*apt*swc[6];
				}

				if (ww[1] != 9898.0) {
					wb[1][8] = (pb[128] * bp[1][1] + pb[129] * bp[3][1]
						+ pb[130] * bp[5][1])*cos(hr*(stl - pb[131] + 6.0))*apt*swc[6];
				}

				if (ww[0] != 9898.0) {
					wc[0][8] = (pc[128] * bp[2][1] + pc[129] * bp[4][1]
						+ pc[130] * bp[6][1])*cos(hr*(stl - pc[131]))*apt*swc[6];
				}

				if (ww[1] != 9898.0) {
					wc[1][8] = -(pc[96] * bt[1][0] + pc[54] * bt[3][0])*apt
						+ (pc[128] * bt[2][1] + pc[129] * bt[4][1]
						+ pc[130] * bt[6][1])*cos(hr*(stl - pc[131] + 6.0))*apt*swc[6];
				}
			}
		}
	}

	wbt[0] = 0.0;
	wbt[1] = 0.0;
	wct[0] = 0.0;
	wct[1] = 0.0;


	//sum winds and change meridional sign to + North
	for (k = 0; k < nsw; ++k) {
		wbt[0] = wbt[0] - fabs(sw[k])*wb[0][k];
		wct[0] = wct[0] - fabs(sw[k])*wc[0][k];
		wbt[1] = wbt[1] + fabs(sw[k])*wb[1][k];
		wct[1] = wct[1] + fabs(sw[k])*wc[1][k];
	}

	if (ww[0] != 9898.0) {
		ww[0] = wbt[0] * sw[23] + wct[0] * sw[24];
	}

	if (ww[1] != 9898.0) {
		ww[1] = wbt[1] * sw[23] + wct[1] * sw[24];
	}

	return;



}

void HWM::tselec(double *sv)
{

	//tselec member function from HWM class
	//To turn on and off particular variations call tselec(sv), where
	//sv is a 25 element array containing:  0.0 for off, 1.0 for on, or
	// 2.0 for main effects off but cross terms on to get current values
	//of sw: call tretrv(sw)

	double sav[25];
	int i;

	for (i = 0; i < 25; i++){
		sav[i] = sv[i];
		sw[i] = fmod(sv[i], 2.0);

		if ((fabs(sv[i]) == 1.0) || (fabs(sv[i]) == 2.0)){
			swc[i] = 1.0;
		}
		else {
			swc[i] = 0.0;
		}
	}

	isw = 64999;

	return;

}


void HWM::gws5(int iyd, double sec, double alt, double glat, double glon, double stl, double f107a, double f107, double *ap, double *w)
{
	//gws5 member function from HWM class
	//Horizontal wind model HWM93 covering all altitude regions, A.E. Hedin (1/25/93) (4/9/93)
	//Calling argument list made similar to GTs5 subroutine for MSIS-86 density model and GWs4 
	//for thermospheric winds.

	const double s = 0.016, zl = 200.0;
	const int nnn = 3, mn1 = 5, mn2 = 14;
	int mn2s = 1, mn2m = 1;
	int i, ii, iz, mn2e, mnn;
	double zn1[5] = { 200.0, 150.0, 130.0, 115.0, 100.0 }, zn2[14] = { 100.0, 90.0, 82.5, 75.0, 67.5, 60.0, 52.5,
		45.0, 37.5, 30.0, 22.5, 15.0, 7.5, 0.0 }, sv[25] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
		1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
	double yrd, windf[2] = { 0.0 }, ww[2] = { 0.0 }, wzl[2] = { 0.0 }, wdzl[2] = { 0.0 }, un1[2][mn1] = { { 0.0 } },
		un2[2][mn2] = { { 0.0 } }, ugn1[2][2] = { { 0.0 } }, ugn2[2][2] = { { 0.0 } };

	if (isw != 64999){
		tselec(sv);
	}

	yrd = iyd;
	ww[0] = w[0];
	ww[1] = w[1];

	if (alt > zn1[mn1 - 1]){

		//Exospheric wind
		glbw5e(yrd, sec, glat, glon, stl, f107a, f107, ap, pwb, pwc, windf);
		windf[0] = sw[15] * windf[0];
		windf[1] = sw[15] * windf[1];

		//Wind at z1
		glbw5m(yrd, sec, glat, glon, stl, f107a, f107, ap, pwbl, pwcl, ww);
		wzl[0] = (pwbl[0] * windf[0] + ww[0])*sw[16] * sw[17];
		wzl[1] = (pwbl[0] * windf[1] + ww[1])*sw[16] * sw[17];
		un1[0][0] = wzl[0];
		un1[1][0] = wzl[1];

		//wind derivative at z1
		ww[0] = 0.0;
		ww[1] = 0.0;
		glbw5m(yrd, sec, glat, glon, stl, f107a, f107, ap, pwbld, pwcld, ww);
		wdzl[0] = (pwbld[0] * windf[0] + ww[0])*sw[18] * sw[17];
		wdzl[1] = (pwbld[0] * windf[1] + ww[1])*sw[18] * sw[17];
		ugn1[0][0] = wdzl[0] * s;
		ugn1[1][0] = wdzl[1] * s;

		if (alt < zl){
			//wind at zn1[1] (150)
			glbw5m(yrd, sec, glat, glon, stl, f107a, f107, ap, pb12, pc12, ww);
			un1[0][1] = (pb12[0] * windf[0] + ww[0])*sw[17];
			un1[1][1] = (pb12[0] * windf[1] + ww[1])*sw[17];

			//wind at zn1[2] (130)
			glbw5m(yrd, sec, glat, glon, stl, f107a, f107, ap, pb13, pc13, ww);
			un1[0][2] = ww[0] * sw[17];
			un1[1][2] = ww[1] * sw[17];

			//wind at zn1[3] (115)
			glbw5m(yrd, sec, glat, glon, stl, f107a, f107, ap, pb14, pc14, ww);
			un1[0][3] = ww[0] * sw[17];
			un1[1][3] = ww[1] * sw[17];

			goto L50;
		}

	}

	else {

	L50:
		mnn = max(1, min(mn2, nnn + 1));

		if (alt >= zn2[mnn - 1]){

			//wind at zn1[4] (100)
			glbw5m(yrd, sec, glat, glon, stl, f107a, f107, ap, pb15, pc15, ww);
			un1[0][4] = ww[0] * sw[17];
			un1[1][4] = ww[1] * sw[17];

			//wind derivative at zn1[4] (100)
			glbw5m(yrd, sec, glat, glon, stl, f107a, f107, ap, pb15d, pc15d, ww);
			ugn1[0][1] = ww[0] * sw[17];
			ugn1[1][1] = ww[1] * sw[17];

			if (alt < zn1[mn1 - 1]){
				ugn2[0][0] = ugn1[0][1];
				ugn2[1][0] = ugn1[1][1];
				un2[0][0] = un1[0][4];
				un2[1][0] = un1[1][4];
			}

			else {
				goto L90;
			}

		}

		else {
			ugn2[0][0] = 1.0e30;
			ugn2[1][0] = 1.0e30;
			un2[0][0] = 0.0;
			un2[1][0] = 0.0;
		}

		for (i = 1; i < mn2 + 1; i++){
			if (alt > zn2[i - 1]){
				break;
			}
		}

		iz = i;

		mn2s = max(1, min(iz - 1, iz - nnn));
		mn2e = min(mn2, max(mn2s + 1, iz - 1 + nnn));

		for (i = mn2s; i < mn2e + 1; i++){
			ii = 2 * (i - 2) + 1;

			if (i > 1){
				glbw5s(iyd, glat, glon, stl, &pwp[ii - 1][0], &pwp[ii][0], ww);
				un2[0][i - 1] = ww[0] * sw[19];
				un2[1][i - 1] = ww[1] * sw[19];
			}
		}

		mn2m = mn2e - mn2s + 1;
		ugn2[0][1] = 1.0e30;
		ugn2[1][1] = 1.0e30;


	}

L90:

	//wind at altitude
	if (w[0] != 9898.0){
		w[0] = wprof(alt, zl, s, windf[0], wzl[0], wdzl[0], mn1, zn1, &un1[0][0], &ugn1[0][0], mn2m, &zn2[mn2s - 1], &un2[0][mn2s - 1], &ugn2[0][0]);
	}

	if (w[1] != 9898.0){
		w[1] = wprof(alt, zl, s, windf[1], wzl[1], wdzl[1], mn1, zn1, &un1[1][0], &ugn1[1][0], mn2m, &zn2[mn2s - 1], &un2[1][mn2s - 1], &ugn2[1][0]);

	}

	return;
}

void HWM::spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
{
	//spline member function from HWM class
	//Calculate 2nd derivative of cubic spline interp function

	const int nmax = 100;
	int i, k;
	double u[nmax], sig, p, qn, un;

	if (yp1 > 0.99e30){
		y2[0] = 0.0;
		u[0] = 0.0;
	}
	else {
		y2[0] = -0.5;
		u[0] = (3.0 / (x[1] - x[0]))*((y[1] - y[0]) / (x[1] - x[0]) - yp1);
	}

	for (i = 1; i < n - 1; i++){
		sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
		p = sig*y2[i - 1] + 2.0;
		y2[i] = (sig - 1.0) / p;
		u[i] = (6.0*((y[i + 1] - y[i]) / (x[i + 1] - x[i])
			- (y[i] - y[i - 1]) / (x[i] - x[i - 1])) / (x[i + 1]
			- x[i - 1]) - sig*u[i - 1]) / p;

	}

	if (ypn > 0.99e30){
		qn = 0.0;
		un = 0.0;
	}
	else {
		qn = 0.5;
		un = (3.0 / (x[n - 1] - x[n - 2]))*(ypn - (y[n - 1] - y[n - 2])
			/ (x[n - 1] - x[n - 2]));
	}

	y2[n - 1] = (un - qn*u[n - 2]) / (qn*y2[n - 2] + 1.0);

	for (k = n - 2; k > -1; k--){
		y2[k] = y2[k] * y2[k + 1] + u[k];
	}

	return;


}



void HWM::splint(double *xa, double *ya, double *y2, int n, double x, double *y)
{

	//splint member function from HWM class
	//Calculate cubic spline interp value adapted from numerical recipes

	int klo, khi, k, count;
	double h, a, b;

	klo = 1;
	khi = n;

	count = 0;
	while (khi - klo > 1){
		count++;
		k = (khi + klo) / 2;
		if (x < xa[k - 1]){
			khi = k;
		}
		else {
			klo = k;
		}
		if (count > 20){
			break;
		}
	}

	h = xa[khi - 1] - xa[klo - 1];

	if (h == 0.0){
		h = 0.1;
	}
	a = (xa[khi - 1] - x) / h;
	b = (x - xa[klo - 1]) / h;
	*y = a*ya[klo - 1] + b*ya[khi - 1] + ((a*a*a - a)*y2[klo - 1] +
		(b*b*b - b)*y2[khi - 1])*h*h / 6.0;

	return;

}

//Getters
double HWM::getCswSw(int index)
{
	return sw[index];
}

int HWM::getCswIsw()
{
	return isw;
}

double HWM::getCswSwc(int index)
{
	return swc[index];
}


void HWM::initw5()
{
	// For wind model GWs

	isw = 0;

	xvl = -999.0;
	lvl = -1;
	mvl = -1;

	tll = -999.0;
	nsvl = -1;

	xll = -999.0;
	ngvl = -1;

	return;

}

void HWM::gwsbk5()
{
	// HWM93    28-Jan-93

	int i;
	int j;

	isdate[0] = "28-J";
	isdate[1] = "AN-9";
	isdate[2] = "3   ";
	istime[0] = "20:3";
	istime[1] = "5:39";
	name[0] = "HWM9";
	name[1] = "3   ";


	double pwb_init[200] =
		// WINF
	{
		0.00000e+00, -1.31640e+01, -1.52352e+01, 1.00718e+02, 3.94962e+00,
		2.19452e-01, 8.03296e+01, -1.02032e+00, -2.02149e-01, 5.67263e+01,
		0.00000e+00, -6.05459e+00, 6.68106e+00, -8.49486e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 8.39399e+01, 0.00000e+00, 9.96285e-02,
		0.00000e+00, -2.66243e-02, 0.00000e+00, -1.32373e+00, 1.39396e-02,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 3.36523e+01, -7.42795e-01, -3.89352e+00, -7.81354e-01,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 3.76631e+00, -1.22024e+00,
		-5.47580e-01, 1.09146e+00, 9.06245e-01, 2.21119e-02, 0.00000e+00,
		7.73919e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		-3.82415e-01, 0.00000e+00, 1.76202e-01, 0.00000e+00, -6.77651e-01,
		1.10357e+00, 2.25732e+00, 0.00000e+00, 1.54237e+04, 0.00000e+00,
		1.27411e-01, -2.84314e-03, 4.62562e-01, -5.34596e+01, -7.23808e+00,
		0.00000e+00, 0.00000e+00, 4.52770e-01, -8.50922e+00, -2.85389e-01,
		2.12000e+01, 6.80171e+02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		-2.72552e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.64109e+03,
		0.00000e+00, 0.00000e+00, 0.00000e+00, -1.47320e+00, -2.98179e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.05412e-02,
		4.93452e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 7.98332e-02, -5.30954e+01, 2.10211e-02, 3.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.79843e-01,
		1.81152e-01, 0.00000e+00, 0.00000e+00, -6.24673e-02, -5.37589e-02,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, -8.94418e-02, 3.70413e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -4.84645e+00,
		4.24178e-01, 0.00000e+00, 0.00000e+00, 1.86494e-01, -9.56931e-02,
		2.08426e+00, 1.53714e+00, -2.87496e-01, 4.06380e-01, -3.59788e-01,
		-1.87814e-01, 0.00000e+00, 0.00000e+00, 2.01362e-01, -1.21604e-01,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 7.86304e+00,
		2.51878e+00, 2.91455e+00, 4.32308e+00, 6.77054e-02, -2.39125e-01,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		1.57976e+00, -5.44598e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		-5.30593e-01, -5.02237e-01, -2.05258e-01, 2.62263e-01, -2.50195e-01,
		4.28151e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
	};

	double pwc_init[200] =
	{
		// WINF
		0.00000e+00, 1.31026e+01, -4.93171e+01, 2.51045e+01, -1.30531e+01,
		6.56421e-01, 2.75633e+01, 4.36433e+00, 1.04638e+00, 5.77365e+01,
		0.00000e+00, -6.27766e+00, 2.33010e+00, -1.41351e+01, 2.49653e-01,
		0.00000e+00, 0.00000e+00, 8.00000e+01, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 1.03817e-02, -1.70950e+01, -1.92295e+00, 0.00000e+00,
		0.00000e+00, -1.17490e+01, -7.14788e-01, 6.72649e+00, 0.00000e+00,
		0.00000e+00, -1.57793e+02, -1.70815e+00, -7.92416e+00, -1.67372e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.87973e-01,
		-1.61602e-01, -1.13832e-01, -7.22447e-01, 2.21119e-02, 0.00000e+00,
		-3.01967e+00, -1.72798e-01, -5.15055e-03, -1.23477e-02, 3.60805e-03,
		-1.36730e+00, 0.00000e+00, 1.24390e-02, 0.00000e+00, -1.36577e+00,
		3.18101e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, -1.39334e+01,
		1.42088e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, -4.72219e+00,
		-7.47970e+00, -4.96528e+00, 0.00000e+00, 1.24712e+00, -2.56833e+01,
		-4.26630e+01, 3.92431e+04, -2.57155e+00, -4.35589e-02, 0.00000e+00,
		0.00000e+00, 2.02425e+00, -1.48131e+00, -7.72242e-01, 2.99008e+04,
		4.50148e-03, 5.29718e-03, -1.26697e-02, 3.20909e-02, 0.00000e+00,
		0.00000e+00, 7.01739e+00, 3.11204e+00, 0.00000e+00, 0.00000e+00,
		-2.13088e+00, 1.32789e+01, 5.07958e+00, 7.26537e-02, 2.87495e-01,
		9.97311e-03, -2.56440e+00, 0.00000e+00, 0.00000e+00, 3.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, -9.90073e-03, -3.27333e-02,
		-4.30379e+01, -2.87643e+01, -5.91793e+00, -1.50460e+02, 0.00000e+00,
		0.00000e+00, 6.55038e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 6.18051e-03, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		1.40484e+00, 5.54554e+00, 0.00000e+00, 0.00000e+00, 7.93810e+00,
		1.57192e+00, 1.03971e+00, 9.88279e-01, -4.37662e-02, -2.15763e-02,
		-2.31583e+00, 4.32633e+00, -1.12716e+00, 3.38459e-01, 4.66956e-01,
		7.18403e-01, 5.80836e-02, 4.12653e-01, 1.04111e-01, -8.30672e-02,
		-5.55541e+00, -4.97473e+00, -2.03007e+01, 0.00000e+00, -6.06235e-01,
		-1.73121e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 9.29850e-02, -6.38131e-02,
		3.93037e-02, 5.21942e-02, 2.26578e-02, 4.13157e-02, 0.00000e+00,
		6.28524e+00, 4.43721e+00, -4.31270e+00, 2.32787e+00, 2.55591e-01,
		1.60098e+00, -1.20649e+00, 3.05042e+00, -1.88944e+00, 5.35561e+00,
		2.02391e-01, 4.62950e-02, 3.39155e-01, 7.94007e-02, 6.30345e-01,
		1.93554e-01, 3.93238e-01, 1.76254e-01, -2.51359e-01, -7.06879e-01
	};

	double pwbl_init[150] =
	{
		// UGN1(1)
		6.22831e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		5.90566e+00, 0.00000e+00, 0.00000e+00, -3.20571e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, -8.30368e-01, 1.39396e-02,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 2.40657e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, -4.80790e+00, -1.62744e+00, 2.21119e-02, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 4.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.10531e-01,
		-8.94829e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
	};

	double pwcl_init[150] =
	{
		// UGN1(1)
		5.45009e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		-3.60304e+00, 0.00000e+00, 0.00000e+00, -5.04071e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 5.62113e-01, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 1.14657e+01, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 4.65483e-01, 1.73636e+00, 2.21119e-02, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 4.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -8.30769e-01,
		7.73649e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
	};

	double pwbld_init[150] =
	{
		// UN1(1)
		6.09940e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.39396e-02,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 2.21119e-02, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 4.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
	};

	double pwcld_init[150] =
	{
		// UN1(1)
		5.46739e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 2.21119e-02, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 4.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
	};

	double pb12_init[150] =
	{
		// UN1(2)
		4.99007e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		2.59994e+00, 0.00000e+00, 0.00000e+00, -1.78418e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, -5.24986e+00, 1.39396e-02,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 2.77918e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 2.21119e-02, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 4.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.68996e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
	};

	double pc12_init[150] =
	{
		// UN1(2)
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		-7.26156e+00, 0.00000e+00, 0.00000e+00, -4.12416e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, -2.88934e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 3.65720e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 2.21119e-02, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 4.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.01835e-01,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
	};

	double pb13_init[150] =
	{
		// UN1(3)
		0.00000e+00, -1.37217e+01, 0.00000e+00, 2.38712e-01, -3.92230e+00,
		6.11035e+00, -1.57794e+00, -5.87709e-01, 1.21178e+01, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 5.23202e+01, 0.00000e+00,
		0.00000e+00, 0.00000e+00, -2.22836e+03, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, -3.94006e+00, 1.39396e-02,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 3.99844e-01, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, -1.38936e+00, 2.22534e+00, 2.21119e-02, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 4.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 4.35518e-01, 8.40051e-01, 0.00000e+00, -8.88181e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 6.81729e-01, 9.67536e-01,
		0.00000e+00, -9.67836e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
	};

	double pc13_init[150] =
	{
		// UN1(3)
		0.00000e+00, -2.75655e+01, -6.61134e+00, 4.85118e+00, 8.15375e-01,
		-2.62856e+00, 2.99508e-02, -2.00532e-01, -9.35618e+00, 1.17196e+01,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, -2.43848e+00, 1.90065e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, -3.37525e-01, 1.76471e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 4.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, -9.23682e-01, -8.84150e-02, 0.00000e+00, -9.88578e-01,
		0.00000e+00, 0.00000e+00, 0.00000e+00, -1.00747e+00, -1.07468e-02,
		0.00000e+00, -3.66376e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
	};

	double pb14_init[150] =
	{
		// UN1(4)
		0.00000e+00, 1.02709e+01, 0.00000e+00, -1.42016e+00, -4.90438e+00,
		-9.11544e+00, -3.80570e+00, -2.09013e+00, 1.32939e+01, -1.28062e+01,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 1.23024e+01, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 3.92126e+02, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.39396e-02,
		0.00000e+00, 0.00000e+00, -5.56532e+00, -1.27046e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, -3.03553e+00, -9.09832e-01, 2.21119e-02, 0.00000e+00,
		8.89965e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 9.19210e-01, 0.00000e+00, 0.00000e+00, 4.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, -2.46693e-01, 7.44650e-02, 3.84661e-01, 9.44052e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, -2.25083e-01, 1.54206e-01,
		4.41303e-01, 8.74742e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
	};

	double pc14_init[150] =
	{
		// UN1(4)
		0.00000e+00, 3.61143e+00, -8.24679e+00, 1.70751e+00, 1.16676e+00,
		6.24821e+00, -5.68968e-01, 8.53046e-01, -6.94168e+00, 1.04152e+01,
		-3.70861e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, -1.23336e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 5.33958e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, -6.43682e-01, -1.00000e+00, 0.00000e+00,
		0.00000e+00, -1.00000e+00, 0.00000e+00, -5.47300e-01, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, -8.58764e-01, 4.72310e-01, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 4.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 3.37325e-01, -3.57698e-02, -6.97393e-01, 1.35387e+01,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 2.78162e-01, -2.33383e-01,
		-7.12994e-01, 1.29234e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
	};

	double pb15_init[150] =
	{
		// UN1(5)
		0.00000e+00, -1.71856e+00, 5.32877e+00, 5.33548e-01, -2.66034e+00,
		6.76192e-01, 2.25618e+00, -5.78954e-01, -2.69685e+00, 1.21933e+00,
		-6.13650e+00, 7.79531e-01, 1.63652e+00, 3.63835e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 7.51539e+00, -5.27337e-01, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 1.06625e-01, 1.39396e-02,
		0.00000e+00, 0.00000e+00, -1.07240e+00, -8.31257e-01, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 7.04016e-01, 0.00000e+00,
		7.56158e-01, -4.21268e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 1.02843e+00, 5.21034e-01, 2.21119e-02, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 4.12504e+00, 1.08459e-01, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		-3.16261e-01, 0.00000e+00, -1.44288e-01, 0.00000e+00, 4.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.36181e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
	};

	double pc15_init[150] =
	{
		// UN1(5)
		0.00000e+00, 3.47155e+00, 1.76102e+01, 2.80371e+00, -2.08556e+00,
		1.10473e+00, 6.74582e+00, -5.75988e-01, 1.02708e+00, -2.23381e+01,
		8.60171e+00, 5.12046e-01, -8.12861e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 9.11036e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 3.89742e+00, 2.01725e-01, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 5.06308e-01, 2.04596e-01, 0.00000e+00,
		4.40377e+00, 0.00000e+00, 0.00000e+00, 2.20760e+00, 0.00000e+00,
		-1.36478e+00, 2.38097e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, -7.08949e-02, -1.61277e-01, 2.21119e-02, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, -2.16898e+00, -5.31596e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		2.53060e+00, 0.00000e+00, -7.17287e-01, 0.00000e+00, 4.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -1.91762e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
	};

	double pb15d_init[150] =
	{
		// UGN1(2)
		0.00000e+00, -7.70936e-01, 1.58158e+00, 3.61790e+00, -1.51748e+00,
		-5.66098e-01, 1.69393e+00, -4.60489e-01, -8.31527e-01, -4.66437e-01,
		-1.21750e+00, 0.00000e+00, 0.00000e+00, 1.56505e+02, 0.00000e+00,
		0.00000e+00, 0.00000e+00, -5.19321e+01, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.39396e-02,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 3.09223e-01, 1.33715e-01, 2.21119e-02, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 4.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
	};

	double pc15d_init[150] =
	{
		// UGN1(2)
		0.00000e+00, 1.72324e-01, 3.08033e-01, 4.55771e-01, 1.46516e-01,
		1.97176e-01, -1.53329e-01, 6.91877e-02, -3.07184e-01, 2.65686e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, -2.24369e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 4.04079e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		4.99627e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, -7.83317e-03, -6.88967e-02, 2.21119e-02, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 4.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
	};

	double pwp_init[26][100] =
	{
		// UN2(2)
		{
			0.00000e+00, -7.99767e-01, -3.24774e-01, 7.70975e-01, 6.71796e-01,
			5.65483e-01, -2.99727e+00, 3.32448e+00, -9.15018e-01, 5.97656e+00,
			0.00000e+00, -1.19515e+00, -8.30457e-01, 3.26074e+00, 0.00000e+00,
			0.00000e+00, -1.58365e+00, 7.44825e-02, 5.91372e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -1.41511e-01, -3.01048e+00,
			2.35960e+01, 0.00000e+00, -1.70352e+00, -2.39746e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 1.30488e+00, 0.00000e+00,
			5.95132e-01, 5.64301e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 5.30317e-01, 5.66569e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 5.72367e+00, 1.58411e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			1.04557e-01, 0.00000e+00, -2.04710e-01, 0.00000e+00, 5.00000e+00
		},

		// UN2(2)
		{
			0.00000e+00, 6.34487e+00, 9.84162e+00, 3.42136e+00, -5.10607e+00,
			-8.58745e-02, 3.11501e+00, 5.34570e-01, 1.18027e+00, 4.28027e+00,
			4.75123e+00, 6.40947e-01, -4.15165e+00, -1.38154e+01, 0.00000e+00,
			0.00000e+00, 1.13145e+01, -5.15954e+00, 0.00000e+00, 0.00000e+00,
			1.35576e+01, 0.00000e+00, -5.78982e+00, -2.22043e+00, 3.36776e+00,
			3.04791e+01, 0.00000e+00, 2.94709e+00, -4.17536e-01, -1.59855e+00,
			-2.18320e+00, 1.68269e+01, 0.00000e+00, 1.00829e+00, 0.00000e+00,
			-6.85096e-01, 2.07822e-01, 3.50168e-01, -3.03662e+01, 0.00000e+00,
			0.00000e+00, -1.65726e-01, -8.97831e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -5.24159e+00, 0.00000e+00, -3.52218e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 5.69093e-01, -7.44918e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			2.10865e+00, 0.00000e+00, 1.76776e-01, 1.54755e+00, 5.00000e+00
		},
		// UN2(3)
		{
			0.00000e+00, 2.28657e+00, 4.96548e-01, 6.99915e+00, -2.31540e+00,
			-1.82163e-01, -5.00779e-01, 3.18199e-01, -6.14645e-01, 6.34816e+00,
			0.00000e+00, 7.94635e-01, -5.55565e-01, 3.85494e+00, 0.00000e+00,
			0.00000e+00, -3.96109e+00, 1.90775e-01, 4.51396e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -5.04618e-01, -4.14385e+00,
			2.30244e+01, 0.00000e+00, 1.00689e+00, 5.75680e-02, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 8.56741e-01, 0.00000e+00,
			9.54921e-02, 5.56659e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 1.38503e-01, 4.50415e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 2.22813e-01, -8.63549e-02, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			1.37970e-01, 0.00000e+00, -3.25612e-01, 0.00000e+00, 5.00000e+00
		},
		// UN2(3)
		{
			0.00000e+00, 5.07608e+00, 3.31479e+00, 3.01548e-01, -1.12100e+00,
			-7.63711e-02, 2.29748e+00, -1.36699e+00, 7.53433e-01, 3.60702e+01,
			-1.55266e+00, 1.47382e+00, -2.53895e+00, -1.47720e+01, 0.00000e+00,
			0.00000e+00, 1.11787e+01, -1.06256e+01, 0.00000e+00, 0.00000e+00,
			7.86391e+00, 0.00000e+00, -8.61020e+00, -1.59313e+00, -5.17013e+00,
			1.20468e+00, 0.00000e+00, 5.76509e-01, 9.96195e-01, -1.45539e+00,
			-1.79950e+01, 8.76957e+00, 0.00000e+00, -1.22863e+00, 0.00000e+00,
			-6.19019e-01, -1.09571e-01, -4.31325e-02, -4.21981e+01, 0.00000e+00,
			0.00000e+00, -1.51519e-01, -1.24067e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -6.39248e+00, 0.00000e+00, 6.64508e-01,
			-7.33184e-01, -9.72031e-03, 1.36789e+00, -8.62311e-01, -3.06395e-03,
			2.53354e-01, -2.40918e-01, -4.06932e-02, -5.82223e-01, 0.00000e+00,
			-8.70285e-01, 7.72318e-01, -6.54213e-01, -2.19231e+01, -1.56509e-01,
			2.71745e-01, 5.93538e-01, 2.27757e-01, -5.98215e-01, 3.96457e-01,
			2.98705e-01, 1.78618e-01, -5.24538e-01, 1.16439e-01, 7.56829e-02,
			-4.26809e-01, 5.77187e-01, 8.65450e-01, -7.53614e-01, 1.38381e-01,
			-1.82265e-01, 2.85263e-01, 4.51322e-01, 1.02775e-01, 3.55731e-01,
			-4.60896e-01, -3.13037e+01, -2.70818e+00, -7.84847e-01, 0.00000e+00,
			-1.03473e-01, -3.87649e-01, -1.22613e-01, 0.00000e+00, 0.00000e+00,
			8.91325e-01, 0.00000e+00, 1.06189e-01, 9.13030e-02, 5.00000e+00
		},
		// UN2(4)
		{
			0.00000e+00, 2.94921e+00, 2.79238e+00, 2.58949e+00, 3.56459e-01,
			3.12952e-01, 3.34337e+00, -2.83209e+00, -1.05979e+00, 3.92313e+00,
			0.00000e+00, 1.73703e-01, -3.23441e-01, 4.15836e+00, 0.00000e+00,
			0.00000e+00, -1.77156e+00, 6.44113e-01, 1.88743e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -4.64778e-01, -4.23560e+00,
			2.27271e+01, 0.00000e+00, -4.89468e-01, 1.82689e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 4.38217e-02, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 8.62449e-02, 4.46041e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -1.40393e-01, 1.01821e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.00000e+00
		},
		// UN2(4)
		{
			0.00000e+00, 6.04465e+00, 4.50924e+00, 3.84425e-02, -8.70772e-01,
			-9.55408e-02, 2.28287e+00, -4.37834e-01, 3.57839e-01, 7.20721e+01,
			-4.41757e+00, -9.13648e-01, -8.71866e-01, -6.26173e+00, 0.00000e+00,
			0.00000e+00, 5.92817e+00, 6.15853e+00, 0.00000e+00, 0.00000e+00,
			-4.89060e+00, 0.00000e+00, -8.30378e+00, 1.07462e-01, 1.08471e+02,
			3.39150e+01, -4.57863e+00, -7.18349e-02, -2.71703e-01, -8.96297e+00,
			-2.37986e+01, 4.11880e+00, 0.00000e+00, -9.95820e-01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -8.91622e+00, -6.85950e+01, 0.00000e+00,
			0.00000e+00, -3.62769e-02, -1.65893e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -2.94563e+00, 0.00000e+00, 1.23581e+00,
			-6.06026e-01, -6.50229e-01, 1.91330e+00, -1.00314e+00, 1.13346e-01,
			4.21885e-01, -3.97688e-01, -2.77437e-01, -6.65893e-01, 0.00000e+00,
			-1.37646e+00, 1.35171e+00, -9.55595e-01, -1.96450e+01, -2.50039e-01,
			5.93389e-01, 9.87131e-01, 5.43559e-01, -1.04322e+00, 6.32546e-01,
			3.73259e-01, 5.22657e-01, -5.81400e-01, -1.26425e-01, -1.29843e-01,
			-5.36598e-01, 8.02402e-01, 9.04347e-01, -1.10799e+00, 1.24800e-01,
			1.62487e-02, 2.84237e-01, -1.68866e+00, 5.07723e-01, 5.14161e-01,
			-4.71292e-01, -3.03487e+01, 4.17455e-01, -1.12591e+00, 0.00000e+00,
			-3.03544e-01, -6.60313e-01, -1.48606e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 1.00607e+01, 5.00000e+00
		},
		// UN2(5)
		{
			0.00000e+00, 2.52207e+00, 3.84550e+00, 1.68023e+00, 7.93489e-01,
			3.93725e-02, -2.79707e+00, -4.76621e-01, -1.19972e-01, 3.20454e-01,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 4.17146e+00, 0.00000e+00,
			0.00000e+00, -5.30699e-01, 9.14373e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -4.84434e-02, 1.85902e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.00000e+00
		},
		// UN2(5)
		{
			0.00000e+00, 1.55386e+01, 4.21418e+00, -9.70151e-01, -8.77326e-01,
			2.65813e-02, 1.40164e+00, -9.03874e-01, 3.17281e-03, 9.26891e+01,
			-4.96004e+00, 0.00000e+00, 0.00000e+00, -4.17851e+00, 0.00000e+00,
			0.00000e+00, -1.14760e+01, 2.67744e+00, 0.00000e+00, 0.00000e+00,
			-1.60056e+01, 0.00000e+00, -7.14647e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -2.89639e+00, 0.00000e+00, 0.00000e+00, -3.88601e+00,
			-1.65784e+01, 8.44796e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -3.75324e+00, -6.24047e+01, 0.00000e+00,
			0.00000e+00, -2.86808e-02, -1.95891e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -3.10534e-01, 0.00000e+00, -3.37448e+00,
			1.63964e-02, -1.45191e+00, 1.85618e+00, -9.61979e-01, 3.93783e-01,
			4.21681e-01, -5.30254e-01, -2.96232e-01, -7.55211e-01, 0.00000e+00,
			-1.85443e+00, 1.88047e+00, -1.07818e+00, -1.35373e+01, -3.05785e-01,
			7.82159e-01, 1.32586e+00, 2.34413e-01, -7.47152e-01, 9.92893e-01,
			-2.80110e-02, 3.61747e-01, -4.16280e-01, -3.46427e-01, -5.76431e-01,
			-2.13906e-01, 9.51184e-01, 3.69403e-01, -1.35563e+00, 6.59534e-02,
			1.39764e-01, 4.50687e-01, -1.22025e+00, 5.73280e-02, 7.49303e-01,
			-8.37947e-01, -3.01332e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-4.36697e-01, -7.76068e-01, -1.41680e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 1.21958e+01, 5.00000e+00
		},
		// UN2(6)
		{
			0.00000e+00, 3.13842e+00, -8.20417e-01, 3.72282e+00, -5.20477e-01,
			-3.61867e-01, -2.92604e+00, 3.13013e-01, -1.38865e-01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 1.30060e+01, 0.00000e+00,
			0.00000e+00, 1.67696e+00, 9.85990e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -8.46922e-02, 5.59429e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.00000e+00
		},
		// UN2(6)
		{
			0.00000e+00, 1.78539e+01, 1.07314e+01, -1.13212e+00, 1.59867e-02,
			1.53736e-01, 2.25710e+00, -9.39080e-01, -9.72620e-02, 9.89789e+01,
			-5.17469e+00, 0.00000e+00, 0.00000e+00, -2.98597e+00, 0.00000e+00,
			0.00000e+00, -2.04707e+01, 4.92899e+00, 0.00000e+00, 0.00000e+00,
			-1.44316e+01, 0.00000e+00, -3.31557e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -6.22743e+00, 0.00000e+00, 0.00000e+00, -4.34344e+00,
			-8.29640e+00, -3.03800e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 2.79387e+00, -5.23752e+01, 0.00000e+00,
			0.00000e+00, -2.59963e-02, -1.73426e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -5.37220e+00, 0.00000e+00, -6.53478e-01,
			3.48181e-01, -1.88980e+00, 1.47787e+00, -7.92670e-01, 6.49224e-01,
			5.96079e-01, -1.04901e+00, -5.24003e-01, -6.77311e-01, 0.00000e+00,
			-2.26873e+00, 2.80910e+00, -9.84994e-01, -6.79661e+00, -3.71975e-01,
			1.13310e+00, 1.57164e+00, 2.15176e-01, -5.58583e-01, 1.16045e+00,
			2.05395e-02, 2.27714e-01, 1.41203e-01, -3.92231e-01, -8.82859e-01,
			4.90400e-01, 1.14013e+00, -2.25250e-01, -1.64930e+00, 5.73434e-02,
			1.89857e-01, 4.31221e-01, -1.35345e+00, -2.94189e-01, 6.87530e-01,
			-7.78284e-01, -2.88975e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-3.98115e-01, -7.40699e-01, -8.28264e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 2.02069e+00, 5.00000e+00
		},
		// UN2(7)
		{
			0.00000e+00, 2.08818e+00, -1.96235e+00, 4.55317e+00, -1.76012e+00,
			-4.75258e-01, -1.44220e+00, -3.28566e-01, -1.41177e-01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 1.49146e+01, 0.00000e+00,
			0.00000e+00, 1.73222e+00, 9.91286e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.35468e-01, 1.91833e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.00000e+00
		},
		// UN2(7)
		{
			0.00000e+00, 1.25645e+01, 2.43937e+01, -4.89691e-01, -5.46437e-01,
			1.22200e-01, 2.89309e+00, -2.85509e-01, -2.27122e-01, 9.54192e+01,
			-4.07394e+00, 0.00000e+00, 0.00000e+00, -3.04354e+00, 0.00000e+00,
			0.00000e+00, -2.36547e+01, 1.04903e+01, 0.00000e+00, 0.00000e+00,
			-8.32274e+00, 0.00000e+00, -3.34712e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -7.95953e+00, 0.00000e+00, 0.00000e+00, -5.83474e+00,
			-1.48074e+00, 1.02268e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 6.19470e+00, -3.90767e+01, 0.00000e+00,
			0.00000e+00, -3.58136e-03, 1.22289e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -8.49787e+00, 0.00000e+00, -3.97498e+00,
			3.79580e-01, -1.93595e+00, 2.89114e+00, -4.73457e-01, 7.67548e-01,
			5.66859e-01, -1.28683e+00, -8.37174e-01, -3.48022e-01, 0.00000e+00,
			-2.62865e+00, 3.50575e+00, -7.93257e-01, -8.10692e-01, -4.99450e-01,
			1.56654e+00, 1.63039e+00, 7.58900e-02, -4.30952e-01, 1.23068e+00,
			1.06404e-01, 4.73870e-02, 5.50559e-01, -4.11375e-01, -9.94162e-01,
			1.35025e+00, 1.26053e+00, -7.34502e-01, -2.01952e+00, 2.05398e-01,
			-4.77248e-02, 2.41549e-01, -9.32522e-01, -5.63663e-01, 5.34833e-01,
			-5.77563e-01, -2.65033e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-2.42317e-01, -7.33679e-01, -7.85537e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 1.56842e-01, 5.00000e+00
		},
		// UN2(8)
		{
			0.00000e+00, 7.00409e-01, -4.17017e-01, 3.24757e+00, -1.28352e+00,
			-4.23875e-01, 1.64346e+00, -1.20855e+00, -7.65316e-01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -3.39417e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 2.68534e+01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.56444e-01, -4.60043e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.00000e+00
		},
		// UN2(8)
		{
			0.00000e+00, 7.30129e+00, 3.14811e+01, -7.06834e-02, -2.96193e-01,
			1.73817e-01, 1.62127e+00, -2.71556e-01, -2.05844e-01, 8.02088e+01,
			-1.86956e-01, 0.00000e+00, 0.00000e+00, -9.43641e-01, -3.24716e+00,
			0.00000e+00, -2.32748e+01, 1.96724e+01, 0.00000e+00, 0.00000e+00,
			-3.95949e+00, 0.00000e+00, 5.44787e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.00161e+01, 0.00000e+00, 0.00000e+00, -4.57422e+00,
			4.31304e+00, 1.49868e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 5.99489e+00, -2.82120e+01, 0.00000e+00,
			0.00000e+00, 4.03624e-02, 1.19463e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -1.39050e+01, 0.00000e+00, -2.65634e+00,
			6.37036e-01, -1.77461e+00, 3.03103e+00, -1.49839e-01, 7.02027e-01,
			6.08841e-01, -9.27289e-01, -8.52362e-01, 5.61723e-01, 0.00000e+00,
			-2.72061e+00, 3.66183e+00, -2.54943e-01, 2.94668e+00, -3.57898e-01,
			1.71858e+00, 1.58782e+00, -2.42995e-01, -3.57783e-01, 1.20157e+00,
			2.58895e-01, -1.05773e-01, 5.79397e-01, -3.30395e-01, -4.03569e-01,
			1.99175e+00, 1.21688e+00, -8.64350e-01, -1.95569e+00, 4.61136e-01,
			-8.61382e-02, 3.38859e-01, 0.00000e+00, -5.78864e-01, 4.46659e-01,
			-4.57428e-01, -1.99920e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-1.19841e-01, -4.56968e-01, 2.00180e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -1.07368e+00, 5.00000e+00
		},
		// UN2(9)
		{
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 1.75863e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 3.18522e+01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.00000e+00
		},
		// UN2(9)
		{
			0.00000e+00, 4.61019e-02, 3.50615e+01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 6.15349e+01,
			4.28634e+00, 0.00000e+00, 0.00000e+00, 6.03982e+00, -4.72305e+00,
			0.00000e+00, -1.43678e+01, 3.62580e+01, 0.00000e+00, 0.00000e+00,
			1.26574e+00, 0.00000e+00, -2.77285e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.14802e+01, 0.00000e+00, 0.00000e+00, -1.11940e+01,
			-1.39535e+00, 2.63070e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -1.53024e+00, -2.14609e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -1.26956e+01, 0.00000e+00, 5.49926e+00,
			9.80142e-01, -1.19016e+00, 2.75110e+00, 4.23423e-01, 5.89893e-01,
			4.94288e-01, -5.25954e-01, -8.51760e-01, 1.62676e+00, 0.00000e+00,
			-1.90027e+00, 3.19950e+00, 4.72739e-01, 7.04179e+00, -1.43685e-03,
			1.43219e+00, 1.32136e+00, -2.92744e-03, -3.43680e-01, 7.75735e-01,
			6.92202e-01, -1.45519e-01, 6.97813e-02, -3.11588e-01, 6.65750e-01,
			2.33809e+00, 1.06694e+00, -5.77590e-01, -1.33717e+00, 8.13367e-01,
			-5.05737e-01, 5.99169e-01, -8.83386e-01, -4.38123e-01, 2.63649e-01,
			-3.03448e-01, -1.28190e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			1.45478e-02, 1.45491e-01, 2.40080e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -3.86910e+00, 5.00000e+00
		},
		// UN2(10)
		{
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 1.10647e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 3.13252e+01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.00000e+00
		},
		// UN2(10)
		{
			0.00000e+00, -3.03260e+00, 3.15488e+01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 4.42798e+01,
			7.08849e+00, 0.00000e+00, 0.00000e+00, 1.64773e+01, -6.86505e+00,
			0.00000e+00, -6.27112e+00, 3.78373e+01, 0.00000e+00, 0.00000e+00,
			2.97763e+00, 0.00000e+00, -3.44134e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.19424e+01, 0.00000e+00, 0.00000e+00, -1.64645e+01,
			-2.27053e+00, 3.82330e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 1.33140e-01, -2.08131e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -7.04687e+00, 0.00000e+00, 6.52184e+00,
			7.31799e-01, -2.75395e-01, 1.92467e+00, 8.71269e-01, 3.72836e-01,
			3.04967e-01, 7.72480e-02, -5.08596e-01, 1.99828e+00, 0.00000e+00,
			-5.51169e-01, 2.12420e+00, 8.96069e-01, 1.12092e+01, -4.30438e-02,
			7.38391e-01, 6.12050e-01, 3.62981e-02, -1.02054e-01, 1.82404e-01,
			3.70643e-01, -1.68899e-01, -1.79628e-01, -1.21117e-01, 1.45823e+00,
			2.04352e+00, 7.83711e-01, -3.42979e-02, -2.31363e-01, 7.11253e-01,
			-3.16353e-01, 6.21069e-01, -1.05676e+00, -4.03488e-01, 4.11595e-01,
			-2.12535e-01, -6.51453e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			1.48238e-01, 6.38716e-01, 2.99311e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -1.01846e+00, 5.00000e+00
		},
		// UN2(11)
		{
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 2.21764e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 6.77475e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.00000e+00
		},
		// UN2(11)
		{
			0.00000e+00, -1.74115e+00, 2.66621e+01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.13017e+01,
			6.86985e+00, 0.00000e+00, 0.00000e+00, 2.08835e+01, -7.86030e+00,
			0.00000e+00, -3.77141e+00, 3.87788e+01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 1.31580e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -9.98927e+00, 0.00000e+00, 0.00000e+00, -1.71002e+01,
			-9.88358e-01, 4.47756e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 5.95029e-01, -2.11313e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -3.84164e+00, 0.00000e+00, 0.00000e+00,
			3.07191e-01, 4.79094e-02, 6.72159e-01, 5.54185e-01, 1.82847e-01,
			-1.23768e-02, 1.91637e-01, -2.89429e-02, 1.18297e+00, 0.00000e+00,
			2.37450e-01, 9.23551e-01, 6.05670e-01, 1.35990e+01, -1.64210e-01,
			5.38355e-03, -4.91246e-02, -1.06966e-01, -2.09635e-01, -3.23023e-02,
			-3.41663e-02, -3.48871e-02, -2.62450e-01, 2.21492e-01, 1.43749e+00,
			1.08677e+00, 3.97778e-01, 3.61526e-01, 5.55950e-01, 3.53058e-01,
			-5.93339e-02, 4.14203e-01, -6.05024e-01, -1.38714e-01, 2.78897e-01,
			-8.92889e-02, -3.59033e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			9.90623e-02, 4.36170e-01, 7.95418e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -1.11426e+00, 5.00000e+00
		},
		// UN2(12)
		{
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 3.07320e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 1.60738e+01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.00000e+00
		},
		// UN2(12)
		{
			0.00000e+00, 1.26217e+01, 2.30787e+01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00029e+01,
			-2.88682e+00, 0.00000e+00, 0.00000e+00, 2.09439e+01, -4.56923e+00,
			0.00000e+00, -2.15929e+00, 3.87149e+01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -7.98039e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -6.63423e+00, 0.00000e+00, 0.00000e+00, -5.84850e+00,
			3.72111e+00, 4.52300e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 3.21872e-01, 0.00000e+00, 0.00000e+00,
			1.09405e-02, -4.35341e-02, 8.00586e-02, 1.48577e-01, 1.01602e-01,
			-1.01104e-01, -1.98993e-02, 3.51174e-02, 2.41112e-01, 0.00000e+00,
			2.76479e-01, 1.97043e-01, 2.68708e-01, 1.39832e+01, -1.56638e-01,
			-2.39101e-01, -1.50605e-01, -2.17139e-01, -2.59057e-01, -4.36362e-01,
			-1.43496e-01, 7.51305e-02, -2.40850e-01, 1.34858e-01, 7.59193e-01,
			3.52708e-01, 1.29922e-01, 3.27957e-01, 5.35491e-01, 1.19120e-01,
			-2.94029e-02, 1.76113e-01, -6.51597e-01, 3.61575e-02, 4.26836e-02,
			-2.29297e-02, -4.27373e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-2.78548e-02, 5.77322e-02, -1.02411e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.00000e+00
		},
		// UN2(13)
		{
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 3.69447e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 2.34073e+01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.00000e+00
		},
		// UN2(13)
		{
			0.00000e+00, 1.22096e+01, 1.92342e+01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 8.13667e+00,
			-6.19078e+00, 0.00000e+00, 0.00000e+00, 2.37009e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -7.87365e+01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -1.12371e+01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -2.76047e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 1.85864e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.00000e+00
		},
		// UN2(14)
		{
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 1.01008e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 2.21469e+01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.00000e+00
		},
		// UN2(14)
		{
			0.00000e+00, -1.40697e+00, 6.88709e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 3.67624e+02, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 1.58312e+01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -2.46486e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.90327e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 1.13248e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.00000e+00
		}
	};


	for (i = 0; i < 200; ++i) {
		pwb[i] = pwb_init[i];
		pwc[i] = pwc_init[i];
	}

	for (i = 0; i < 150; ++i) {
		pwbl[i] = pwbl_init[i];
		pwcl[i] = pwcl_init[i];
		pwbld[i] = pwbld_init[i];
		pwcld[i] = pwcld_init[i];
		pb12[i] = pb12_init[i];
		pc12[i] = pc12_init[i];
		pb13[i] = pb13_init[i];
		pc13[i] = pc13_init[i];
		pb14[i] = pb14_init[i];
		pc14[i] = pc14_init[i];
		pb15[i] = pb15_init[i];
		pc15[i] = pc15_init[i];
		pb15d[i] = pb15d_init[i];
		pc15d[i] = pc15d_init[i];
	}

	for (i = 0; i < 26; ++i) {
		for (j = 0; j < 100; ++j) {
			pwp[i][j] = pwp_init[i][j];
		}
	}


	return;


}

