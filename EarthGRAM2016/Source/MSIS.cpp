//P. White
//MSIS thermosphere model Class
//Calculates atmosphere variables and constituents in thermosphere region

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include "MSIS.h"
using namespace std;

MSIS::MSIS()
{
	intializeMemberVariables();

}

MSIS::~MSIS()
{

	for (int i = 0; i < 25; i++){

		delete[] pdl[i];

	}

	delete[] pdl;
}

void MSIS::intializeMemberVariables()
{

	//Initialize member variables for MSIS class

	pi = 3.1415926535897931;
	pi180 = pi / 180.0;
	for (int i = 0; i < 9; i++){
		for (int j = 0; j < 4; j++){
			plg[i][j] = 0.0;
		}
	}
	// gts3c
	dd = 0.0, ctloc = 0.0, stloc = 0.0, c2tloc = 0.0, s2tloc = 0.0;
	c3tloc = 0.0, s3tloc = 0.0, day = 0.0, dfa = 0.0, apdf = 0.0, xlong = 0.0;
	gsurf = 0.0, re = 0.0;

	for (int i = 0; i < 4; ++i) {
		apt[i] = 0.0;
	}

	// meso7
	for (int i = 0; i < 5; ++i) {
		tn1[i] = 0.0;       // K Lower thermosphere temperature variations

		if (i < 4) {
			tn2[i] = 0.0;   // K Lower mesosphere/upper stratosphere temperature variations
		}

		tn3[i] = 0.0;       // K Lower stratosphere and troposphere temperature variations
	}

	for (int i = 0; i < 2; ++i) {
		tgn1[i] = 0.0;      // K Lower thermosphere temperature gradient
		tgn2[i] = 0.0;      // K Lower mesospherer/upper stratosphere temperature gradient
		tgn3[i] = 0.0;      // K Lower stratosphere and troposphere temperature gradient
	}


	pdl = new float*[25];
	for (int i = 0; i < 25; ++i){
		pdl[i] = new float[2];
	}


	for (int i = 0; i < 25; ++i){
		for (int j = 0; j < 2; ++j){
			pdl[i][j] = 0.0;
		}
	}

	//Initialize tabulated MSIS data
	gtd7bk();

	dtz = 0.0;
	re = 0.0;
	mn = 0;
}

double MSIS::densu(double alt, double dlb, double tinf, double tlb, double xm, double alpha, double *tz, double zlb, double s2, int mn1, double *zn1, double *tn1, double *tgn1)
{
	//densu member function from MSIS class
	//Calculate temperature and density profiles for MSIS models
	//New lower thermo polynomial 10/30/89

	const double rgas = 831.4;

	double densu_out, za, z, zg2, tt, dta, z1, z2, t1, t2, zg, zgdiff;
	double xs[5], ys[5], y2out[5], x, y, glb, ta, yd1, yd2, expll, densa, gamma1;
	double gamm, yi;

	int mn, k;

	densu_out = 1.0;

	//Joining altitudes of Bates and spline
	za = zn1[0];
	z = fmax(alt, za);

	//Geopotantial altitude difference from zlb
	zg2 = zeta(z, zlb, re);

	//Bates temperature
	tt = tinf - (tinf - tlb)*exp(-s2*zg2);
	ta = tt;
	*tz = tt;
	densu_out = *tz;

	if (alt < za){
		//Calculate temperature below za
		//Temperature gradient at za from Bates profile
		dta = (tinf - ta)*s2*pow(((re + zlb) / (re + za)), 2.0);
		tgn1[0] = dta;
		tn1[0] = ta;
		z = fmax(alt, zn1[mn1 - 1]);
		mn = mn1;
		z1 = zn1[0];
		z2 = zn1[mn - 1];
		t1 = tn1[0];
		t2 = tn1[mn - 1];
		//Geopotential difference from z1
		zg = zeta(z, z1, re);
		zgdiff = zeta(z2, z1, re);

		//Set up spline nodes
		for (k = 0; k < mn; k++){
			xs[k] = zeta(zn1[k], z1, re) / zgdiff;
			ys[k] = 1.0 / tn1[k];
		}

		//End node derivatives
		yd1 = -tgn1[0] / (t1*t1)*zgdiff;
		yd2 = -tgn1[1] / (t2*t2)*zgdiff*pow((re + z2) / (re + z1), 2.0);

		//Calculate spline coefficients
		jb2008.hwm.spline(xs, ys, mn, yd1, yd2, y2out);
		x = zg / zgdiff;

		jb2008.hwm.splint(xs, ys, y2out, mn, x, &y);

		//Temperature at altitude
		*tz = 1.0 / y;

		densu_out = *tz;
	}

	if (xm != 0.0){
		//Calculate density above za
		glb = gsurf / pow(1.0 + zlb / re, 2.0);

		gamma1 = xm*glb / (s2*rgas*tinf);
		expll = exp(-s2*gamma1*zg2);

		if ((expll > 50.0) || (0.0 >= tt)){
			expll = 50.0;
		}

		//Density at altitude
		densa = dlb*pow(tlb / tt, 1.0 + alpha + gamma1)*expll;
		densu_out = densa;

		if (alt < za){
			//Calculate density below za
			glb = gsurf / pow(1.0 + z1 / re, 2.0);

			gamm = xm*glb*zgdiff / rgas;

			//Integrate spline temperatures
			splini(xs, ys, y2out, mn, x, &yi);
			expll = gamm*yi;

			if ((expll > 50.0) || (*tz <= 0.0)){
				expll = 50.0;
			}

			densu_out = densu_out*pow(t1 / (*tz), 1.0 + alpha)*exp(-expll);
		}
	}


	return densu_out;
}

double MSIS::densm(double alt, double d0, double xm, double *tz, int mn3, double *zn3,
	double *tn3, double tgn3[2], int mn2, double *zn2, double *tn2, double tgn2[2])
{
	//densm member function from MSIS class
	//Calculate temperature and density profiles for lower atmos

	const double rgas = 831.4;
	double densm_out, xs[10], ys[10], y2out[10], z, z1, z2, t1, zg, zgdif, yd1, yd2,
		t2, x, y, glb, gamm, yi, expll;
	int mn, k;

	densm_out = d0;

	if (alt <= zn2[0]){
		//Stratosphere/mesosphere temperature
		z = fmax(alt, zn2[mn2 - 1]);
		mn = mn2;
		z1 = zn2[0];
		z2 = zn2[mn - 1];
		t1 = tn2[0];
		t2 = tn2[mn - 1];
		zg = zeta(z, z1, re);
		zgdif = zeta(z2, z1, re);

		//Set up slpine nodes
		for (k = 0; k < mn; k++){
			xs[k] = zeta(zn2[k], z1, re) / zgdif;
			ys[k] = 1.0 / tn2[k];
		}

		yd1 = -tgn2[0] / (t1*t1)*zgdif;
		yd2 = -tgn2[1] / (t2*t2)*zgdif*pow((re + z2) / (re + z1), 2);

		//Calculate spline coefficients
		jb2008.hwm.spline(xs, ys, mn, yd1, yd2, y2out);
		x = zg / zgdif;
		jb2008.hwm.splint(xs, ys, y2out, mn, x, &y);

		//Temperature at altitude
		*tz = 1.0 / y;

		if (xm != 0.0){
			//Calculate stratosphere/mesosphere density
			glb = gsurf / pow(1.0 + z1 / re, 2);
			gamm = xm*glb*zgdif / rgas;

			//Integrate temperature profile
			splini(xs, ys, y2out, mn, x, &yi);
			expll = gamm*yi;

			if (expll > 50.0){
				expll = 50.0;
			}

			//Density at altitude
			densm_out = densm_out*(t1 / (*tz))*exp(-expll);
		}

		if (alt <= zn3[0]){
			//Troposphere/stratosphere temperature
			z = alt;
			mn = mn3;
			z1 = zn3[0];
			z2 = zn3[mn - 1];
			t1 = tn3[0];
			t2 = tn3[mn - 1];
			zg = zeta(z, z1, re);
			zgdif = zeta(z2, z1, re);

			//Set up spline nodes
			for (k = 0; k < mn; k++){
				xs[k] = zeta(zn3[k], z1, re) / zgdif;
				ys[k] = 1.0 / tn3[k];
			}

			yd1 = -tgn3[0] / (t1*t1)*zgdif;
			yd2 = -tgn3[1] / (t2*t2)*zgdif*pow((re + z2) / (re + z1), 2);

			//Calculate spline coefficients
			jb2008.hwm.spline(xs, ys, mn, yd1, yd2, y2out);
			x = zg / zgdif;
			jb2008.hwm.splint(xs, ys, y2out, mn, x, &y);

			//Temperature at altitude
			*tz = 1.0 / y;

			if (xm != 0.0){
				//Calculate tropospheric/stratosphere density
				glb = gsurf / pow(1.0 + z1 / re, 2);
				gamm = xm*glb*zgdif / rgas;

				//Integrate temperature profile
				splini(xs, ys, y2out, mn, x, &yi);
				expll = gamm*yi;

				if (expll > 50){
					expll = 50.0;
				}

				densm_out = densm_out*(t1 / (*tz))*exp(-expll);
			}
		}
	}

	if (xm == 0.0){
		densm_out = *tz;
	}

	return densm_out;
}

double MSIS::dnet(double dd, double dm, double zhm, double xmm, double xm)
{
	//dnet member function from MSIS class
	//Turbopause correction for MSIS model
	//Root mean density
	//8/20/80
	//dd - diffusive density
	//dm - full mixed density
	//zhm - transition scale length
	//xmm - full mixed molecular weight
	//xm - species molecular weight
	//dnet_out - combined density

	double dnet_out, a, ylog;

	a = zhm / (xmm - xm);

	if ((dm <= 0.0) && (dd <= 0.0)){

		if ((dd != 0.0) && (dm != 0.0)){
			dd = 1.0;
		}

		if (dm == 0.0){
			dnet_out = dd;
			return dnet_out;
		}

		if (dd == 0.0){
			dnet_out = dm;
			return dnet_out;
		}
	}

	ylog = a*log(dm / dd);

	if (ylog < -10.0){
		dnet_out = dd;
		return dnet_out;
	}
	else if (ylog > 10.0){
		dnet_out = dm;
		return dnet_out;
	}

	dnet_out = dd*pow(1.0 + exp(ylog), 1 / a);

	return dnet_out;
}

double MSIS::ccor(double alt, double r, double h1, double zh)
{
	//ccor member function from MSIS class
	//Chemistry/dissociation correction for MSIS models
	//alt - altitude
	//r - target ratio
	//h1 - transition scale length
	//zh - altitude of 1/2 r

	double ccor_out, e, ex;

	e = (alt - zh) / h1;

	if (e > 70.0){
		ccor_out = 0.0;
	}
	else if (e < -70.0){
		ccor_out = r;
	}
	else {
		ex = exp(e);
		ccor_out = r / (1.0 + ex);
	}

	ccor_out = exp(ccor_out);

	return ccor_out;
}

double MSIS::ccor2(double alt, double r, double h1, double zh, double h2)
{
	//ccor2 member function from MSIS class
	//O and O2 chemistry dissociation correction for MSIS models

	double ccor2_out, e1, e2, ex1, ex2;

	e1 = (alt - zh) / h1;
	e2 = (alt - zh) / h2;

	if ((e1 > 70.0) || (e2 > 70.0)){
		ccor2_out = 0.0;
	}
	else if ((e1 < -70.0) && (e2 < -70.0)){
		ccor2_out = r;
	}
	else {
		ex1 = exp(e1);
		ex2 = exp(e2);
		ccor2_out = r / (1.0 + 0.5*(ex1 + ex2));
	}

	ccor2_out = exp(ccor2_out);

	return ccor2_out;

}

void MSIS::meters(bool meter)
{
	//meters member function from MSIS class
	//Convert outputs to kg and meters if meter true

	imr = 0;

	if (meter) {
		imr = 1;
	}

	return;
}

void MSIS::splini(double *xa, double *ya, double *y2a, int n, double x, double *yi)
{
	//splini member function MSIS class
	//Integrate cubic spline function from xa[0] to x
	//xa, ya: arrays of tabulated function in ascending order by x
	//y2a: array of second derivatives
	//n: size of arrays xa, ya, y2a
	//x: abscissa endpoint for integration
	//yi: output value

	int klo, khi;
	double xx, h, a, b, a2, b2;

	*yi = 0;
	klo = 1;
	khi = 2;

	while ((x > xa[klo - 1]) && (khi <= n)){
		xx = x;
		if (khi < n){
			xx = fmin(x, xa[khi - 1]);
		}

		h = xa[khi - 1] - xa[klo - 1];
		a = (xa[khi - 1] - xx) / h;
		b = (xx - xa[klo - 1]) / h;
		a2 = a*a;
		b2 = b*b;
		*yi = *yi + ((1.0 - a2)*ya[klo - 1] / 2.0 + b2*ya[khi - 1] / 2.0 +
			((-(1.0 + a2*a2) / 4.0 + a2 / 2.0)*y2a[klo - 1] + (b2*b2 / 4.0 -
			b2 / 2.0)*y2a[khi - 1])*h*h / 6.0)*h;

		klo++;
		khi++;

		if (klo > 1000){
			break;
		}
	}

	return;
}


void MSIS::glatf(double lat, double *gv, double *refff)
{
	//glatf member function MSIS class
	//Calculates latitude variable gravity (gv) and effective radius (refff)

	double c2;
	c2 = cos(2.0*pi180*lat);
	*gv = 980.616*(1.0 - 0.0026373*c2);
	*refff = 2.0**gv / (3.085462e-6 + 2.27e-9*c2)*1e-5;
	return;

}


void MSIS::ghp7(int iyd, double sec, double *alt, double glat, double glong, double stl, double f107a, double f107, double *ap, double *d, double *t, double press)
{
	//ghp7 member function from MSIS class
	//Find altitude of pressure surface (press) from gtd7
	//Input:
	//	iyd - year and day as yyddd
	//	sec - ut(sec)
	//	glat - geodetic latitude(deg)
	//	glong - geodetic longitude(deg)
	//	stl - local apparent solar time(hrs)
	//	f107a - 3 month average of f10.7 flux
	//	f107 - daily f10.7 flux for previous day
	//	ap - magnetic index(daily) or when sw(9) = -1. :
	//		array containing:
	//		(1) daily ap
	//		(2) 3 hr ap index for current time
	//		(3) 3 hr ap index for 3 hrs before current time
	//		(4) 3 hr ap index for 6 hrs before current time
	//		(5) 3 hr ap index for 9 hrs before current time
	//		(6) average of eight 3 hr ap indicies from 12 to 33 hrs
	//		(7) average of eight 3 hr ap indicies from 36 to 59 hrs
	//	press - pressure level(mb)
	//Output:
	//	alt - altitude(km)
	//	d(1) - he number density(cm-3)
	//	d(2) - o number density(cm-3)
	//	d(3) - n2 number density(cm-3)
	//	d(4) - o2 number density(cm-3)
	//	d(5) - ar number density(cm-3)
	//	d(6) - total mass density(gm/cm3)
	//	d(7) - h number density(cm-3)
	//	d(8) - n numner density(cm-3)
	//	d(9) - hot o number density(cm-3)
	//	t(1) - exospheric temperature
	//	t(2) - temperature at alt

	const double bm = 1.3806e-19, rgas = 831.4, test = 0.00043;
	const int ltest = 12;
	double p, pl, zi, cl, cl2, cd, z, ca, xn, diff, xm, g, sh;
	int iday, l;

	pl = log10(press);

	//Initial altitude estimate
	if (pl >= -5.0){
		if (pl > 2.5){
			zi = 18.06*(3.0 - pl);
		}
		if (pl > 0.75 && pl <= 2.5){
			zi = 14.98*(3.08 - pl);
		}
		if (pl > -1.0 && pl <= 0.75){
			zi = 17.8*(2.72 - pl);
		}
		if (pl > -2.0 && pl <= -1.0){
			zi = 14.28*(3.64 - pl);
		}
		if (pl > -4.0 && pl <= -2.0){
			zi = 12.72*(4.32 - pl);
		}
		if (pl <= -4.0){
			zi = 25.3*(0.11 - pl);
		}

		iday = iyd % 1000;
		cl = glat / 90.0;
		cl2 = cl*cl;

		if (iday < 182){
			cd = 1.0 - iday / 91.25;
		}

		if (iday >= 182){
			cd = iday / 91.25 - 3.0;
		}

		ca = 0.0;

		if (pl > -1.11 && pl <= -0.23){
			ca = 1.0;
		}

		if (pl > -0.23){
			ca = (2.79 - pl) / (2.79 + 0.23);
		}

		if (pl <= -1.11 && pl > -3.0){
			ca = (-2.93 - pl) / (-2.93 + 1.11);
		}

		z = zi - 4.87*cl*cd*ca - 1.64*cl2*ca + 0.31*ca*cl;
	}

	if (pl < -5.0){
		z = 22.0*pow(pl + 4.0, 2) + 110.0;
	}

	//Iteration Loop
	l = 0;

L10:

	l++;

	gtd7(iyd, sec, z, glat, glong, stl, f107a, f107, *ap, 48, d, t);

	xn = d[0] + d[1] + d[2] + d[3] + d[4] + d[5] + d[6] + d[7];

	p = bm*xn*t[1];

	if (imr == 1){
		p = p*1.0e-06;
	}

	diff = pl - log10(p);

	if (fabs(diff) < test || l == ltest){
		goto L20;
	}

	xm = d[5] / xn / 1.66e-24;

	if (imr == 1){
		xm = xm*1.0e03;
	}

	g = gsurf / pow(1.0 + z / re, 2);
	sh = rgas*t[1] / (xm*g);

	//New altitude estimate using scale height
	if (l < 6){
		z = z - sh*diff*2.302;
	}
	else {
		z = z - sh*diff;
	}


	goto L10;

L20:

	*alt = z;

	return;
}


void MSIS::msishwm(int iyr, int mn, int ida, int ihr, double xmin, double z,
	double xlat, double xlon, double f10b, double f10, double *ap, double *ph,
	double *dh, double *th, double *n2nd, double *o2nd, double *ond,
	double *arnd, double *hend, double *hnd, double *nnd, double *wtmol,
	double *tex, double *uh, double *vh)
{

	//msishwm member function from MSIS class

	//Molecular weights
	const double ei[9] = { 4.0026, 15.9994, 28.0134, 31.9988, 39.948, 0.0,
		1.00797, 14.0067, 15.9994 };
	const double r0 = 8314.32;
	double pi, rlat, sda, sha, dd, dy, sra, raloc, dut, utsec, dlst, xlst, d[9],
		t[2], sumn, summn, xmjd, w[2] = { 0.0 };
	int iday, i;

	pi = 4.0*atan(1.0);

	rlat = xlat;

	maps.perts.iperts.mets.tme(mn, ida, iyr, ihr, xmin, xlon, &rlat, &sda, &sha, &dd, &dy, &sra, &raloc, &xmjd);

	iday = 1000 * iyr + int(dd + 1.0);
	dut = 86400 * (dd - int(dd));
	utsec = dut;
	dlst = 12.0*(1.0 + sha / pi);
	xlst = fmod(dlst, 24.0);

	//Call MSIS thermosphere subroutine

	gtd7(iday, utsec, z, xlat, xlon, xlst, f10b, f10, *ap, 48, d, t);

	// MSIS output variables:
	//  d[0] - He number density                (cm-3)
	//  d[1] - O number density                 (cm-3)
	//  d[2] - N2 number density                (cm-3)
	//  d[3] - O2 number density                (cm-3)
	//  d[4] - Ar number density                (cm-3)
	//  d[5] - Total mass density               (GM/cm3)
	//  d[6] - H number density                 (cm-3)
	//  d[7] - N number density                 (cm-3)
	//  d[8] - Anomalous oxygen number density  (cm-3)
	//  t[0] - exospheric temperature
	//  t[1] - temperature at alt

	sumn = 0.0;
	summn = 0.0;

	//Convert MSIS outputs to SI units
	for (i = 0; i < 9; i++){
		d[i] = 1.0e06*d[i];

		if (i != 5){
			sumn = sumn + d[i];
		}

		summn = summn + d[i] * ei[i];
	}

	//Change notation for outputs
	*dh = d[5] / 1000.0;
	*th = t[1];
	*tex = t[0];
	*wtmol = summn / sumn;
	*ph = *dh*r0*(*th) / (*wtmol);
	*hend = d[0];
	*ond = d[1] + d[8];
	*n2nd = d[2];
	*o2nd = d[3];
	*arnd = d[4];
	*hnd = d[6];
	*nnd = d[7];

	//Call MSIS Harmonic Wind Model (HWM) subroutine
	jb2008.hwm.gws5(iday, utsec, z, xlat, xlon, xlst, f10b, f10, ap, w);

	//Change notation for outputs
	*uh = w[1];
	*vh = w[0];

	return;

}

void MSIS::msismod(double h, double phi, double thet, double elt, double ri,
	double g, double *pmj, double *dmj, double *tmj, double *umj, double *vmj,
	double *wmj, double *n2nd, double *o2nd, double *ond, double *arnd,
	double *hend, double *hnd, double *wtmol, double *dmdz, double *nnd)
{
	//msismod member function from MSIS class
	//MSIS and HWM model driver routine to evaluate mean p, d, t, u, v, w

	double dy5, dx5, xmin, tex, phn, dhn, thn, d1, d2, d3, d4, d5, d6, d7, d8,
		d9, d10, d11, phe, dhe, the, dtx, dty, pb, db, tb, cp, dphi;

	//Factor for radians to degrees
	const double pi = 3.14159265, fac = 180 / pi;

	imin = mino + int(elt) / 60;
	sec = seco + elt;
	sec = int(sec) % 60;
	ihr = ihro + imin / 60;
	imin = imin % 60;

	//Distances for 5 degrees of latitude, longitude
	dy5 = 5000.0*ri / fac;
	dx5 = dy5*cos(phi / fac);

	if (dx5 < 2000.0){
		dx5 = 2000.0;
	}

	//Following is the pure thermosphere height range section
	//MSIS and HWM values at current position
	double	ap[7] = { 20.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

	xmin = imin + sec / 60.0;
	
	msishwm(iyr, mn, ida, ihr, xmin, h, phi, thet, f10b, f10, ap, pmj, dmj,
		tmj, n2nd, o2nd, ond, arnd, hend, hnd, nnd, wtmol, &tex, umj, vmj);
	
	//Set latitude increment for temperature gradients 
	dphi = 5.0;

	if (phi >= 85.0){
		dphi = -5.0;
	}

	//MSIS temperature at current position latitude increment
	msishwm(iyr, mn, ida, ihr, xmin, h, phi + dphi, thet, f10b, f10,
		ap, &phn, &dhn, &thn, &d1, &d2, &d3, &d4, &d5, &d6, &d7, &d8, &d9,
		&d10, &d11);

	//MSIS temperature at current position+5 degrees lon
	msishwm(iyr, mn, ida, ihr, xmin, h, phi, thet + 5.0, f10b, f10,
		ap, &phe, &dhe, &the, &d1, &d2, &d3, &d4, &d5, &d6, &d7, &d8, &d9,
		&d10, &d11);

	//dt/dx, dt/dy, and dt/dz for vertical wind
	dtx = the - *tmj;
	dty = thn - *tmj;

	if (dphi < 0.0){
		dty = -dty;
	}

	//MSIS temperature and molecular weight 1 km higher
	msishwm(iyr, mn, ida, ihr, xmin, h + 1.0, phi, thet, f10b, f10,
		ap, &pb, &db, &tb, &d1, &d2, &d3, &d4, &d5, &d6, &d7, &d8, &d9,
		&d10, &d11);

	//Gradients for temperature and molecular weight
	dtz = (tb - *tmj) / 1000.0;
	*dmdz = (d8 - *wtmol) / 1000.0;

	//Compute vertical mean wind
	//Specific heat
	cp = 7.0*(*pmj) / (2.0*(*dmj)*(*tmj));

	//Mean vertical wind from Montgomery stream function
	*wmj = -cp*(*umj*dtx / dx5 + *vmj*dty / dy5) / (g + cp*dtz);

	return;

}

void MSIS::gtd7(int iyd, double sec, double alt, double glat, double glong, double stl, double f107a, double f107, double ap, double mass, double *d, double *t)
{
	//gtd7 member function from MSIS class
	//NRLMSISE-00
	//Neutral Atmosphere Empirical Model fro the surface to lower exosphere


	double alast = 99999.0, zn2[4] = { 72.5, 55.0, 45.0, 32.5 }, zn3[5] = { 32.5,
		20.0, 15.0, 10.0, 0.0 }, sv[25] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
		1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
		1.0, 1.0 };
	double v1, xlat, xmm, altt, ds[9], ts[2], dm28m, dmc, dz28, dmr, tz;
	int mssl = -999, j, mss;
	const int mn2 = 4, mn3 = 5;
	const double zmix = 62.5;

	//Test for changed input
	v1 = vtst7(iyd, sec, glat, glong, stl, f107a, f107, &ap, 1);

	//Latitude variation of gravity (none for hwm.getCswSw(1) = 0.0
	xlat = glat;

	if (jb2008.hwm.getCswSw(1) == 0.0){
		xlat = 45.0;
	}

	glatf(xlat, &gsurf, &re);

	xmm = pdm[4][2];

	//Thermosphere/mesosphere (above zn2[0])
	altt = fmax(alt, zn2[0]);
	mss = mass;

	//Only calculate thermosphere if input parameters changed or altitude 
	//above zn2[0] in mesosphere
	if (v1 == 1.0 || alt > zn2[0] || alast > zn2[0] || mss != mssl){
		gts7(iyd, sec, altt, glat, glong, stl, f107a, f107, &ap, mss, &dm28m, ds, ts);

		//Metric adjustment
		if (imr == 1){
			dm28m = dm28m*1.0e06;
		}

		mssl = mss;
	}

	t[0] = ts[0];
	t[1] = ts[1];

	if (alt >= zn2[0]){

		for (j = 0; j < 9; j++){
			d[j] = ds[j];
		}

		goto L10;
	}

	//The rest of the code in this method never gets called.
	//Only calculate nodes if input changed
	if (v1 == 1.0 || alast >= zn2[0]) {
		tgn2[0] = tgn1[1];
		tn2[0] = tn1[4];
		tn2[1] = pma[0][0] * pavgm[0] / (1.0 - jb2008.hwm.getCswSw(19)*glob7s(&pma[0][0]));
		tn2[2] = pma[1][0] * pavgm[1] / (1.0 - jb2008.hwm.getCswSw(19)*glob7s(&pma[1][0]));
		tn2[3] = pma[2][0] * pavgm[2] / (1.0 - jb2008.hwm.getCswSw(19)*jb2008.hwm.getCswSw(21)
			*glob7s(&pma[2][0]));
		tgn2[1] = pavgm[8] * pma[9][0] * (1.0 + jb2008.hwm.getCswSw(19)*jb2008.hwm.getCswSw(21)
			*glob7s(&pma[9][0]))*tn2[3] * tn2[3] / pow(pma[2][0] * pavgm[2], 2);
		tn3[0] = tn2[3];
	}

	if (alt >= zn3[0]) {
		goto L6;
	}

	//LOWER STRATOSPHERE AND TROPOSPHERE [below ZN3(1)]
	//Temperature at nodes and gradients at end nodes
	//Inverse temperature a linear function of spherical harmonics
	//Only calculate nodes if input changed
	if (v1 == 1.0 || alast >= zn3[0]) {
		tgn3[0] = tgn2[1];
		tn3[1] = pma[3][0] * pavgm[3] / (1.0 - jb2008.hwm.getCswSw(21)*glob7s(&pma[3][0]));
		tn3[2] = pma[4][0] * pavgm[4] / (1.0 - jb2008.hwm.getCswSw(21)*glob7s(&pma[4][0]));
		tn3[3] = pma[5][0] * pavgm[5] / (1.0 - jb2008.hwm.getCswSw(21)*glob7s(&pma[5][0]));
		tn3[4] = pma[6][0] * pavgm[6] / (1.0 - jb2008.hwm.getCswSw(21)*glob7s(&pma[6][0]));
		tgn3[1] = pma[7][0] * pavgm[7] * (1.0 + jb2008.hwm.getCswSw(21)*glob7s(&pma[7][0]))
			*tn3[4] * tn3[4] / pow(pma[6][0] * pavgm[6], 2);
	}

L6:
	if (mass == 0) {
		goto L50;
	}

	// Linear transition to full mixing below zn2[0]
	dmc = 0.0;

	if (alt > zmix) {
		dmc = 1.0 - (zn2[0] - alt) / (zn2[0] - zmix);
	}

	dz28 = ds[2];

	// N2 Density
	dmr = ds[2] / dm28m - 1.0;
	d[2] = densm(alt, dm28m, xmm, &tz, mn3, zn3, tn3, tgn3, mn2, zn2, tn2, tgn2);
	d[2] = d[2] * (1.0 + dmr*dmc);

	// He Density
	d[0] = 0.0;

	if (mass != 4 && mass != 48) {
		goto L204;
	}

	dmr = ds[0] / (dz28*pdm[1][0]) - 1.0;
	d[0] = d[2] * pdm[1][0] * (1.0 + dmr*dmc);

L204:
	// O Density
	d[1] = 0.0;
	d[8] = 0.0;

	// O2 Density
	d[3] = 0.0;

	if (mass != 32 && mass != 48) {
		goto L232;
	}

	dmr = ds[3] / (dz28*pdm[1][3]) - 1.0;
	d[3] = d[2] * pdm[1][3] * (1.0 + dmr*dmc);

L232:
	// Ar Density
	d[4] = 0.0;

	if (mass != 40 && mass != 48) {
		goto L240;
	}

	dmr = ds[4] / (dz28*pdm[1][4]) - 1.0;
	d[4] = d[2] * pdm[1][4] * (1.0 + dmr*dmc);

L240:
	// H Density
	d[6] = 0.0;

	// N Density
	d[7] = 0.0;

	// Total Mass Density
	if (mass == 48) {
		d[5] = 1.66e-24*(4.0*d[0] + 16.0*d[1] + 28.0*d[2]
			+ 32.0*d[3] + 40.0*d[4] + d[6] + 14.0*d[7]);
	}

	if (imr == 1) {
		d[5] = d[5] / 1000.0;
	}

	t[1] = tz;

L10:
	goto L90;

L50:
	dd = densm(alt, 1.0, 0.0, &tz, mn3, zn3, tn3, tgn3, mn2, zn2, tn2, tgn2);

	t[1] = tz;

L90:
	alast = alt;

	return;
}

void MSIS::gts7(int iyd, double sec, double alt, double glat, double glong, double stl, double f107a, double f107, double *ap, int mass, double *dm28, double *d, double *t)
{

	//gts7 member function from MSIS class
	//Thermospheric portion of NRLMSISE-00
	//See GTD7 for more extensive comments

	double alast = -999.0, zn1[5] = { 120.0, 110.0, 100.0, 90.0, 72.5 }, v2, tinf,
		g28, day, zhf, z, zh28, zhm28, b28, xmd, tz, xmm, g4, zh04, b04, zhm04, zc04,
		hc04, g16, zh16, b16, zhm16, hc16, zc16, hc216, hcc16, zcc16, rc16, g32, zh32,
		b32, zhm32, hc32, zc32, hcc32, hcc232, zcc32, rc32, g40, zh40, b40, zhm40,
		hc40, zc40, g1, zh01, b01, zhm01, hc01, zc01, hcc01, zcc01, rc01, g14, zh14,
		b14, zhm14, hc14, zc14, hcc14, zcc14, rc14, g16h, db16h, tho, t2, zsht, zhmo,
		zsho, ddum, dm04, dm16, dm32, dm40, dm01, dm14, tlb, s, db04, db16, db28, db32,
		db40, db48, db01, za, t0, z0, g0, rl, db14, tr12;
	int i, j, iyr, lenyr, yrd;
	const double altl[8] = { 200.0, 300.0, 160.0, 250.0, 240.0, 450.0, 320.0, 450.0 },
		dgtr = 1.74533e-02, dr = 1.72142e-02, alpha[9] = { -0.38, 0.0, 0.0, 0.0, 0.17,
		0.0, -0.38, 0.0, 0.0 };
	const int mn1 = 5, mt[11] = { 48, 0, 4, 16, 28, 32, 40, 1, 49, 14, 17 };


	glatf(glat, &gsurf, &re);

	//Test for changed input
	v2 = vtst7(iyd, sec, glat, glong, stl, f107a, f107, ap, 2);

	yrd = iyd;

	za = pdl[15][1];
	zn1[0] = za;

	for (j = 0; j < 9; j++){
		d[j] = 0.0;
	}

	//Tinf variations not important below za or zn1[0]
	if (alt > zn1[0]){

		if ((v2 == 1.0) || (alast <= zn1[0])){
			tinf = ptm[0] * pt[0] * (1.0 + jb2008.hwm.getCswSw(15)*globe7(yrd, sec, glat, glong,
				stl, f107a, f107, ap, pt));
		}
	}
	else {
		tinf = ptm[0] * pt[0];
	}

	t[0] = tinf;

	//Gradient variations not impontant below zn1[4]
	if (alt > zn1[4]){
		if ((v2 == 1.0) || (alast <= zn1[4])){
			g0 = ptm[3] * ps[0] * (1.0 + jb2008.hwm.getCswSw(18)
				*globe7(yrd, sec, glat, glong, stl, f107a, f107, ap, ps));
		}
	}
	else {
		g0 = ptm[3] * ps[0];
	}

	//Calculate these temperatures only if input changed
	if ((v2 == 1.0) || (alt < 300.0)) {
		tlb = ptm[1] * (1.0 + jb2008.hwm.getCswSw(16)
			*globe7(yrd, sec, glat, glong, stl, f107a, f107, ap, &pd[3][0]))*pd[3][0];
	}
	s = g0 / (tinf - tlb);

	//Lower themosphere temp variations not significant for density above 300 km
	if (alt < 300.0){
		if (v2 == 1.0 || alast >= 300.0) {
			tn1[1] = ptm[6] * ptl[0][0] / (1.0 - jb2008.hwm.getCswSw(17)*glob7s(&ptl[0][0]));
			tn1[2] = ptm[2] * ptl[1][0] / (1.0 - jb2008.hwm.getCswSw(17)*glob7s(&ptl[1][0]));
			tn1[3] = ptm[7] * ptl[2][0] / (1.0 - jb2008.hwm.getCswSw(17)*glob7s(&ptl[2][0]));
			tn1[4] = ptm[4] * ptl[3][0] / (1.0 - jb2008.hwm.getCswSw(17)*jb2008.hwm.getCswSw(19)
				*glob7s(&ptl[3][0]));
			tgn1[1] = ptm[8] * pma[8][0] * (1.0 + jb2008.hwm.getCswSw(17)*jb2008.hwm.getCswSw(19)
				*glob7s(&pma[8][0]))*tn1[4] * tn1[4] / pow(ptm[4] * ptl[3][0], 2);
		}
	}
	else {
		tn1[1] = ptm[6] * ptl[0][0];
		tn1[2] = ptm[2] * ptl[1][0];
		tn1[3] = ptm[7] * ptl[2][0];
		tn1[4] = ptm[4] * ptl[3][0];
		tgn1[1] = ptm[8] * pma[8][0] * tn1[4] * tn1[4] / pow(ptm[4] * ptl[3][0], 2);
	}

	z0 = zn1[3];
	t0 = tn1[3];
	tr12 = 1.0;

	if (mass == 0){
		goto L50;
	}

	//N2 variation factor at Zlb
	g28 = jb2008.hwm.getCswSw(20)*globe7(yrd, sec, glat, glong, stl, f107a, f107, ap, &pd[2][0]);
	day = fmod(yrd, 1000.0);
	iyr = int(yrd) / 1000;

	//Adjust day-of-year dependency for leap year vs non-leap year
	lenyr = 365;

	if (iyr % 4 == 0){
		lenyr = 366;
	}

	day = (day - 1.0 + sec / 86400.0)*(365.0 / lenyr);

	//Variation of turbopause height
	zhf = pdl[24][1] * (1.0 + jb2008.hwm.getCswSw(4)*pdl[24][0] * sin(dgtr*glat)*
		cos(dr*(day - pt[13])));
	yrd = iyd;
	t[0] = tinf;
	xmm = pdm[4][2];
	z = alt;

	for (j = 0; j < 11; j++){

		if (mass == mt[j]) {
			goto L15;
		}

	}

	goto L90;

L15:
	if ((z > altl[5]) && (mass != 28) && (mass != 48)){
		goto L17;
	}

	//N2 Density
	//Diffusive density at zlb
	db28 = pdm[0][2] * exp(g28)*pd[2][0];

	//Diffusive density at alt
	d[2] = densu(z, db28, tinf, tlb, 28.0, alpha[2], &t[1], ptm[5], s, mn1,
		zn1, tn1, tgn1);

	dd = d[2];

	//Turbopause
	zh28 = pdm[2][2] * zhf;
	zhm28 = pdm[3][2] * pdl[5][1];
	xmd = 28.0 - xmm;

	//Mixed density at zlb
	b28 = densu(zh28, db28, tinf, tlb, xmd, alpha[2] - 1.0, &tz,
		ptm[5], s, mn1, zn1, tn1, tgn1);

	if ((z > altl[2]) || (jb2008.hwm.getCswSw(14) == 0.0)){
		goto L17;
	}

	//Mixed density at alt
	*dm28 = densu(z, b28, tinf, tlb, xmm, alpha[2], &tz, ptm[5], s, mn1, zn1,
		tn1, tgn1);

	//Net density at alt
	d[2] = dnet(d[2], *dm28, zhm28, xmm, 28.0);

L17:
	switch (j) {
	case 0:
		goto L20;
	case 1:
		goto L50;
	case 2:
		goto L20;
	case 3:
		goto L25;
	case 4:
		goto L90;
	case 5:
		goto L35;
	case 6:
		goto L40;
	case 7:
		goto L45;
	case 8:
		goto L25;
	case 9:
		goto L48;
	case 10:
		goto L46;
	default:
		break;
	}

L20:
	//He Density
	//Density variation factor at zlb
	g4 = jb2008.hwm.getCswSw(20)*globe7(yrd, sec, glat, glong, stl, f107a, f107, ap, &pd[0][0]);

	//Diffusive density at zlb
	db04 = pdm[0][0] * exp(g4)*pd[0][0];

	//Diffusive density at alt
	d[0] = densu(z, db04, tinf, tlb, 4.0, alpha[0], &t[1], ptm[5], s, mn1, zn1, tn1,
		tgn1);
	dd = d[0];

	if ((z > altl[0]) || (jb2008.hwm.getCswSw(14) == 0.0)){
		goto L24;
	}

	//Turbopause
	zh04 = pdm[2][0];

	//Mixed density at zlb
	b04 = densu(zh04, db04, tinf, tlb, 4.0 - xmm, alpha[0] - 1.0, &t[1], ptm[5], s, mn1,
		zn1, tn1, tgn1);

	//Mixed density at alt
	dm04 = densu(z, b04, tinf, tlb, xmm, 0.0, &t[1], ptm[5], s, mn1, zn1, tn1, tgn1);
	zhm04 = zhm28;
	//Net density at alt
	d[0] = dnet(d[0], dm04, zhm04, xmm, 4.0);

	//Correction to specified mixing ratio at ground
	rl = log(b28*pdm[1][0] / b04);
	zc04 = pdm[4][0] * pdl[0][1];
	hc04 = pdm[5][0] * pdl[1][1];

	//Net density corrected at alt
	d[0] = d[0] * ccor(z, rl, hc04, zc04);

L24:
	if (mass != 48){
		goto L90;
	}

L25:
	//O Density

	g16 = jb2008.hwm.getCswSw(20)*globe7(yrd, sec, glat, glong, stl, f107a, f107, ap, &pd[1][0]);

	//Diffusive density at zlb
	db16 = pdm[0][1] * exp(g16)*pd[1][0];

	//Diffusive density at Alt
	d[1] = densu(z, db16, tinf, tlb, 16.0, alpha[1], &t[1], ptm[5], s, mn1, zn1, tn1,
		tgn1);

	dd = d[1];

	if (z > altl[1] || jb2008.hwm.getCswSw(14) == 0.0){
		goto L34;
	}

	//Corrected from low.pdm[2][0] to low.pdm[2][1] 12/2/85
	//Turbopause
	zh16 = pdm[2][1];

	//Mixed density at zlb
	b16 = densu(zh16, db16, tinf, tlb, 16.0 - xmm, alpha[1] - 1.0, &t[1], ptm[5], s,
		mn1, zn1, tn1, tgn1);

	//Mixed density at alt
	dm16 = densu(z, b16, tinf, tlb, xmm, 0.0, &t[1], ptm[5], s, mn1, zn1, tn1, tgn1);
	zhm16 = zhm28;

	//Net density at alt
	d[1] = dnet(d[1], dm16, zhm16, xmm, 16.0);

	// 3/16/99 change form to match O2 departure from diff equil near 150 km 
	//and add dependence on F10.7
	rl = pdm[1][1] * pdl[16][1] * (1.0 + jb2008.hwm.getCswSw(0)*pdl[23][0] *
		(f107a - 150.0));
	hc16 = pdm[5][1] * pdl[3][1];
	zc16 = pdm[4][1] * pdl[2][1];
	hc216 = pdm[5][1] * pdl[4][1];
	d[1] = d[1] * ccor2(z, rl, hc16, zc16, hc216);

	//Chemistry correction
	hcc16 = pdm[7][1] * pdl[13][1];
	zcc16 = pdm[6][1] * pdl[12][1];
	rc16 = pdm[3][1] * pdl[14][1];

	//Net density corrected at alt
	d[1] = d[1] * ccor(z, rc16, hcc16, zcc16);

L34:
	if ((mass != 48) && (mass != 49)){
		goto L90;
	}

L35:
	//O2 Density
	//Density variation factor at zlb
	g32 = jb2008.hwm.getCswSw(20)*globe7(yrd, sec, glat, glong, stl, f107a,
		f107, ap, &pd[4][0]);

	//Diffusive density at zlb
	db32 = pdm[0][3] * exp(g32)*pd[4][0];

	//Diffusive density at alt
	d[3] = densu(z, db32, tinf, tlb, 32.0, alpha[3], &t[1], ptm[5], s, mn1,
		zn1, tn1, tgn1);

	if (mass == 49){
		dd = dd + 2.0*d[3];
	}
	else{
		dd = d[3];
	}

	if (jb2008.hwm.getCswSw(14) == 0.0){
		goto L39;
	}

	if (z > altl[3]){
		goto L38;
	}

	//Turbopause
	zh32 = pdm[2][3];

	//Mixed density at zlb
	b32 = densu(zh32, db32, tinf, tlb, 32.0 - xmm, alpha[3] - 1.0, &t[1],
		ptm[5], s, mn1, zn1, tn1, tgn1);

	//Mixed density at alt
	dm32 = densu(z, b32, tinf, tlb, xmm, 0.0, &t[1], ptm[5], s, mn1, zn1, tn1,
		tgn1);

	zhm32 = zhm28;

	//Net density at Alt
	d[3] = dnet(d[3], dm32, zhm32, xmm, 32.0);

	//Correction to specified mixing ratio at ground
	rl = log(b28*pdm[1][3] / b32);
	hc32 = pdm[5][3] * pdl[7][1];
	zc32 = pdm[4][3] * pdl[6][1];
	d[3] = d[3] * ccor(z, rl, hc32, zc32);

L38:
	//Correction for general deprature from diffusive equilibrium above zlb
	hcc32 = pdm[7][3] * pdl[22][1];
	hcc232 = pdm[7][3] * pdl[22][0];
	zcc32 = pdm[6][3] * pdl[21][1];
	rc32 = pdm[3][3] * pdl[23][1] * (1.0 + jb2008.hwm.getCswSw(0)*pdl[23][0] *
		(f107a - 150.0));

	//Net density corrected at alt
	d[3] = d[3] * ccor2(z, rc32, hcc32, zcc32, hcc232);

L39:
	if (mass != 48){
		goto L90;
	}

L40:
	//AR Density
	//Density variation factor at zlb
	g40 = jb2008.hwm.getCswSw(20)*globe7(yrd, sec, glat, glong, stl, f107a,
		f107, ap, &pd[5][0]);

	//Diffusive density at zlb
	db40 = pdm[0][4] * exp(g40)*pd[5][0];

	//Diffusive density at alt
	d[4] = densu(z, db40, tinf, tlb, 40.0, alpha[4], &t[1], ptm[5], s, mn1,
		zn1, tn1, tgn1);
	dd = d[4];

	if ((z > altl[4]) || (jb2008.hwm.getCswSw(14) == 0.0)){
		goto L44;
	}

	//Turbopause
	zh40 = pdm[2][4];

	//Mixed density at zlb
	b40 = densu(zh40, db40, tinf, tlb, 40.0 - xmm, alpha[4] - 1.0, &t[1], ptm[5], s,
		mn1, zn1, tn1, tgn1);

	//Mixed density at alt
	dm40 = densu(z, b40, tinf, tlb, xmm, 0.0, &t[1], ptm[5], s, mn1, zn1, tn1, tgn1);
	zhm40 = zhm28;

	//Net density at alt
	d[4] = dnet(d[4], dm40, zhm40, xmm, 40.0);

	//Correction to specified mixing ratio at ground
	rl = log(b28*pdm[1][4] / b40);
	hc40 = pdm[5][4] * pdl[9][1];
	zc40 = pdm[4][4] * pdl[8][1];

	//Net density corrected at alt
	d[4] = d[4] * ccor(z, rl, hc40, zc40);

L44:
	if (mass != 48){
		goto L90;
	}

L45:
	//Hydrogen Density
	//Density variation factor at zlb

	g1 = jb2008.hwm.getCswSw(20)*globe7(yrd, sec, glat, glong, stl, f107a,
		f107, ap, &pd[6][0]);

	//Diffusive density at zlb
	db01 = pdm[0][5] * exp(g1)*pd[6][0];

	//Diffusive density at alt
	d[6] = densu(z, db01, tinf, tlb, 1.0, alpha[6], &t[1], ptm[5], s, mn1,
		zn1, tn1, tgn1);
	dd = d[6];

	if ((z > altl[6]) || (jb2008.hwm.getCswSw(14) == 0.0)){
		goto L47;
	}

	//Turbopause
	zh01 = pdm[2][5];

	//Mixed density at zlb
	b01 = densu(zh01, db01, tinf, tlb, 1.0 - xmm, alpha[6] - 1.0, &t[1],
		ptm[5], s, mn1, zn1, tn1, tgn1);

	//Mixed density at alt
	dm01 = densu(z, b01, tinf, tlb, xmm, 0.0, &t[1], ptm[5], s, mn1, zn1, tn1,
		tgn1);
	zhm01 = zhm28;
	//Net density at alt
	d[6] = dnet(d[6], dm01, zhm01, xmm, 1.0);

	//Correction to specified mixing ratio at ground
	rl = log(b28*pdm[1][5] * fabs(pdl[17][1]) / b01);
	hc01 = pdm[5][5] * pdl[11][1];
	zc01 = pdm[4][5] * pdl[10][1];
	d[6] = d[6] * ccor(z, rl, hc01, zc01);

	//Chemistry correction
	hcc01 = pdm[7][5] * pdl[19][1];
	zcc01 = pdm[6][5] * pdl[18][1];
	rc01 = pdm[3][5] * pdl[20][1];

	//Net density corrected at alt
	d[6] = d[6] * ccor(z, rc01, hcc01, zcc01);

L47:
	if (mass != 48){
		goto L90;
	}

L48:

	//Atomic Nitrogen Density
	//Density variation at zlb
	g14 = jb2008.hwm.getCswSw(20)*globe7(yrd, sec, glat, glong, stl, f107a,
		f107, ap, &pd[7][0]);

	//Diffusive density at zlb
	db14 = pdm[0][6] * exp(g14)*pd[7][0];

	//Diffusive density at alt
	d[7] = densu(z, db14, tinf, tlb, 14.0, alpha[7], &t[1], ptm[5], s, mn1,
		zn1, tn1, tgn1);
	dd = d[7];

	if ((z > altl[7]) || (jb2008.hwm.getCswSw(14) == 0.0)){
		goto L49;
	}

	//Turbopause
	zh14 = pdm[2][6];

	//Mixed density at zlb
	b14 = densu(zh14, db14, tinf, tlb, 14.0 - xmm, alpha[7] - 1.0, &t[1],
		ptm[5], s, mn1, zn1, tn1, tgn1);

	//Mixed density at alt
	dm14 = densu(z, b14, tinf, tlb, xmm, 0.0, &t[1], ptm[5], s, mn1, zn1, tn1,
		tgn1);
	zhm14 = zhm28;
	//Net density at alt
	d[7] = dnet(d[7], dm14, zhm14, xmm, 14.0);

	//Correction to specified mixing ratio at ground
	rl = log(b28*pdm[1][6] * fabs(pdl[2][0]) / b14);
	hc14 = pdm[5][6] * pdl[1][0];
	zc14 = pdm[4][6] * pdl[0][0];
	d[7] = d[7] * ccor(z, rl, hc14, zc14);

	//Chemistry correction
	hcc14 = pdm[7][6] * pdl[4][0];
	zcc14 = pdm[6][6] * pdl[3][0];
	rc14 = pdm[3][6] * pdl[5][0];

	//Net density corrected at alt
	d[7] = d[7] * ccor(z, rc14, hcc14, zcc14);

L49:
	if (mass != 48){
		goto L90;
	}

L46:

	//Anomalous Oxygen Density
	g16h = jb2008.hwm.getCswSw(20)*globe7(yrd, sec, glat, glong, stl, f107a,
		f107, ap, &pd[8][0]);
	db16h = pdm[0][7] * exp(g16h)*pd[8][0];
	tho = pdm[9][7] * pdl[6][0];
	dd = densu(z, db16h, tho, tho, 16.0, alpha[8], &t2, ptm[5], s, mn1, zn1,
		tn1, tgn1);
	zsht = pdm[5][7];
	zhmo = pdm[4][7];
	zsho = scalh(zhmo, 16.0, tho);
	d[8] = dd*exp(-zsht / zsho*(exp(-(z - zhmo) / zsht) - 1.0));

	if (mass != 48){
		goto L90;
	}

	//Total mass Density
	d[5] = 1.66e-24*(4.0*d[0] + 16.0*d[1] + 28.0*d[2] + 32.0*d[3] + 40.0*d[4] + d[6]
		+ 14.0*d[7]);

	db48 = 1.66e-24*(4.0*db04 + 16.0*db16 + 28.0*db28 + 32.0*db32 + 40.0*db40 + db01
		+ 14.0*db14);

	goto L90;

	//Temperature at altitude
L50:
	z = fabs(alt);
	ddum = densu(z, 1.0, tinf, tlb, 0.0, 0.0, &t[1], ptm[5], s, mn1, zn1, tn1, tgn1);

L90:
	//Adjust densities from cgs to kgm
	if (imr == 1){
		for (i = 0; i < 9; i++){
			d[i] = d[i] * 1.0e+06;
		}
		d[5] = d[5] / 1000.0;
	}

	alast = alt;

	return;
}

double MSIS::scalh(double alt, double xm, double temp)
{
	//scalh member function from MSIS class
	//Calculate scale height (km)

	const double rgas = 831.4;
	double g, scalh_out;

	g = gsurf / pow(1.0 + alt / re, 2);

	scalh_out = rgas*temp / (g*xm);

	return scalh_out;

}


void MSIS::namelist(string namef)
{

	//namelist member function from MSIS class
	//Open and read namelist input file for setting model parameters.

	ifstream namelist;
	string dummy, NCEPpath, NCEPmn, atmpath, trapath, prtpath, nprpath,
		conpath, rndpath, rrapath, rralist, profile, NCEPpath1;
	double h1, phi1, thet1, s10, s10b, xm10, xm10b, y10,
		y10b, dstdtc, dphi, dthet, dhgt, rpscale, ruscale, rwscale,
		rdinit, rtinit, ruinit, rvinit, rwinit, sitelim, sitenear, patchy,
		z0in;
	int iopt, NCEPyr, NCEPhr, nr1, iurra, iyrrra, initpert, itherm,
		ibltest, nmax, iaux, mc;


	string mn1 = "12";
	NCEPmn = "Nb9008" + mn1 + ".bin";


	namelist.open(namef.c_str());

	namelist >> dummy;

	namelist >> dummy >> dummy >> atmpath;

	namelist >> dummy >> dummy >> NCEPpath;

	namelist >> dummy >> dummy >> trapath;

	namelist >> dummy >> dummy >> prtpath;

	namelist >> dummy >> dummy >> nprpath;

	namelist >> dummy >> dummy >> conpath;

	namelist >> dummy >> dummy >> rrapath;

	namelist >> dummy >> dummy >> rralist;

	namelist >> dummy >> dummy >> profile;

	namelist >> dummy >> dummy >> h1;

	namelist >> dummy >> dummy >> phi1;

	namelist >> dummy >> dummy >> thet1;

	namelist >> dummy >> dummy >> f10;

	namelist >> dummy >> dummy >> f10b;

	namelist >> dummy >> dummy >> ap;

	namelist >> dummy >> dummy >> s10;

	namelist >> dummy >> dummy >> s10b;

	namelist >> dummy >> dummy >> xm10;

	namelist >> dummy >> dummy >> xm10b;

	namelist >> dummy >> dummy >> y10;

	namelist >> dummy >> dummy >> y10b;

	namelist >> dummy >> dummy >> dstdtc;

	namelist >> dummy >> dummy >> mn;

	namelist >> dummy >> dummy >> ida;

	namelist >> dummy >> dummy >> iyr;

	namelist >> dummy >> dummy >> ihro;

	namelist >> dummy >> dummy >> mino;

	namelist >> dummy >> dummy >> seco;

	namelist >> dummy >> dummy >> dphi;

	namelist >> dummy >> dummy >> dthet;

	namelist >> dummy >> dummy >> dhgt;

	namelist >> dummy >> dummy >> nmax;

	namelist >> dummy >> dummy >> delt;

	namelist >> dummy >> dummy >> iopt;

	namelist >> dummy >> dummy >> iaux;

	namelist >> dummy >> dummy >> NCEPyr;

	namelist >> dummy >> dummy >> NCEPhr;

	namelist >> dummy >> dummy >> nr1;

	namelist >> dummy >> dummy >> mc;

	namelist >> dummy >> dummy >> rpscale;

	namelist >> dummy >> dummy >> ruscale;

	namelist >> dummy >> dummy >> rwscale;

	namelist >> dummy >> dummy >> iurra;

	namelist >> dummy >> dummy >> iyrrra;

	namelist >> dummy >> dummy >> sitelim;

	namelist >> dummy >> dummy >> sitenear;

	namelist >> dummy >> dummy >> initpert;

	namelist >> dummy >> dummy >> rdinit;

	namelist >> dummy >> dummy >> rtinit;

	namelist >> dummy >> dummy >> ruinit;

	namelist >> dummy >> dummy >> rvinit;

	namelist >> dummy >> dummy >> rwinit;

	namelist >> dummy >> dummy >> patchy;

	namelist >> dummy >> dummy >> itherm;

	namelist >> dummy >> dummy >> z0in;

	namelist >> dummy >> dummy >> ibltest;

	namelist.close();


	NCEPpath1 = NCEPpath + NCEPmn;

	imin = mino;
	ihr = ihro;
	sec = seco;

}

double MSIS::globe7(double yrd, double sec, double lat, double lon, double tloc, double f107a, double f107, double *ap, double *p)
{

	//globe7 member function from MSIS class
	//Calculate g(l) function
	//Upper Thermosphere Parameters

	double globe7_out, c, c2, c4, s, s2, cd14, cd18, cd32, cd39, f1, f2, t71, t72,
		t81, t82, p44, p45, exp1, xl = 1000.0, tll = 1000.0, sw9 = 1.0, dayl = -1.0,
		p14 = -1000.0, p18 = -1000.0, p32 = -1000.0, p39 = -1000.0;
	const double dgtr = 1.74533e-02, dr = 1.72142e-02, hr = 0.2618, sr = 7.2722e-05;
	double sv[25] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
		1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
	double apd, df, t[15], tinf;
	const int nsw = 14;
	int i, j, lenyr, iyr;



	if (jb2008.hwm.getCswIsw() != 64999){
		jb2008.hwm.tselec(sv);
	}

	for (j = 0; j < 14; j++){
		t[j] = 0.0;
	}

	if (jb2008.hwm.getCswSw(8) > 0){
		sw9 = 1.0;
	}

	if (jb2008.hwm.getCswSw(8) < 0){
		sw9 = -1.0;
	}

	iyr = int(yrd / 1000.0);
	day = yrd - iyr*1000.0;

	//Adjust day-of-year dependency of coefficients for leap year vs non-leap year
	lenyr = 365;
	if (iyr % 4 == 0){
		lenyr = 366;
	}

	day = (day - 1.0 + sec / 86400.0)*(365.0 / lenyr);

	xlong = lon;

	//Eq. A22 (remainder of code)
	if (xl != lat){
		//Calculate Legendre polynomials
		c = sin(lat*dgtr);
		s = cos(lat*dgtr);
		c2 = c*c;
		c4 = c2*c2;
		s2 = s*s;
		plg[1][0] = c;
		plg[2][0] = 0.5*(3.0*c2 - 1.0);
		plg[3][0] = 0.5*(5.0*c*c2 - 3.0*c);
		plg[4][0] = (35.0*c4 - 30.0*c2 + 3.0) / 8.0;
		plg[5][0] = (63.0*c2*c2*c - 70.0*c2*c + 15.0*c) / 8.0;
		plg[6][0] = (11.0*c*plg[5][0] - 5.0*plg[4][0]) / 6.0;
		plg[1][1] = s;
		plg[2][1] = 3.0*c*s;
		plg[3][1] = 1.5*(5.0*c2 - 1.0)*s;
		plg[4][1] = 2.5*(7.0*c2*c - 3.0*c)*s;
		plg[5][1] = 1.875*(21.0*c4 - 14.0*c2 + 1.0)*s;
		plg[6][1] = (11.0*c*plg[5][1] - 6.0*plg[4][1]) / 5.0;
		plg[2][2] = 3.0*s2;
		plg[3][2] = 15.0*s2*c;
		plg[4][2] = 7.5*(7.0*c2 - 1.0)*s2;
		plg[5][2] = 3.0*c*plg[4][2] - 2.0*plg[3][2];
		plg[6][2] = (11.0*c*plg[5][2] - 7.0*plg[4][2]) / 4.0;
		plg[7][2] = (13.0*c*plg[6][2] - 8.0*plg[5][2]) / 5.0;
		plg[3][3] = 15.0*s2*s;
		plg[4][3] = 105.0*s2*s*c;
		plg[5][3] = (9.0*c*plg[4][3] - 7.0*plg[3][3]) / 2.0;
		plg[6][3] = (11.0*c*plg[5][3] - 8.0*plg[4][3]) / 3.0;
		xl = lat;
	}

	if ((tll != tloc) || (jb2008.hwm.getCswSw(6) != 0.0) && (jb2008.hwm.getCswSw(7) != 0.0) &&
		(jb2008.hwm.getCswSw(13) != 0.0)){
		stloc = sin(hr*tloc);
		ctloc = cos(hr*tloc);
		s2tloc = sin(2.0*hr*tloc);
		c2tloc = cos(2.0*hr*tloc);
		s3tloc = sin(3.0*hr*tloc);
		c3tloc = cos(3.0*hr*tloc);
		tll = tloc;
	}

	if ((day != dayl) || (*(p + 13) != p14)){
		cd14 = cos(dr*(day - *(p + 13)));
	}

	if ((day != dayl) || (*(p + 17) != p18)){
		cd18 = cos(2.0*dr*(day - *(p + 17)));
	}

	if ((day != dayl) || (*(p + 31) != p32)){
		cd32 = cos(dr*(day - *(p + 31)));
	}

	if ((day != dayl) || (*(p + 38) != p39)){
		cd39 = cos(2.0*dr*(day - *(p + 38)));
	}

	dayl = day;
	p14 = *(p + 13);
	p18 = *(p + 17);
	p32 = *(p + 31);
	p39 = *(p + 38);

	//f10.7 effect
	df = f107 - f107a;
	dfa = f107a - 150.0;
	t[0] = *(p + 19)*df*(1.0 + *(p + 59)*dfa) + *(p + 20)*df*df + *(p + 21)*dfa
		+ *(p + 29)*pow(dfa, 2.0);
	f1 = 1.0 + (*(p + 47)*dfa + *(p + 19)*df + *(p + 20)*df*df)*jb2008.hwm.getCswSwc(0);
	f2 = 1.0 + (*(p + 49)*dfa + *(p + 19)*df + *(p + 20)*df*df)*jb2008.hwm.getCswSwc(0);

	//Time independent
	t[1] = (*(p + 1)*plg[2][0] + *(p + 2)*plg[4][0] + *(p + 22)*plg[6][0])
		+ (*(p + 14)*plg[2][0])*dfa*jb2008.hwm.getCswSwc(0) + *(p + 26)*plg[1][0];

	//Symmetrical annual
	t[2] = (*(p + 18))*cd32;

	//Symmetrical semiannual
	t[3] = (*(p + 15) + *(p + 16)*plg[2][0])*cd18;


	//Asymmetrical annual
	t[4] = f1*(*(p + 9)*plg[1][0] + *(p + 10)*plg[3][0])*cd14;

	//Asymmetrical semiannual
	t[5] = *(p + 37)*plg[1][0] * cd39;

	//Diurnal
	if (jb2008.hwm.getCswSw(6) != 0.0) {
		t71 = (*(p + 11)*plg[2][1])*cd14*jb2008.hwm.getCswSwc(4);
		t72 = (*(p + 12)*plg[2][1])*cd14*jb2008.hwm.getCswSwc(4);
		t[6] = f2*((*(p + 3)*plg[1][1] + *(p + 4)*plg[3][1] + *(p + 27)*plg[5][1] + t71)*(ctloc)
			+(*(p + 6)*plg[1][1] + *(p + 7)*plg[3][1] + *(p + 28)*plg[5][1] + t72)*(stloc));
	}

	//Semidiurnal
	if (jb2008.hwm.getCswSw(7) != 0.0){
		t81 = (*(p + 23)*plg[3][2] + *(p + 35)*plg[5][2])*cd14*jb2008.hwm.getCswSwc(4);
		t82 = (*(p + 33)*plg[3][2] + *(p + 36)*plg[5][2])*cd14*jb2008.hwm.getCswSwc(4);

		t[7] = f2*((*(p + 5)*plg[2][2] + *(p + 41)*plg[4][2] + t81)*(c2tloc)
			+(*(p + 8)*plg[2][2] + *(p + 42)*plg[4][2] + t82)*(s2tloc));
	}

	//Terdiurnal
	if (jb2008.hwm.getCswSw(13) != 0.0){
		t[13] = f2*((*(p + 39)*plg[3][3] + (*(p + 93)*plg[4][3]
			+ *(p + 46)*plg[6][3])*cd14*jb2008.hwm.getCswSwc(4))*s3tloc
			+ (*(p + 40)*plg[3][3] + (*(p + 94)*plg[4][3]
			+ *(p + 48)*plg[6][3])*cd14*jb2008.hwm.getCswSwc(4))*c3tloc);
	}

	//Magnetic activity based on daily ap
	if (sw9 != -1.0){
		apd = ap[0] - 4.0;
		p44 = *(p + 43);
		p45 = *(p + 44);

		if (p44 < 0.0){
			p44 = 1.0e-05;
		}

		apdf = apd + (p45 - 1.0)*(apd + (exp(-p44*apd) - 1.0) / p44);

		if (0.0 != jb2008.hwm.getCswSw(8)) {
			t[8] = apdf*(*(p + 32) + *(p + 45)*plg[2][0]
				+ *(p + 34)*plg[4][0] + (*(p + 100)*plg[1][0] + *(p + 101)*plg[3][0]
				+ *(p + 102)*plg[5][0])*cd14*jb2008.hwm.getCswSwc(4)
				+ (*(p + 121)*plg[1][1] + *(p + 122)*plg[3][1]
				+ *(p + 123)*plg[5][1])*jb2008.hwm.getCswSwc(6)*cos(hr*(tloc - *(p + 124))));
		}
	}

	else {
		if (*(p + 51) != 0.0) {

			exp1 = exp(-10800.0*fabs(*(p + 51)) / (1.0 + *(p + 138)*(45.0 - fabs(lat))));

			if (0.99999 < exp1) {
				exp1 = 0.99999;
			}

			if (1.0e-04 > *(p + 24)) {
				*(p + 24) = 1.0e-04;
			}

			apt[0] = sg0(exp1, ap, p);


			if (jb2008.hwm.getCswSw(8) != 0.0) {
				t[8] = apt[0] * (*(p + 50) + *(p + 96)*plg[2][0]
					+ *(p + 54)*plg[4][0] + (*(p + 125)*plg[1][0] + *(p + 126)*plg[3][0]
					+ *(p + 127)*plg[5][0])*cd14*jb2008.hwm.getCswSwc(4)
					+ (*(p + 128)*plg[1][1] + *(p + 129)*plg[3][1]
					+ *(p + 130)*plg[5][1])*jb2008.hwm.getCswSwc(6)*cos(hr*(tloc - *(p + 131))));
			}
		}
	}

	if ((jb2008.hwm.getCswSw(9) != 0.0) || (lon > -1000.0)) {

		//Longitudinal
		if (0.0 != jb2008.hwm.getCswSw(10)) {
			t[10] = (1.0 + *(p + 80)*dfa*jb2008.hwm.getCswSwc(0))*((*(p + 64)*plg[2][1]
				+ *(p + 65)*plg[4][1] + *(p + 66)*plg[6][1] + *(p + 103)*plg[1][1]
				+ *(p + 104)*plg[3][1] + *(p + 105)*plg[5][1] + jb2008.hwm.getCswSwc(4)*(*(p + 109)*plg[1][1]
				+ *(p + 110)*plg[3][1] + *(p + 111)*plg[5][1])*cd14)*cos(dgtr*lon)
				+ (*(p + 90)*plg[2][1] + *(p + 91)*plg[4][1] + *(p + 92)*plg[6][1]
				+ *(p + 106)*plg[1][1] + *(p + 107)*plg[3][1] + *(p + 108)*plg[5][1]
				+ jb2008.hwm.getCswSwc(4)*(*(p + 112)*plg[1][1] + *(p + 113)*plg[3][1]
				+ *(p + 114)*plg[5][1])*cd14)*sin(dgtr*lon));
		}

		//Ut and mixed ut, longitude
		if (jb2008.hwm.getCswSw(11) != 0.0) {
			t[11] = (1.0 + *(p + 95)*plg[1][0])*(1.0 + *(p + 81)*dfa*jb2008.hwm.getCswSwc(0))*(1.0
				+ *(p + 119)*plg[1][0] * jb2008.hwm.getCswSwc(4)*cd14)*((*(p + 68)*plg[1][0]
				+ *(p + 69)*plg[3][0] + *(p + 70)*plg[5][0])*cos(sr*(sec - *(p + 71))));



			t[11] = t[11] + jb2008.hwm.getCswSwc(10)*(*(p + 76)*plg[3][2] + *(p + 77)*plg[5][2]
				+ *(p + 78)*plg[7][2])*cos(sr*(sec - *(p + 79)) + 2.0*dgtr*lon)*(1.0
				+ *(p + 137)*dfa*jb2008.hwm.getCswSwc(0));
		}

		//Ut, longitude magnetic activity
		if (jb2008.hwm.getCswSw(12) != 0.0) {

			if (sw9 != -1.0) {
				t[12] = apdf*jb2008.hwm.getCswSwc(10)*(1.0 + *(p + 120)*plg[1][0])*((*(p + 60)*plg[2][1]
					+ *(p + 61)*plg[4][1] + *(p + 62)*plg[6][1])*cos(dgtr*(lon - *(p + 63))))
					+ apdf*jb2008.hwm.getCswSwc(10)*jb2008.hwm.getCswSwc(4)*(*(p + 115)*plg[1][1]
					+ *(p + 116)*plg[3][1] + *(p + 117)*plg[5][1])*cd14*cos(dgtr*(lon - *(p + 118)))
					+ apdf*jb2008.hwm.getCswSwc(11)*(*(p + 83)*plg[1][0] + *(p + 84)*plg[3][0]
					+ *(p + 85)*plg[5][0])*cos(sr*(sec - *(p + 75)));
			}


			else if (*(p + 51) != 0.0) {
				t[12] = apt[0] * jb2008.hwm.getCswSwc(10)*(1.0 + *(p + 132)*plg[1][0])*((*(p + 52)*plg[2][1]
					+ *(p + 98)*plg[4][1] + *(p + 67)*plg[6][1])*cos(dgtr*(lon - *(p + 97))))
					+ apt[0] * jb2008.hwm.getCswSwc(10)*jb2008.hwm.getCswSwc(4)*(*(p + 133)*plg[1][1]
					+ *(p + 134)*plg[3][1] + *(p + 135)*plg[5][1])*cd14*cos(dgtr*(lon - *(p + 136)))
					+ apt[0] * jb2008.hwm.getCswSwc(11)*(*(p + 55)*plg[1][0] + *(p + 56)*plg[3][0]
					+ *(p + 57)*plg[5][0])*cos(sr*(sec - *(p + 58)));
			}
		}
	}

	//Parms not used: 83, 90, 100, 140-150
	tinf = *(p + 30);

	for (i = 0; i < nsw; ++i) {
		tinf = tinf + fabs(jb2008.hwm.getCswSw(i))*t[i];
	}

	globe7_out = tinf;

	return globe7_out;

}

double MSIS::glob7s(double *p)
{
	//glob7s member function from MSIS class
	//Version of globe for lower atmosphere

	double t[14], glob7s_out, cd14, cd18, cd32, cd39, t71, t72, t81, t82, tt,
		dayl = -1.0, p14 = -1000.0, p18 = -1000.0, p32 = -1000.0, p39 = -1000.0;
	int i, j;
	const double dr = 1.72142e-02, dgtr = 1.74533e-02, pset = 2.0;



	if (*(p + 99) == 0.0){
		*(p + 99) = pset;
	}

	if (pset != *(p + 99)){
		cout << "Wrong parameter set for glob7s" << '\n';
		system("pause");
		exit(1);
	}

	for (j = 0; j < 14; j++){
		t[j] = 0.0;
	}

	if ((day != dayl) || (*(p + 13) != p14)) {
		cd14 = cos(dr*(day - *(p + 13)));
	}
	if ((day != dayl) || (*(p + 17) != p18)) {
		cd18 = cos(2.0*dr*(day - *(p + 17)));
	}
	if ((day != dayl) || (*(p + 31) != p32)) {
		cd32 = cos(dr*(day - *(p + 31)));
	}
	if ((day != dayl) || (*(p + 38) != p39)) {
		cd39 = cos(2.0*dr*(day - *(p + 38)));
	}

	dayl = day;
	p14 = *(p + 13);
	p18 = *(p + 17);
	p32 = *(p + 31);
	p39 = *(p + 38);

	//F10.7
	t[0] = *(p + 21)*(dfa);

	//Time independent
	t[1] = *(p + 1)*plg[2][0] + *(p + 2)*plg[4][0] + *(p + 22)*plg[6][0] + *(p + 26)*plg[1][0]
		+ *(p + 14)*plg[3][0] + *(p + 59)*plg[5][0];

	//Symmetrical annual
	t[2] = (*(p + 18) + *(p + 47)*plg[2][0] + *(p + 29)*plg[4][0])*cd32;

	//Symmetrical semiannual
	t[3] = (*(p + 15) + *(p + 16)*plg[2][0] + *(p + 30)*plg[4][0])*cd18;

	//Asymmetrical annual
	t[4] = (*(p + 9)*plg[1][0] + *(p + 10)*plg[3][0] + *(p + 20)*plg[5][0])*cd14;

	//Asymmetrical semiannual
	t[5] = (*(p + 37)*plg[1][0])*cd39;

	//Diurnal
	if (jb2008.hwm.getCswSw(6) != 0.0) {
		t71 = *(p + 11)*plg[2][1] * cd14*jb2008.hwm.getCswSwc(4);
		t72 = *(p + 12)*plg[2][1] * cd14*jb2008.hwm.getCswSwc(4);
		t[6] = ((*(p + 3)*plg[1][1] + *(p + 4)*plg[3][1] + t71)*ctloc + (*(p + 6)*plg[1][1]
			+ *(p + 7)*plg[3][1] + t72)*stloc);
	}

	//Semidiurnal
	if (jb2008.hwm.getCswSw(7) != 0.0) {
		t81 = (*(p + 23)*plg[3][2] + *(p + 35)*plg[5][2])*cd14*jb2008.hwm.getCswSwc(4);
		t82 = (*(p + 33)*plg[3][2] + *(p + 36)*plg[5][2])*cd14*jb2008.hwm.getCswSwc(4);
		t[7] = ((*(p + 5)*plg[2][2] + *(p + 41)*plg[4][2] + t81)*c2tloc + (*(p + 8)*plg[2][2]
			+ *(p + 42)*plg[4][2] + t82)*s2tloc);
	}

	//Terdiurnal
	if (jb2008.hwm.getCswSw(13) != 0.0) {
		t[13] = *(p + 39)*plg[3][3] * s3tloc + *(p + 40)*plg[3][3] * c3tloc;
	}

	//Magnetic activity
	if (jb2008.hwm.getCswSw(8) != 0.0) {
		if (1.0 == jb2008.hwm.getCswSw(8)) {
			t[8] = apdf*(*(p + 32) + *(p + 45)*plg[2][0] * jb2008.hwm.getCswSwc(1));
		}
		else if (jb2008.hwm.getCswSw(8) == -1.0) {
			t[8] = (*(p + 50)*apt[0] + *(p + 96)*plg[2][0] * apt[0] * jb2008.hwm.getCswSwc(1));
		}
	}

	if ((jb2008.hwm.getCswSw(9) != 0.0) || (jb2008.hwm.getCswSw(10) != 0.0) || (xlong > -1000.0)){
		// Longitudinal
		t[10] = (1.0 + plg[1][0] * (*(p + 80)*jb2008.hwm.getCswSwc(4)*cos(dr*(day - *(p + 81)))
			+ *(p + 85)*jb2008.hwm.getCswSwc(5)*cos(2.0*dr*(day - *(p + 86))))
			+ *(p + 83)*jb2008.hwm.getCswSwc(2)*cos(dr*(day - *(p + 84)))
			+ *(p + 87)*jb2008.hwm.getCswSwc(3)*cos(2.0*dr*(day - *(p + 88))))*((*(p + 64)*plg[2][1]
			+ *(p + 65)*plg[4][1] + *(p + 66)*plg[6][1] + *(p + 74)*plg[1][1] + *(p + 75)*plg[3][1]
			+ *(p + 76)*plg[5][1])*cos(dgtr*xlong) + (*(p + 90)*plg[2][1] + *(p + 91)*plg[4][1]
			+ *(p + 92)*plg[6][1] + *(p + 77)*plg[1][1] + *(p + 78)*plg[3][1]
			+ *(p + 79)*plg[5][1])*sin(dgtr*xlong));
	}

	tt = 0.0;

	for (i = 0; i < 14; ++i) {
		tt = tt + fabs(jb2008.hwm.getCswSw(i))*t[i];
	}

	glob7s_out = tt;

	return glob7s_out;

}

double MSIS::vtst7(int iyd, double sec, double glat, double glong, double stl, double f107a, double f107, double *ap, int ic)
{
	//vtst7 member function from MSIS class
	//Test if geophysical variables or switches changed and save
	//Return 0 if unchanged and 1 if changed

	int iydl[2] = { -999, -999 };
	double secl[2] = { -999.0, -999.0 }, glatl[2] = { -999.0, -999.0 }, gll[2] = { -999.0, -999.0 },
		stll[2] = { -999.0, -999.0 }, fal[2] = { -999.0, -999.0 }, fl[2] = { -999.0, -999.0 },
		apl[7][2], swl[25][2], swcl[25][2];


	for (int i = 0; i < 25; i++){
		for (int j = 0; j < 2; j++){
			swl[i][j] = -999.0;
			swcl[i][j] = -999.0;

			if (i < 7){
				apl[i][j] = -999.0;
			}
		}
	}

	double vtst7_out;
	int i, j;

	ic = ic - 1;

	vtst7_out = 0;

	for (i = 0; i < 7; ++i) {
		for (j = 0; j < 25; ++j) {
			if ((iyd != iydl[ic]) ||
				(sec != secl[ic]) ||
				(glat != glatl[ic]) ||
				(glong != gll[ic]) ||
				(stl != stll[ic]) ||
				(f107a != fal[ic]) ||
				(f107 != fl[ic]) ||
				(ap[i] != apl[i][ic]) ||
				(jb2008.hwm.getCswSw(j) != swl[j][ic]) ||
				(swcl[j][ic] != jb2008.hwm.getCswSwc(j))) {

				vtst7_out = 1;
				iydl[ic] = iyd;
				secl[ic] = sec;
				glatl[ic] = glat;
				gll[ic] = glong;
				stll[ic] = stl;
				fal[ic] = f107a;
				fl[ic] = f107;

				for (i = 0; i < 7; ++i) {
					apl[i][ic] = ap[i];
				}

				for (i = 0; i < 25; ++i) {
					swl[i][ic] = jb2008.hwm.getCswSw(i);
					swcl[i][ic] = jb2008.hwm.getCswSwc(i);
				}

				return vtst7_out;
			}
		}
	}

	return vtst7_out;

}

void MSIS::gtd7bk()
{
	// MSISE-00 01-Feb-02
	//Tabulated MSIS data

	int i;
	int j;

	imr = 0;

	isdate[0] = "01-F";
	isdate[1] = "EB-0";
	isdate[2] = "2   ";
	istime[0] = "15:4";
	istime[1] = "9:27";
	name[0] = "MSIS";
	name[1] = "E-00";

	double pavgm_init[10] =
	{
		// Middle Atmosphere Averages
		2.61000e+02, 2.64000e+02, 2.29000e+02, 2.17000e+02, 2.17000e+02,
		2.23000e+02, 2.86760e+02, -2.93940e+00, 2.50000e+00, 0.00000e+00
	};

	double ptm_init[10] =
	{
		// Lower Boundary
		1.04130e+03, 3.86000e+02, 1.95000e+02, 1.66728e+01, 2.13000e+02,
		1.20000e+02, 2.40000e+02, 1.87000e+02, -2.00000e+00, 0.00000e+00
	};

	double pdm_init[8][10] =
	{
		{
			2.45600e+07, 6.71072e-06, 1.00000e+02, 0.00000e+00, 1.10000e+02,
			1.00000e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
		},
		{
			8.59400e+10, 1.00000e+00, 1.05000e+02, -8.00000e+00, 1.10000e+02,
			1.00000e+01, 9.00000e+01, 2.00000e+00, 0.00000e+00, 0.00000e+00,
		},
		{
			2.81000e+11, 0.00000e+00, 1.05000e+02, 2.80000e+01, 2.89500e+01,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
		},
		{
			3.30000e+10, 2.68270e-01, 1.05000e+02, 1.00000e+00, 1.10000e+02,
			1.00000e+01, 1.10000e+02, -1.00000e+01, 0.00000e+00, 0.00000e+00
		},
		{
			1.33000e+09, 1.19615e-02, 1.05000e+02, 0.00000e+00, 1.10000e+02,
			1.00000e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
		},
		{
			1.76100e+05, 1.00000e+00, 9.50000e+01, -8.00000e+00, 1.10000e+02,
			1.00000e+01, 9.00000e+01, 2.00000e+00, 0.00000e+00, 0.00000e+00
		},
		{
			1.00000e+07, 1.00000e+00, 1.05000e+02, -8.00000e+00, 1.10000e+02,
			1.00000e+01, 9.00000e+01, 2.00000e+00, 0.00000e+00, 0.00000e+00
		},
		{
			1.00000e+06, 1.00000e+00, 1.05000e+02, -8.00000e+00, 5.50000e+02,
			7.60000e+01, 9.00000e+01, 2.00000e+00, 0.00000e+00, 4.00000e+03
		}
	};

	double pt_init[150] =
	{
		//        TEMPERATURE

		9.86573e-01, 1.62228e-02, 1.55270e-02, -1.04323e-01, -3.75801e-03,
		-1.18538e-03, -1.24043e-01, 4.56820e-03, 8.76018e-03, -1.36235e-01,
		-3.52427e-02, 8.84181e-03, -5.92127e-03, -8.61650e+00, 0.00000e+00,
		1.28492e-02, 0.00000e+00, 1.30096e+02, 1.04567e-02, 1.65686e-03,
		-5.53887e-06, 2.97810e-03, 0.00000e+00, 5.13122e-03, 8.66784e-02,
		1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, -7.27026e-06,
		0.00000e+00, 6.74494e+00, 4.93933e-03, 2.21656e-03, 2.50802e-03,
		0.00000e+00, 0.00000e+00, -2.08841e-02, -1.79873e+00, 1.45103e-03,
		2.81769e-04, -1.44703e-03, -5.16394e-05, 8.47001e-02, 1.70147e-01,
		5.72562e-03, 5.07493e-05, 4.36148e-03, 1.17863e-04, 4.74364e-03,

		6.61278e-03, 4.34292e-05, 1.44373e-03, 2.41470e-05, 2.84426e-03,
		8.56560e-04, 2.04028e-03, 0.00000e+00, -3.15994e+03, -2.46423e-03,
		1.13843e-03, 4.20512e-04, 0.00000e+00, -9.77214e+01, 6.77794e-03,
		5.27499e-03, 1.14936e-03, 0.00000e+00, -6.61311e-03, -1.84255e-02,
		-1.96259e-02, 2.98618e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		6.44574e+02, 8.84668e-04, 5.05066e-04, 0.00000e+00, 4.02881e+03,
		-1.89503e-03, 0.00000e+00, 0.00000e+00, 8.21407e-04, 2.06780e-03,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		-1.20410e-02, -3.63963e-03, 9.92070e-05, -1.15284e-04, -6.33059e-05,
		-6.05545e-01, 8.34218e-03, -9.13036e+01, 3.71042e-04, 0.00000e+00,

		4.19000e-04, 2.70928e-03, 3.31507e-03, -4.44508e-03, -4.96334e-03,
		-1.60449e-03, 3.95119e-03, 2.48924e-03, 5.09815e-04, 4.05302e-03,
		2.24076e-03, 0.00000e+00, 6.84256e-03, 4.66354e-04, 0.00000e+00,
		-3.68328e-04, 0.00000e+00, 0.00000e+00, -1.46870e+02, 0.00000e+00,
		0.00000e+00, 1.09501e-03, 4.65156e-04, 5.62583e-04, 3.21596e+00,
		6.43168e-04, 3.14860e-03, 3.40738e-03, 1.78481e-03, 9.62532e-04,
		5.58171e-04, 3.43731e+00, -2.33195e-01, 5.10289e-04, 0.00000e+00,
		0.00000e+00, -9.25347e+04, 0.00000e+00, -1.99639e-03, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
	};

	double pd_init[9][150] =
	{
		//         HE DENSITY
		{
			1.09979e+00, -4.88060e-02, -1.97501e-01, -9.10280e-02, -6.96558e-03,
			2.42136e-02, 3.91333e-01, -7.20068e-03, -3.22718e-02, 1.41508e+00,
			1.68194e-01, 1.85282e-02, 1.09384e-01, -7.24282e+00, 0.00000e+00,
			2.96377e-01, -4.97210e-02, 1.04114e+02, -8.61108e-02, -7.29177e-04,
			1.48998e-06, 1.08629e-03, 0.00000e+00, 0.00000e+00, 8.31090e-02,
			1.12818e-01, -5.75005e-02, -1.29919e-02, -1.78849e-02, -2.86343e-06,
			0.00000e+00, -1.51187e+02, -6.65902e-03, 0.00000e+00, -2.02069e-03,
			0.00000e+00, 0.00000e+00, 4.32264e-02, -2.80444e+01, -3.26789e-03,
			2.47461e-03, 0.00000e+00, 0.00000e+00, 9.82100e-02, 1.22714e-01,
			-3.96450e-02, 0.00000e+00, -2.76489e-03, 0.00000e+00, 1.87723e-03,

			-8.09813e-03, 4.34428e-05, -7.70932e-03, 0.00000e+00, -2.28894e-03,
			-5.69070e-03, -5.22193e-03, 6.00692e-03, -7.80434e+03, -3.48336e-03,
			-6.38362e-03, -1.82190e-03, 0.00000e+00, -7.58976e+01, -2.17875e-02,
			-1.72524e-02, -9.06287e-03, 0.00000e+00, 2.44725e-02, 8.66040e-02,
			1.05712e-01, 3.02543e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-6.01364e+03, -5.64668e-03, -2.54157e-03, 0.00000e+00, 3.15611e+02,
			-5.69158e-03, 0.00000e+00, 0.00000e+00, -4.47216e-03, -4.49523e-03,
			4.64428e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			4.51236e-02, 2.46520e-02, 6.17794e-03, 0.00000e+00, 0.00000e+00,
			-3.62944e-01, -4.80022e-02, -7.57230e+01, -1.99656e-03, 0.00000e+00,

			-5.18780e-03, -1.73990e-02, -9.03485e-03, 7.48465e-03, 1.53267e-02,
			1.06296e-02, 1.18655e-02, 2.55569e-03, 1.69020e-03, 3.51936e-02,
			-1.81242e-02, 0.00000e+00, -1.00529e-01, -5.10574e-03, 0.00000e+00,
			2.10228e-03, 0.00000e+00, 0.00000e+00, -1.73255e+02, 5.07833e-01,
			-2.41408e-01, 8.75414e-03, 2.77527e-03, -8.90353e-05, -5.25148e+00,
			-5.83899e-03, -2.09122e-02, -9.63530e-03, 9.77164e-03, 4.07051e-03,
			2.53555e-04, -5.52875e+00, -3.55993e-01, -2.49231e-03, 0.00000e+00,
			0.00000e+00, 2.86026e+01, 0.00000e+00, 3.42722e-04, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
		},
		{
			//         O DENSITY

			1.02315e+00, -1.59710e-01, -1.06630e-01, -1.77074e-02, -4.42726e-03,
			3.44803e-02, 4.45613e-02, -3.33751e-02, -5.73598e-02, 3.50360e-01,
			6.33053e-02, 2.16221e-02, 5.42577e-02, -5.74193e+00, 0.00000e+00,
			1.90891e-01, -1.39194e-02, 1.01102e+02, 8.16363e-02, 1.33717e-04,
			6.54403e-06, 3.10295e-03, 0.00000e+00, 0.00000e+00, 5.38205e-02,
			1.23910e-01, -1.39831e-02, 0.00000e+00, 0.00000e+00, -3.95915e-06,
			0.00000e+00, -7.14651e-01, -5.01027e-03, 0.00000e+00, -3.24756e-03,
			0.00000e+00, 0.00000e+00, 4.42173e-02, -1.31598e+01, -3.15626e-03,
			1.24574e-03, -1.47626e-03, -1.55461e-03, 6.40682e-02, 1.34898e-01,
			-2.42415e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 6.13666e-04,

			-5.40373e-03, 2.61635e-05, -3.33012e-03, 0.00000e+00, -3.08101e-03,
			-2.42679e-03, -3.36086e-03, 0.00000e+00, -1.18979e+03, -5.04738e-02,
			-2.61547e-03, -1.03132e-03, 1.91583e-04, -8.38132e+01, -1.40517e-02,
			-1.14167e-02, -4.08012e-03, 1.73522e-04, -1.39644e-02, -6.64128e-02,
			-6.85152e-02, -1.34414e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			6.07916e+02, -4.12220e-03, -2.20996e-03, 0.00000e+00, 1.70277e+03,
			-4.63015e-03, 0.00000e+00, 0.00000e+00, -2.25360e-03, -2.96204e-03,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			3.92786e-02, 1.31186e-02, -1.78086e-03, 0.00000e+00, 0.00000e+00,
			-3.90083e-01, -2.84741e-02, -7.78400e+01, -1.02601e-03, 0.00000e+00,

			-7.26485e-04, -5.42181e-03, -5.59305e-03, 1.22825e-02, 1.23868e-02,
			6.68835e-03, -1.03303e-02, -9.51903e-03, 2.70021e-04, -2.57084e-02,
			-1.32430e-02, 0.00000e+00, -3.81000e-02, -3.16810e-03, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -9.05762e-04, -2.14590e-03, -1.17824e-03, 3.66732e+00,
			-3.79729e-04, -6.13966e-03, -5.09082e-03, -1.96332e-03, -3.08280e-03,
			-9.75222e-04, 4.03315e+00, -2.52710e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
		},
		{
			//         N2 DENSITY

			1.16112e+00, 0.00000e+00, 0.00000e+00, 3.33725e-02, 0.00000e+00,
			3.48637e-02, -5.44368e-03, 0.00000e+00, -6.73940e-02, 1.74754e-01,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 1.74712e+02, 0.00000e+00,
			1.26733e-01, 0.00000e+00, 1.03154e+02, 5.52075e-02, 0.00000e+00,
			0.00000e+00, 8.13525e-04, 0.00000e+00, 0.00000e+00, 8.66784e-02,
			1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -2.50482e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.48894e-03,
			6.16053e-04, -5.79716e-04, 2.95482e-03, 8.47001e-02, 1.70147e-01,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,

			0.00000e+00, 2.47425e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
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
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
		},
		{
			//         TLB

			9.44846e-01, 0.00000e+00, 0.00000e+00, -3.08617e-02, 0.00000e+00,
			-2.44019e-02, 6.48607e-03, 0.00000e+00, 3.08181e-02, 4.59392e-02,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 1.74712e+02, 0.00000e+00,
			2.13260e-02, 0.00000e+00, -3.56958e+02, 0.00000e+00, 1.82278e-04,
			0.00000e+00, 3.07472e-04, 0.00000e+00, 0.00000e+00, 8.66784e-02,
			1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 3.83054e-03, 0.00000e+00, 0.00000e+00,
			-1.93065e-03, -1.45090e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.23493e-03, 1.36736e-03, 8.47001e-02, 1.70147e-01,
			3.71469e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,

			5.10250e-03, 2.47425e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 3.68756e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,

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
		},
		{
			//         O2 DENSITY

			1.35580e+00, 1.44816e-01, 0.00000e+00, 6.07767e-02, 0.00000e+00,
			2.94777e-02, 7.46900e-02, 0.00000e+00, -9.23822e-02, 8.57342e-02,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 2.38636e+01, 0.00000e+00,
			7.71653e-02, 0.00000e+00, 8.18751e+01, 1.87736e-02, 0.00000e+00,
			0.00000e+00, 1.49667e-02, 0.00000e+00, 0.00000e+00, 8.66784e-02,
			1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -3.67874e+02, 5.48158e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 8.47001e-02, 1.70147e-01,
			1.22631e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,

			8.17187e-03, 3.71617e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.10826e-03,
			-3.13640e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-7.35742e-02, -5.00266e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 1.94965e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,

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
		},
		{
			//         AR DENSITY

			1.04761e+00, 2.00165e-01, 2.37697e-01, 3.68552e-02, 0.00000e+00,
			3.57202e-02, -2.14075e-01, 0.00000e+00, -1.08018e-01, -3.73981e-01,
			0.00000e+00, 3.10022e-02, -1.16305e-03, -2.07596e+01, 0.00000e+00,
			8.64502e-02, 0.00000e+00, 9.74908e+01, 5.16707e-02, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 8.66784e-02,
			1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 3.46193e+02, 1.34297e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -3.48509e-03,
			-1.54689e-04, 0.00000e+00, 0.00000e+00, 8.47001e-02, 1.70147e-01,
			1.47753e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,

			1.89320e-02, 3.68181e-05, 1.32570e-02, 0.00000e+00, 0.00000e+00,
			3.59719e-03, 7.44328e-03, -1.00023e-03, -6.50528e+03, 0.00000e+00,
			1.03485e-02, -1.00983e-03, -4.06916e-03, -6.60864e+01, -1.71533e-02,
			1.10605e-02, 1.20300e-02, -5.20034e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-2.62769e+03, 7.13755e-03, 4.17999e-03, 0.00000e+00, 1.25910e+04,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -2.23595e-03, 4.60217e-03,
			5.71794e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-3.18353e-02, -2.35526e-02, -1.36189e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 2.03522e-02, -6.67837e+01, -1.09724e-03, 0.00000e+00,

			-1.38821e-02, 1.60468e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.51574e-02,
			-5.44470e-04, 0.00000e+00, 7.28224e-02, 6.59413e-02, 0.00000e+00,
			-5.15692e-03, 0.00000e+00, 0.00000e+00, -3.70367e+03, 0.00000e+00,
			0.00000e+00, 1.36131e-02, 5.38153e-03, 0.00000e+00, 4.76285e+00,
			-1.75677e-02, 2.26301e-02, 0.00000e+00, 1.76631e-02, 4.77162e-03,
			0.00000e+00, 5.39354e+00, 0.00000e+00, -7.51710e-03, 0.00000e+00,
			0.00000e+00, -8.82736e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
		},
		{
			//          H DENSITY

			1.26376e+00, -2.14304e-01, -1.49984e-01, 2.30404e-01, 2.98237e-02,
			2.68673e-02, 2.96228e-01, 2.21900e-02, -2.07655e-02, 4.52506e-01,
			1.20105e-01, 3.24420e-02, 4.24816e-02, -9.14313e+00, 0.00000e+00,
			2.47178e-02, -2.88229e-02, 8.12805e+01, 5.10380e-02, -5.80611e-03,
			2.51236e-05, -1.24083e-02, 0.00000e+00, 0.00000e+00, 8.66784e-02,
			1.58727e-01, -3.48190e-02, 0.00000e+00, 0.00000e+00, 2.89885e-05,
			0.00000e+00, 1.53595e+02, -1.68604e-02, 0.00000e+00, 1.01015e-02,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.84552e-04,
			-1.22181e-03, 0.00000e+00, 0.00000e+00, 8.47001e-02, 1.70147e-01,
			-1.04927e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, -5.91313e-03,

			-2.30501e-02, 3.14758e-05, 0.00000e+00, 0.00000e+00, 1.26956e-02,
			8.35489e-03, 3.10513e-04, 0.00000e+00, 3.42119e+03, -2.45017e-03,
			-4.27154e-04, 5.45152e-04, 1.89896e-03, 2.89121e+01, -6.49973e-03,
			-1.93855e-02, -1.48492e-02, 0.00000e+00, -5.10576e-02, 7.87306e-02,
			9.51981e-02, -1.49422e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			2.65503e+02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 6.37110e-03, 3.24789e-04,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			6.14274e-02, 1.00376e-02, -8.41083e-04, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.27099e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,

			-3.94077e-03, -1.28601e-02, -7.97616e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -6.71465e-03, -1.69799e-03, 1.93772e-03, 3.81140e+00,
			-7.79290e-03, -1.82589e-02, -1.25860e-02, -1.04311e-02, -3.02465e-03,
			2.43063e-03, 3.63237e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
		},
		{
			//          N DENSITY

			7.09557e+01, -3.26740e-01, 0.00000e+00, -5.16829e-01, -1.71664e-03,
			9.09310e-02, -6.71500e-01, -1.47771e-01, -9.27471e-02, -2.30862e-01,
			-1.56410e-01, 1.34455e-02, -1.19717e-01, 2.52151e+00, 0.00000e+00,
			-2.41582e-01, 5.92939e-02, 4.39756e+00, 9.15280e-02, 4.41292e-03,
			0.00000e+00, 8.66807e-03, 0.00000e+00, 0.00000e+00, 8.66784e-02,
			1.58727e-01, 9.74701e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 6.70217e+01, -1.31660e-03, 0.00000e+00, -1.65317e-02,
			0.00000e+00, 0.00000e+00, 8.50247e-02, 2.77428e+01, 4.98658e-03,
			6.15115e-03, 9.50156e-03, -2.12723e-02, 8.47001e-02, 1.70147e-01,
			-2.38645e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.37380e-03,

			-8.41918e-03, 2.80145e-05, 7.12383e-03, 0.00000e+00, -1.66209e-02,
			1.03533e-04, -1.68898e-02, 0.00000e+00, 3.64526e+03, 0.00000e+00,
			6.54077e-03, 3.69130e-04, 9.94419e-04, 8.42803e+01, -1.16124e-02,
			-7.74414e-03, -1.68844e-03, 1.42809e-03, -1.92955e-03, 1.17225e-01,
			-2.41512e-02, 1.50521e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			1.60261e+03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -3.54403e-04, -1.87270e-02,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			2.76439e-02, 6.43207e-03, -3.54300e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -2.80221e-02, 8.11228e+01, -6.75255e-04, 0.00000e+00,

			-1.05162e-02, -3.48292e-03, -6.97321e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.45546e-03, -1.31970e-02, -3.57751e-03, -1.09021e+00,
			-1.50181e-02, -7.12841e-03, -6.64590e-03, -3.52610e-03, -1.87773e-02,
			-2.22432e-03, -3.93895e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
		},
		{
			//        HOT O DENSITY

			6.04050e-02, 1.57034e+00, 2.99387e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -1.51018e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -8.61650e+00, 1.26454e-02,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 5.50878e-03, 0.00000e+00, 0.00000e+00, 8.66784e-02,
			1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 6.23881e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 8.47001e-02, 1.70147e-01,
			-9.45934e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,

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
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
		}
	};

	double ps_init[150] =
	{
		//          S PARAM

		9.56827e-01, 6.20637e-02, 3.18433e-02, 0.00000e+00, 0.00000e+00,
		3.94900e-02, 0.00000e+00, 0.00000e+00, -9.24882e-03, -7.94023e-03,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 1.74712e+02, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 2.74677e-03, 0.00000e+00, 1.54951e-02, 8.66784e-02,
		1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, -6.99007e-04, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 1.24362e-02, -5.28756e-03, 8.47001e-02, 1.70147e-01,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,

		0.00000e+00, 2.47425e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
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
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
	};

	double pdl_init[2][25] =
	{
		//          TURBO
		{
			1.09930e+00, 3.90631e+00, 3.07165e+00, 9.86161e-01, 1.63536e+01,
			4.63830e+00, 1.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 1.28840e+00, 3.10302e-02, 1.18339e-01
		},
		{
			1.00000e+00, 7.00000e-01, 1.15020e+00, 3.44689e+00, 1.28840e+00,
			1.00000e+00, 1.08738e+00, 1.22947e+00, 1.10016e+00, 7.34129e-01,
			1.15241e+00, 2.22784e+00, 7.95046e-01, 4.01612e+00, 4.47749e+00,
			1.23435e+02, -7.60535e-02, 1.68986e-06, 7.44294e-01, 1.03604e+00,
			1.72783e+02, 1.15020e+00, 3.44689e+00, -7.46230e-01, 9.49154e-01
		}
	};

	double ptl_init[4][100] = {
		//         TN1(2)
		{
			1.00858e+00, 4.56011e-02, -2.22972e-02, -5.44388e-02, 5.23136e-04,
			-1.88849e-02, 5.23707e-02, -9.43646e-03, 6.31707e-03, -7.80460e-02,
			-4.88430e-02, 0.00000e+00, 0.00000e+00, -7.60250e+00, 0.00000e+00,
			-1.44635e-02, -1.76843e-02, -1.21517e+02, 2.85647e-02, 0.00000e+00,
			0.00000e+00, 6.31792e-04, 0.00000e+00, 5.77197e-03, 8.66784e-02,
			1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -8.90272e+03, 3.30611e-03, 3.02172e-03, 0.00000e+00,
			-2.13673e-03, -3.20910e-04, 0.00000e+00, 0.00000e+00, 2.76034e-03,
			2.82487e-03, -2.97592e-04, -4.21534e-03, 8.47001e-02, 1.70147e-01,
			8.96456e-03, 0.00000e+00, -1.08596e-02, 0.00000e+00, 0.00000e+00,

			5.57917e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 9.65405e-03, 0.00000e+00, 0.00000e+00, 2.00000e+00
		},
		{
			//         TN1(3)

			9.39664e-01, 8.56514e-02, -6.79989e-03, 2.65929e-02, -4.74283e-03,
			1.21855e-02, -2.14905e-02, 6.49651e-03, -2.05477e-02, -4.24952e-02,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 1.19148e+01, 0.00000e+00,
			1.18777e-02, -7.28230e-02, -8.15965e+01, 1.73887e-02, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -1.44691e-02, 2.80259e-04, 8.66784e-02,
			1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 2.16584e+02, 3.18713e-03, 7.37479e-03, 0.00000e+00,
			-2.55018e-03, -3.92806e-03, 0.00000e+00, 0.00000e+00, -2.89757e-03,
			-1.33549e-03, 1.02661e-03, 3.53775e-04, 8.47001e-02, 1.70147e-01,
			-9.17497e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,

			3.56082e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.00902e-02, 0.00000e+00, 0.00000e+00, 2.00000e+00
		},
		{
			//         TN1(4)

			9.85982e-01, -4.55435e-02, 1.21106e-02, 2.04127e-02, -2.40836e-03,
			1.11383e-02, -4.51926e-02, 1.35074e-02, -6.54139e-03, 1.15275e-01,
			1.28247e-01, 0.00000e+00, 0.00000e+00, -5.30705e+00, 0.00000e+00,
			-3.79332e-02, -6.24741e-02, 7.71062e-01, 2.96315e-02, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 6.81051e-03, -4.34767e-03, 8.66784e-02,
			1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 1.07003e+01, -2.76907e-03, 4.32474e-04, 0.00000e+00,
			1.31497e-03, -6.47517e-04, 0.00000e+00, -2.20621e+01, -1.10804e-03,
			-8.09338e-04, 4.18184e-04, 4.29650e-03, 8.47001e-02, 1.70147e-01,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,

			-4.04337e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -9.52550e-04,
			8.56253e-04, 4.33114e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.21223e-03,
			2.38694e-04, 9.15245e-04, 1.28385e-03, 8.67668e-04, -5.61425e-06,
			1.04445e+00, 3.41112e+01, 0.00000e+00, -8.40704e-01, -2.39639e+02,
			7.06668e-01, -2.05873e+01, -3.63696e-01, 2.39245e+01, 0.00000e+00,
			-1.06657e-03, -7.67292e-04, 1.54534e-04, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
		},
		{
			//         TN1(5) TN2(1)

			1.00320e+00, 3.83501e-02, -2.38983e-03, 2.83950e-03, 4.20956e-03,
			5.86619e-04, 2.19054e-02, -1.00946e-02, -3.50259e-03, 4.17392e-02,
			-8.44404e-03, 0.00000e+00, 0.00000e+00, 4.96949e+00, 0.00000e+00,
			-7.06478e-03, -1.46494e-02, 3.13258e+01, -1.86493e-03, 0.00000e+00,
			-1.67499e-02, 0.00000e+00, 0.00000e+00, 5.12686e-04, 8.66784e-02,
			1.58727e-01, -4.64167e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			4.37353e-03, -1.99069e+02, 0.00000e+00, -5.34884e-03, 0.00000e+00,
			1.62458e-03, 2.93016e-03, 2.67926e-03, 5.90449e+02, 0.00000e+00,
			0.00000e+00, -1.17266e-03, -3.58890e-04, 8.47001e-02, 1.70147e-01,
			0.00000e+00, 0.00000e+00, 1.38673e-02, 0.00000e+00, 0.00000e+00,

			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.60571e-03,
			6.28078e-04, 5.05469e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -1.57829e-03,
			-4.00855e-04, 5.04077e-05, -1.39001e-03, -2.33406e-03, -4.81197e-04,
			1.46758e+00, 6.20332e+00, 0.00000e+00, 3.66476e-01, -6.19760e+01,
			3.09198e-01, -1.98999e+01, 0.00000e+00, -3.29933e+02, 0.00000e+00,
			-1.10080e-03, -9.39310e-05, 1.39638e-04, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
		}
	};

	double pma_init[10][100] = {
		//          TN2(2)
		{
			9.81637e-01, -1.41317e-03, 3.87323e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -3.58707e-02,
			-8.63658e-03, 0.00000e+00, 0.00000e+00, -2.02226e+00, 0.00000e+00,
			-8.69424e-03, -1.91397e-02, 8.76779e+01, 4.52188e-03, 0.00000e+00,
			2.23760e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -7.07572e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-4.11210e-03, 3.50060e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -8.36657e-03, 1.61347e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -1.45130e-02, 0.00000e+00, 0.00000e+00,

			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.24152e-03,
			6.43365e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.33255e-03,
			2.42657e-03, 1.60666e-03, -1.85728e-03, -1.46874e-03, -4.79163e-06,
			1.22464e+00, 3.53510e+01, 0.00000e+00, 4.49223e-01, -4.77466e+01,
			4.70681e-01, 8.41861e+00, -2.88198e-01, 1.67854e+02, 0.00000e+00,
			7.11493e-04, 6.05601e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
		},
		{
			//          TN2(3)

			1.00422e+00, -7.11212e-03, 5.24480e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -5.28914e-02,
			-2.41301e-02, 0.00000e+00, 0.00000e+00, -2.12219e+01, -1.03830e-02,
			-3.28077e-03, 1.65727e-02, 1.68564e+00, -6.68154e-03, 0.00000e+00,
			1.45155e-02, 0.00000e+00, 8.42365e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -4.34645e-03, 0.00000e+00, 0.00000e+00, 2.16780e-02,
			0.00000e+00, -1.38459e+02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 7.04573e-03, -4.73204e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 1.08767e-02, 0.00000e+00, 0.00000e+00,

			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -8.08279e-03,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.21769e-04,
			-2.27387e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.26769e-03,
			3.16901e-03, 4.60316e-04, -1.01431e-04, 1.02131e-03, 9.96601e-04,
			1.25707e+00, 2.50114e+01, 0.00000e+00, 4.24472e-01, -2.77655e+01,
			3.44625e-01, 2.75412e+01, 0.00000e+00, 7.94251e+02, 0.00000e+00,
			2.45835e-03, 1.38871e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
		},
		{
			//          TN2(4) TN3(1)

			1.01890e+00, -2.46603e-02, 1.00078e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -6.70977e-02,
			-4.02286e-02, 0.00000e+00, 0.00000e+00, -2.29466e+01, -7.47019e-03,
			2.26580e-03, 2.63931e-02, 3.72625e+01, -6.39041e-03, 0.00000e+00,
			9.58383e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.85291e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 1.39717e+02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 9.19771e-03, -3.69121e+02, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -1.57067e-02, 0.00000e+00, 0.00000e+00,

			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -7.07265e-03,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.92953e-03,
			-2.77739e-03, -4.40092e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.47280e-03,
			2.95035e-04, -1.81246e-03, 2.81945e-03, 4.27296e-03, 9.78863e-04,
			1.40545e+00, -6.19173e+00, 0.00000e+00, 0.00000e+00, -7.93632e+01,
			4.44643e-01, -4.03085e+02, 0.00000e+00, 1.15603e+01, 0.00000e+00,
			2.25068e-03, 8.48557e-04, -2.98493e-04, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
		},
		{
			//          TN3(2)

			9.75801e-01, 3.80680e-02, -3.05198e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.85575e-02,
			5.04057e-02, 0.00000e+00, 0.00000e+00, -1.76046e+02, 1.44594e-02,
			-1.48297e-03, -3.68560e-03, 3.02185e+01, -3.23338e-03, 0.00000e+00,
			1.53569e-02, 0.00000e+00, -1.15558e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 4.89620e-03, 0.00000e+00, 0.00000e+00, -1.00616e-02,
			-8.21324e-03, -1.57757e+02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 6.63564e-03, 4.58410e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -2.51280e-02, 0.00000e+00, 0.00000e+00,

			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 9.91215e-03,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -8.73148e-04,
			-1.29648e-03, -7.32026e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -4.68110e-03,
			-4.66003e-03, -1.31567e-03, -7.39390e-04, 6.32499e-04, -4.65588e-04,
			-1.29785e+00, -1.57139e+02, 0.00000e+00, 2.58350e-01, -3.69453e+01,
			4.10672e-01, 9.78196e+00, -1.52064e-01, -3.85084e+03, 0.00000e+00,
			-8.52706e-04, -1.40945e-03, -7.26786e-04, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
		},
		{
			//          TN3(3)

			9.60722e-01, 7.03757e-02, -3.00266e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.22671e-02,
			4.10423e-02, 0.00000e+00, 0.00000e+00, -1.63070e+02, 1.06073e-02,
			5.40747e-04, 7.79481e-03, 1.44908e+02, 1.51484e-04, 0.00000e+00,
			1.97547e-02, 0.00000e+00, -1.41844e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 5.77884e-03, 0.00000e+00, 0.00000e+00, 9.74319e-03,
			0.00000e+00, -2.88015e+03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -4.44902e-03, -2.92760e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 2.34419e-02, 0.00000e+00, 0.00000e+00,

			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.36685e-03,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -4.65325e-04,
			-5.50628e-04, 3.31465e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.06179e-03,
			-3.08575e-03, -7.93589e-04, -1.08629e-04, 5.95511e-04, -9.05050e-04,
			1.18997e+00, 4.15924e+01, 0.00000e+00, -4.72064e-01, -9.47150e+02,
			3.98723e-01, 1.98304e+01, 0.00000e+00, 3.73219e+03, 0.00000e+00,
			-1.50040e-03, -1.14933e-03, -1.56769e-04, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
		},
		{
			//          TN3(4)

			1.03123e+00, -7.05124e-02, 8.71615e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -3.82621e-02,
			-9.80975e-03, 0.00000e+00, 0.00000e+00, 2.89286e+01, 9.57341e-03,
			0.00000e+00, 0.00000e+00, 8.66153e+01, 7.91938e-04, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 4.68917e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 7.86638e-03, 0.00000e+00, 0.00000e+00, 9.90827e-03,
			0.00000e+00, 6.55573e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -4.00200e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 7.07457e-03, 0.00000e+00, 0.00000e+00,

			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.72268e-03,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.04970e-04,
			1.21560e-03, -8.05579e-06, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.49941e-03,
			-4.57256e-04, -1.59311e-04, 2.96481e-04, -1.77318e-03, -6.37918e-04,
			1.02395e+00, 1.28172e+01, 0.00000e+00, 1.49903e-01, -2.63818e+01,
			0.00000e+00, 4.70628e+01, -2.22139e-01, 4.82292e-02, 0.00000e+00,
			-8.67075e-04, -5.86479e-04, 5.32462e-04, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
		},
		{
			//          TN3(5) SURFACE TEMP TSL

			1.00828e+00, -9.10404e-02, -2.26549e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.32420e-02,
			-9.08925e-03, 0.00000e+00, 0.00000e+00, 3.36105e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -1.24957e+01, -5.87939e-03, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 2.79765e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 2.01237e+03, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -1.75553e-02, 0.00000e+00, 0.00000e+00,

			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.29699e-03,
			1.26659e-03, 2.68402e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.17894e-03,
			1.48746e-03, 1.06478e-04, 1.34743e-04, -2.20939e-03, -6.23523e-04,
			6.36539e-01, 1.13621e+01, 0.00000e+00, -3.93777e-01, 2.38687e+03,
			0.00000e+00, 6.61865e+02, -1.21434e-01, 9.27608e+00, 0.00000e+00,
			1.68478e-04, 1.24892e-03, 1.71345e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
		},
		{
			//          TGN3(2) SURFACE GRAD TSLG

			1.57293e+00, -6.78400e-01, 6.47500e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -7.62974e-02,
			-3.60423e-01, 0.00000e+00, 0.00000e+00, 1.28358e+02, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 4.68038e+01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.67898e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 2.90994e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 3.15706e+01, 0.00000e+00,
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
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
		},
		{
			//          TGN2(1) TGN1(2)

			8.60028e-01, 3.77052e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -1.17570e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 7.77757e-03, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 1.01024e+02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 6.54251e+02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,

			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -1.56959e-02,
			1.91001e-02, 3.15971e-02, 1.00982e-02, -6.71565e-03, 2.57693e-03,
			1.38692e+00, 2.82132e-01, 0.00000e+00, 0.00000e+00, 3.81511e+02,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
		},
		{
			//          TGN3(1) TGN2(2)

			1.06029e+00, -5.25231e-02, 3.73034e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.31072e-02,
			-3.88409e-01, 0.00000e+00, 0.00000e+00, -1.65295e+02, -2.13801e-01,
			-4.38916e-02, -3.22716e-01, -8.82393e+01, 1.18458e-01, 0.00000e+00,
			-4.35863e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.19782e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 2.62229e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -5.37443e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -4.55788e-01, 0.00000e+00, 0.00000e+00,

			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.84009e-02,
			3.96733e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.05494e-02,
			7.39617e-02, 1.92200e-02, -8.46151e-03, -1.34244e-02, 1.96338e-02,
			1.50421e+00, 1.88368e+01, 0.00000e+00, 0.00000e+00, -5.13114e+01,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			5.11923e-02, 3.61225e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
		}
	};

	double sam_init[100] = {
		//          SEMIANNUAL MULT SAM

		1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
		1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
		1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
		1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
		1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
		1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
		1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
		1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
		1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
		1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,

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

	for (i = 0; i < 10; ++i) {
		pavgm[i] = pavgm_init[i];
		ptm[i] = ptm_init[i];

		for (j = 0; j < 8; ++j) {
			pdm[i][j] = pdm_init[j][i];
		}
	}

	for (i = 0; i < 25; ++i) {

		for (j = 0; j < 2; ++j) {
			pdl[i][j] = float(pdl_init[j][i]);
		}
	}

	for (i = 0; i < 100; ++i) {
		sam[i] = sam_init[i];
	}

	for (i = 0; i < 150; ++i) {
		pt[i] = pt_init[i];
		ps[i] = ps_init[i];
	}

	for (i = 0; i < 9; ++i) {
		for (j = 0; j < 150; ++j) {
			pd[i][j] = pd_init[i][j];
		}
	}

	for (i = 0; i < 10; ++i) {
		for (j = 0; j < 100; ++j) {
			pma[i][j] = pma_init[i][j];
		}
	}


	for (i = 0; i < 4; ++i) {
		for (j = 0; j < 100; ++j) {
			ptl[i][j] = ptl_init[i][j];
		}
	}

	return;
}

/*int main() 
{

	MSIS ms;
	Map map;
	ofstream msis;

	msis.open("C:\\EarthGRAM2010V4.0\\My Test\\MSIS.txt");
	double h1 = 0.0, phi1 = 33.31, thet1 = -87.55, elt1 = 0.0, ri1, g1, pmj1, dmj1, tmj1, umj1, vmj1, wmj1, n2nd, o2nd,
		ond, arnd, hend, hnd, wtmol, dmdz, nnd;
	double pi = 3.1415926535897931;
	double phir, pi180, g0, re, r0;
	double dh = 2.0;
	pi180 = pi / 180.0;
	phir = phi1 * pi180;

	ms.namelist("C:\\EarthGRAM2010V4.0\\My Test\\NameRef.txt");
	
	for (int i = 0; i < 71; i++){

		map.perts.rig(h1, phir, &g0, &g1, &re, &r0, &ri1);


		ms.msismod(h1, phi1, thet1, elt1, ri1, g1, &pmj1, &dmj1, &tmj1, &umj1, &vmj1, &wmj1, &n2nd, &o2nd, &ond,
			&arnd, &hend, &hnd, &wtmol, &dmdz, &nnd);

		//cout << h1 << '\n';
		//cout << tmj1 << '\n';
		//cout << umj1 << '\n';
		//cout << vmj1 << '\n';

		msis << h1 << " " << phi1 << " " << thet1 << " " << pmj1 << " " << dmj1 << " " <<
			tmj1 << " " << umj1 << " " << vmj1 << '\n';

		h1 = h1 + dh;

	}


	msis.close();
	//void rig(double ch, double phir, double *g0, double *g, double *re, double *r0, double *ri);
	void msismod(double h, double phi, double thet, double elt, double ri, double g, double *pmj,
		double *dmj, double *tmj, double *umj, double *vmj, double *wmj, double *n2nd,
		double *o2nd, double *ond, double *arnd, double *hend, double *hnd,
		double *wtmol, double *dmdz, double *nnd);
		

}
*/