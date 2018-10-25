//P. White
//Marshall Engineering Thermosphere (MET) model class

#include <iostream>
#include <fstream>

#include <sstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include "MET.h"


using namespace std;

Met::Met()
{
	initializeMemberVariables();
}


void Met::initializeMemberVariables()
{

	//Initialize member variables for Met class

	pi = 3.1415926535897931;
	pi180 = pi / 180;
	sec = 0.0;
	dtz = 0.0;
	ida = 0;
	iyr = 0;
	mn = 0;
	wtmol1 = 0.0;
	dmdz = 0.0;
}

double Met::molwt(double a)
{

	//molwt member function from Met class
	//molwt calculates the molecular weight for altitudes between 90 and 105
	//km according to equation (1) of SAO report 313.  Otherwise, molwt is 
	//set to unity.

	int i;
	double u, molwt_out;
	double b[7] = { 28.15204, -0.085586, 1.284e-04, -1.0056e-05,
		-1.021e-05, 1.5044e-06, 9.9826e-08 };

	if (a > 105.0) molwt_out = 1.0;
	else {
		u = a - 100.0;
		molwt_out = b[0];
		for (i = 1; i < 7; i++){
			molwt_out = molwt_out + b[i] * pow(u, i);

		}
	}

	return molwt_out;
}

double Met::temp(double alt, double tx, double t1, double t3, double t4, double a2)
{

	//temp member function from Met class
	//temp calculates the temperature at altitude alt using equation (10) for 
	//altitudes between 90 and 125 km and equation (13) for altitudes greater
	//than 125 km, from SAO Report 313.

	double u, temp_e10, bb;
	bb = 4.5e-6;
	u = alt - 125.0;
	if (u > 0.0){
		temp_e10 = tx + a2*atan(t1*u*(1.0 + bb*pow(u, 2.5)) / a2);
	}
	else {
		temp_e10 = tx + t1*u + t3*pow(u, 3.0) + t4*pow(u, 4.0);
	}

	return temp_e10;
}

void Met::gauss(int nmin, double z2, double tx, double t1, double t3, double t4, double a2, double gphi, double rphi, double *r)
{

	//gauss member function from Met class
	//Subdivide total integration-altitude range into intervals suitable for 
	//applying Gaussian Quadrature, set the number of points for integration
	//for each sub-interval, and then perform Gaussian Quadrature.

	int ngauss, i, j, k;
	double z, del, rr, d, a;
	int ng[8] = { 4, 5, 6, 6, 6, 6, 6, 6 };
	double altmin[9] = { 90.0, 105.0, 125.0, 160.0, 200.0, 300.0, 500.0, 1500.0, 2500.0 };

	//coefficients for gaussian quadrature
	const double c[8][6] = {
		{ 0.5555556, 0.3478548, 0.2369269, 0.1713245, 0.1294850, 0.1012285 },
		{ 0.8888889, 0.6521452, 0.4786287, 0.3607616, 0.2797054, 0.2223810 },
		{ 0.5555556, 0.6521452, 0.5688889, 0.4679139, 0.3818301, 0.3137067 },
		{ 0.0000000, 0.3478548, 0.4786287, 0.4679139, 0.4179592, 0.3626838 },
		{ 0.0000000, 0.0000000, 0.2369269, 0.3607616, 0.3818301, 0.3626838 },
		{ 0.0000000, 0.0000000, 0.0000000, 0.1713245, 0.2797054, 0.3137067 },
		{ 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.1294850, 0.2223810 },
		{ 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.1012285 }
	};

	//abscissas for gaussian quadrature
	const double x[8][6] =
	{
		{ -0.7745967, -0.8611363, -0.9061798, -0.9324695, -0.9491079, -0.9602899 },
		{ 0.0000000, -0.3399810, -0.5384693, -0.6612094, -0.7415312, -0.7966665 },
		{ 0.7745967, 0.3399810, 0.0000000, -0.2386192, -0.4058452, -0.5255324 },
		{ 0.0000000, 0.8611363, 0.5384693, 0.2386192, 0.0000000, -0.1834346 },
		{ 0.0000000, 0.0000000, 0.9061798, 0.6612094, 0.4058452, 0.1834346 },
		{ 0.0000000, 0.0000000, 0.0000000, 0.9324695, 0.7415312, 0.5255324 },
		{ 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.9491079, 0.7966665 },
		{ 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.9602899 }
	};

	*r = 0.0;

	for (k = nmin - 1; k < 8; k++){
		ngauss = ng[k];
		a = altmin[k];
		d = fmin(z2, altmin[k + 1]);
		rr = 0.0;
		del = 0.5*(d - a);
		j = ngauss - 3;

		for (i = 0; i < ngauss; i++){
			z = del*(x[i][j] + 1.0) + a;
			//calculation includes using molwt function and temp function
			rr = rr + c[i][j] * molwt(z) *(gphi / (pow((1.0 + z / rphi), 2.0))
				/ temp(z, tx, t1, t3, t4, a2));

		}

		rr = del*rr;

		*r = *r + rr;

		if (d == z2) {
			return;
		}

	}

	return;
}

void Met::met_07_jac(double z, double t, double *tz, double *an, double *ao2, double *ao, double *aa, double *ahe, double *ah, double *em, double *dens, double *dl, double gphi, double rphi)
{

	//met_07_jac member function from Met class
	//jac calculates the temperature tz, the total density dens, and its 
	//logarithm dl, the mean molecular weight em, the individual specie number
	//densities for n, o2, o, a, he, h at altitude z given the exospheric 
	//temperature t.  This member function uses the member function gauss,
	//temp, and molwt.

	double av = 6.02214179e+23, qn = 0.7811, qo2 = 0.20955, qa = 0.009343, qhe = 1.289e-05, rgas = 8.314472, t0 = 183.0;
	int i;
	double s, a1, par, factor, td, d, t4, t3, t1, tx_t0, a2, tx, r;

	double alpha[6] = { 0.0, 0.0, 0.0, 0.0, -0.380, 0.0 };
	double ei[6] = { 28.0134, 31.9988, 15.9994, 39.948, 4.0026, 1.00797 };
	double di[6], dit[6];



	tx = 444.3807 + 0.02385*t - 392.8292*exp(-0.0021357*t);
	a2 = 2.0*(t - tx) / pi;
	tx_t0 = tx - t0;
	t1 = 1.9*tx_t0 / 35.0;
	t3 = -1.7*tx_t0 / pow(35.0, 3);
	t4 = -0.8*tx_t0 / pow(35.0, 4);
	*tz = temp(z, tx, t1, t3, t4, a2);

	//Section 1

	d = min(z, 105.0);

	//Integrate gM/T from 90 to minimum of Z or 105 km

	gauss(1, d, tx, t1, t3, t4, a2, gphi, rphi, &r);

	//The number 2.1926e-08 = density x temperature / mean molecular weight 
	//at 90 km.

	*em = molwt(d);
	td = temp(d, tx, t1, t3, t4, a2);
	*dens = 2.1926e-08**em*exp(-r / rgas) / td;
	factor = av**dens;
	par = factor / (*em);
	factor = factor / 28.96;

	//for altitudes below and at 105 calculate the individual specie number 
	//densities from the mean molecular weight and total density

	if (z <= 105.0){

		*dl = log10(*dens);
		*an = log10(qn*factor);
		*aa = log10(qa*factor);
		*ahe = log10(qhe*factor);
		*ao = log10(2.0*par*(1.0 - *em / 28.96));
		*ao2 = log10(par*(*em*((1.0 + qo2) / 28.96) - 1.0));
		*ah = 0.0;

		return;
	}

	//Section 2: This section is only for altitudes above 105 km

	//Note that having reached this section means that D in section 1 is 
	//105 km.

	//Calculate individual specie number densities from the total and mean 
	//molecular weight at 105 km.

	di[0] = qn*factor;
	di[1] = par*(*em*(1.0 + qo2) / 28.96 - 1.0);
	di[2] = 2.0*par*(1 - *em / 28.96);
	di[3] = qa*factor;
	di[4] = qhe*factor;


	//Integrate g/T from 105 km z km:

	gauss(2, z, tx, t1, t3, t4, a2, gphi, rphi, &r);

	for (i = 0; i < 5; i++){
		dit[i] = di[i] * pow((td / (*tz)), (1.0 + alpha[i])) * exp((-ei[i] * r / rgas));

		if (dit[i] <= 0.0){
			dit[i] = 1.0e-06;
		}
	}

	//This section calculates atomic hydrogen densities above 500 km altitude.
	//Below this altitude, H densities are set to 10**6.

	// Section 3

	if (z > 500.0){

		a1 = 500.0;
		s = temp(a1, tx, t1, t3, t4, a2);

		di[5] = pow(10.0, (73.13 - 39.4*log10(s) + 5.5*log10(s)*log10(s)));

		gauss(7, z, tx, t1, t3, t4, a2, gphi, rphi, &r);

		dit[5] = di[5] * (s / (*tz))*exp((-ei[5] * r / rgas));

	}

	else {
		dit[5] = 1.0;
	}

	//For altitudes greater than 105 km, calculate total density and mean 
	//molecular weight from individual specie number densities.
	*dens = 0;
	for (i = 0; i < 6; i++){

		*dens = *dens + ei[i] * dit[i] / av;
	}

	*em = *dens*av / (dit[0] + dit[1] + dit[2] + dit[3] + dit[4] + dit[5]);
	*dl = log10(*dens);
	*an = log10(dit[0]);
	*ao2 = log10(dit[1]);
	*ao = log10(dit[2]);
	*aa = log10(dit[3]);
	*ahe = log10(dit[4]);
	*ah = log10(dit[5]);


}


void Met::caltojul(int iy, int im, int id, int ihour, int imin, double sec, double *xjd)
{
	//caltojul member function from Met class
	//Compute Julian day by method of Meeus, Astronomical Algorithms, 
	//2nd Edition, 1998, page 61. Inputs are year iy, month im, day of month
	//id, and time of day in minutes hours and seconds. Output is Julian
	//day, xjd.

	int y, m, a, b;
	double d;

	y = iy;
	m = im;

	//Consider Jan or Feb as if months 13 and 14 of previous year
	if (im <= 2){
		y = iy - 1;
		m = im + 12;
	}

	//Compute day of month plus fractional part
	d = id + (ihour / 24) + (imin / 1440) + (sec / 86400);
	a = int(0.01*y);
	b = 2 - a + int(0.25*a);

	//Compute julian day with fractional part
	*xjd = int(365.25*(y + 4716)) + int(30.6001*(m + 1)) + d + b - 1524.5;

	return;
}


void Met::slv(double alt, double xlat, double day, double *den)
{

	//slv member function from Met class
	//slv computes the seasonal-latitudinal variation of density in the lower
	//thermosphere in accordance with L. Jacchia, SAO 332, 1971.  This affects
	//the densities between 90 and 170 km.  This member function need not be 
	//called for densities above 170 km, because no effect is observed.

	//The variation should be computed after the calculation of density due 
	//to temperature variations and the density must be in the form of a 
	//base 10 log.  No adjustments are made to the temperature or consituent
	//number densities in the region affected by this variation.

	double d, s, sp, p, y, x, z;


	*den = 0.0;

	//check if altitude exceeds 170 km
	if (alt > 170.0) return;

	//Compute density change in lower thermosphere
	z = alt - 90.0;
	x = -0.0013*z*z;
	y = 0.0172*day + 1.72;
	p = sin(y);
	sp = pow(sin(xlat), 2.0);
	s = 0.014*z*exp(x);
	d = s*p*sp;

	//check to compute absolute value of xlat
	if (xlat < 0.0) d = -d;

	*den = d;

	return;
}

void Met::slvh(double xlat, double sda, double *den, double *denhe)
{

	//slvh member function from Met class
	//slvh computes the seasonal-latitudinal variation of the helium number
	//density according to L. Jacchia, SAO 332, 1971.  This correction is not
	//important below about 500 km.

	double drho, rho, del, d1, dhe, y, x, b, a, d0;

	d0 = pow(10.0, *denhe);
	a = abs(0.65*(sda / 0.40909079));
	b = 0.5*xlat;

	//Check to compute absolute value of 'b'

	if (sda < 0.0) b = -b;

	//Compute x, y, dhe, and denhe

	x = 0.7854 - b;
	y = pow(sin(x), 3);
	dhe = a*(y - 0.35356);
	*denhe = *denhe + dhe;

	//Compute helium number density change

	d1 = pow(10.0, *denhe);
	del = d1 - d0;
	rho = pow(10.0, *den);
	drho = 6.646e-24 * del;
	rho = rho + drho;
	*den = log10(rho);
}

void Met::tme(int mn, int ida, int iyr, int ihr, double xmin, double xlng, double *xlat, double *sda, double *sha, double *dd, double *dy, double *sra, double *raloc, double *xmjd)
{

	//tme member function from MET class
	//'tme' performs the calculations of the solar declination and solar
	//hour angle.

	int i, id;
	double degs, gmt, dp2k, sunlon, g, eclon, oblq, alpha, temp, year, xhr, eqt;

	int iday[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };



	*xlat = *xlat*pi180;
	degs = 360.0;
	year = 365.0;

	//Compute day of year, dd
	//Correct year 2000 leap year error in original MET code
	if (iyr % 4 == 0)
	{
		iday[1] = 29;
		if (iyr % 100 == 0 && iyr % 400 != 0) iday[1] = 28;
	}
	if (iday[1] == 29) year = 366.0;
	id = 0;

	if (mn > 1)
	{

		for (i = 0; i < mn - 1; i++){
			id += iday[i];
		}
	}


	id += ida;
	*dd = id + ihr / 24.0 + xmin / 1440.0 - 1.0;
	*dy = *dd / year;

	//Compute mean Julian date at 00 UT
	caltojul(iyr, mn, ida, 0, 0, 0.0, xmjd);

	//Other calculations
	xhr = ihr;
	gmt = 60.0*xhr + xmin;
	dp2k = *xmjd - 2451545.0 + gmt / 1440.0;

	//Compute mean longitude of Sun, degrees
	sunlon = fmod(280.466 + 0.9856474*dp2k, degs);
	if (sunlon < 0.0) sunlon += degs;

	//Compute mean anamoly, radians
	g = pi180*fmod(357.5280 + 0.9856003*dp2k, degs);
	if (g < 0.0) g += 2.0*pi;

	//Compute ecliptic longitude, radians
	eclon = pi180*(sunlon + 1.915*sin(g) + 0.02*sin(2.0*g));

	//Compute obliquity, radians
	oblq = pi180*(23.44 - 0.0000004*dp2k);

	//Compute right ascension, radians
	alpha = atan(cos(oblq)*tan(eclon));

	//Put alpha in same quadrant as eclon
	alpha = abs(alpha);
	temp = abs(eclon);

	if (temp <= pi && temp > (pi / 2.0)){
		alpha = pi - alpha;
	}
	else if (temp <= 1.5*pi && temp > pi){
		alpha = pi + alpha;
	}
	else if (temp > 1.5*pi){
		alpha = 2.0*pi - alpha;
	}
	if (alpha < 0.0) alpha = -alpha;

	//Compute declination of sun, radians

	*sda = asin(sin(oblq)*sin(eclon));

	//Compute equation of time, radians

	eqt = pi180*sunlon - alpha;

	//Compute solar hour angle - sha (in rad)
	*sha = fmod(pi180*(gmt / 4.0 + xlng) + eqt + pi, (2 * pi));

	if (*sha > pi){
		*sha = *sha - (2 * pi);
	}
	if (*sha < -pi) {
		*sha = *sha + (2 * pi);
	}

	//rename Solar RA for output
	*sra = alpha;

	//Compute local RA
	*raloc = alpha + *sha;
	if (*raloc < 0.0) *raloc = *raloc + (2 * pi);
	if (*raloc >(2 * pi)) *raloc = *raloc + (2 * pi);
}



void Met::met07(double *indata, double *outdata)
{

	//met07 member function from Met class
	//met07 is the driving program for the MET model
	//The atmospheric model is a modified Jacchia 1970 model and is given
	//in the subroutine J70.  All of the other subroutines were designed to 
	//allow flexible use of the model so that various input parameters could
	//be varied within a driving program with very little software development.
	//This, for example, driving routines can be written quite easily to 
	//facilitate the plotting of output as line or contour plots.

	double rgas = 8.314472e+03, bfh = 440.0, av = 6.02214179e+26, r0 = 8314.472;
	int i, i1, ida, mn;
	double p, em, sumn, summn, tz, te, dens, fdhel, dlg2, dlg1, dhel2, dhel1, dl, den, dummy;
	double gphi, rphi, denlg, gi, f10b, f10, fdlg, xlng, xlat, z, sda, sha, dy;
	double xmin, avx, dd, sra, raloc, xmjd;

	double a[6];

	//Molecular weights
	double ei[6] = { 28.0134, 31.9988, 15.9994, 39.948, 4.0026, 1.00797 };

	//Set parameters to INDATA values

	z = indata[0];
	xlat = indata[1];

	xlng = indata[2];
	iyr = int(indata[3]);
	mn = int(indata[4]);
	ida = int(indata[5]);
	ihr = int(indata[6]);
	xmin = indata[7];
	i1 = int(indata[8]);
	f10 = indata[9];
	f10b = indata[10];
	gi = indata[11];
	gphi = indata[12];
	rphi = indata[13];



	tme(mn, ida, iyr, ihr, xmin, xlng, &xlat, &sda, &sha, &dd, &dy, &sra, &raloc, &xmjd);

	tinf(i1, f10, f10b, gi, xlat, sda, sha, dy, &te);

	met_07_jac(z, te, &tz, &a[0], &a[1], &a[2], &a[3], &a[4], &a[5], &em, &dens, &dl, gphi, rphi);

	denlg = 0.0;
	dummy = dl;
	den = dl;

	if (z <= 170.0){
		slv(z, xlat, dd, &dummy);

		denlg = dummy;

	}

	//Fair helium number density between base fairing height (bfh) and 500 km

	if (z >= 500.0){

		slvh(xlat, sda, &den, &a[4]);

		dl = den;
	}
	else if (z > bfh){
		dhel1 = a[4];
		dhel2 = a[4];
		dlg1 = dl;
		dlg2 = dl;
		slvh(xlat, sda, &dlg2, &dhel2);

		fair5(dhel1, dhel2, dlg1, dlg2, z, &fdhel, &fdlg);

		dl = fdlg;
		a[4] = fdhel;

	}

	dl = dl + denlg;


	dens = pow(10, dl);
	xlat = xlat / pi180;

	//Fill out outdata array
	outdata[0] = te;
	outdata[1] = tz;
	for (i = 0; i < 6; i++){
		outdata[i + 2] = 1.0e+6*pow(10, a[i]);
	}
	outdata[8] = em;
	outdata[9] = dens*1000.0;
	outdata[10] = dl + 3.0;
	p = outdata[9] * rgas*tz / em;
	outdata[11] = p;

	//Compute sumn = sum of number densities and sumn = sum of product of 
	//number density and molecular weight
	sumn = 0.0;
	summn = 0.0;
	for (i = 0; i < 6; i++){
		sumn += outdata[i + 2];
		summn += ei[i] * outdata[i + 2];
	}

	if (z < 170.0){
		//Adjust number densities to be consistent with mass density
		for (i = 0; i < 6; i++){
			avx = av*outdata[9] / summn;
			outdata[i + 2] *= avx;
		}
	}

	else if (z > 440.0){
		//Adjust molecular weight to be consistent with number densities
		outdata[8] = summn / sumn;
		//Re-compute pressure with revised molecular weight
		outdata[11] = outdata[9] * r0*outdata[1] / outdata[8];
	}

}


void Met::tinf(int i1, double f10, double f10b, double gi, double xlat, double sda, double sha, double dy, double *te)
{

	//tinf member function from Met class
	//tinf calculates the exospheric temperature according to L. Jacchia SAO
	//313, 1970.

	double e4, e6, e10, ts, g2, g1, tau1, g3, tg, tl, tv, b2, b1, a3, a2, a1, tau, theta, eta, tc, twopi;
	double xm = 2.5, xnn = 3.0, c1 = 383.0, c2 = 3.32, c3 = 1.8, d1 = 28.0, d2 = 0.03, d3 = 1.0, d4 = 100.0, d5 = -0.08, e1 = 2.41, e2 = 0.349, e3 = 0.206, e5 = 3.9531708, e7 = 4.3214352, e8 = 0.1145, e9 = 0.5, e11 = 5.974262, e12 = 2.16, beta = -0.6457718, gamma = 0.7504916, p = 0.1047198, re = 0.31;

	twopi = 2.0 * pi;

	//Ei are semiannual variation variables
	e6 = 2.0*twopi;
	e4 = twopi;
	e10 = twopi;

	//solar activity variation
	tc = c1 + c2*f10b + c3*(f10 - f10b);

	//diurnal variation
	eta = 0.5*abs(xlat - sda);
	theta = 0.5*abs(xlat + sda);
	tau = sha + beta + p*sin(sha + gamma);

	if (tau > pi) tau = tau - twopi;
	if (tau < (-pi)) tau = tau + twopi;

	a1 = pow(sin(theta), xm);
	a2 = pow(cos(eta), xm);
	a3 = pow(cos(tau / 2.0), xnn);
	b1 = 1.0 + re*a1;
	b2 = (a2 - a1) / b1;
	tv = b1*(1.0 + re*b2*a3);
	tl = tc*tv;

	//geomagnetic variation

	if (i1 == 1){
		tg = d1*gi + d2*exp(gi);
	}
	else {
		tg = d3*gi + d4*(1.0 - exp(d5*gi));
	}

	//semiannual variation
	g3 = 0.5*(1.0 + sin(e10*dy + e11));
	g3 = pow(g3, e12);
	tau1 = dy + e8*(g3 - e9);

	g1 = e2 + e3*sin(e4*tau1 + e5);
	g2 = sin(e6*tau1 + e7);
	ts = e1 + f10b*g1*g2;

	//exospheric temperature
	*te = tl + tg + ts;

	return;
}

void Met::jacch(int iyr, int m, int ida, int ihr, double xmin, double z, double phi, double thet, double f10, double f10b, double ap, double *ph, double *dh, double *th, double *n2nd, double *o2nd, double *ond, double *arnd, double *hend, double *hnd, double *wtmol, double *tex)
{

	//jacch member function from Met class

	double gphi, rphi, b, a, ri, r0, g, phir;
	double indata[14];
	double outdata[12];


	indata[0] = z;
	indata[1] = phi;
	indata[2] = thet;
	indata[3] = double(iyr);
	indata[4] = double(m);
	indata[5] = double(ida);
	indata[6] = double(ihr);
	indata[7] = xmin;
	indata[8] = 2.0;
	indata[9] = f10;
	indata[10] = f10b;
	indata[11] = ap;
	phir = phi*pi180;
	rig(z, phir, &gphi, &g, &rphi, &r0, &ri, &a, &b);
	indata[12] = gphi;
	indata[13] = rphi;

	met07(indata, outdata);

	*ph = outdata[11];
	*dh = outdata[9];
	*th = outdata[1];
	*n2nd = outdata[2];
	*o2nd = outdata[3];
	*ond = outdata[4];
	*arnd = outdata[5];
	*hend = outdata[6];
	*hnd = outdata[7];
	*wtmol = outdata[8];
	*tex = outdata[0];

}


void Met::jacmod(double h, double phi, double thet, double elt, double *pmj, double *dmj, double *tmj, double *umj, double *vmj, double *wmj, double *arnda, double *henda, double *hnda, double *o2nda, double *n2nda, double *onda, double *wtmola)
{

	//jacmod member function from Met class
	//MET (Jacchia) model driver routine to evaluate mean p,d,t,u,v,w

	double th, dh, ph;
	double xmin, tjx, pjx, djx, dphij, phn, phe, the, thn, tb, d7, uplus,
		vplus, wplus, uminus, vminus, wminus, factint, d5, d3, db, pb, dhe,
		texp, d6, d4, d2, d1, dhn;
	double  phir;
	double gphi, rphi, r0, ri, a, b, g;
	double dy5, dx5;


	double  dpy, dpx, dty, dtx;



	imin = mino + int(elt) / 60;
	sec = seco + elt;

	sec = fmod(sec, 60.0);

	ihr = ihro + imin / 60;
	imin = imin % 60;


	phir = phi*pi180;
	double n2nd1, o2nd1, ond1, arnd1, hend1, hnd1, wtmol1, tex1, flat = 15,
		ugh, vgh, wgh;

	//Distances for 5 degrees of geocentric latitude, longitude
	rig(h, phir, &gphi, &g, &rphi, &r0, &ri, &a, &b);

	dy5 = 5000.0*ri*pi180;
	dx5 = dy5*cos(phi*pi180);
	dx5 = max(2000.0, dx5);

	//Following is the pure Jacchia height range section
	//Jacchia values at current position

	xmin = imin + sec / 60.0;

	jacch(iyr, mn, ida, ihr, xmin, h, phi, thet, f10, f10b, ap, pmj, dmj, tmj,
		&n2nd1, &o2nd1, &ond1, &arnd1, &hend1, &hnd1, &wtmol1, &tex1);

	th = *tmj;
	ph = *pmj;
	dh = *dmj;


	dphij = 5.0;
	if (phi >= 85.0) {
		dphij = -5.0;
	}
	if (abs(phi) > flat){
		//Jacchia values at current position+5 degrees lat
		jacch(iyr, mn, ida, ihr, xmin, h, (phi + dphij), thet, f10, f10b, ap, &phn, &dhn, &thn, &d1, &d2, &d3, &d4, &d5, &d6, &d7, &texp);
		//Jacchia values at current position+5 degrees lon
		jacch(iyr, mn, ida, ihr, xmin, h, phi, thet + 5.0, f10, f10b, ap, &phe, &dhe, &the, &d1, &d2, &d3, &d4, &d5, &d6, &d7, &texp);

		//dp/dy and dp/dx for geostrophic wind
		dpy = phn - *pmj;
		dpx = phe - *pmj;

		//dt/dx, dt/dy and dt/dz for vertical wind
		dtx = the - *tmj;
		dty = thn - *tmj;
		if (dphij < 0.0){
			dpy = -dpy;
			dty = -dty;
		}
		jacch(iyr, mn, ida, ihr, xmin, h + 1.0, phi, thet, f10, f10b, ap, &pb, &db, &tb, &d1, &d2, &d3, &d4, &d5, &d6, &d7, &texp);
		dtz = (tb - *tmj) / 1000.0;
		dmdz = (d7 - wtmol1) / 1000.0;
		wind(ph, dh, th, h, g, phi, ri, dpx, dpy, dtx, dty, dtz, &ugh, &vgh, &wgh);
		//change notation for output
		*umj = ugh;
		*vmj = vgh;
		*wmj = wgh;

	}
	else
	{
		//Jacchia values at +flat and +flat+5
		jacch(iyr, mn, ida, ihr, xmin, h, flat, thet, f10, f10b, ap, &pjx, &djx, &tjx, &d1, &d2, &d3, &d4, &d5, &d6, &d7, &texp);
		jacch(iyr, mn, ida, ihr, xmin, h, flat + dphij, thet, f10, f10b, ap, &phn, &dhn, &thn, &d1, &d2, &d3, &d4, &d5, &d6, &d7, &texp);
		//Jacchia values at +flat and lon+5
		jacch(iyr, mn, ida, ihr, xmin, h, flat, thet + 5.0, f10, f10b, ap, &phe, &dhe, &the, &d1, &d2, &d3, &d4, &d5, &d6, &d7, &texp);
		//dp/py and dp/dx for geostrophic wind
		dpy = phn - pjx;
		dpx = phe - pjx;
		//dt/dx, dt/dy and dt/dz for vertical wind
		dtx = the - tjx;
		dty = thn - tjx;
		if (dphij < 0.0){
			dpy = -dpy;
			dty = -dty;
		}
		jacch(iyr, mn, ida, ihr, xmin, h + 1.0, phi, thet, f10, f10b, ap, &pb, &db, &tb, &d1, &d2, &d3, &d4, &d5, &d6, &d7, &texp);
		dtz = (tb - th) / 1000.0;
		dmdz = (d7 - wtmol1) / 1000.0;
		wind(ph, dh, th, h, g, flat, ri, dpx, dpy, dtx, dty, dtz, &ugh, &vgh, &wgh);

		uplus = ugh;
		vplus = vgh;
		wplus = wgh;

		//Jacchia values at -flat and -flat+5
		jacch(iyr, mn, ida, ihr, xmin, h, (-flat), thet, f10, f10b, ap, &pjx, &djx, &tjx, &d1, &d2, &d3, &d4, &d5, &d6, &d7, &texp);
		jacch(iyr, mn, ida, ihr, xmin, h, (-flat) + dphij, thet, f10, f10b, ap, &phn, &dhn, &thn, &d1, &d2, &d3, &d4, &d5, &d6, &d7, &texp);
		//Jacchia values at -fla and lon+5
		jacch(iyr, mn, ida, ihr, xmin, h, (-flat), thet + 5.0, f10, f10b, ap, &phe, &dhe, &the, &d1, &d2, &d3, &d4, &d5, &d6, &d7, &texp);
		//dp/dy and dp/dx for geostrophic wind
		dpy = phn - pjx;
		dpx = phe - pjx;
		//dt/dx, dt/dy and dt/dz for vertical wind
		dtx = the - tjx;
		dty = thn - tjx;
		if (dphij < 0.0){
			dpy = -dpy;
			dty = -dty;
		}
		jacch(iyr, mn, ida, ihr, xmin, h + 1.0, (-flat), thet, f10, f10b, ap, &pb, &db, &tb, &d1, &d2, &d3, &d4, &d5, &d6, &d7, &texp);
		dtz = (tb - th) / 1000.0;
		wind(ph, dh, th, h, g, (-flat), ri, dpx, dpy, dtx, dty, dtz, &ugh, &vgh, &wgh);

		//Change notation for output
		uminus = ugh;
		vminus = vgh;
		wminus = wgh;

		//Interpolate wind to latitude
		factint = 0.5*(phi / flat + 1.0);
		*umj = uminus + factint*(uplus - uminus);
		*vmj = vminus + factint*(vplus - vminus);
		*wmj = wminus + factint*(wplus - wminus);
	}

	//Atmospheric Constituents
	*arnda = arnd1;
	*henda = hend1;
	*hnda = hnd1;
	*o2nda = o2nd1;
	*n2nda = n2nd1;
	*onda = ond1;
	*wtmola = wtmol1;

}

void Met::wind(double ph, double dh, double th, double h, double g, double phid, double ri, double dpx, double dpy, double dtx, double dty, double dtz, double *ugh, double *vgh, double *wgh)
{

	//wind member function from Met class
	//Computes winds for MET (Jacchia) section (geostrophic winds with 
	//viscosity modifications at high altitude)

	double beta = 1.458e-6, sval = 110.4;
	double coriol, cp, denom, dpdx, dpdy, visc, viscfac, vls, splim, dx, dy, fac;

	fac = 180 / pi;

	//Avoid cases with 0 temperature or density
	if (dh <= 0.0 || th <= 0.0) {
		*ugh = 0.0;
		*vgh = 0.0;
		*wgh = 0.0;
	}
	else {
		//distances for 5 degrees separation points
		dy = 5000.0*ri / fac;
		dx = dy*cos(pi180*phid);

		//Coriolis factor
		coriol = pi180*sin(pi180*phid) / 120.0;

		//Viscosity coefficients
		visc = beta*pow(th, 1.5) / (th + sval);
		vls = 5.3 + 0.0622*pow(h, 1.5);
		vls = min(h, vls);
		viscfac = visc / (1.0e6*dh*pow(vls, 2));

		//Horizontal pressure gradients
		dpdx = dpx / (dx*dh);
		dpdy = dpy / (dy*dh);
		denom = pow(coriol, 2) + pow(viscfac, 2);

		//Viscosity-modified geostrophic winds (horizontal components)
		*ugh = ((-coriol*dpdy) - (viscfac*dpdx)) / denom;
		*vgh = ((coriol*dpdx) - (viscfac*dpdy)) / denom;

		//Insure wind component magnitudes < 0.7 times sound speed
		splim = 0.7*sqrt(1.4*ph / dh);
		if (fabs(*ugh) > splim) {
			*ugh = copysign(splim, *ugh);
		}

		if (fabs(*vgh) > splim) {
			*vgh = copysign(splim, *vgh);
		}

		cp = 7.0*ph / (2.0*dh*th);

		//Mean vertical wind from Montgomery stream function
		*wgh = -cp*(*ugh*dtx / dx + *vgh*dty / dy) / (g + cp*dtz);
	}


}

void Met::rig(double ch, double phir, double *g0, double *g, double *re, double *r0, double *ri, double *a, double *b)
{

	//rig member function from Met class
	//Computes surface gravity g0 (m/s^2) and effective Earth radius re (km) from 
	//input geocentric latitude phir (radians).  Also computes gravity g (m/s^2)
	//and total radius ri (km) at input height ch (km).

	double g0a = 9.80616, g0b = 0.0026373, g0c = 0.0000059, rea = 3.085462e-3, reb = 2.27e-6, rec = 2.0e-9;
	double eps, c2phi, cphi2, c4phi;

	//Parameters for computing effective Earth radius

	//Set up cosines:  cphi2 = [cos(phir)]**2, ...
	cphi2 = pow(cos(phir), 2);
	c2phi = 2.0*cphi2 - 1.0;
	c4phi = 8.0*cphi2*(cphi2 - 1.0) + 1.0;

	*a = 6378.137;
	*b = 6356.752314;
	eps = 1.0 - pow(*b / *a, 2);

	//Compute Earth radius, r0
	*r0 = *b / sqrt(1.0 - eps*cphi2);

	//Compute surface gravity
	*g0 = g0a*(1.0 - g0b*c2phi + g0c*pow(c2phi, 2));

	//Compute effective Earth radius, re
	*re = 2.0**g0 / (rea + reb*c2phi - rec*c4phi);

	//Compute g at height ch
	*g = *g0 / pow((1.0 + ch / (*re)), 2);

	//Compute radius at height ch
	*ri = *r0 + ch;

	return;
}

void Met::fair5(double dhel1, double dhel2, double dlg1, double dlg2, double h, double *fdhel, double *fdlg)
{

	//fair5 member function from Met class
	//This subroutine fairs between the region above 500 km, which invokes
	//the seasonal-latitudinal variation of the helium number density (slvh),
	//and the region below, which does not invoke any seasonal-latitudinal 
	//variation at all.


	double bfh = 440.0;
	double szi, czi, hi;


	//height fairing factor
	hi = 1.50*pi180*(h - bfh);

	//Non-slvh fairing coefficient
	czi = pow(cos(hi), 2);

	//slvh fairing factor
	szi = 1.0 - czi;

	//Faired density
	*fdlg = dlg1*czi + dlg2*szi;

	//Faired helium number density
	*fdhel = dhel1*czi + dhel2*szi;


}


void Met::namelist(string namef)
{

	//namelist member function from Met class
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

	//Set initial time variables
	imin = mino;
	ihr = ihro;
	sec = seco;



}
