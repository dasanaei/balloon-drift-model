//P. White
//JB2008 class

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include "JB2008.h"
#include <cstdlib>
#include <cmath>

using namespace std;


JB2008::JB2008()
{

	intializeMemberVariables();

}

void JB2008::intializeMemberVariables()
{

	//Initialize member variables from JB2008 class

	mn = 0, ida = 0, iyr = 0, ihr = 0;
	dtz = 0.0;
}

double JB2008::xambar(double z)
{

	//xambar member function from JB2008 class

	//Evalutes equation (1) 
	const double c[7] = { 28.15204, -8.5586e-2, 1.284e-4, -1.0056e-5, -1.021e-5, 1.5044e-6, 9.9826e-8 };
	double dz, amb, xambar_out;
	int i, j;

	dz = z - 100.0;

	amb = c[6];

	for (i = 0; i < 6; i++){
		j = 5 - i;
		amb = dz*amb + c[j];

	}
	xambar_out = amb;

	return xambar_out;
}

void JB2008::semian08(int iyr, double day, double ht, double f10b, double s10b, double xm10b, double *fzz, double *gtz, double *drlog)
{
	//semian08 member function from JB2008 class
	//Compute semiannual variation (delta log rho)
	//Input iyr, day, height, f10b, s10b, m10b fsmb
	//Output fz, gt, and del log rho value

	double const twopi = 6.2831853072, gtm[10] = { -0.3633, 0.8506e-01, 0.2401, -0.1897, -0.2554, -0.1790e-01, 0.565e-03, -0.6407e-03, -0.3418e-02, -0.1252e-02 },
		fzm[5] = { 0.2689, -0.1176e-01, 0.2782e-01, -0.2782e-01, 0.3470e-03 };
	double fsmb, htz, yrln, tau, sin1p, cos1p, sin2p, cos2p;

	//Compute new 81-day centered solar index for fz
	fsmb = 1.00*f10b - 0.7*s10b - 0.04*xm10b;

	htz = ht / 1000.0;

	*fzz = fzm[0] + fzm[1] * fsmb + fzm[2] * fsmb*htz + fzm[3] * fsmb*pow(htz, 2) +
		fzm[4] * pow(fsmb, 2)*htz;

	//Compute daily 81-day centered solar index for gt
	fsmb = 1.00*f10b - 0.75*s10b - 0.37*xm10b;

	yrln = 365.0;
	if (iyr % 4 == 0)yrln = 366.0;

	tau = (day - 1.0) / yrln;
	sin1p = sin(twopi*tau);
	cos1p = cos(twopi*tau);
	sin2p = sin(2.0*twopi*tau);
	cos2p = cos(2.0*twopi*tau);

	*gtz = gtm[0] + gtm[1] * sin1p + gtm[2] * cos1p + gtm[3] * sin2p +
		gtm[4] * cos2p + gtm[5] * fsmb + gtm[6] * fsmb*sin1p + gtm[7] * fsmb*cos1p
		+ gtm[8] * fsmb*sin2p + gtm[9] * fsmb*cos2p;

	if (*fzz < 1.0e-6) *fzz = 1.0e-06;

	*drlog = *fzz**gtz;

	return;

}

void JB2008::dtsub(double f10, double xlst, double xlat, double zht, double *dtc)
{
	//dtsub member function from JB2008 class
	//Compute dTc correction for Jacchia-Bowman model

	double const b[19] = { -0.457512297e+01, -0.512114909e+01, -0.693003609e+02, 0.203716701e+03,
		0.703316291e+03, -0.194349234e+04, 0.110651308e+04, -0.174378996e+03, 0.188594601e+04,
		-0.709371517e+04, 0.922454523e+04, -0.384508073e+04, -0.645841789e+01, 0.409703319e+02,
		-0.482006560e+03, 0.181870931e+04, -0.237389204e+04, 0.996703815e+03, 0.361416936e+02 };

	double const c[23] = { -0.155986211e+02, -0.512114909e+01, -0.693003609e+02, 0.203716701e+03,
		0.703316291e+03, -0.194349234e+04, 0.110651308e+04, -0.220835117e+03, 0.143256989e+04,
		-0.318481844e+04, 0.328981513e+04, -0.135332119e+04, 0.199956489e+02, -0.127093998e+02,
		0.212825156e+02, -0.275555432e+01, 0.110234982e+02, 0.148881951e+03, -0.751640284e+03,
		0.637876542e+03, 0.127093998e+02, -0.212825156e+02, 0.275555432e+01 };

	double tx, ycs, f, h, dtc200, dtc200dz, sum, cc, dd, zp, aa, bb, dtc300, dtc300dz, hp;

	*dtc = 0.0;
	tx = xlst / 24.0;
	ycs = cos(xlat);
	f = (f10 - 100.0) / 100.0;

	//Calculate dTc

	if ((zht >= 120.0) & (zht <= 200.0)){
		h = (zht - 200.0) / 50.0;
		dtc200 = c[16] + c[17] * tx*ycs + c[18] * pow(tx, 2)*ycs +
			c[19] * pow(tx, 3)*ycs + c[20] * f*ycs + c[21] * tx*f*ycs +
			c[22] * pow(tx, 2)*f*ycs;
		sum = c[0] + b[1] * f + c[2] * tx*f + c[3] * pow(tx, 2)*f +
			c[4] * pow(tx, 3)*f + c[5] * pow(tx, 4)*f + c[6] * pow(tx, 5)*f +
			c[7] * tx*ycs + c[8] * pow(tx, 2)*ycs + c[9] * pow(tx, 3)*ycs +
			c[10] * pow(tx, 4)*ycs + c[11] * pow(tx, 5)*ycs + c[12] * ycs +
			c[13] * f*ycs + c[14] * tx*f*ycs + c[15] * pow(tx, 2)*f*ycs;
		dtc200dz = sum;
		cc = 3.0*dtc200 - dtc200dz;
		dd = dtc200 - cc;
		zp = (zht - 120.0) / 80.0;
		*dtc = cc*zp*zp + dd*zp*zp*zp;

	}

	if ((zht > 200.0) & (zht <= 240.0)){
		h = (zht - 200.0) / 50.0;
		sum = c[0] * h + b[1] * f*h + c[2] * tx*f*h + c[3] * pow(tx, 2)*f*h +
			c[4] * pow(tx, 3)*f*h + c[5] * pow(tx, 4)*f*h + c[6] * pow(tx, 5)*f*h +
			c[7] * tx*ycs*h + c[8] * pow(tx, 2)*ycs*h + c[9] * pow(tx, 3)*ycs*h +
			c[10] * pow(tx, 4)*ycs*h + c[11] * pow(tx, 5)*ycs*h + c[12] * ycs*h +
			c[13] * f*ycs*h + c[14] * tx*f*ycs*h + c[15] * pow(tx, 2)*f*ycs*h +
			c[16] + c[17] * tx*ycs + c[18] * pow(tx, 2)*ycs + c[19] * pow(tx, 3)*ycs +
			c[20] * f*ycs + c[21] * tx*f*ycs + c[22] * pow(tx, 2)*f*ycs;
		*dtc = sum;

	}

	if ((zht > 240.0) & (zht <= 300.0)){
		h = (40.0) / 50.0;
		sum = c[0] * h + b[1] * f*h + c[2] * tx*f*h + c[3] * pow(tx, 2)*f*h +
			c[4] * pow(tx, 3)*f*h + c[5] * pow(tx, 4)*f*h + c[6] * pow(tx, 5)*f*h +
			c[7] * tx*ycs*h + c[8] * pow(tx, 2)*ycs*h + c[9] * pow(tx, 3)*ycs*h +
			c[10] * pow(tx, 4)*ycs*h + c[11] * pow(tx, 5)*ycs*h + c[12] * ycs*h +
			c[13] * f*ycs*h + c[14] * tx*f*ycs*h + c[15] * pow(tx, 2)*f*ycs*h +
			c[16] + c[17] * tx*ycs + c[18] * pow(tx, 2)*ycs + c[19] * pow(tx, 3)*ycs +
			c[20] * f*ycs + c[21] * tx*f*ycs + c[22] * pow(tx, 2)*f*ycs;
		aa = sum;
		bb = c[0] + b[1] * f + c[2] * tx*f + c[3] * pow(tx, 2)*f +
			c[4] * pow(tx, 3)*f + c[5] * pow(tx, 4)*f + c[6] * pow(tx, 5)*f +
			c[7] * tx*ycs + c[8] * pow(tx, 2)*ycs + c[9] * pow(tx, 3)*ycs +
			c[10] * pow(tx, 4)*ycs + c[11] * pow(tx, 5)*ycs + c[12] * ycs +
			c[13] * f*ycs + c[14] * tx*f*ycs + c[15] * pow(tx, 2)*f*ycs;
		h = 300.0 / 100.0;
		sum = b[0] + b[1] * f + b[2] * tx*f + b[3] * pow(tx, 2)*f +
			b[4] * pow(tx, 3)*f + b[5] * pow(tx, 4)*f + b[6] * pow(tx, 5)*f +
			b[7] * tx*ycs + b[8] * pow(tx, 2)*ycs + b[9] * pow(tx, 3)*ycs +
			b[10] * pow(tx, 4)*ycs + b[11] * pow(tx, 5)*ycs + b[12] * h*ycs +
			b[13] * tx*h*ycs + b[14] * pow(tx, 2)*h*ycs + b[15] * pow(tx, 3)*h*ycs +
			b[16] * pow(tx, 4)*h*ycs + b[17] * pow(tx, 5)*h*ycs + b[18] * ycs;
		dtc300 = sum;
		sum = b[12] * ycs + b[13] * tx*ycs + b[14] * pow(tx, 2)*ycs + b[15] * pow(tx, 3)*ycs +
			b[16] * pow(tx, 4)*ycs + b[17] * pow(tx, 5)*ycs;
		dtc300dz = sum;
		cc = 3.0*dtc300 - dtc300dz - 3.0*aa - 2.0*bb;
		dd = dtc300 - aa - bb - cc;
		zp = (zht - 240.0) / 60.0;
		*dtc = aa + bb*zp + cc*zp*zp + dd*zp*zp*zp;
	}

	if ((zht > 300) & (zht <= 600.0)){
		h = zht / 100.0;
		sum = b[0] + b[1] * f + b[2] * tx*f + b[3] * pow(tx, 2)*f +
			b[4] * pow(tx, 3)*f + b[5] * pow(tx, 4)*f + b[6] * pow(tx, 5)*f +
			b[7] * tx*ycs + b[8] * pow(tx, 2)*ycs + b[9] * pow(tx, 3)*ycs +
			b[10] * pow(tx, 4)*ycs + b[11] * pow(tx, 5)*ycs + b[12] * h*ycs +
			b[13] * tx*h*ycs + b[14] * pow(tx, 2)*h*ycs + b[15] * pow(tx, 3)*h*ycs +
			b[16] * pow(tx, 4)*h*ycs + b[17] * pow(tx, 5)*h*ycs + b[18] * ycs;
		*dtc = sum;
	}

	if ((zht > 600.0) & (zht <= 800.0)){
		zp = (zht - 600.0) / 100.0;
		hp = 600.0 / 100.0;
		aa = b[0] + b[1] * f + b[2] * tx*f + b[3] * pow(tx, 2)*f +
			b[4] * pow(tx, 3)*f + b[5] * pow(tx, 4)*f + b[6] * pow(tx, 5)*f +
			b[7] * tx*ycs + b[8] * pow(tx, 2)*ycs + b[9] * pow(tx, 3)*ycs +
			b[10] * pow(tx, 4)*ycs + b[11] * pow(tx, 5)*ycs + b[12] * hp*ycs +
			b[13] * tx*hp*ycs + b[14] * pow(tx, 2)*hp*ycs + b[15] * pow(tx, 3)*hp*ycs +
			b[16] * pow(tx, 4)*hp*ycs + b[17] * pow(tx, 5)*hp*ycs + b[18] * ycs;
		bb = b[12] * ycs + b[13] * tx*ycs + b[14] * pow(tx, 2)*ycs + b[15] * pow(tx, 3)*ycs +
			b[16] * pow(tx, 4)*ycs + b[17] * pow(tx, 5)*ycs;
		cc = -(3.0*aa + 4.0*bb) / 4.0;
		dd = (aa + bb) / 4.0;
		*dtc = aa + bb*zp + cc*zp*zp + dd*zp*zp*zp;
	}

	return;
}

void JB2008::tmoutd(double d1950, int *iyr, double *day)
{

	//tmoutd member function from JB2008 class
	//Compute day and year from the d1950 (days since 1950)

	int iyday, itemp;
	double frac;

	iyday = d1950;
	frac = d1950 - iyday;
	iyday = iyday + 364;
	itemp = iyday / 1461;
	iyday = iyday - itemp * 1461;
	*iyr = 1949 + 4 * itemp;
	itemp = iyday / 365;
	if (itemp >= 3) itemp = 3;
	*iyr = *iyr + itemp;
	iyday = iyday - 365 * itemp + 1;
	*iyr = *iyr - 1900;
	*day = iyday + frac;
	if (*iyr >= 100)*iyr = *iyr - 100;

	return;
}

double JB2008::xlocal(double z, double tc[])
{

	//xlocal member function from JB2008 class
	//Evaluates equation (10) or equation (13), depending on z

	double dz, xlocal_out;

	dz = z - 125.0;
	if (dz > 0.0) {
		xlocal_out = tc[0] + tc[2] * atan(tc[3] * dz*(1.0 + 4.5e-06*pow(dz, 2.5)));
	}

	else {

		xlocal_out = ((-9.8204695e-06*dz - 7.3039742e-04)*pow(dz, 2) + 1)*dz*tc[1] + tc[0];
	}


	return xlocal_out;

}

double JB2008::xgrav(double z)
{
	//xgrav member function from JB2008 class
	//Evaluates equation (8)

	double xgrav_out;

	xgrav_out = 9.80665 / pow((1.0 + z / 6356.766), 2);

	return xgrav_out;
}


void JB2008::JB08(double amjd, double sun[], double sat[], double f10,
	double f10b, double s10, double s10b, double xm10, double xm10b,
	double y10, double y10b, double dstdtc, double temp[], double *rho,
	double *pres, double *avgmw, double and1[], double *sumn)
{
	//JB08 member function from JB2008 class
	//Jacchia-Bowman 2008 Model Atmosphere
	//This is the CIRA "Integration Form" of a Jacchia Model.  There are no tabular
	//values of density.  Instead, the barometric equation and diffusion equation are 
	//integrated numerically using the Newton-Coates method to produce the density
	//profile up to the input position.

	//The alpha are the thermal diffusion coefficients in Eq. (6)
	const double alpha[5] = { 0.0, 0.0, 0.0, 0.0, -0.38 };

	//al10 is log(10.0)
	const double al10 = 2.3025851;

	//The AMW are the molecular weights in order:  N2, O2, O, Ar, He, and H
	const double amw[6] = { 28.0134, 31.9988, 15.9994, 39.948, 4.0026, 1.00797 };

	//avogad is Avogadro's number in mks units (molecules/kmol)
	const double avogad = 6.02257e+26;

	const double twopi = 6.2831853, pi = 3.1415927, piov2 = 1.5707963;

	//the frac are assumed sea-level volume fractions in order:  N2, O2, Ar, and He
	const double frac[4] = { 0.7811, 0.20955, 9.34e-03, 1.289e-05 };

	//rstar is the universal gas-constant in mks units (joules/K/kmol)
	const double rstar = 8314.32;

	//The R# are values used to establish height step sizes in the regimes 90km
	//to 105km, 105km to 500km and 500km upward.
	const double r1 = 0.01, r2 = 0.025, r3 = 0.075;

	//The wt are weights for the Newton-Coates five-point quad. formula
	const double wt[5] = { 0.311111111111111, 1.422222222222222, 0.533333333333333,
		1.422222222222222, 0.311111111111111 };

	//The cht are coefficients for high altitude density correction
	const double cht[4] = { 0.22, -0.20e-02, 0.115e-02, -0.211e-05 };

	const double degrad = pi / 180.0;

	double fn, fsb, tsubc, eta, theta, h, tau, glat, zht, glst, glsthr, c, s, df, tsubl,
		dtclst, tinf, tsubx, gsubx, tc[4], z1, z2, al, oor, zr, ambar1, zend, sum2,
		ain, tloc1, z, dz, ambar2, tloc2, gravl, sum1, fact1, anm, an, fact2, aln[6], z3,
		tloc3, z4, r, sum3, tloc4, altr, hsign, alnh5, trash, capphi, dlrsl, al10t5,
		dlrsa, d1950, yrday, fzz, gtz, dlr, sumnm, al10n[6], fex, zeta, zeta2, zeta3,
		f15c, f15c_zeta, fex2, fex3;

	int n, i, j, iyr;

	//Equation (14)
	fn = pow((f10b / 240.0), (1.0 / 4.0));
	if (fn > 1.0) fn = 1.0;
	fsb = f10b*fn + s10b*(1.0 - fn);
	tsubc = 392.4 + 3.227*fsb + 0.298*(f10 - f10b) + 2.259*(s10 - s10b) +
		0.312*(xm10 - xm10b) + 0.178*(y10 - y10b);

	//Equation (15)
	eta = 0.5*abs(sat[1] - sun[1]);
	theta = 0.5*abs(sat[1] + sun[1]);

	//Equation (16)
	h = sat[0] - sun[0];
	tau = h - 0.64577182 + 0.10471976 * sin(h + 0.75049158);
	glat = sat[1];
	zht = sat[2];
	glst = h + pi;
	glsthr = (glst / degrad)*(24.0 / 360.0);
	if (glsthr >= 24.0) glsthr = glsthr - 24.0;
	if (glsthr < 0.0) glsthr = glsthr + 24.0;

	//Equation (17)
	c = pow(cos(eta), 2.5);
	s = pow(sin(theta), 2.5);
	df = s + (c - s)*pow(abs(cos(0.5*tau)), 3);
	tsubl = tsubc*(1.0 + 0.31 * df);

	//Compute correction to dTc for local solar time and lat correction
	dtsub(f10, glsthr, glat, zht, &dtclst);

	//Compute the local exospheric temperature.
	//Add geomagnetic storm effect from input dTc value

	temp[0] = tsubl + dstdtc;
	tinf = tsubl + dstdtc + dtclst;

	//Equation (9)
	tsubx = 444.3807 + 0.02385 * tinf - 392.8292 * exp(-0.0021357*tinf);

	//Equation (11)
	gsubx = 0.054285714*(tsubx - 183.0);

	//The tc array will be an argument in the call to xlocal, which evaluates 
	//Equation (10) and Equation (13)
	tc[0] = tsubx;
	tc[1] = gsubx;

	//a and gsubx/a of equation (13)
	tc[2] = (tinf - tsubx) / piov2;
	tc[3] = gsubx / tc[2];

	//Equation (5)
	z1 = 90.0;

	z2 = min(sat[2], 105.0);

	al = log(z2 / z1);

	oor = 1.0 / r1;

	n = int(al*oor) + 1;

	zr = exp(al / double(n));
	ambar1 = xambar(z1);
	tloc1 = xlocal(z1, tc);
	zend = z1;
	sum2 = 0.0;
	ain = ambar1 * xgrav(z1) / tloc1;

	for (i = 0; i < n; i++){
		z = zend;

		zend = zr * z;
		dz = 0.25*(zend - z);
		sum1 = wt[0] * ain;
		for (j = 1; j < 5; j++){
			z = z + dz;

			ambar2 = xambar(z);

			tloc2 = xlocal(z, tc);

			gravl = xgrav(z);
			ain = ambar2*gravl / tloc2;
			sum1 = sum1 + wt[j] * ain;
		}
		sum2 = sum2 + dz*sum1;
	}

	fact1 = 1000.0 / rstar;
	*rho = 3.46e-6*ambar2*tloc1*exp(-fact1*sum2) / ambar1 / tloc2;

	//Equation (2)
	anm = avogad * *rho;
	an = anm / ambar2;

	//Equation (3)
	fact2 = anm / 28.96;
	aln[0] = log(frac[0] * fact2);

	aln[3] = log(frac[2] * fact2);
	aln[4] = log(frac[3] * fact2);

	//Equation (4)
	aln[1] = log(fact2 * (1.0 + frac[1]) - an);
	aln[2] = log(2.0 * (an - fact2));

	if (sat[2] > 105.0) goto three;
	temp[1] = tloc2;

	//Put in negligible hydrogen for use in loop 13
	aln[5] = aln[4] - 25.0;
	goto eleven;

	//Equation (6)
three:
	z3 = min(sat[2], 500.0);
	al = log(z3 / z);
	oor = 1.0 / r2;
	n = int(al*oor) + 1;
	zr = exp(al / double(n));
	sum2 = 0.0;

	ain = gravl / tloc2;

	for (i = 0; i < n; i++){
		z = zend;
		zend = zr * z;
		dz = 0.25*(zend - z);
		sum1 = wt[0] * ain;
		for (j = 1; j < 5; j++){
			z = z + dz;
			tloc3 = xlocal(z, tc);
			gravl = xgrav(z);
			ain = gravl / tloc3;
			sum1 = sum1 + wt[j] * ain;
		}
		sum2 = sum2 + dz*sum1;
	}

	z4 = max(sat[2], 500.0);
	al = log(z4 / z);
	r = r2;
	if (sat[2] > 500.0) r = r3;
	oor = 1.0 / r;
	n = int(al*oor) + 1;
	zr = exp(al / double(n));
	sum3 = 0.0;

	for (i = 0; i < n; i++){
		z = zend;
		zend = zr*z;
		dz = 0.25*(zend - z);
		sum1 = wt[0] * ain;
		for (j = 1; j < 5; j++){
			z = z + dz;
			tloc4 = xlocal(z, tc);
			gravl = xgrav(z);
			ain = gravl / tloc4;
			sum1 = sum1 + wt[j] * ain;
		}
		sum3 = sum3 + dz * sum1;
	}

	if (sat[2] > 500.0) goto eight;

	temp[1] = tloc3;
	altr = log(tloc3 / tloc2);
	fact2 = fact1 * sum2;
	hsign = 1.0;
	goto nine;
eight:
	temp[1] = tloc4;
	altr = log(tloc4 / tloc2);
	fact2 = fact1*(sum2 + sum3);

	hsign = -1.0;

nine:
	for (i = 0; i < 5; i++){
		aln[i] = aln[i] - (1.0 + alpha[i]) * altr - fact2*amw[i];
	}



	//Equation (7) - Note that in CIRA72, al10t5 = log10(t500)
	al10t5 = log10(tinf);
	alnh5 = (5.5*al10t5 - 39.40) * al10t5 + 73.13;
	aln[5] = al10*(alnh5 + 6.0) + hsign * (log(tloc4 / tloc3) + fact1 * sum3 * amw[5]);



eleven:

	//Equation (24) - J70 Seasonal-Latitudinal Variation
	trash = (amjd - 36204.0) / 365.2422;
	capphi = fmod(trash, 1.0);

	dlrsl = 0.02*(sat[2] - 90.0)*exp(-0.045*(sat[2] - 90.0))*copysign(1.0, sat[1])*
		sin(twopi*capphi + 1.72)*pow(sin(sat[1]), 2);

	//Equation (23) - Computes the semiannual variational
	dlrsa = 0.0;
	if (z < 2000.0){
		d1950 = amjd - 33281.0;
		tmoutd(d1950, &iyr, &yrday);
		//Use new semiannual model
		semian08(iyr, yrday, zht, f10b, s10b, xm10b, &fzz, &gtz, &dlrsa);
		if (fzz < 0.0)dlrsa = 0.0;
	}


	//Sum the delta-log-rhos and apply to the number densities.
	dlr = al10 * (dlrsl + dlrsa);

	for (i = 0; i < 6; i++){
		aln[i] = aln[i] + dlr;
	}


	//Compute mass-density and mean-molecular-weight and convert number density
	//logs from natural to common.
	*sumn = 0.0;
	sumnm = 0.0;

	for (i = 0; i < 6; i++){
		an = exp(aln[i]);
		and1[i] = an;
		*sumn = *sumn + an;
		sumnm = sumnm + an*amw[i];
		al10n[i] = aln[i] / al10;
	}

	*rho = sumnm / avogad;
	*avgmw = sumnm / *sumn;


	//Compute the high altitude exospheric density correction factor
	fex = 1.0;
	if ((zht >= 1000.0) & (zht < 1500.0)){
		zeta = (zht - 1000.0)*0.002;
		zeta2 = zeta*zeta;
		zeta3 = zeta*zeta2;
		f15c = cht[0] + cht[1] * f10b + cht[2] * 1500.0 + cht[3] * f10b * 1500;
		f15c_zeta = (cht[2] + cht[3] * f10b) *500.0;
		fex2 = 3.0*f15c - 2.0*f15c + 2.0;
		fex3 = f15c_zeta - 2.0*f15c + 2.0;
		fex = 1.0 + fex2*zeta2 + fex3*zeta3;
	}

	if (zht >= 1500.0){
		fex = cht[0] + cht[1] * f10b + cht[2] * zht + cht[3] * f10b*zht;
	}

	//Apply the exospheric density correction factor.
	*rho = fex* *rho;

	*pres = *rho*rstar*temp[1] / *avgmw;
	for (i = 0; i < 6; i++){
		and1[i] = and1[i] * fex;
	}

	*sumn = *sumn*fex;

	return;

}

void JB2008::tme(int mn, int ida, int iyr, int ihr, double xmin, double xlng, double *xlat, double *sda, double *sha, double *dd, double *dy, double *sra, double *raloc, double *xmjd)
{

	//tme member function from the JB2008 class
	//'tme' performs the calculations of the solar declination and solar
	//hour angle.


	int i, id;
	double degs, gmt, dp2k, sunlon, g, eclon, oblq, alpha, temp, year, xhr,
		eqt;
	double pi, pi180;
	int iday[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };


	pi = 3.1415926535897931;
	pi180 = pi / 180;
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

	//Compute days past J2000

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

	return;
}



void JB2008::caltojul(int iy, int im, int id, int ihour, int imin, double sec, double *xjd)
{

	//caltojul member function from JB2008 class
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

}

void JB2008::JB08HWM(double xmin, double z, double xlat, double xlon,
	double *ph, double *dh, double *th, double *n2nd, double *o2nd,
	double *ond, double *arnd, double *hend, double *hnd, double *nnd,
	double *wtmol, double *tex, double *uh, double *vh)
{
	//JB08HWM member function from the JB2008 class


	int iday;
	double rlat, sda, utsec, sha, xlst, dy, w[2], aph[2], dd, sra, raloc, xmjd, sun[2],
		sat[3], dlst, temp[2], rho, pres, avgmw, and1[6], sumn, dut, pi = 3.1415926535897931, pi180;

	rlat = xlat;
	pi180 = pi / 180.0;
	tme(mn, ida, iyr, ihr, xmin, xlon, &rlat, &sda, &sha, &dd, &dy, &sra, &raloc, &xmjd);


	//JB2008 thermosphere subroutine
	sun[0] = sra;
	sun[1] = sda;
	sat[0] = raloc;
	sat[1] = xlat*pi180;
	sat[2] = z;

	//Compute modified Julian date from mean Julian date plus time of day
	xmjd = xmjd - 2400000.5 + (60.0*ihr + xmin) / 1440.0;

	JB08(xmjd, sun, sat, f10, f10b, s10, s10b, xm10, xm10b, y10, y10b, dstdtc, temp, &rho, &pres, &avgmw, and1, &sumn);

	//Change of notation for outputs
	*dh = rho;
	*th = temp[1];
	*tex = temp[0];
	*wtmol = avgmw;
	*ph = pres;
	*n2nd = and1[0];
	*o2nd = and1[1];
	*ond = and1[2];
	*arnd = and1[3];
	*hend = and1[4];
	*hnd = and1[5];
	*nnd = 0.0;


	//Compute inputs needed for HWM
	iday = 1000 * iyr + int(dd + 1.0);
	dut = 86400.0*(dd - int(dd));
	utsec = dut;
	dlst = 12.0*(1.0 + sha / pi);
	xlst = fmod(dlst, 24.0);

	//Call harmonic wind model (HWM) subroutine
	aph[0] = ap;
	aph[1] = 0.0;

	hwm.gws5(iday, utsec, z, xlat, xlon, xlst, f10b, f10, aph, w);

	//change notation for outputs
	*uh = w[1];
	*vh = w[0];

	return;
}

void JB2008::JB08mod(double h, double phi, double thet, double elt, double ri, double g, double *pmj, double *dmj, double *tmj, double *umj, double *vmj, double *wmj, double *n2nd,
	double *o2nd, double *ond, double *arnd, double *hend, double *hnd, double *wtmol,
	double *dmdz, double *nnd)
{
	//JB08mod member function from JB2008 class
	//JB2008 and HWM model driver routine to evaluate mean p, d, t, u, v, w

	double dy5, dx5, xmin, dphij, dtx, dty, cp, the, thn, tb, d8, d10, d9, d11,
		tex, dhn, d1, d2, d4, d6, dhe, pb, db, d3, d5, d7, phe, phn, pi, pi180;

	imin = mino + int(elt) / 60;
	sec = seco + elt;
	sec = int(sec) % 60;
	ihr = ihro + imin / 60;
	imin = imin % 60;

	pi = 3.1415926535897931;
	pi180 = pi / 180.0;

	//Distances for 5 degrees of geocentric latitude, longitude
	dy5 = 5000.0*ri*pi180;
	dx5 = dy5*cos(phi*pi180);
	dx5 = max(2000.0, dx5);

	//Following is the pure thermospheric height range section
	//JB2008 and HWM values at current position
	xmin = imin + sec / 60.0;

	JB08HWM(xmin, h, phi, thet, pmj, dmj, tmj, n2nd, o2nd, ond, arnd, hend, hnd, nnd, wtmol, &tex, umj, vmj);

	//Set geocentric latitude increment for temperature gradients
	dphij = 5.0;
	if (phi > 85.0) dphij = -5.0;

	//JB2008 temperature at current position plus geocentric latitude increment
	JB08HWM(xmin, h, phi + dphij, thet, &phn, &dhn, &thn, &d1, &d2, &d3, &d4, &d5, &d6, &d7, &d8, &d9, &d10, &d11);

	//JB2008 temperature at current position plus 5 degrees lon
	JB08HWM(xmin, h, phi, thet + 5.0, &phe, &dhe, &the, &d1, &d2, &d3, &d4, &d5, &d6, &d7, &d8, &d9, &d10, &d11);

	//dt/dx, dt/dy and dt/dz for vertical wind
	dtx = the - *tmj;
	dty = thn - *tmj;
	if (dphij < 0.0) dty = -dty;


	//JB2008 temperature and molecular weight at current lat-lon, 1 km higher
	JB08HWM(xmin, h + 1.0, phi, thet, &pb, &db, &tb, &d1, &d2, &d3, &d4, &d5, &d6, &d7, &d8, &d9, &d10, &d11);

	//Gradients for temperature and molecular weight
	dtz = (tb - *tmj) / 1000.0;
	*dmdz = (d8 - *wtmol) / 1000.0;

	//Compute vertical mean wind
	//Specific heat
	cp = 7.0**pmj / (2.0**dmj**tmj);

	//Mean vertical wind from Montgomery stream function
	*wmj = -cp*(*umj*dtx / dx5 + *vmj*dty / dy5) / (g + cp*dtz);

	return;

}


void JB2008::namelist(string namef)
{

	//namelist member function from JB2008 class
	//Open and read namelist input file for setting model parameters.

	ifstream namelist;
	string dummy, NCEPpath, NCEPmn, atmpath, trapath, prtpath, nprpath,
		conpath, rndpath, rrapath, rralist, profile, NCEPpath1;
	double h1, phi1, thet1, dphi, dthet, dhgt, rpscale, ruscale, rwscale,
		rdinit, rtinit, ruinit, rvinit, rwinit, sitelim, sitenear, patchy,
		z0in;
	int iopt, NCEPyr, NCEPhr, nr1, iurra, iyrrra, initpert, itherm, ibltest,
		nmax, iaux, mc;


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


	if (itherm == 3){
		if (s10 <= 0.0 || s10b == 0.0){
			s10 = f10;
			s10b = f10b;
		}

		if (xm10 <= 0.0 || xm10b == 0.0){
			xm10 = f10;
			xm10b = f10b;
		}

		if (y10 <= 0.0 || y10b == 0.0){
			y10 = f10;
			y10b = f10b;
		}

		if (dstdtc <= 0.0){
			dstdtc = 0.0;
		}
	}

	imin = mino;
	ihr = ihro;
	sec = seco;
}

