//P. White
//Auxilliary Profile Class


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <string>
#include <cmath>
#include "AuxProf.h"
using namespace std;


AuxProf::AuxProf()
{


	initializeMemberVariables();

}


void AuxProf::initializeMemberVariables()
{

	nprof = 0;
	profnear = 0.0, proffar = 0.0;
	pi = 4.0*atan(1.0);
	pi180 = pi / 180.0;

}



void AuxProf::rdprof(double sitenear, double sitelim, string home1, string profile1)
{
	//rdprof member function from AuxProf class
	//Reads alternate profile data file.  The rdprof member function reads 
	//height (zhgt), latitude (xlat), longitude (xlon), temperature (t),
	//pressure (p), density (d), east-west wind (u), north-south wind (v),
	//standard deviation of temperature (st), pressure (sp), density (sd),
	//east-west wind (su), and north-south wind (sv).

	ifstream altprof;

	int n;
	double zhgt, xlat, xlon, t, p, d, u, v, st, sp, sd, su, sv,
		g0, g, re, r0, ri, a, b, fac;
	string dummy[13];
	string altp = home1 + profile1;
	
	fac = atan(1.0) / 45;

	profnear = sitenear;
	proffar = sitelim;

	//Open profile data file
	altprof.open(altp.c_str());

	//File open error check
	if (altprof.is_open()){
	}
	else {
		cout << "File Open Error!  " << home1 + profile1 << '\n';
		system("pause");
		exit(1);
	}

	//Read and ingnore header line
	for (int i = 0; i < 13; i++){
		altprof >> dummy[i];
	}
	n = 0;

	//Read auxiliary profile while not eof
	while (!altprof.eof()){

		altprof >> zhgt >> xlat >> xlon >> t >> p >> d >> u >> v >>
			st >> sp >> sd >> su >> sv;
		
		//Interpret input zhgt as radius if zhgt > 6000 km
		if (zhgt > 6000.0){
			//Get local earth radius r0
			rig(0.0, xlat*fac, &g0, &g, &re, &r0, &ri, &a, &b);
			zhgt = zhgt - r0;
		}
		
		//Convert negative longitudes
		if (xlon < 0.0){
			xlon += 360.0;
		}

		//Store profile height (zght) data in vector (pght)
		phgt.push_back(zhgt);

		//Stop if two successive heights are the same
		if ((n > 0) && (phgt[n] == phgt[n - 1])){
			//	cout << "Successive heights the same." << '\n';
			goto L40;
		}

		//If temperature <= 0, set standard deviaiton to zero
		if (t <= 0.0){
			st = 0.0;
		}
		//Check for large temperature standard deviation
		else if (st > 0.3*t){
			st = 0.3;
		}
		else {
			st = st / t;
		}
		//If pressure <= 0, set standard deviation to zero
		if (p <= 0.0){
			sp = 0.0;
		}
		//Check for large pressure standard deviation
		else if (sp > 0.3*p){
			sp = 0.3;
		}
		else {
			sp = sp / p;
		}
		//If density <= 0, set standard deviation to zero
		if (d <= 0.0){
			sd = 0.0;
		}
		//Check for large density standard deviation 
		else if (sd > 0.3*d){
			sd = 0.3;
		}
		else {
			sd = sd / d;
		}

		//Store profile data in vectors
		plat.push_back(xlat);
		plon.push_back(xlon);
		ptmp.push_back(t);
		pprs.push_back(p);
		pden.push_back(d);
		puwn.push_back(u);
		pvwn.push_back(v);
		pstmp.push_back(st);
		psprs.push_back(sp);
		psden.push_back(sd);
		psuwn.push_back(su);
		psvwn.push_back(sv);

		//If previous standard > 0.0 and current standard deviation <= 0.0,
		//set current standard deviation = previous standard deviation.
		if (n > 0){
			if (pstmp[n] <= 0.0 && pstmp[n - 1] > 0.0) pstmp[n] = pstmp[n - 1];
			if (psprs[n] <= 0.0 && psprs[n - 1] > 0.0) psprs[n] = psprs[n - 1];
			if (psden[n] <= 0.0 && psden[n - 1] > 0.0) psden[n] = psden[n - 1];
			if (psuwn[n] <= 0.0 && psuwn[n - 1] > 0.0) psuwn[n] = psuwn[n - 1];
			if (psvwn[n] <= 0.0 && psvwn[n - 1] > 0.0) psvwn[n] = psvwn[n - 1];
		}

		n++;
		nprof = n;
	}

L40:

	//Close profile input file when end-of-file encountered
	altprof.close();
	
	return;
}


void AuxProf::profterp(double chgt, double clat, double clon, double tin,
	double pin, double din, double uin, double vin, double *ptemp,
	double *ppres, double *pdens, double *puwin, double *pvwin,
	double *profwgt)
{
	//profterp member function from AuxProf class
	//Interpolates profile data to current position (chgt, clat, clon)
	//and weights results (with factor profwgt) with input values
	//(tin,pin,din,uin,vin), yielding weighted (ptemp,ppres,pdens,puwin,
	//pvwin).  Input profnear is lat-lon radius over which profile is 
	//weighted with 1.0; proffar is lat-lon radius beyond which profile
	//is given zero weight.

	int ia[2] = { 0, 0 };
	double adll[2] = { 0.0, 0.0 };

	int i1, i2, ni, i;
	double piby2, radius1, radius2, factor, pilat, dplon, pilon, facthgt,
		factll, radius, tpdwgt, uvwgt;


	//Return input values if profnear = 0.0
	if (profnear <= 0.0){
		*ptemp = tin;
		*ppres = pin;
		*pdens = din;
		*puwin = uin;
		*pvwin = vin;

		return;
	}

	//Calculate pi/2
	piby2 = 1.5707963267948966;

	i1 = 0;
	i2 = 0;
	ni = 0;

	*profwgt = 0.0;

	//Find nearest pair of points above and below current height
	for (i = 0; i < nprof - 1; i++){

		if (((chgt - phgt[i])*(chgt - phgt[i + 1]) < 0.0) ||
			(chgt == phgt[i])){
			ni++;
			if (ni > 2){
				cout << "Too many height pairs in profile." << '\n';
				system("pause");
				exit(1);
			}

			ia[ni - 1] = i + 1;

			//Lat-lon radius from position of points i1 and i2
			radius1 = radll(clat, clon, plat[i], plon[i]);
			radius2 = radll(clat, clon, plat[i + 1], plon[i + 1]);
			adll[ni - 1] = (radius1 + radius2) / 2.0;
		}
	}

	if (ni > 0){

		i1 = ia[0];

		if ((ni == 2) && (adll[1] < adll[0])){
			i1 = ia[1];
		}

		i2 = i1 + 1;
	}

	if (i1 == 0) {
		*pdens = 0.0;
		*puwin = 0.0;
		*pvwin = 0.0;
		goto L50;
	}

	//Compute factor for linear height interpolation
	factor = (chgt - phgt[i1 - 1]) / (phgt[i2 - 1] - phgt[i1 - 1]);

	//Linear height intepolation for lat, lon, temperature, winds
	pilat = plat[i1 - 1] + factor*(plat[i2 - 1] - plat[i1 - 1]);
	dplon = plon[i2-1] - plon[i1-1];

	if (dplon > 180.0){
		dplon -= 360.0;
	}

	if (dplon < -180.0){
		dplon += 360.0;
	}

	pilon = plon[i1] + factor*dplon;
	*ptemp = ptmp[i1 - 1] + factor*(ptmp[i2 - 1] - ptmp[i1 - 1]);
	*puwin = puwn[i1 - 1] + factor*(puwn[i2 - 1] - puwn[i1 - 1]);
	*pvwin = pvwn[i1 - 1] + factor*(pvwn[i2 - 1] - pvwn[i1 - 1]);

	//Power-law interpolation for pressure (unless profile density
	//is zero, for which zero weight will be used)
	*ppres = 0.0;

	if (pprs[i1] > 0.0){
		*ppres = pprs[i1 - 1] * pow(pprs[i2 - 1] / pprs[i1 - 1], factor);
	}

	// Power-law interpolation for density (unless profile density 
	//is zero, for which zero weight will be used)
	*pdens = 0.0;

	if (pden[i1] > 0.0) {
		*pdens = pden[i1 - 1] * pow(pden[i2 - 1] / pden[i1 - 1], factor);
	}


	//Initialize weighting factor components for height and lat-lon
	facthgt = 1.0;
	factll = 1.0;

	if (i1 == 1){

		//Sine-squared variation of height weighting from 0 at 1st point
		//to 1 at 2nd point
		facthgt = abs((chgt - phgt[0]) / (phgt[1] - phgt[0]));
		facthgt = pow(sin(piby2*facthgt), 2.0);

	}
	else if (i2 == nprof){

		//Sine-squared variation of height weighting from 0 at 
		//next-to-last point to 1 at last point
		facthgt = (chgt - phgt[nprof]) / (phgt[nprof - 2] - phgt[nprof]);
		facthgt = pow(sin(piby2*facthgt), 2.0);
	}

	//Lat-lon radius of current position from profile lat-lon
	radius = radll(clat, clon, pilat, pilon);

	//Use weight = 0 if radius > proffar, weight = 1 if radius < profnear,
	//with sine squared variation between proffar and profnear
	if (radius >= proffar){
		factll = 0.0;
	}
	else if (radius <= profnear){
		factll = 1.0;
	}
	else {
		factll = (proffar - radius) / (proffar - profnear);
		factll = pow(sin(piby2*factll), 2.0);
	}

	//Total weight = product of weights for lat-lon and height
	*profwgt = factll*facthgt;

L50:
	tpdwgt = *profwgt;
	uvwgt = *profwgt;

	//Set profile weight to zero for p, d, & t if profile values are 0
	if (*ptemp*(*ppres)*(*pdens) == 0.0){
		tpdwgt = 0.0;
	}

	//Set profile weight to zero u & v if profile values are 0
	if (fabs(*puwin) + fabs(*pvwin) == 0.0){
		uvwgt = 0.0;
	}

	//Apply weighted averaging of profile values with input values
	*ptemp = tpdwgt*(*ptemp) + (1.0 - tpdwgt)*tin;
	*ppres = tpdwgt*(*ppres) + (1.0 - tpdwgt)*pin;
	*pdens = tpdwgt*(*pdens) + (1.0 - tpdwgt)*din;
	*puwin = uvwgt*(*puwin) + (1.0 - uvwgt)*uin;
	*pvwin = uvwgt*(*pvwin) + (1.0 - uvwgt)*vin;

	return;

}



void AuxProf::profsigs(double chgt, double clat, double clon, double tin,
	double pin, double din, double uin, double vin, double *ptemp, double *ppres,
	double *pdens, double *puwin, double *pvwin, double *profwgt)
{

	//profsigs member function from AuxProf class
	//Interpolates profile data to current position (chgt,clat,clon)
	//and weights results (with factor profwgt) with input values
	//(tin,pin,din,uin,vin), yielding weighted average (ptemp,ppres,
	//pdens,puwin,pvwin).  Input profnear is lat-lon radius over which
	//profile is weighted with 1.0; proffar is lat-lon radius beyond
	//which profile is given zero weight.

	const int mpmax = 100000;
	int ia[2] = { 0, 0 };
	double adll[2] = { 0.0, 0.0 };

	int i1, i2, ni, i;
	double piby2, radius1, radius2, factor, pilat, dplon, pilon, facthgt,
		factll, radius, tpdwgt, uvwgt;

	piby2 = 1.5707963267948966;


	//Return input values if profnear = 0
	if (profnear <= 0.0){
		*ptemp = tin;
		*ppres = pin;
		*pdens = din;
		*puwin = uin;
		*pvwin = vin;
		*profwgt = 0.0;
		return;
	}

	i1 = 0;
	i2 = 0;
	ni = 0;
	*profwgt = 0.0;

	//Find nearest pair of points above and below current height
	for (i = 0; i < nprof - 1; i++){

		if (((chgt - phgt[i])*(chgt - phgt[i + 1]) < 0.0)
			|| (chgt == phgt[i])) {

			++ni;

			if (ni > 2) {
				cout << "Too many height pairs in profile." << '\n';
				system("pause");
				exit(1);
			}

			ia[ni - 1] = i + 1;

			//Lat-lon radius from position of points i1 and i2
			radius1 = radll(clat, clon, plat[i], plon[i]);
			radius2 = radll(clat, clon, plat[i + 1], plon[i + 1]);

			adll[ni - 1] = (radius1 + radius2) / 2.0;
		}
	}

	if (ni > 0){

		i1 = ia[0];

		if ((ni == 2) && (adll[1] < adll[0])){
			i1 = ia[1];
		}

		i2 = i1 + 1;
	}

	if (i1 == 0){
		*pdens = 0.0;
		*puwin = 0.0;
		*pvwin = 0.0;
		goto L50;
	}

	//Compute factor for linear height interpolation
	factor = (chgt - phgt[i1 - 1]) / (phgt[i2 - 1] - phgt[i1 - 1]);

	//Linear height interpolation for lat, lon, temperature, winds
	pilat = plat[i1-1] + factor*(plat[i2-1] - plat[i1-1]);
	dplon = plon[i2-1] - plon[i1-1];

	if (dplon > 180.0){
		dplon -= 360.0;
	}

	if (dplon < -180.0){
		dplon += 360.0;
	}

	pilon = plon[i1] + factor*dplon;
	*ptemp = pstmp[i1 - 1] + factor*(pstmp[i2 - 1] - pstmp[i1 - 1]);
	*puwin = psuwn[i1 - 1] + factor*(psuwn[i2 - 1] - psuwn[i1 - 1]);
	*pvwin = psvwn[i1 - 1] + factor*(psvwn[i2 - 1] - psvwn[i1 - 1]);
	*ppres = psprs[i1 - 1] + factor*(psprs[i2 - 1] - psprs[i1 - 1]);
	*pdens = psden[i1 - 1] + factor*(psden[i2 - 1] - psden[i1 - 1]);

	//Initialize weighting factor components for height and lat-lon
	facthgt = 1.0;
	factll = 1.0;

	if (i1 == 1){
		//Sine-squared variation of height weighting from 0 at 1st 
		//point to 1 at 2nd point
		facthgt = (chgt - phgt[0]) / (phgt[1] - phgt[0]);
		facthgt = pow(sin(piby2*facthgt), 2.0);
	}
	else if (i2 == nprof){
		//Sine-squared variation of height weighting from 0 to 
		//netx-to-last point to 1 at last point
		facthgt = (chgt - phgt[nprof]) / (phgt[nprof - 2] - phgt[nprof]);
		facthgt = pow(sin(piby2*facthgt), 2.0);

	}

	//Lat-lon radius of current position from profile lat-lon
	radius = radll(clat, clon, pilat, pilon);

	//Use weight=0 if radius>proffar, weight=1 if radius<profnear,
	//with sine-squared variation between proffar and profnear
	if (radius >= proffar){
		factll = 0.0;
	}
	else if (radius <= profnear){
		factll = 1.0;
	}
	else {
		factll = (proffar - radius) / (proffar - profnear);
		factll = pow(sin(piby2*factll), 2.0);
	}

	//Total weight = product of weights for lat-lon and height
	*profwgt = factll*facthgt;
L50:
	tpdwgt = *profwgt;
	uvwgt = *profwgt;

	//Set profile weight to zero for p, d, & t if profile values are 0
	if (*ptemp*(*ppres)*(*pdens) == 0.0){
		tpdwgt = 0.0;
	}

	//Set profile weight to zero for u & v if profile values are 0
	if (fabs(*puwin) + fabs(*pvwin) == 0.0){
		uvwgt = 0.0;
	}

	//Apply weighted averaging of profile values with input values
	*ptemp = tpdwgt*(*ptemp) + (1.0 - tpdwgt)*tin;
	*ppres = tpdwgt*(*ppres) + (1.0 - tpdwgt)*pin;
	*pdens = tpdwgt*(*pdens) + (1.0 - tpdwgt)*din;
	*puwin = uvwgt*(*puwin) + (1.0 - uvwgt)*uin;
	*pvwin = uvwgt*(*pvwin) + (1.0 - uvwgt)*vin;

	return;

}

double AuxProf::radll(double phi1, double thet1, double phi2, double thet2)
{

	//radll member function from AuxProf class
	//Returns great-circle distance (degrees) between two input lat-lon positions (in degrees)

	double r1[3], r2[3], r1xr2[3], r1xr2mag;
	double radll_out;

	//Get components of unit-magnitude vector toward 1st lat-lon
	r1[0] = cos(pi180*phi1)*cos(pi180*thet1);
	r1[1] = cos(pi180*phi1)*sin(pi180*thet1);
	r1[2] = sin(pi180*phi1);

	//Get components of unit-magnitude vector toward 2nd lat-lon
	r2[0] = cos(pi180*phi2)*cos(pi180*thet2);
	r2[1] = cos(pi180*phi2)*sin(pi180*thet2);
	r2[2] = sin(pi180*phi2);

	//Get cross product vector components from these two unit vectors
	r1xr2[0] = (r1[1] * r2[2]) - (r1[2] * r2[1]);
	r1xr2[1] = (r1[2] * r2[0]) - (r1[0] * r2[2]);
	r1xr2[2] = (r1[0] * r2[1]) - (r1[1] * r2[0]);

	//Get magnitude of cross product vector (Sine of great-circle distance)
	r1xr2mag = sqrt(r1xr2[0] * r1xr2[0] + r1xr2[1] * r1xr2[1] + r1xr2[2] * r1xr2[2]);

	//Get great-circle distance from cross-product magnitude
	if (r1xr2mag >= 1.0){
		radll_out = 90.0;
	}
	else {
		radll_out = asin(r1xr2mag) / pi180;
	}

	if ((r1[0] * r2[0] + r1[1] * r2[1] + r1[2] * r2[2]) < 0.0){
		radll_out = 180.0 - radll_out;
	}

	return radll_out;


}

void AuxProf::rig(double ch, double phir, double *g0, double *g, double *re, double *r0, double *ri, double *a, double *b)
{
	//rig member function from AuxProf class
	//Computes surface gravity g0 (m/s^2) and effective Earth radius re (km) from 
	//input geocentric latitude phir (radians).  Also computes gravity g (m/s^2)
	//and total radius ri (km) at input height ch (km).

	double g0a = 9.80616, g0b = 0.0026373, g0c = 0.0000059, rea = 3.085462e-3, reb = 2.27e-6, rec = 2.0e-9;
	double eps, c2phi, cphi2, c4phi;



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

}


