//P. White
//Map class for calculating atmospheric values in the middle atmosphere region. 
//As well as calculating concentrations data from atmosdat file.


#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include "Map.h"


using namespace std;


Map::Map()
{


}

void Map::gterp(int ih, double phi, double *p, double *d, double *t, double *u, double *dpy, double *dty){

	//gterp member function from the Map class
	//Interpolates pg, dg, tg, ug arrays of middle atmosphere values to p,
	//d, t, u, at height ih (integer km) and geocentric latitude phi
	//(degrees).  Also computes N-S pressure and temperature gradients dp/dy
	//dt/dy.


	int i, j, jp, mn1 = 2;
	double dl, r, r1, r2, tl, phif, chk;


	//Interpolates zonal mean data to height ih and geocentric latitude phi
	//Height index
	i = (ih - 20) / 5;
	i = int(fmin(20, i));

	//Lower latitude index
	j = int((phi + 100.0) / 10.0) - 1;
	j = int(fmax(0, j));
	j = int(fmin(17, j));
	//Upper latitude index
	jp = j + 1;

	//Check for density of temperature leq 0
	chk = perts.iperts.inits1.dg[i][j] * perts.iperts.inits1.tg[i][j] * perts.iperts.inits1.dg[i][jp] * perts.iperts.inits1.tg[i][jp];
	if (chk <= 0.0){
		*p = perts.iperts.inits1.pg[i][j];
		*d = perts.iperts.inits1.dg[i][j];
		*t = perts.iperts.inits1.tg[i][j];
	}

	//Geocentric latitude deviation from zonal mean position
	phif = (phi + 100.0 - 10.0 * (j + 1)) / 10.0;

	//Latitude interpolation
	tl = perts.iperts.inits1.tg[i][j] + (perts.iperts.inits1.tg[i][jp] - perts.iperts.inits1.tg[i][j])*phif;
	dl = perts.iperts.inits1.dg[i][j] + (perts.iperts.inits1.dg[i][jp] - perts.iperts.inits1.dg[i][j])*phif;
	*u = perts.iperts.inits1.ug[i][j] + (perts.iperts.inits1.ug[i][jp] - perts.iperts.inits1.ug[i][j])*phif;
	r1 = perts.iperts.inits1.pg[i][j] / (perts.iperts.inits1.dg[i][j] * perts.iperts.inits1.tg[i][j]);
	r2 = perts.iperts.inits1.pg[i][jp] / (perts.iperts.inits1.dg[i][jp] * perts.iperts.inits1.tg[i][jp]);


	//Interpolated gas constant
	r = r1 + (r2 - r1)*phif;
	//Pressure computed from interpolated gas constant
	*p = dl*r*tl;
	*d = dl;
	*t = tl;
	//dt/dy for vertical wind
	*dty = (perts.iperts.inits1.tg[i][jp] - perts.iperts.inits1.tg[i][j])*0.5;
	//dp/py for geostrophic winds
	*dpy = (perts.iperts.inits1.pg[i][jp] - perts.iperts.inits1.pg[i][j])*0.5;

	return;
}

void Map::pdtuv(int ih, double clat, double clon, double *ps, double *ds, double *ts, double *us, double *vs, double *dpx, double *dpy, double *dtx, double *dty)
{
	//pdtuv member function from the Map class
	//Interpolates stationary perturbations on geocentric latitude and 
	//longitude at height ih.

	double  dlon, dlat, xlon, dpxa, dpya, dtxa, dtya;
	int i, ip, k, j, jp;

	//Height index k
	k = (ih - 20) / 5;
	xlon = clon;

	//dlon - relative longitude deviation from corner reference point
	if (xlon >= 180.0) xlon = xlon - 360.0;
	if (xlon < (-180.0)) xlon = xlon + 360.0;
	//Lower longitude index j
	j = int((xlon + 200.0) / 20.0) - 1;
	dlon = (xlon + 200.0 - 20.0*(j + 1)) / 20.0;
	//Upper longitude index jp
	jp = j + 1;
	if (jp > 17) jp = 0;
	//Lower latitude index i
	i = int((clat + 100.0) / 10.0) - 1;

	//Upper latitude index ip
	ip = i + 1;
	ip = int(fmin(18, ip));

	//dlat - relative latitude deviation from corner reference location
	dlat = (clat - 10.0*(i + 1) + 100.0) / 10.0;
	//Pressure lat-lon interpolation
	*ps = double(perts.iperts.inits1.psp[k][i][j] + (perts.iperts.inits1.psp[k][ip][j] - perts.iperts.inits1.psp[k][i][j])*dlat + (perts.iperts.inits1.psp[k][i][jp] - perts.iperts.inits1.psp[k][i][j])*dlon + (perts.iperts.inits1.psp[k][ip][jp] - perts.iperts.inits1.psp[k][i][jp] - perts.iperts.inits1.psp[k][ip][j] + perts.iperts.inits1.psp[k][i][j])*dlat*dlon);
	//Density lat-lon interpolation
	*ds = double(perts.iperts.inits1.dsp[k][i][j] + (perts.iperts.inits1.dsp[k][ip][j] - perts.iperts.inits1.dsp[k][i][j])*dlat + (perts.iperts.inits1.dsp[k][i][jp] - perts.iperts.inits1.dsp[k][i][j])*dlon + (perts.iperts.inits1.dsp[k][ip][jp] - perts.iperts.inits1.dsp[k][i][jp] - perts.iperts.inits1.dsp[k][ip][j] + perts.iperts.inits1.dsp[k][i][j])*dlat*dlon);
	//Temperature lat-lon interpolation
	*ts = double(perts.iperts.inits1.tsp[k][i][j] + (perts.iperts.inits1.tsp[k][ip][j] - perts.iperts.inits1.tsp[k][i][j])*dlat + (perts.iperts.inits1.tsp[k][i][jp] - perts.iperts.inits1.tsp[k][i][j])*dlon + (perts.iperts.inits1.tsp[k][ip][jp] - perts.iperts.inits1.tsp[k][i][jp] - perts.iperts.inits1.tsp[k][ip][j] + perts.iperts.inits1.tsp[k][i][j])*dlat*dlon);
	//Zonal wind lat-lon interpolation
	*us = double(perts.iperts.inits1.usp[k][i][j] + (perts.iperts.inits1.usp[k][ip][j] - perts.iperts.inits1.usp[k][i][j])*dlat + (perts.iperts.inits1.usp[k][i][jp] - perts.iperts.inits1.usp[k][i][j])*dlon + (perts.iperts.inits1.usp[k][ip][jp] - perts.iperts.inits1.usp[k][i][jp] - perts.iperts.inits1.usp[k][ip][j] + perts.iperts.inits1.usp[k][i][j])*dlat*dlon);
	//Meridional wind lat-lon interpolation
	*vs = double(perts.iperts.inits1.vsp[k][i][j] + (perts.iperts.inits1.vsp[k][ip][j] - perts.iperts.inits1.vsp[k][i][j])*dlat + (perts.iperts.inits1.vsp[k][i][jp] - perts.iperts.inits1.vsp[k][i][j])*dlon + (perts.iperts.inits1.vsp[k][ip][jp] - perts.iperts.inits1.vsp[k][i][jp] - perts.iperts.inits1.vsp[k][ip][j] + perts.iperts.inits1.vsp[k][i][j])*dlat*dlon);
	//dpx - dp/dx for geostrophic winds
	dpxa = double((perts.iperts.inits1.psp[k][i][jp] - perts.iperts.inits1.psp[k][i][j]) / 4.0);
	*dpx = double(dpxa + ((perts.iperts.inits1.psp[k][ip][jp] - perts.iperts.inits1.psp[k][ip][j]) / 4.0 - dpxa)*dlat);
	//dpy - dp/dy for geostrophic winds
	dpya = double((perts.iperts.inits1.psp[k][ip][j] - perts.iperts.inits1.psp[k][i][j]) / 2.0);
	*dpy = double(dpya + ((perts.iperts.inits1.psp[k][ip][jp] - perts.iperts.inits1.psp[k][i][jp]) / 2.0 - dpya)*dlon);
	//dtx - dt/dx for vertical wind
	dtxa = double((perts.iperts.inits1.tsp[k][i][jp] - perts.iperts.inits1.tsp[k][i][j]) / 4.0);
	*dtx = double(dtxa + ((perts.iperts.inits1.tsp[k][ip][jp] - perts.iperts.inits1.tsp[k][ip][j]) / 4.0 - dtxa)*dlat);
	//dty - dt/dy for vertical wind
	dtya = double((perts.iperts.inits1.tsp[k][ip][j] - perts.iperts.inits1.tsp[k][i][j]) / 2.0);
	*dty = double(dtya + ((perts.iperts.inits1.tsp[k][ip][jp] - perts.iperts.inits1.tsp[k][i][jp]) / 2.0 - dtya)*dlon);

	return;
}

void Map::mapmod(double h, double phi, double thet, double ri, double g, double *pmm, double *dmm, double *tmm, double *umm, double *vmm, double *wmm, double *dtz)
{
	//mapmod member function from the Map class
	//Computes middle-atmosphere values of p, d, t, u, v and w from input
	//and arrays read from atmosdat file.

	int ihgb, ihga, ihsb, ihsa;
	double pi180, pi, dy5, dx5, pgb, dgb, tgb, ugb, dpygb, dtygb, pga, dga, tga, uga, dpyga, dtyga, psb, dsb, tsb, dpxsb, dpysb, dtxsb, dtysb, usb, vsb, psa, dsa, usa, vsa, tsa, dpxsa, dpysa, dtxsa, dtysa, umb, dtyg, spu, spv, ush, vsh, pmb, dmb, tmb, psh, dsh, tsh, dtxs, dtys, dtx, dty, cp;
	double hgb, hga, hsb, hsa, hs;
	using namespace std;
	pi = 3.1415926535897931;
	pi180 = pi / 180.0;

	//Distances for 5 degrees of lat, lon
	dy5 = 5000.0*ri*pi180;
	dx5 = dy5*cos(phi * pi180);
	dx5 = fmax(2000.0, dx5);

	//The following section is for zonal mean or mixed zonal mean
	//Upper height index
	ihgb = 5 * (int(h) / 5) + 5;
	//Upper height
	hgb = double(ihgb);

	//Zonal mean at upper height
	gterp(ihgb, phi, &pgb, &dgb, &tgb, &ugb, &dpygb, &dtygb);
	ihsb = 5 * (int(h) / 5) + 5;
	ihsb = int(fmin(90, ihsb));
	//Upper stationary perturbation height
	hsb = double(ihsb);

	//Stationary perturbations at upper height
	pdtuv(ihsb, phi, thet, &psb, &dsb, &tsb, &usb, &vsb, &dpxsb, &dpysb, &dtxsb, &dtysb);

	//Lower height index
	ihga = ihgb - 5;
	hga = double(ihga);
	//Zonal mean at lower height
	gterp(ihga, phi, &pga, &dga, &tga, &uga, &dpyga, &dtyga);
	ihsa = ihsb - 5;

	//Lower stationary perturbation height
	hsa = double(ihsa);
	hs = h;
	hs = fmin(90.0, hs);

	//Stationary perturbations at lower height
	pdtuv(ihsa, phi, thet, &psa, &dsa, &tsa, &usa, &vsa, &dpxsa, &dpysa, &dtxsa, &dtysa);

	//Interpolate for winds
	perts.iperts.interw(uga, dtyga, hga, ugb, dtygb, hgb, &umb, &dtyg, h);
	perts.iperts.interw(usa, vsa, hsa, usb, vsb, hsb, &spu, &spv, hs);
	ush = umb + spu;
	vsh = spv;

	//Zonal mean values height interpolation
	inter2(pga, dga, tga, ihga, pgb, dgb, tgb, ihgb, &pmb, &dmb, &tmb, h);

	//Stationary perturbations height interpolation
	perts.iperts.interz(psa, dsa, tsa, hsa, psb, dsb, tsb, hsb, &psh, &dsh, &tsh, hs);

	//Height interpolation of stationary perturbation, dt/dx and dt/dy
	perts.iperts.interw(dtxsa, dtysa, hsa, dtxsb, dtysb, hsb, &dtxs, &dtys, h);

	//Unperturbed (monthly mean) values for output
	*tmm = tmb*(1.0 + tsh);
	*pmm = pmb*(1.0 + psh);
	*dmm = dmb*(1.0 + dsh);
	*umm = ush;
	*vmm = vsh;

	//Total dt/dx
	dtx = dtxs**tmm;
	//Total dt/dy
	dty = *tmm*dtys + dtyg*(1.0 + tsh + dtys);
	cp = 7.0**pmm / (2.0**dmm**tmm);
	*dtz = (tgb*(1.0 + tsb) - tga*(1.0 + tsa)) / 5000.0;

	//Vertical mean wind
	*wmm = -cp*(*umm*dtx / dx5 + *vmm*dty / dy5) / (g + cp**dtz);

	return;

}

void Map::inter2(double p1, double d1, double t1, double z1, double p2, double d2, double t2, double z2, double *p, double *d, double *t, double z)
{

	//inter2 member function from Map class
	//Interpolates between p1, d1, t1 at height z1 and p2, d2, t2 at 
	//height z2 to output values of p, d, t at height z.
	//Checks for t1, d1, t2, d2 product = 0, for gas interpolation.

	double a, tz, b, pb, pz, r1, r2, r;

	//Set p=d=t=0 if some input values are negative or zero
	if (p1*d1*t1*p2*d2*t2 <= 0.0 || p2*p1 <= 0.0){
		*p = 0.0;
		*d = 0.0;
		*t = 0.0;
		return;
	}

	//Sets p,d,t = p1,d1,t1 if z1=z2
	if (abs(z1 - z2) <= 0.001){
		*p = p1;
		*d = d1;
		*t = t1;
		return;
	}

	a = (z - z1) / (z2 - z1);
	//Linear interpolation on t
	tz = t1 + a*(t2 - t1);
	r1 = p1 / (d1*t1);
	r2 = p2 / (d2*t2);
	//Linear interpolation on gas constant r
	r = (r2 - r1)*a + r1;
	//Logarithmic interpolation if temperature gradient = 0
	if (abs(t2 - t1) <= 0.001){
		b = (z2 - z1) / log(p1 / p2);
		*p = p1*exp((z1 - z) / b);
		*d = *p / (r*tz);
		*t = tz;
	}
	else {
		//Power law interpolation on pressure
		b = log(p2 / p1) / log(t1 / t2);
		pb = (t1 / tz);
		pz = p1*pow(pb, b);
		//Density from perfect gas law
		*d = pz / (r*tz);
		*p = pz;
		*t = tz;
	}

	return;
}

double Map::dedt(double t, double p)
{
	//dedt member function from Map class
	//Wexler formulation for the first derivative of saturation vapor 
	//pressure (Pa/K) as function of temperature T (kelvin), and pressure
	//(Pa), as given by Flatau et al., J. Appl. Meteorol., 31(12), 1507, 
	//Dec., 1992.

	const double g0 = -0.29912729E4, g1 = -0.60170128E4, g3 = -0.028354721,
		g4 = 0.17838301E-4, g5 = -0.84150417E-9, g6 = 0.44412543E-12,
		g7 = 2.858487;
	double esat, dedt_out;

	esat = wexler(t, p);

	if (esat <= 0.0) {
		dedt_out = 0.0;
		return dedt_out;
	}
	dedt_out = esat*(((g7 + (g3 + (2.0*g4 + (3.0*g5 +
		4.0*g6*t)*t)*t)*t)*t - g1)*t - 2.0*g0) / pow(t, 3);

	return dedt_out;

}


double Map::d2edt2(double t, double p)
{
	//d2edt2 member function from Map class
	//Wexler formulation for the 2nd derivative of saturation vapor pressure
	//(Pa/K^2) as a function of temperature T (kelvin), and pressure (Pa)
	//as given by the method of Flatau et al., J. Appl. Meteorol., 31(12),
	//1507, Dec., 1992.

	const double g0 = -0.29912729e4, g1 = -0.60170128e4, g4 = 0.17838301e-4,
		g5 = -0.84150417e-9, g6 = 0.44412543e-12, g7 = 2.858487;
	double esat, fde, d2edt2_out;

	esat = wexler(t, p);
	fde = dedt(t, p);

	if (esat <= 0.0){
		d2edt2_out = 0.0;
		return d2edt2_out;
	}

	d2edt2_out = (esat*(((-g7 + ((2.0*g4 + (6.0*g5 + 12.0*g6*t)*t)*t)*t)*t
		+ 2.0*g1)*t + 6.0*g0) / pow(t, 4)) + (pow(fde, 2) / esat);

	return d2edt2_out;

}


double Map::valint(int nlat, int nz, int ilat, int iz, double alpha,
	double beta, float **value)
{
	//valint member function from Map class
	//Latitude-height interpolation of concentrations that are dependent on
	//both variables.  Interpolation is logarithmic.

	double alphap, betap, val0, val1, val2, val3, vali, valint_out;

	alphap = 1.0 - alpha;
	betap = 1.0 - beta;
	val0 = log(value[ilat - 1][iz - 1]);
	val1 = log(value[ilat][iz - 1]);
	val2 = log(value[ilat - 1][iz]);
	val3 = log(value[ilat][iz]);
	vali = alphap*betap*val0 + alpha*betap*val1 + alphap*beta*val2 +
		alpha*beta*val3;

	valint_out = exp(vali);

	return valint_out;


}


double Map::valz(int nz, int iz, double beta, double *value)
{
	//valz member function from Map class
	//Height interpolation of concentrations that are altitude dependent only.
	//Interpolation is logarithmic.

	double betap, val0, val1, vali, valz_out;

	betap = 1.0 - beta;
	val0 = log(value[iz - 1]);
	val1 = log(value[iz]);
	vali = betap*val0 + beta*val1;
	valz_out = exp(vali);

	return valz_out;

}


double Map::wexler(double t, double p)
{

	//wexler member function from Map class
	//Wexler formulation for saturation vapor pressure (Pa) as function of
	//temperature T (kelvin) and pressure p (Pa), as given by Flatau et al.,
	//J. Appl Meteorol., 31(12), 1507, Dec., 1992, with correction factor (f3)
	//as given by Buck, J. Appl. Meteorol., 20(12), 1527, Dec. 1981.  
	//Wexler (T dry bulb) gives saturation vapor pressure.
	//Wexler (T dew point) gives actual vapor pressure.
	//Relative humidity (0-1) is wexler(Td)/wexler(T)

	const double g0 = -0.29912729E4, g1 = -0.60170128E4, g2 = 18.87643854,
		g3 = -0.028354721, g4 = 0.17838301E-4, g5 = -0.84150417E-9,
		g6 = 0.44412543E-12, g7 = 2.858487;
	double a = 1.0007, b = 3.46E-8, wexler_out;


	if (t <= 75.0){
		wexler_out = 1.0E-23;
	}
	else {
		wexler_out = exp((g0 + (g1 + (g2 + g7*log(t) + (g3 + (g4 + (g5 +
			g6*t)*t)*t)*t)*t)*t) / pow(t, 2))*(a + b*p);
	}

	return wexler_out;

}

void Map::afglconc(double z, double phi, double *watera)
{
	//afglconc member function from Map class
	//Evaluate the AFGL concentration data for a given location.

	const double xlat[8] = { -90.0, -60.0, -45.0, -15.0, 15.0, 45.0, 60.0,
		90.0 };


	int i, ilat, iz;
	double alpha, beta;

	if ((z < perts.iperts.inits1.hgta[0]) ||
		(z > perts.iperts.inits1.hgta[49])){
		*watera = 0.0;
		perts.iperts.inits1.ozone = 0.0;
		perts.iperts.inits1.nitrous = 0.0;
		perts.iperts.inits1.carbmon = 0.0;
		perts.iperts.inits1.methane = 0.0;
		perts.iperts.inits1.carbdiox = 0.0;
		perts.iperts.inits1.moloxy = 0.0;
		perts.iperts.inits1.nitrogen = 0.0;
		return;
	}

	//Find height index and height interpolation coefficient
	for (i = 0; i < 50; i++){
		if ((z >= perts.iperts.inits1.hgta[i]) &&
			(z <= perts.iperts.inits1.hgta[i + 1])){
			iz = i + 1;
			beta = (z - perts.iperts.inits1.hgta[i]) /
				(perts.iperts.inits1.hgta[i + 1] -
				perts.iperts.inits1.hgta[i]);
		}
	}

	//Find latitude index and latitude interpolation coefficient
	for (i = 0; i < 7; i++){
		if ((phi >= xlat[i]) && (phi <= xlat[i + 1])){
			ilat = i + 1;
			alpha = (phi - xlat[i]) / (xlat[i + 1] - xlat[i]);
		}
	}


	//Latitude-height interpolation for H2O, O3, N2O, CO, and CH4
	*watera = valint(8, 50, ilat, iz, alpha, beta, perts.iperts.inits1.h2oa);
	perts.iperts.inits1.ozone = valint(8, 50, ilat, iz, alpha, beta,
		perts.iperts.inits1.o3);
	perts.iperts.inits1.nitrous = valint(8, 50, ilat, iz, alpha, beta,
		perts.iperts.inits1.n2o);
	perts.iperts.inits1.carbmon = valint(8, 50, ilat, iz, alpha, beta,
		perts.iperts.inits1.co);
	perts.iperts.inits1.methane = valint(8, 50, ilat, iz, alpha, beta,
		perts.iperts.inits1.ch4);

	//Height interpolation only for CO2, O2, N2
	perts.iperts.inits1.carbdiox = valz(50, iz, beta,
		perts.iperts.inits1.co2);
	perts.iperts.inits1.moloxy = valz(50, iz, beta,
		perts.iperts.inits1.o2);
	perts.iperts.inits1.nitrogen = valz(50, iz, beta,
		perts.iperts.inits1.n2);

	return;
}

void Map::mapconc(double z, double pgh, double phi, double *o3m, double *h2om,
	double *sh2om, double *n2om, double *ch4m, double *oxym)
{
	//mapconc member function from Map class
	//Evaluates the MAP concentration data at a given point

	const double xlat[8] = { -90.0, -60.0, -45.0, -15.0, 15.0, 45.0, 60.0,
		90.0 };

	int i, ilat, iz, ilatm;
	double alpha, alpham, beta, alp;

	*o3m = 0.0;
	*h2om = 0.0;
	*sh2om = 0.0;
	*n2om = 0.0;
	*ch4m = 0.0;
	*oxym = 0.0;

	//Find latitude index and latitude interpolation coefficient
	//Index and coefficient for H2O
	for (i = 0; i < 7; i++){
		if ((phi >= xlat[i]) && (phi <= xlat[i + 1])){
			ilat = i + 1;
			alpha = (phi - xlat[i]) / (xlat[i + 1] - xlat[i]);
		}
	}

	//Index and coefficient for all other MAP species
	ilatm = int((phi + 90.0) / 10.0) + 1;

	if (ilatm > 18){
		ilatm = 18;
	}
	alpham = (phi + 100.0 - 10.0*ilatm) / 10.0;

	//Convert pressure to log scale
	alp = log(pgh);

	//Find ozone pressure level index and interpolation coefficient
	iz = 0;
	for (i = 0; i < 23; i++){
		if ((alp >= perts.iperts.inits1.po3map[i]) &&
			(alp <= perts.iperts.inits1.po3map[i + 1])){
			iz = i + 1;
			beta = (alp - perts.iperts.inits1.po3map[i]) /
				(perts.iperts.inits1.po3map[i + 1] -
				perts.iperts.inits1.po3map[i]);
		}
	}

	//Find MAP ozone value
	if (iz > 0){
		*o3m = valint(19, 24, ilatm, iz, alpham, beta,
			perts.iperts.inits1.o3map);
	}

	//Find H2O pressure level index and interpolation coefficient
	iz = 0;
	int count = 0;
	for (i = 0; i < 18; i++){
		if ((alp < perts.iperts.inits1.ph2omap[i]) ||
			(alp > perts.iperts.inits1.ph2omap[i + 1])){
			continue;
		}
		iz = i + 1;
		beta = (alp - perts.iperts.inits1.ph2omap[i]) /
			(perts.iperts.inits1.ph2omap[i + 1] -
			perts.iperts.inits1.ph2omap[i]);

		count++;
	}

	//Find MAP water vapor value and standard deviation
	if (iz > 0) {
		*h2om = valint(8, 19, ilat, iz, alpha, beta,
			perts.iperts.inits1.h2omap);
		*sh2om = valint(8, 19, ilat, iz, alpha, beta,
			perts.iperts.inits1.sh2omap);
	}

	//Find N2O and CH4 level index and interpolation coefficient
	iz = 0;
	for (i = 0; i < 16; i++){
		if ((alp >= perts.iperts.inits1.pmap31[i]) &&
			(alp <= perts.iperts.inits1.pmap31[i + 1])){
			iz = i + 1;
			beta = (alp - perts.iperts.inits1.pmap31[i]) /
				(perts.iperts.inits1.pmap31[i + 1] -
				perts.iperts.inits1.pmap31[i]);
		}
	}

	//Find MAP N2O and CH4 values (convert N2O to ppm)
	if (iz > 0){
		*n2om = valint(19, 17, ilatm, iz, alpham, beta,
			perts.iperts.inits1.n2omap)
			/ 1000.0;
		*ch4m = valint(19, 17, ilatm, iz, alpham, beta,
			perts.iperts.inits1.ch4map);
	}

	//Find atmoic oxygen height level and interpolation coefficient
	iz = 0;
	for (i = 0; i < 18; i++){
		if ((z <= perts.iperts.inits1.mapzox[i]) &&
			(z >= perts.iperts.inits1.mapzox[i + 1])){
			iz = i + 1;
			beta = (z - perts.iperts.inits1.mapzox[i]) /
				(perts.iperts.inits1.mapzox[i + 1] -
				perts.iperts.inits1.mapzox[i]);
		}
	}

	//Find MAP atomic oxygen values (convert to #/m^3)
	if (iz > 0){
		*oxym = valint(19, 19, ilatm, iz, alpham, beta,
			perts.iperts.inits1.oxymap)*1.0e06;
	}

	return;
}

void Map::larcwat(double z, double phi, double *waterl)
{
	//larcwat member function from Map class
	//Water vapor concentration from NASA Langley climatology, as a function
	//of height z (km) and latitude phi (degrees)

	int iz, ilat;
	double xlat, alpha, alphap, beta, betap, val0, val1, val2, val3, vali;

	if ((z < perts.iperts.inits1.hgtl[0]) ||
		(z > perts.iperts.inits1.hgtl[34])){
		*waterl = 0.0;
		return;
	}

	iz = 1 + (int)(z - perts.iperts.inits1.hgtl[0]);

	if (iz >= 35) {
		iz = 34;
	}

	xlat = (phi + 90.0) / 20.0;
	ilat = 1 + (int)xlat;

	if (ilat >= 10){
		ilat = 9;
	}

	alpha = xlat + 1.0 - ilat;
	alphap = 1.0 - alpha;
	beta = z - perts.iperts.inits1.hgtl[iz - 1];
	betap = 1.0 - beta;
	val0 = log(perts.iperts.inits1.h2ol[ilat - 1][iz - 1]);
	val1 = log(perts.iperts.inits1.h2ol[ilat][iz - 1]);
	val2 = log(perts.iperts.inits1.h2ol[ilat - 1][iz]);
	val3 = log(perts.iperts.inits1.h2ol[ilat][iz]);
	vali = alphap*betap*val0 + alpha*betap*val1 + alphap*beta*val2 +
		alpha*beta*val3;

	*waterl = exp(vali);

	return;
}


void Map::concvals(double z, double phi, int iyr, double pgh, double waterg,
	double *oxygen)
{
	//concvals member function from Map class
	//The main driver for the species concentration data.

	int iypt0;
	double waterl, ro3, watera, al0p3, al0p4, al10, al12p5, al2000, al1500,
		alp, al1, al1p5, alpg1, alpg2, o3m, h2om, n2om, ch4m, oxym;

	//Fractional rate of change (per year) for constituents
	const double rco2 = 0.5 / 100.0, rch4 = 0.9 / 100.0, rn2o = 0.3 / 100.0,
		rco = 0.7 / 100.0, ro30 = 0.23 / 100.0, ro340 = -0.5 / 100.0,
		pg1 = 3.0e04, pg2 = 2.5e04;

	al0p3 = log(0.3);
	al0p4 = log(0.4);
	al10 = log(10.0);
	al12p5 = log(12.5);
	al1500 = log(1500.0);
	al2000 = log(2000.0);
	al1 = log(1.0);
	al1p5 = log(1.5);
	alpg1 = log(pg1);
	alpg2 = log(pg2);
	iypt0 = iyr - 1976;

	//Evaluate concentrations from MAP Handbook vol 31
	mapconc(z, pgh, phi, &o3m, &h2om, &perts.iperts.inits1.sigwater,
		&n2om, &ch4m, &oxym);

	//Evaluate LaRC and AFGL concentrations unless z > 120 km
	if (z > perts.iperts.inits1.hgta[49]){
		return;
	}

	larcwat(z, phi, &waterl);
	afglconc(z, phi, &watera);

	//Update concentrations from 1976 to year of output:
	//CO2, CH4, CO, and N2O annual rates based on Table 14.5, page 306
	//of Graedel and Crutzen, "Atmospheric Change", 1993.
	perts.iperts.inits1.nitrous = perts.iperts.inits1.nitrous*
		pow(1.0 + rn2o, iypt0);
	perts.iperts.inits1.carbmon = perts.iperts.inits1.carbmon*
		pow(1.0 + rco, iypt0);
	perts.iperts.inits1.methane = perts.iperts.inits1.methane*
		pow(1.0 + rch4, iypt0);
	perts.iperts.inits1.carbdiox = perts.iperts.inits1.carbdiox*
		pow(1.0 + rco2, iypt0);

	//Height-dependent annual rate for ozone change based on 
	//Figure 17.1, page 373 of Graedel and Crutzen.
	if (z < 15.0){
		ro3 = ro30*((15.0 - z) / 15.0);
	}
	else if (z < 30.0){
		ro3 = 0.0;
	}
	else if (z < 40.0){
		ro3 = ro340*((z - 30.0) / 10.0);
	}
	else {
		ro3 = ro340*((120.0 - z) / 80.0);
	}

	perts.iperts.inits1.ozone = perts.iperts.inits1.ozone*
		pow(1.0 + ro3, iypt0);

	//Rates of change for water vapor, oxygen and nitrogen are assumed
	//to be zero.  Update concetrations to year of output for MAP data.
	iypt0 = iyr - 1981;
	n2om = n2om*pow(1.0 + rn2o, iypt0);
	ch4m = ch4m*pow(1.0 + rch4, iypt0);
	o3m = o3m*pow(1.0 + ro3, iypt0);

	//Get log of p for interpolations
	alp = log(pgh);

	//For water vapor use the following:
	//(a) AFGL value (watera) above the 0.01 mb height level
	//(b) Fair between AFGL and MAP vol 31 (h2om) for 0.015-0.01 bm
	//(c) MAP vol 31 value between 40.5 km and 0.015 mb
	//(d) Fair between MAP vol 31 and LaRC (waterl) for 39.5-40.5 km
	//(e) LaRC value between 39.5 km and 250 mb level
	//(f) Fair between LaRC and NCEP (waterg) values 300-250 mb
	//(g) NCEP (or RH-extended) (waterg) for surface to 300 mb
	if (pgh <= 1.5){
		if (pgh <= 1.0){
			perts.iperts.inits1.water = watera;
		}
		else {
			perts.iperts.inits1.water = ((alp - al1)*h2om +
				(al1p5 - alp)*watera) / (al1p5 - al1);
		}
	}
	else if (z >= 39.5){
		if (z >= 40.5){
			perts.iperts.inits1.water = h2om;

		}
		else {
			perts.iperts.inits1.water = (z - 39.5)*h2om +
				(40.5 - z)*waterl;
		}
	}
	else if (pgh <= pg1) {
		if (pgh <= pg2) {
			perts.iperts.inits1.water = waterl;
		}
		else {
			perts.iperts.inits1.water = ((alp - alpg2)*waterg +
				(alpg1 - alp)*waterl) / (alpg1 - alpg2);
		}
	}
	else {
		perts.iperts.inits1.water = waterg;
	}

	//Use AFGL O3, N2O, and CH4 if pressure > 2000 Pa; fair 1500 -
	//2000 Pa
	//Use AFGL O3 if p < 0.3 Pa; fair 0.3 - 0.4 Pa
	//Use AFGL N2O and CH4 if p < 10 Pa; fair 10 - 12.5 Pa
	//Otherwise use MAP concentrations
	if (pgh <= 2000.0){
		if (pgh >= 1500.0){
			perts.iperts.inits1.nitrous = ((alp - al1500)*
				perts.iperts.inits1.nitrous + (al2000 - alp)*n2om) /
				(al2000 - al1500);
			perts.iperts.inits1.methane = ((alp - al1500)*
				perts.iperts.inits1.methane + (al2000 - alp)*ch4m) /
				(al2000 - al1500);
			perts.iperts.inits1.ozone = ((alp - al1500)*
				perts.iperts.inits1.ozone + (al2000 - alp)*o3m) /
				(al2000 - al1500);
		}
		else if (pgh <= 12.5) {
			if (pgh >= 10.0) {
				perts.iperts.inits1.nitrous = ((alp - al10)*n2om +
					(al12p5 - alp)*perts.iperts.inits1.nitrous) /
					(al12p5 - al10);
				perts.iperts.inits1.methane = ((alp - al10)*ch4m +
					(al12p5 - alp)*perts.iperts.inits1.methane) /
					(al12p5 - al10);
			}

			if (pgh <= 0.4) {
				if (pgh >= 0.3) {
					perts.iperts.inits1.ozone = ((alp - al0p3)*o3m +
						(al0p4 - alp)*perts.iperts.inits1.ozone) /
						(al0p4 - al0p3);
				}
			}
			else {
				perts.iperts.inits1.ozone = o3m;
			}
		}
		else {
			perts.iperts.inits1.nitrous = n2om;
			perts.iperts.inits1.methane = ch4m;
			perts.iperts.inits1.ozone = o3m;
		}
	}

	//Use MET oxygen z > 100; MAP oxygen z < 90; fair otherwise
	if (z < 100.0){
		if (z < 90.0){
			*oxygen = oxym;
		}
		else {
			*oxygen = ((z - 90.0)*(*oxygen) + (100.0 - z)*oxym) / 10.0;
		}
	}

	return;

}


