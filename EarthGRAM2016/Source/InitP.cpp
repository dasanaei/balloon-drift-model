//Class for calculating Initial Perturbations
//P. White
#include "InitP.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

InitPert::InitPert()
{
	initializeMemberVariables();
}


InitPert::~InitPert()
{
}

void InitPert::initializeMemberVariables()
{

	//Initialize the member variables for the InitPert class.
	xbar1 = 0.0;
	zbar1 = 0.0;
	xl1 = 0.0;
	zl1 = 0.0;
	sxl1 = 0.0;
	szl1 = 0.0;
	prh1 = 0.0;
	drh1 = 0.0;
	trh1 = 0.0;
	urh1 = 0.0;
	vrh1 = 0.0;
	sp1l = 0.0;
	sp1s = 0.0;
	sd1l = 0.0;
	sd1s = 0.0;
	st1l = 0.0;
	st1s = 0.0;
	su1s = 0.0;
	su1l = 0.0;
	sv1s = 0.0;
	sv1l = 0.0;
	uds1 = 0.0;
	vds1 = 0.0;
	rp1s = 0.0;
	rd1s = 0.0;
	rt1s = 0.0;
	ru1s = 0.0;
	rv1s = 0.0;
	sd1 = 0.0;
	sp1 = 0.0;
	st1 = 0.0;
	su1 = 0.0;
	sv1 = 0.0;

	swh = 0.0, sw1 = 0.0, udlsrf = 0.0;
	rw1 = 0.0;
	rp1l = 0.0, rd1l = 0.0, rt1l = 0.0, ru1l = 0.0, rv1l = 0.0;

	waverand = 0.0;
	phidens = 0.0;
	ampfact = 0.0;

	for (int i = 0; i < 24; i++){
		mseeds[i] = 0.0;
	}

	mcarry = 0.0;
	mi24 = 0;
	mj24 = 0;
	hsrf1 = 0.0, usrf1 = 0.0, vsrf1 = 0.0, susrf1 = 0.0, svsrf1 = 0.0;
	tsrf1 = 0.0, uvtsrf1 = 0.0;
	pi = 4.0*atan(1.0);
	pi180 = pi / 180.0;
	twopi = 2.0 * pi, uvth = 0.0;



}


void InitPert::initpert(double h, double phi, double thet, double elt, int nr1)
{

	//Member function initpert from the class InitPert.  This member function calculates
	//the initial perturbations for pressure, density, temperature, east-west wind,
	//north-south wind, vertical wind.

	double pz1, rhoz1, tz1, uz1, vz1, wz1, tdz1, spz1, srhoz1, stz1,
		suz1, svz1, stdz1, rhn1, srhn1, vpn1, svpn1, spdavz1, spdsdz1,
		uvcorrz1;
	double udlph, sp, sd, st, plp1, dlp1, plmax = 0.85, ulp1, vlp1, tlp1;
	double rpd, rdtlim, rpdlim, rvec[7], suhl, svhl;
	double rpdl, stsr, r0;
	double sqrt2, dphid, av, hwn, vll, vlls, dxdlat, dxdlon, alpha, phidenz;
	double dphiu, dphiv, polefac;
	double sph = 0.0, sdh = 0.0, sdu, sdv, sth = 0.0, suh = 0.0, svh = 0.0, rrawt;
	string rra_id;
	int l, mc = 0, i = 0;
	double uvtsrf, ulphsrf, vlphsrf, rpinit;
	double st2, sp2, sd2, su2, sv2, profwgt, st3 = 0, sp3 = 0,
		sd3 = 0, su3 = 0, sv3 = 0, sitenear = inits1.sitenear, sitelim = inits1.sitelim;
	string home1 = inits1.home;
	string profile1 = inits1.profile;
	int iurra1 = inits1.iurra;


	//Calculate surface height from 1-degree topography map
	zinterp(phi, thet);
	//Calculate u-density and u-v correlations at surface height
	intr25(inits1.udl, inits1.uvt, hsrf1, phi, &udlsrf, &uvtsrf);


	//Get u-density and u-v cross correlations from global data
	intr25(inits1.udl, inits1.uvt, h, phi, &udlph, &uvth);

	if (h > 25.0){
		//Get interpolated standard deviations of pressure, density, temperature
		rterp(h, phi, &sp, &sd, &st);

		//Get interpolated standard deviations of east-west and north-south wind
		intruv(inits1.ur, inits1.vr, h, phi, &sdu, &sdv);

		sdu = sqrt(sdu);
		sdv = sqrt(sdv);

		sp1 = sp;
		sd1 = sd;
		st1 = st;
		su1 = sdu;
		sv1 = sdv;
	}

	else {


		//Call ncepmd from NCEPmods class to calculate NCEP data atmospheric variables
		//at surface height
		ncps.ncepmd(hsrf1, phi, thet, &pz1, &rhoz1, &tz1, &uz1, &vz1, &wz1, &tdz1, &spz1,
			&srhoz1, &stz1, &suz1, &svz1, &stdz1, &rhn1, &srhn1, &vpn1, &svpn1, &spdavz1,
			&spdsdz1, &uvcorrz1);

		//Atmospheric variables needed to calculate vertical wind sigma 
		usrf1 = uz1;
		vsrf1 = vz1;
		susrf1 = suz1;
		svsrf1 = svz1;
		tsrf1 = tz1;
		uvtsrf1 = uvcorrz1;
		stsrf1 = stz1;
		spdavsrf1 = spdavz1;
		spdsdsrf1 = spdsdz1;
		//Call ncepmd from NCEPmods class to calculate NCEP data atmospheric variables
		//at current height
		ncps.ncepmd(h, phi, thet, &pz1, &rhoz1, &tz1, &uz1, &vz1, &wz1, &tdz1, &spz1,
			&srhoz1, &stz1, &suz1, &svz1, &stdz1, &rhn1, &srhn1, &vpn1, &svpn1, &spdavz1,
			&spdsdz1, &uvcorrz1);

		//Standard deviations and u-v correlation needed to calculate initial 
		//perturbations.
		sp1 = spz1;
		sd1 = srhoz1;
		st1 = stz1;
		su1 = suz1;
		sv1 = svz1;
		uvth = uvcorrz1;

	}



	stsrf1 = tsrf1*stsrf1;

	//Use RRA data if conditions are met.
	if ((h <= 70.0) & (iurra1 > 0)){
		r0 = 8314.427;
		sph = sqrt(sp1);
		sdh = sqrt(sd1);
		sth = sqrt(st1);
		suh = su1;
		svh = sv1;

		//Call rrasigs from RRA class to calculate standard deviaitions needed for
		//calculating perturbations
		rras.rrasigs(h, phi, thet, inits1.mn, hsrf1, usrf1, vsrf1, tsrf1, uvtsrf1,
			susrf1, svsrf1, stsrf1, spdavsrf1, spdsdsrf1, r0, uvth, &sph, &sdh, &sth,
			&suh, &svh, rra_id, &rrawt);

		if (rras.z1[rras.num1 - 1] && rras.rrawt1 > 0.0){

			su1 = suh;
			sv1 = svh;
			uvth = rras.uvrra;
			if (rrawt > 0.0){
				hsrf1 = rras.rrsrf;
				usrf1 = rras.usrf;
				vsrf1 = rras.vsrf;
				susrf1 = rras.susrf;
				svsrf1 = rras.svsrf;
				uvtsrf1 = rras.uvtsrf;
				spdavsrf1 = rras.spdavsrf;
				spdsdsrf1 = rras.spdsdsrf;
			}
		}

		if (rras.z2[rras.num2 - 1] && rras.rrawt2 > 0.0){

			sp1 = pow(sph, 2.0);
			sd1 = pow(sdh, 2.0);
			st1 = pow(sth, 2.0);
			if (rrawt > 0.0){
				stsrf1 = rras.stsrf;
				tsrf1 = rras.tsrf;
			}
		}
	}


	//If iaux > 0, use auxiliary profile data
	if (inits1.iaux > 0){

		sp1 = sqrt(sp1);
		sd1 = sqrt(sd1);
		st1 = sqrt(st1);

		//Call rdprof and profsigs from AuxProf class to calculate standard deviations
		//needed to calculate initial perturbations.
		auxs.profsigs(h, phi, thet, st1, sp1, sd1, su1, sv1, &st2, &sp2,
			&sd2, &su2, &sv2, &profwgt);

		st1 = pow(st2, 2.0);
		sp1 = pow(sp2, 2.0);
		sd1 = pow(sd2, 2.0);
		su1 = su2;
		sv1 = sv2;



	}

	//Compute large scale fractional variances
	intr25(inits1.plp, inits1.dlp, h, phi, &plp1, &dlp1);
	intr25(inits1.ulp, inits1.vlp, h, phi, &ulp1, &vlp1);

	plp1 = min(0.99, plp1);
	dlp1 = min(plmax, dlp1);
	tlp1 = dlp1;
	vlp1 = min(plmax, vlp1);
	ulp1 = vlp1;


	//Compute large scale (l) and small scale (s) standard deviations
	//from the total standard deviations
	sp1 = sqrt(abs(sp1));
	sd1 = sqrt(abs(sd1));
	st1 = sqrt(abs(st1));
	sp1 = sp1*inits1.rpscale;
	sd1 = sd1*inits1.rpscale;
	st1 = st1*inits1.rpscale;
	su1 = su1*inits1.ruscale;
	sv1 = sv1*inits1.ruscale;
	sp1l = sqrt(plp1)*abs(sp1);
	sp1s = sqrt(1.0 - plp1)*abs(sp1);
	sd1l = sqrt(dlp1)*abs(sd1);
	sd1s = sqrt(1.0 - dlp1)*abs(sd1);
	st1l = sqrt(tlp1)*abs(st1);
	st1s = sqrt(1.0 - tlp1)*abs(st1);
	su1l = sqrt(ulp1)*su1;
	su1s = sqrt(1.0 - ulp1)*su1;
	sv1l = sqrt(vlp1)*sv1;
	sv1s = sqrt(1.0 - vlp1)*sv1;

	intr25(inits1.uds, inits1.vds, h, phi, &uds1, &vds1);

	//Get p-d cross correlation
	rpd = (sp1*sp1 + sd1*sd1 - st1*st1) / (2.0*sp1*sd1);

	//Get upper limit on rdt (which should be negative
	if (h < 10.0){
		rdtlim = -0.9 + 0.05*h;
	}
	else if (h < 15.0){
		rdtlim = -0.4 - 0.08*(h - 10.0);
	}
	else if (h < 30.0){
		rdtlim = -0.8 + 0.04*(h - 15.0);
	}
	else{
		rdtlim = -0.2;
	}

	//Get upper limit on rpd from rdt limit
	rpdlim = (sd1 + rdtlim*st1) / sp1;
	if (abs(rpdlim) > 0.99)rpdlim = copysign(0.99, rpdlim);
	if (rpd > rpdlim) rpd = rpdlim;

	//Initialize random number generator with nr1 seed input
	rcarin(nr1);

	//Get 7 random numbers
	rcarry(rvec, 7);

	//Evaluate initial small-scale perturbations values from the initial standard 
	//deviations
	rd1s = sd1s*ppnd(rvec[0], &l);

	//Adjust temperature perturbation, assuming rpdl = rpds
	rpdl = rpd*sp1*sd1 / (sp1l*sd1l + sp1s*sd1s);
	if (abs(rpdl) > 0.99)rpdl = copysign(0.99, rpdl);
	rp1s = sp1s*(rpdl*rd1s / sd1s + sqrt(1 - pow(rpdl, 2))*ppnd(rvec[1], &l));
	rt1s = rp1s - rd1s;
	stsr = sqrt(abs(pow(sp1s, 2) + pow(sd1s, 2) - 2 * rpdl*sp1s*sd1s));
	rt1s = rt1s*st1s / stsr;

	//Get initial small-scale velocity perturbations
	ru1s = su1s*ppnd(rvec[2], &l);
	rv1s = sv1s*ppnd(rvec[3], &l);


	//Compute initial large-scale wave amplitude factor and wave phase
	sqrt2 = sqrt(2.0);
	ampfact = 0.4808 + 0.96*rvec[6];
	phidens = 2.0*pi*rvec[4];

	//Compute initial vertical scale vll and wave numbers dxdlat, dxdlon
	waverand = rvec[5];
	dphid = ppnd(waverand, &l);
	if (abs(dphid) > 3.0) dphid = copysign(3.0, dphid);
	av = 25.0 + 3.0*dphid;
	hwn = double(int(round(4.0 + 0.833*dphid)));
	vll = av + 0.045*sqrt(pow(abs(h), 3));
	vll = min(200.0, vll);
	vlls = av + 0.045*sqrt(pow(abs(hsrf1), 3));
	dxdlat = hwn*pi / 180.0;
	dxdlon = hwn*pi / 180.0;
	alpha = dxdlon*thet + dxdlat*phi;
	dphiu = acos(rpdl);
	phidenz = phidens + 2.0*pi*h / vll;
	rd1l = sqrt2*sd1l*cos(alpha + phidenz)*ampfact;
	rp1l = sqrt2*sp1l*cos(alpha + phidenz + dphiu)*ampfact;

	//Get phase angles for large-scale U and V components
	dphiu = acos(udlph);
	dphiv = uvth / (0.02 + 0.98*vlp1);
	if (abs(dphiv) > 1.0)dphiv = copysign(1.0, dphiv);
	dphiv = dphiu + copysign(acos(dphiv), waverand - 0.5);

	ru1l = sqrt2*su1l*cos(alpha + phidenz + dphiu)*ampfact;
	rv1l = sqrt2*sv1l*cos(alpha + phidenz + dphiv)*ampfact;


	//Get large-scale temperature from pressure and density
	rt1l = rp1l - rd1l;
	//Adjust large-scale T perturbation using rpdl correlation
	stsr = sqrt(abs(pow(sp1l, 2) + pow(sd1l, 2) - 2.0*rpdl*sp1l*sd1l));
	rt1l = rt1l*st1l / stsr;
	//Adjust large-scale pressure using 2nd-order gas law
	rp1l = rd1l + rd1s + rt1l + rt1s + (rd1l + rd1s)*(rt1l + rt1s) - rp1s;

	//Insure that large scale perturbations approach zero at poles
	if (abs(phi) >= 85.0){
		polefac = cos(pi*(abs(phi) - 85.0) / 10.0);
		rd1l = rd1l*polefac;
		rp1l = rp1l*polefac;
		ru1l = ru1l*polefac;
		rv1l = rv1l*polefac;
	}

	//Substitute user-selected initial perturbations if initpert != 0
	if (inits1.initpert != 0){
		rpinit = inits1.rdinit + inits1.rtinit;
		rp1s = rpinit / 100.0 - rp1l;
		rd1s = inits1.rdinit / 100.0 - rd1l;
		rt1s = inits1.rtinit / 100.0 - rt1l;
		ru1s = inits1.ruinit - ru1l;
		rv1s = inits1.rvinit - rv1l;
	}


	//Get 3 random numbers
	rcarry(rvec, 3);

	//Compute initial vertical wind
	intrw(inits1.wr, h, &sw1);
	suh = su1;
	svh = sv1;
	swh = sw1;
	suhl = su1l;
	svhl = sv1l;

	//Get phase angles for large-scale U and V components at surface
	intr25(inits1.udl, inits1.uvt, hsrf1, phi, &udlsrf, &uvtsrf);
	dphiu = acos(udlsrf);

	intr25(inits1.ulp, inits1.vlp, hsrf1, phi, &ulphsrf, &vlphsrf);
	dphiv = uvtsrf1 / (0.02 + 0.98*vlphsrf);

	if (abs(dphiv) > 1.0) dphiv = copysign(1.0, dphiv);
	dphiv = dphiu + copysign(acos(dphiv), waverand - 0.5);

	//Calculate sigma-w from boundary layer (BL) model
	getsigw(h, phi, thet, alpha, alpha, vlls, dphiu, dphiv, elt);

	sw1 = swh;

	rw1 = sw1*ppnd(rvec[0], &l);
	//Substitute user-selected vertical wind perturbation if initpert not 0
	if (inits1.initpert != 0){
		rw1 = inits1.rwinit;
	}

	//Compute variable small-scale lengths xl1 = horizontal, zl1 = vertical
	intrw(inits1.xsigl, h, &sxl1);
	intrw(inits1.zsigl, h, &szl1);
	intrw(inits1.xlbar, h, &xbar1);
	intrw(inits1.zlbar, h, &zbar1);
	xl1 = xbar1 + sxl1*ppnd(rvec[1], &l);
	zl1 = zbar1 + szl1*ppnd(rvec[2], &l);

	//Total perturbation
	prh1 = rp1s + rp1l;
	drh1 = rd1s + rd1l;
	trh1 = rt1s + rt1l;
	urh1 = ru1s + ru1l;
	vrh1 = rv1s + rv1l;

}

void InitPert::intr25(float **xarray, float **yarray, double h, double phi, double *xint, double *yint)
{

	//intr25 member function from InitPert class
	//Finds large-scale fractional variance at height h (km), geocentric latitude phi (degrees), 
	//from ur and vr arrays.

	int i, ip, j, jp;
	double phi1, phi2, u1, u2, v1, v2, z1, z2;

	if (h < 95.0){
		//i - lower height index
		i = int(h) / 5;
	}
	else{
		i = 18 + (int(h) - 80) / 20;
	}

	i = min(24, i);
	//Upper height index
	ip = i + 1;
	ip = min(24, ip);
	//Lower latitude index
	j = int(phi + 110.0) / 20 - 1;
	//Upper latitude index
	jp = j + 1;
	jp = min(9, jp);
	//phi1 - lower latitude index for ur and vr arrays
	phi1 = (-110.0) + 20.0*(j + 1);
	//phi2 - upper latitude index for ur and vr arrays
	phi2 = (-110.0) + 20.0*(jp + 1);

	if (i > 18){
		//lower height for ur and vr array values
		z1 = 20.0*(i - 14);
	}
	else{
		z1 = 5.0*i;
	}
	if (ip > 18){
		//Uppert height for ur and vr array values
		z2 = 20.0*(ip - 14);
	}
	else{
		z2 = 5.0*ip;
	}


	//Interpolate on geocentric latitude at lower height
	//Interpolate with member function interw from Initpert Class
	interw(xarray[i][j], yarray[i][j], phi1, xarray[i][jp], yarray[i][jp], phi2, &u1, &v1, phi);

	//Interpolate on geocentric latitude at upper height
	interw(xarray[ip][j], yarray[ip][j], phi1, xarray[ip][jp], yarray[ip][jp], phi2, &u2, &v2, phi);

	//Interpolate on height
	interw(u1, v1, z1, u2, v2, z2, xint, yint, h);


}

void InitPert::interw(double u1, double v1, double z1, double u2, double v2, double z2, double *u, double *v, double z)
{
	//interw member function from the InitPert class
	//Linear interpolation between u1, v1 at height z1 and u2, v2 at 
	//height z2.  Output is u, v at height z.

	double a;
	//Sets u, v = u1, v1 if z1 = z2
	if (fabs(z1 - z2) <= 0.001) {
		*u = u1;
		*v = v1;
	}
	else{
		a = (z - z1) / (z2 - z1);
		*u = u1 + (u2 - u1)*a;
		*v = v1 + (v2 - v1)*a;
	}

	return;
}

void InitPert::rterp(double h, double phi, double *x, double *y, double *z)
{


	//rterp member function from InitPert class
	//Computes random perturbation standard deviations pressure, density, temperture 
	//at height h (km), geocentric latitude phi (degrees) from sigma arrays pr, dr, tr.

	double d1, d2, p1, p2, phi1, phi2, t1, t2, z1, z2;
	int i, ip, j, jp;


	if (h < 120.0){
		//i = lower height index
		i = int(0.2*(h + 5.0)) - 1;
	}
	else{
		i = 24 + int(0.05*(h - 120.0));
	}
	i = max(0, i);
	i = min(28, i);
	ip = i + 1;
	ip = min(28, ip);
	//Lower latitude index
	j = int((phi + 100.0) / 10.0) - 1;
	jp = j + 1;
	jp = min(18, jp);
	if (i > 24){
		//Lower height for pr,tr,dr arrays
		z1 = 120.0 + 20.0*(i - 24);
	}
	else
	{
		z1 = 5.0*(i + 1) - 5.0;
	}
	if (ip > 24){
		//Upper height for pr,tr,dr arrays
		z2 = 120.0 + 20.0*(ip - 24);
	}
	else
	{
		z2 = 5.0*(ip + 1) - 5.0;
	}

	phi1 = (-100.0) + 10.0*(j + 1);
	phi2 = (-100.0) + 10.0*(jp + 1);

	//Interpolate on latitude at lower height
	//interz member function from InitPert class
	interz(inits1.pr[i][j], inits1.dr[i][j], inits1.tr[i][j], phi1, inits1.pr[i][jp], inits1.dr[i][jp], inits1.tr[i][jp], phi2, &p1, &d1, &t1, phi);

	//Interpolate on latitude at upper height
	interz(inits1.pr[ip][j], inits1.dr[ip][j], inits1.tr[ip][j], phi1, inits1.pr[ip][jp], inits1.dr[ip][jp], inits1.tr[ip][jp], phi2, &p2, &d2, &t2, phi);

	//Interpolate on height using latitude interpolated values
	interz(p1, d1, t1, z1, p2, d2, t2, z2, x, y, z, h);

}

void InitPert::interz(double p1, double d1, double t1, double z1, double p2, double d2, double t2, double z2, double *p, double *d, double *t, double z)
{

	//interz member function from InitPert class
	//Linear interpolation between p1, d1, t1 at height z1 and p2, d2, t2
	//at height z2 to output values of p, d, t at height z.

	double a;

	if (abs(z1 - z2) <= 0.001){
		//Sets p, d, t = p1,d1,t1, if z1 = z2
		*p = p1;
		*d = d1;
		*t = t1;

	}
	else{
		a = (z - z1) / (z2 - z1);

		*t = t1 + (t2 - t1)*a;
		*d = d1 + (d2 - d1)*a;
		*p = p1 + (p2 - p1)*a;

	}

	return;
}

void InitPert::intruv(float **uarr, float **varr, double h, double phi, double *su, double *sv)
{

	//intruv member function from InitPert class
	//Finds random wind standard deviation at height h (km), geocentric latitude phi (degrees),
	//from ur and vr arrays.  

	int i, ip, j, jp;
	double phi1, phi2, u1, u2, v1, v2, z1, z2;

	//i - lower height index
	i = int(h) / 5;
	if (h >= 125.0){
		i = 24 + (int(h) - 120) / 20;
	}
	i = min(28, i);
	//Upper height index
	ip = i + 1;
	ip = min(28, ip);

	//Lower latitude index
	j = int(phi + 100.0) / 10 - 1;
	//Upper latitude index
	jp = j + 1;
	jp = min(18, jp);
	//phi1 - lower latitude for ur and vr array values
	phi1 = (-100.0) + 10.0*(j + 1);
	//phi2 - upper latitude for ur and vr array values
	phi2 = (-100.0) + 10.0*(jp + 1);
	if (i > 24){
		//Lower height for ur and vr arrays
		z1 = 20.0*(i - 18);
	}
	else{
		z1 = 5.0*i;
	}
	if (ip > 24){
		//Upper height for ur and vr array values
		z2 = 20.0*(ip - 18);
	}
	else{
		z2 = 5.0*ip;
	}

	//Interpolate on geocentric latitude at lower height
	//interw is member function from InitPert class
	interw(uarr[i][j], varr[i][j], phi1, uarr[i][jp], varr[i][jp], phi2, &u1, &v1, phi);

	//Interpolate on geocentric latitude at upper height
	interw(uarr[ip][j], varr[ip][j], phi1, uarr[ip][jp], varr[ip][jp], phi2, &u2, &v2, phi);

	//Interpolate on height
	interw(u1, v1, z1, u2, v2, z2, su, sv, h);

	return;
}

void InitPert::intrw(double sarr[], double h, double *sarrw)
{

	//intrw member function from InitPert class
	//Finds random wind standard deviation at height h (km), from wr arrays.

	int i, ip;
	double a, z1, z2;

	//i - lower height index
	i = int(h) / 5;
	if (h >= 125.0) i = 24 + (int(h) - 120) / 20;
	i = min(28, i);

	//Upper height index
	ip = i + 1;
	ip = min(28, ip);
	if (i == ip){
		*sarrw = sarr[i];
		return;
	}

	if (i > 24){
		//Lower height for wr array values
		z1 = 20 * (i - 18);
	}
	else {
		z1 = 5.0 * i;
	}

	if (ip > 24){
		//Upper height for wr array values
		z2 = 20.0*(ip - 18);
	}
	else {
		z2 = 5.0*ip;
	}

	//Interpolate on height
	a = (h - z1) / (z2 - z1);

	*sarrw = sarr[i] + (sarr[ip] - sarr[i])*a;



	return;
}

void InitPert::rcarin(int ijkl)
{

	//rcarin member function from InitPert class
	//Initializing routine for RCARRY.  Must be called before generating any pseudorandom
	//numbers with RCARRY.  The input value ijkl should be in the range: 0 <= ijkl <=
	//900,000,000.  Adapted from RMARIN subroutine of James, Comp. Phys. Comm. 60, 
	//329-344 (1990).

	int ii, ij, kl, i, j, k, l, m, jj;
	double s, t;

	ij = ijkl / 30082;
	kl = ijkl - 30082 * ij;
	i = (ij / 177) % 177 + 2;
	j = ij % 177 + 2;
	k = (kl / 169) % 178 + 1;
	l = kl % 169;
	for (ii = 0; ii < 24; ii++){
		s = 0.0;
		t = 0.5;
		for (jj = 0; jj < 24; jj++){
			m = (((i*j) % 179)*k) % 179;
			i = j;
			j = k;
			k = m;
			l = (53 * l + 1) % 169;
			if ((l*m) % 64 >= 32){
				s = s + t;
			}
			t = 0.5*t;
		}

		mseeds[ii] = s;
	}

	mi24 = 24;
	mj24 = 10;
	mcarry = 0.0;


}


void InitPert::rcarry(double *rvec, int lenv)
{

	//rcarry member function from InitPert class
	//Portable Pseudorandom number generator with period of about (1/48)*(2^24)^24 =
	//2^570 = 10^171.  Author F James, CERN, 1989.  Algorithm due to:  G. Marsgalia
	//and A. Zaman.  Code as in James, Comp. Phys. Comm., 60, 329-344 (1990).

	double twop24 = 16777216.0;
	double twom24 = 1.0 / twop24;
	int ivec;
	double uni;

	for (ivec = 0; ivec < lenv; ivec++){
		uni = mseeds[mi24 - 1] - mseeds[mj24 - 1] - mcarry;
		if (uni < 0.0){
			uni = uni + 1.0;
			mcarry = twom24;
		}
		else {
			mcarry = 0.0;
		}
		mseeds[mi24 - 1] = uni;
		mi24 = mi24 - 1;
		if (mi24 == 0){
			mi24 = 24;
		}
		mj24 = mj24 - 1;
		if (mj24 == 0){
			mj24 = 24;
		}
		//Avoid random number of exactly zero (see James, p.344)
		if (uni == 0.0){
			uni = mseeds[mi24 - 1] * twom24;

			if (uni == 0.0){
				uni = pow(2.0, -48);
			}
		}

		rvec[ivec] = uni;

	}


}

double InitPert::correl(double x)
{
	//correl member function from InitPert class
	//Correlation function for relative seperation distance x (r/L)

	double correl_out;

	//Exponential function, with defaults for small and large x
	if (abs(x) < 1.0e-7) {
		correl_out = 0.9999999;
		return correl_out;
	}
	else if (abs(x) > 23.0) {
		correl_out = exp(-23.0);
		return correl_out;
	}
	else {
		correl_out = exp(-abs(x));
		return correl_out;
	}

}

double InitPert::ptail(double x1)
{
	//ptail member function from InitPert class
	//Tail probability from Normal distribution (Integral x to infinity of Gaussian),
	//where x is normalized deviate (deviation from mean divided by standard deviation
	//about mean).  Note Ptail = (1/2)*(1 - erf(x/sqrt(2)), where erf is the error function,
	//or Ptail = (1/2)erfcomp(x/sqrt(2)), where erfcomp is the complementary error function.
	//This implementation of erfcomp is based on Chebyshev fitting.  Coefficient values
	//are from "Numerical Recipes in Fortran", 2nd Edition, page 214.

	double ptail_out, y, q, t, erfcomp;

	y = x1 / sqrt(2.0);

	q = abs(y);

	t = 1.0 / (1.0 + q / 2.0);

	erfcomp = t*exp((-q*q) - 1.26551223 + t*(1.00002368 + t*(0.37409196 +
		t*(0.09678418 + t*(-0.18628806 + t*(0.27886807 + t*(-1.13520398 +
		t*(1.48851587 + t*(-0.82215223 + t*(0.17087277))))))))));

	if (x1 < 0.0){
		erfcomp = 2.0 - erfcomp;
	}

	ptail_out = erfcomp / 2.0;

	return ptail_out;
}

double InitPert::ppnd(double p, int *ifault)
{
	//ppnd member function from InitPert class
	//Algorithm AS 111 Appl. Statist. (1977) vol. 26, p. 118
	//Produces normal deviate corresponding to lower tail area of p.  Returns ifault = 1
	//in input p >= 1 or <= 0, ifault = 0 otherwise.  If ifault = 1, ppnd value is set
	//to 0.  Single precision version with error epsilon = 2^(-31).  The hash sums are
	//the sums of the moduli of the coefficients.  They have no inherent meanings, but
	//are included for use in checking transpositions.


	double zero = 0.0, split = 0.42, half = 0.5, one = 1.0;
	double q, r, ppnd_out;

	double	a0 = 2.50662823884, a1 = -18.61500062529, a2 = 41.39119773534, a3 = -25.44106049637;
	double b1 = -8.47351093090, b2 = 23.08336743743, b3 = -21.06224101826, b4 = 3.13082909833;

	//Hash sum for a and b = 143.70383558076

	double c0 = -2.78718931138, c1 = -2.29796479134, c2 = 4.85014127135, c3 = 2.32121276858,
		d1 = 3.54388924762, d2 = 1.63706781897;


	//Hash sum for c and d = 17.43746520924

	*ifault = 0;
	q = p - half;

	if (fabs(q) <= split){
		r = q*q;
		ppnd_out = q*(((a3*r + a2)*r + a1)*r + a0) / ((((b4*r + b3)*r + b2)*r + b1)*r + one);

		return ppnd_out;
	}

	r = p;
	if (q > zero) {
		r = one - p;
	}

	if (r > zero){
		r = sqrt(-log(r));
		ppnd_out = (((c3*r + c2)*r + c1)*r + c0) / ((d2*r + d1)*r + one);

		if (q < zero){
			ppnd_out = -ppnd_out;
		}

		return ppnd_out;

	}

	*ifault = 1;

	ppnd_out = zero;

	return ppnd_out;


}


void InitPert::getsigw(double h, double phi, double thet, double alpha, double alphav, double vlls, double dphiu, double dphiv, double elt)
{

	//getsigw member function from InitPert class
	//Uses boundary layer (BL) model to get sigma-w as a function of height, lat, lon,
	//wind speed at 10-m, surface roughness (z0), and atmospheric stability (which
	//depends on time of day).

	int iwb, i, j, m, intday, lc;

	ofstream bltest;

	double sqrt2, phidenz, sul2b, svl2b, ul2b, vl2b, hbl, hspbl, chb = 0.0, spdsrf,
		swb = 0.0, xmin, xlatm, sha, sda, dd, dy, sra, raloc, xmjd, sdec, cdec, slat, clat,
		el, coriol, elmn, elmd, blfact, lst, nri, nrimn, f, s, aofs, ool, psi10, denom,
		vonk, a, b, bvfsq, theta, hn, ha = 0.0, hb = 0.0, r = 0.0, sigrat = 0.0, siglim = 0.0, third, gamma, fact = 0.0,
		dz0m = 0.0, rday = 0.0, ulphsrf, vlphsrf, ch = 0.0, z0, sec;

	int imin, ihr;

	double zruf[14] = { 3.2e-04, 0.6, 0.48, 0.42, 0.0056, 0.45, 0.12,
		0.046, 0.015, 0.042, 0.065, 0.45, 1.0e-04, 3.2e-04 };

	if (elt == 0.0){
		imin = inits1.mino;
		sec = inits1.seco;
		ihr = inits1.ihro;
	}
	else{
		imin = inits1.mino + int(elt) / 60;
		sec = inits1.seco + elt;
		sec = int(sec) % 60;
		ihr = inits1.ihro + imin / 60;
		imin = imin % 60;
	}

	sqrt2 = sqrt(2.0);
	vonk = 0.4;
	third = 1.0 / 3.0;

	//Default landcode and z0
	lc = 99;
	z0 = 9.99999;

	//Use sigm-w from atmosdat table if height > 10 km
	if (h > 10.0) {
		intrw(inits1.wr, h, &swh);
		return;
	}
	//Get large-scale wind perturbations at surface for boundary-layer model
	phidenz = phidens + twopi*hsrf1 / vlls;
	intr25(inits1.ulp, inits1.vlp, hsrf1, phi, &ulphsrf, &vlphsrf);
	sul2b = inits1.ruscale*susrf1*sqrt(vlphsrf);
	svl2b = inits1.ruscale*svsrf1*sqrt(vlphsrf);
	ul2b = sqrt2*sul2b*cos(alpha + phidenz + dphiu)*ampfact;
	vl2b = sqrt2*svl2b*cos(alphav + phidenz + dphiv)*ampfact;

	//Compute surface (10 m) wind speed
	spdsrf = sqrt(pow((usrf1 + ul2b), 2.0) + pow((vsrf1 + vl2b), 2.0));
	if (spdsrf < 0.1) spdsrf = 0.1;

	//Index for first sigma-w above boundary layer
	iwb = 2;
	if (hsrf1 > 2.0) iwb = 3;

	//Get solar hour angle (sha) and solar declination (sda)
	xmin = imin + sec / 60.0;
	xlatm = phi;

	//Note: MET thermosphere tme changes input latitude in degrees to 
	//output in radians
	mets.tme(inits1.mn, inits1.ida, inits1.iyr, ihr, xmin, thet, &xlatm, &sda, &sha, &dd, &dy, &sra, &raloc,
		&xmjd);

	//Get some sines and cosines
	slat = sin(phi*pi180);
	clat = cos(phi*pi180);
	sdec = sin(sda);
	cdec = cos(sda);

	//Coriolis parameter (absolute value) with minimum near 20 degrees
	//latitude
	coriol = 1.4584e-04*abs(slat);
	if (coriol < 5.0e-05) coriol = 5.0e-05;

	//Get solar elevation angle (deg.)
	el = asin(sdec*slat + cdec*clat*cos(sha)) / pi180;

	//Convert solar hour angle to degrees
	sha = sha / pi180;

	//Solar elevation at midnight (deg.)
	elmn = asin(sdec*slat - cdec*clat) / pi180;
	if (elmn == 0.0) elmn = 1.0e-03;

	//Solar elevation at midday (deg.)
	elmd = asin(sdec*slat + cdec*clat) / pi180;
	if (elmd == 0.0) elmd = 1.0e-03;

	//Height factor for unstable B.L. during early daytime
	blfact = 1.0;
	if ((elmn < 0.0) & (el > 0.0) & (el <= elmd) & (sha < 0.0)){
		blfact = 0.3 + 0.7*el / elmd;
	}

	//Local solar time (hours)
	lst = 12.0 + sha / 15.0;
	if (lst < 0.0) lst = lst + 24.0;
	if (lst > 24.0) lst = lst - 24.0;

	//Get (modified) Net Radiation Index (nri) from simplified versions
	//of Table 4-7 and 4-8 of Justus (1978)

	//Case when dark all 24 hours
	if (sdec*slat + cdec*clat < 0.0){
		nri = -3.5;
	}
	else {
		if ((sha < 0.0) & (el > elmn) & (el < 0.0)){
			//Special interpolated nri from midnight to dawn
			nrimn = 0.5 + 6.1154e-02*elmn - 2.6390e-06*pow(elmn, 3);
			nri = nrimn + (3.5 + nrimn)*(el / elmn - 1.0);
		}
		else {
			//Usual nri versus elevation (-90 to +90)
			nri = 0.5 + 6.1154e-02*el - 2.6390e-06*pow(el, 3);
		}
	}

	if (nri < -3.5) nri = -3.5;
	if (nri > 4.5) nri = 4.5;

	//Convert declination to degrees
	sda = sda / pi180;

	//Use input z0 value (z0in) or compute from model
	if (inits1.z0in > 0.0){
		z0 = inits1.z0in;
		lc = 99;
	}
	else {
		if (inits1.z0in == 0.0){
			lc = 0;
		}
		else {
			//Get surface type from landcode array
			i = int(thet + 180.0);
			if (i > 359)i = 0;
			j = int(phi + 90.0);
			if (j > 179)j = 179;
			lc = inits1.landcd[i][j];
		}


		//Get surface type from landcode value
		if (lc == 0){
			//use iterative solution if water type
			z0 = zruf[0];

			for (i = 0; i < 4; i++){
				z0 = 2.612e-04*pow((spdsrf / log(10.0 / z0)), 2);
			}
		}
		else {
			z0 = zruf[lc];
		}

		if (z0 < 1.0e-05) z0 = 1.0e-05;

		//Model to increase z0 when in high-altitude (mountainous) terrain
		if (hsrf1 < 1.5){
			dz0m = 0.0;
		}
		else {
			dz0m = hsrf1 - 1.5;
		}
		z0 = z0 + dz0m;
		if (z0 > 3.0) z0 = 3.0;

	}

	//Get stability category S and inverse Monin-Obukhov length 1/L = 
	//ool, versus wind speed u, and Net Radiation Index (nri).
	//Wind speed factor for stability category [for simplified version
	//of Table 4-7 of Justus (1978)]
	if (spdsrf <= 6.0){
		f = (1.0 - spdsrf / 7.5);
	}
	else {
		f = 0.2*exp(12.0 - 2.0*spdsrf);
	}

	//Stability category for wind speed and net radiation index
	//[simplified version of Justus (1978) Pages 57-58, neglecting
	//cloud cover and ceiling information, unavailable in GRAM]
	s = 4.229 - nri*f;
	if (s < 0.5) s = 0.5;
	if (s > 7.5) s = 7.5;

	//Inverse fo Monin-Obukhov length (ool=1/L, in m**-1) vs stability 
	//and z0 from Figure 4-9 of Justus (1978)
	aofs = -0.2161 + 0.0511*s;
	ool = 0.25*aofs*log10(10.0 / z0);
	if (ool >= 0.0){
		//Stable case psi(10/L) [psi(z/L) at z = 10 m]
		psi10 = -50.0*ool;
	}
	else {
		//Simplified version of Paulson unstable relation for psi(10/L),
		//from equation (5) of Hsu et al. (1999) 
		psi10 = 1.0496*pow((-10.0*ool), 0.4591);
	}

	//Friction velocity ustar (m/s) versus z0 and psi(10/L)
	denom = log(10.0 / z0) - psi10;
	if (denom < 0.1) denom = 0.1;
	ustar = vonk*spdsrf / denom;
	if (ustar < 0.04 / denom) ustar = 0.04 / denom;
	if (ustar > 4.0) ustar = 4.0;

	//Calculate hbl = height of boundary layer
	//Method for stable-to-neutral, Section 2.1, Sugiyama and Nasstrom,
	//1999, with neutral BL height changed from hN = 0.2*ustar/coriol to 
	//hN = ustar*(4.0*B/(coriol*BVfsq))**third.  As expressed in Seibert 
	//(2000), with time dependent term d(hbl)/dt replaced by hbl*coriol/2.
	a = 0.4;
	b = 20.0;
	theta = tsrf1;

	gamma = 3.3e-03;
	bvfsq = (9.8 / theta)*gamma;
	hn = ustar*pow(((4.0*b) / (coriol*bvfsq)), third);
	if (ool >= 0.0){
		hbl = 2.0*hn / (1.0 + sqrt(1.0 + 4.0*hn*ool));
	}
	else {
		ha = hn;
		r = -(1.0 + 2.0*a)*ool / (2.0*b*vonk);
		for (m = 0; m < 7; m++){
			hb = hn*pow((1.0 + r*ha), third);
			if (hb < hb) hb = hn;
			if (abs(hb - ha) <= 5.0) goto Fourty;
			ha = hb;
		}

		//BL height hbl in meters
	Fourty: hbl = hb;
	}

	//Factor blfact for unstable conditions, dawn to mid-morning,
	//to lessen sudden increase in B.L. height at dawn
	if (hbl > 3000.0) hbl = 3000.0;
	if (ool < -0.0001) hbl = blfact*hbl;
	if (hbl < 200.0) hbl = 200.0;

	//Height (above MSL) for top of BL (in km)
	hspbl = hsrf1 + (hbl / 1000.0);

	if (h < 5.0*(iwb - 1.0)){
		//Use height (AGL) at top of BL (in meters) if current height 
		//above BL
		chb = hbl;
		//Use actual height (AGL in meters) if current height within BL
		if (h < hspbl) chb = (h - hsrf1)*1000.0;
		if (chb < 0.0) chb = 0.0;
		//Revised BL sigma-w from z0 and ustar
		//Ratio sigma-w/ustar
		if (ool >= 0.0){
			//Eq. 1.33, Kaimel and Finnigan (1994); larger range from
			//Pahlow et al. (2001)
			sigrat = 1.25*(1.0 + 0.2*chb*ool);
			siglim = 3.0;
			if (sigrat > siglim) sigrat = siglim;
		}
		else {
			//Limit sigma-w/ustar with convective velocity wstar: Eq. 1.48,
			//1.51 at z/zi = 0.5, and Fig 1.10 of Kaimal and Finnigan (1994)
			siglim = 0.62*pow((-hbl*ool / vonk), third);
			if (siglim < 1.25) siglim = 1.25;
			//Equation (2), Page 161, Panofsky and Dutton (1984)
			sigrat = 1.25*pow((1.0 - 3.0*chb*ool), third);
			if (sigrat > siglim) sigrat = siglim;
		}

		//sigma-w values from ustar and ratios
		swb = ustar*sigrat;
		if (swb < 0.1) swb = 0.1;
		if (swb > 5.0) swb = 5.0;
		//If height above top of BL, interpolate between sigma-w at top
		//of BL and next highest level
		if (h > hspbl){
			fact = (h - hspbl) / (5.0*(iwb - 1.0) - hspbl);
			swh = swb + (inits1.wr[iwb - 1] - swb)*fact;

		}
		else {
			fact = 0.0;
			swh = swb;
		}

		//Get integer day (local solar time from start)
		rday = inits1.ida + (60.0*(ihr - inits1.ihro) + xmin - inits1.mino) / 1440.0 - lst / 24.0;
		intday = 1 + int(rday);
		rday = intday + lst / 24.0;

		//Output boudary layer parameters if ibltest > 1
		if (inits1.ibltest > 1){




			inits1.ibl << fixed << setw(7) << setprecision(3) << rday << setw(6) <<
				setprecision(2) << lst << setw(7) << setprecision(3) << h << setw(8) <<
				setprecision(3) << phi << setw(9) << setprecision(3) << thet << setw(6) <<
				setprecision(3) << hsrf1 << setw(7) << setprecision(2) << spdsrf <<
				setw(3) << lc << setw(8) << setprecision(5) << z0 << setw(6) << setprecision(1) <<
				elmn << setw(6) << setprecision(1) << el << setw(6) << setprecision(1) <<
				elmd << setw(8) << setprecision(2) << sha << setw(6) << setprecision(3) << blfact <<
				setw(6) << setprecision(2) << nri << setw(5) << setprecision(2) << s << setw(9) <<
				setprecision(5) << ool << setw(6) << setprecision(3) << ustar << setw(9) <<
				setprecision(5) << bvfsq << setw(6) << setprecision(0) << hn << setw(6) <<
				setprecision(0) << hbl << setw(6) << setprecision(0) << chb << setw(7) << setprecision(3) << sigrat <<
				setw(6) << setprecision(2) << swb << setw(6) << setprecision(2) << swh <<
				setw(9) << setprecision(2) << spdavsrf1 << setw(9) << setprecision(2) << spdsdsrf1 <<
				setw(7) << setprecision(2) << tsrf1 << setw(7) << setprecision(2) << stsrf1 <<

				'\n';
		}


	}

	else {
		intrw(inits1.wr, h, &swh);
	}







}


void InitPert::zinterp(double clat, double clon)
{
	//zinterp member function form InitPert class
	//Interpolates surface height array ztopo to geocentric latitude and 
	//longitude.

	using namespace std;

	double dlat, dlon, xlon, xloni, xlonp, xlat, xlati, xlatp, dlatp,
		dlonp;
	int ip, jp, i, j;



	xlon = clon;
	//Lower longitude index, i
	if (xlon > 180.0) xlon = xlon - 360.0;
	if (xlon < -180.0) xlon = xlon + 360.0;
	i = int(xlon + 180.5) - 1;
	if (i < 0) i = i + 360;
	if (i > 359) i = i - 360;

	//Upper longitude index, ip
	ip = i + 1;
	if (ip > 359) ip = ip - 360;
	xloni = i - 179.5;
	xlonp = ip - 179.5;

	//dlon - relative longitude deviation from corner reference point
	dlon = xlon - xloni;
	if (dlon < 0.0) dlon = dlon + 360.0;
	dlonp = 1.0 - dlon;

	//Lowe latitude index, j
	xlat = clat;
	if (xlat < -90.0) xlat = -90.0;
	if (xlat > 90.0) xlat = 90.0;
	j = int(xlat + 89.5);
	if (j < 0) j = 0;
	if (j > 178) j = 178;

	//Upper latitude index, jp
	jp = j + 1;
	xlati = j - 89.5;
	xlatp = jp - 89.5;

	//dlat - relative latitude deviation from corner reference location
	dlat = xlat - xlati;
	if (dlat < 0.0) dlat = 0.0;
	if (dlat > 1.0) dlat = 1.0;
	dlatp = 1.0 - dlat;

	//lat-lon interpolation of topographic height
	hsrf1 = (inits1.ztopo[i][j] / 1000.0) * dlonp*dlatp + (inits1.ztopo[ip][j] / 1000.0)
		* dlon*dlatp + (inits1.ztopo[i][jp] / 1000.0) * dlonp*dlat +
		(inits1.ztopo[ip][jp] / 1000.0) * dlat*dlon;

}






