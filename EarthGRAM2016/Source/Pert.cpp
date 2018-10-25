//Pert class for calculating perturbations for pressure, density, temperature,
//east-west wind, north-south wind, vertical wind.
//P. White

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "Pert.h"

using namespace std;


Pert::Pert()
{
	initializeMemberVariables();
}


Pert::~Pert()
{
}

void Pert::initializeMemberVariables()
{

	//Initialize member variables for Pert class

	dz = 0.0;
	dx = 0.0;
	sd1s = 0.0;
	rp1s = 0.0, rd1s = 0.0, rt1s = 0.0, ru1s = 0.0, rv1s = 0.0;
	rp1l = 0.0, rd1l = 0.0, rt1l = 0.0, ru1l = 0.0, rv1l = 0.0;
	sp1s = 0.0;
	st1s = 0.0;
	su1s = 0.0;
	sv1s = 0.0;
	sd1l = 0.0;
	phi = 0.0;
	thet = 0.0;
	sp1l = 0.0;
	uds1 = 0.0;
	vds1 = 0.0;
	prh = 0.0;
	drh = 0.0;
	trh = 0.0;
	urh = 0.0;
	vrh = 0.0;
	phi1 = 0.0;
	thet1 = 0.0;
	dphi = 0.0;
	dthet = 0.0;
	sph = 0.0;
	sdh = 0.0;
	sth = 0.0;
	suh = 0.0;
	svh = 0.0;
	wrh = 0.0;
	xlh = 0.0, zlh = 0.0;
	pi = 4.0*atan(1.0);
	pi180 = pi / 180.0;
	elt = 0.0;
	isev = 0;
	prhl = 0.0, drhl = 0.0, trhl = 0.0, urhl = 0.0, vrhl = 0.0, prhs = 0.0, drhs = 0.0,
		trhs = 0.0, urhs = 0.0, vrhs = 0.0, sphs = 0.0, sphs = 0.0, sths = 0.0, suhs = 0.0,
		svhs = 0.0, sphl = 0.0, sdhl = 0.0, sthl = 0.0, suhl = 0.0, svhl = 0.0, uvt2 = 0.0;
}

void Pert::initsigs()
{

	//initsigs member function for Pert class
	//Sets initial standard deviations and perturbations used for calculating
	//perturbations.

	sd1s = iperts.sd1s;
	sp1s = iperts.sp1s;
	st1s = iperts.st1s;
	su1s = iperts.su1s;
	sv1s = iperts.sv1s;
	uds1 = iperts.uds1;
	vds1 = iperts.vds1;
	sd1l = iperts.sd1l;
	sp1l = iperts.sp1l;
	st1l = iperts.st1l;
	su1l = iperts.su1l;
	sv1l = iperts.sv1l;
	rd1s = iperts.rd1s;
	rp1s = iperts.rp1s;
	rt1s = iperts.rt1s;
	ru1s = iperts.ru1s;
	rv1s = iperts.rv1s;
	rp1l = iperts.rp1l;
	rd1l = iperts.rd1l;
	rt1l = iperts.rt1l;
	ru1l = iperts.ru1l;
	rv1l = iperts.rv1l;



}

void Pert::pert1(double h, double phi, double thet, int iupdate, int nm)
{

	//pert1 member function from Pert class
	//Large-scale and small-scale perturbation values in pressure, density,
	//temperature, horizontal and vertical wind components.


	double xbarh, zbarh, sxlh, szlh, xmin, zmin, xlen, zlen, rvec[7];
	double tpers, as, bs, alph, beta, gamm, rcor, denom;
	double fbig, xtail, ztail, ptailx, scalfact, fsmall, xmax, zmax, zfactor, vdsm, vts, vus, hls, densfact;
	double twopi, pi, rds, rps, rvs, rpd, rdtlim, sp, sd, st, su, sv;
	int l;
	double plf = 3.0, bv = 0.045, slim = 1.0e-5;
	double plph, dlph, vlph, ulph, udl2, uds2, vds2, plmax = 0.85, tlph, rlim = 0.99;
	double rpdlim, rpdl, cs, ds, es, fs, gs, hs, ais, ajs, aks, z2, zd, zt;
	double spsr, stsC, umax, vmax;
	double wper, wperv, dphid, av, hwn, vll, vlls, dxdlat, dxdlon;
	double pi180, alpha, alphav, phidenz, dphiu, dphiv;
	double dsqrt2, polefac, phir, g0, g, re, r0, ri, ulphsrf, vlphsrf;
	double aw, bw, hg1, hg2, uvt2a;
	double pz1, rhoz1, tz1, uz1, vz1, wz1, tdz1, spz1, srhoz1, stz1,
		suz1, svz1, stdz1, rhn1, srhn1, vpn1, svpn1, spdavz1, spdsdz1,
		uvcorrz1, swh, spj, sdj, stj, suj, svj, swj, czi;
	double rrawt;
	double uvtsrf;
	double st2, sp2, sd2, su2, sv2, profwgt, sitenear = iperts.inits1.sitenear,
		sitelim = iperts.inits1.sitelim;
	double a, b, c, d, e, f;
	int iurra1 = iperts.inits1.iurra;
	string home1 = iperts.inits1.home, rra_id;

	string profile1 = iperts.inits1.profile;
	dsqrt2 = sqrt(2);
	pi = 3.1415926535897931;
	pi180 = pi / 180;
	twopi = 2 * pi;

	//call zinterp from InitPert class to calculate suface height from 
	//topographic map
	iperts.zinterp(phi, thet);

	//call rig to calculate earth radius at current height
	phir = phi * pi180;
	rig(h, phir, &g0, &g, &re, &r0, &ri);



	//Compute parameters for variable small scale lengths xlen, zlen
	iperts.intrw(iperts.inits1.xlbar, h, &xbarh);
	iperts.intrw(iperts.inits1.zlbar, h, &zbarh);
	iperts.intrw(iperts.inits1.xsigl, h, &sxlh);
	iperts.intrw(iperts.inits1.zsigl, h, &szlh);
	iperts.intrw(iperts.inits1.xlmin, h, &xmin);
	iperts.intrw(iperts.inits1.zlmin, h, &zmin);
	iperts.intrw(iperts.inits1.xscale, h, &xlen);
	iperts.intrw(iperts.inits1.zscale, h, &zlen);

	if (nm == 0){
		xlh = xbarh + (iperts.xl1 - iperts.xbar1) / iperts.sxl1;
		zlh = zbarh + (iperts.zl1 - iperts.zbar1) / iperts.szl1;
		goto RS;
	}

	//if using trajectory file data use deltas from trajectory calculation
	if (iperts.inits1.iopt > 0){
		dz = iperts.inits1.dz1;
		dphi = iperts.inits1.dphi1;
		dthet = iperts.inits1.dthet1;
		dtime = iperts.inits1.delt1;
		elt = iperts.inits1.time1;
		phi1 = phi - dphi;
		thet1 = thet - dthet;
	}

	//Use deltas from namelist file
	else {
		dz = iperts.inits1.dhgt;
		dtime = iperts.inits1.delt;
		dphi = iperts.inits1.dphi, dthet = iperts.inits1.dthet;
		phi1 = phi - dphi;
		thet1 = thet - dthet;
		elt = elt + dtime;
	}
	

	//horizontal great circle distance between successive positions
	dx = pi180*ri*radll(phi, thet, phi1, thet1);

	//If certain conditions are met new perturbations are not calculated
	if (dx < slim && dz == 0.0 && dtime == 0.0 || iupdate < 0 || nm == 0){
		prhl = rp1l;
		drhl = rd1l;
		trhl = rt1l;
		urhl = ru1l;
		vrhl = rv1l;
		prhs = rp1s;
		drhs = rd1s;
		trhs = rt1s;
		urhs = ru1s;
		vrhs = rv1s;
		prh = prhs + prhl;
		drh = drhs + drhl;
		trh = trhs + trhl;
		urh = urhs + urhl;
		vrh = vrhs + vrhl;
		wrh = iperts.rw1;
		xlh = xbarh + (iperts.xl1 - iperts.xbar1) / iperts.sxl1;
		zlh = zbarh + (iperts.zl1 - iperts.zbar1) / iperts.szl1;
		goto RS;
	}

	//Get 7 random numbers
	iperts.rcarry(rvec, 7);

	//Compute coefficients as, bs, alph, beta, gamm for small scale length model
	//Time scale (sec) for small-scale perturbations
	tpers = 86400.0*0.735*pow(abs(h), 0.116);
	tpers = max(10800.0, tpers);
	as = iperts.correl(abs(dx / xlen) + abs(dz / zlen) + abs(dtime / tpers));
	bs = sqrt(1.0 - as*as);

	//Compute small horizontal scale (xlh) from as and bs
	xlh = xbarh + (as*(iperts.xl1 - iperts.xbar1) / iperts.sxl1 +
		bs*iperts.ppnd(rvec[5], &l))*sxlh;

	//Evaluate Rcor = correlation between horizontal and vertical small scales
	if (h >= 200.0){
		rcor = 0.9;
	}
	else {
		rcor = 0.5 + 0.002*h;
	}
	denom = 1.0 - pow((rcor*as), 2);
	alph = as*(1.0 - pow(rcor, 2)) / denom;
	beta = rcor*(1.0 - pow(as, 2)) / denom;
	gamm = sqrt(abs(1.0 - pow(alph, 2) - pow(beta, 2) - 2.0*alph*beta*rcor*as));

	//Compute small vertical scale (zlh) from alph, beta and gamm
	zlh = zbarh + (alph*(iperts.zl1 - iperts.zbar1) / iperts.szl1 +
		beta*(xlh - xbarh) / sxlh + gamm*iperts.ppnd(rvec[6], &l))*szlh;

	//Compute fbig factor for severe perturbed standard deviations
	if (h <= 10.0) {
		fbig = 6.0 + 0.6*h;
	}
	else if (h >= 16.0){
		fbig = 6.0;
	}
	else{
		fbig = 22.0 - h;
	}

	//Compute probability for severe perturbed conditions
	xtail = (xbarh - xmin) / sxlh;
	ztail = (zbarh - zmin) / szlh;
	xtail = min(ztail, xtail);
	ptailx = iperts.ptail(xtail);


	//Compute factor for perturbed length scales
	scalfact = ptailx*xtail / (1.0 - ptailx);

	//Compute fsmall factor for normal perturbation conditions
	fsmall = (1.0 - fbig*ptailx) / (1.0 - ptailx);

	if (fsmall < 0.0005){
		fsmall = 0.0005;
		fbig = (0.9995 + 0.0005*ptailx) / ptailx;
	}

	//Compute scale sizes for normal perturbed conditions
	xmax = xbarh + scalfact*sxlh;

	zmax = zbarh + scalfact*szlh;

	//Small-scale multiplier for better compaison with KSC shears
	zfactor = 1.0;
	if (h <= 30.0) zfactor = exp(1.4*sin(twopi*h / 30.0));

	//If severe perturbed conditions, set scales to min and standard deviation to big
	if ((iperts.inits1.patchy > 0.0) && ((zlh < zmin) || (xlh < xmin))){
		vdsm = zmin*zfactor;
		vts = zmin*zfactor;
		vus = zmin*zfactor;
		hls = xmin*zfactor*2.0;
		densfact = sqrt(fbig);

	}

	//If normal perturbed conditions set scales to max and standard deviation to small
	else {
		vdsm = zmax*zfactor / 2.0;
		vts = zmax*zfactor / 2.0;
		vus = zmax*zfactor / 2.0;
		hls = xmax*zfactor*2.0;
		densfact = sqrt(fsmall);
	}

	//Compute relative separations (x/Lx and z/Lz) for small and large-scale models, and compute the correlations across the separation distance with the correl function
	//Include small-scale time effect with explicit time scale tpers
	hls = abs(dx / hls) + abs(dtime / tpers);
	rds = hls + abs(dz / vdsm);
	rds = iperts.correl(rds);
	rps = hls + abs(dz / vts);
	rps = iperts.correl(rps / plf);
	rvs = hls + abs(dz / vus);
	rvs = iperts.correl(rvs);

	//Fractional variances
	iperts.intr25(iperts.inits1.plp, iperts.inits1.dlp, h, phi, &plph, &dlph);
	iperts.intr25(iperts.inits1.ulp, iperts.inits1.vlp, h, phi, &ulph, &vlph);
	plph = min(0.99, plph);
	dlph = min(plmax, dlph);
	tlph = dlph;
	vlph = min(plmax, vlph);
	ulph = vlph;

	iperts.intr25(iperts.inits1.udl, iperts.inits1.uvt, h, phi, &udl2, &uvt2a);
	iperts.intr25(iperts.inits1.uds, iperts.inits1.vds, h, phi, &uds2, &vds2);

	//Calculate standard deviations for calculating perturbations
	//Call the gethgs member function from the NCEPmods class to get the height 
	//at the 10 mb (hg2) and 20 mb (hg1) level
	iperts.ncps.gethgs(phi, thet, &hg1, &hg2);
	//Call ncepmd from NCEPmods class to get NCEP data at surface height
	iperts.ncps.ncepmd(iperts.hsrf1, phi, thet, &pz1, &rhoz1, &tz1, &uz1, &vz1, &wz1, &tdz1, &spz1,
		&srhoz1, &stz1, &suz1, &svz1, &stdz1, &rhn1, &srhn1, &vpn1, &svpn1, &spdavz1,
		&spdsdz1, &uvcorrz1);

	iperts.usrf1 = uz1;
	iperts.vsrf1 = vz1;
	iperts.susrf1 = suz1;
	iperts.svsrf1 = svz1;
	iperts.tsrf1 = tz1;
	iperts.uvtsrf1 = uvcorrz1;
	iperts.stsrf1 = stz1;
	iperts.spdavsrf1 = spdavz1;
	iperts.spdsdsrf1 = spdsdz1;

	//Call ncepmd from NCEPmods class to get NCEP data at current height
	iperts.ncps.ncepmd(h, phi, thet, &pz1, &rhoz1, &tz1, &uz1, &vz1, &wz1, &tdz1, &spz1,
		&srhoz1, &stz1, &suz1, &svz1, &stdz1, &rhn1, &srhn1, &vpn1, &svpn1, &spdavz1,
		&spdsdz1, &uvcorrz1);

	sph = spz1;
	sdh = srhoz1;
	sth = stz1;
	suh = suz1;
	svh = svz1;

	if (h > hg1){
		//Call intrw and intruv from InitPert class to get interpolated standard 
		//deviations of u, v, w from atmosdat data
		iperts.intrw(iperts.inits1.wr, h, &swh);
		iperts.intruv(iperts.inits1.ur, iperts.inits1.vr, h, phi, &su, &sv);
		suh = sqrt(su);
		svh = sqrt(sv);
		//Call rterp from InitPert class to get interpolated standard deviations 
		//of p, d, t
		iperts.rterp(h, phi, &sp, &sd, &st);
		sph = sp;
		sdh = sd;
		sth = st;

		//Fair between NCEP and MAP data
		if (h < hg2){
			fair(2, hg1, spz1, srhoz1, stz1, suz1, svz1, swh, hg2, sph, sdh, sth, suh, svh, swh,
				h, &spj, &sdj, &stj, &suj, &svj, &swj, &czi);
			sph = spj;
			sdh = sdj;
			sth = stj;
			suh = suj;
			svh = svj;
			swh = swj;
		}

	}
	else {
		sph = spz1;
		sdh = srhoz1;
		sth = stz1;
		suh = suz1;
		svh = svz1;
	}

	sph = sqrt(abs(sph));
	sdh = sqrt(abs(sdh));
	sth = sqrt(abs(sth));



	if (h <= hg2){
		uvt2 = uvcorrz1;
	}
	else {
		uvt2 = uvt2a;
	}

	iperts.stsrf1 = iperts.tsrf1*sqrt(iperts.stsrf1);

	//If conditions are met use RRA data
	if ((h <= 70.0) & (iurra1 > 0)){
		r0 = 8314.427;
		sph = sph;
		sdh = sdh;
		sth = sth;
		suh = suh;
		svh = svh;

		//Call rrasigs member function from RRA class to get standard deviations
		iperts.rras.rrasigs(h, phi, thet, iperts.inits1.mn, iperts.hsrf1, iperts.usrf1, iperts.vsrf1,
			iperts.tsrf1, iperts.uvtsrf1, iperts.susrf1, iperts.svsrf1, iperts.stsrf1,
			iperts.spdavsrf1, iperts.spdsdsrf1, r0, uvt2, &sph, &sdh, &sth, &suh, &svh,
			rra_id, &rrawt);

		if (iperts.rras.z1[iperts.rras.num1 - 1] && iperts.rras.rrawt1 > 0.0){
			suh = suh;
			svh = svh;
			uvt2 = iperts.rras.uvrra;
			if (rrawt > 0.0){
				iperts.hsrf1 = iperts.rras.rrsrf;
				iperts.usrf1 = iperts.rras.usrf;
				iperts.vsrf1 = iperts.rras.vsrf;
				iperts.susrf1 = iperts.rras.susrf;
				iperts.svsrf1 = iperts.rras.svsrf;
				iperts.uvtsrf1 = iperts.rras.uvtsrf;
				iperts.spdavsrf1 = iperts.rras.spdavsrf;
				iperts.spdsdsrf1 = iperts.rras.spdsdsrf;
			}
		}

		if (iperts.rras.z2[iperts.rras.num2 - 1] && iperts.rras.rrawt2 > 0.0){
			sph = sph;
			sdh = sdh;
			sth = sth;
			if (rrawt > 0.0){
				iperts.tsrf1 = iperts.rras.tsrf;
				iperts.stsrf1 = iperts.rras.stsrf;
			}

		}
	}

	//Use auxiliary profile data if conditions are met
	if (iperts.inits1.iaux > 0){
		//Call member functions rdprof and profsigs from AuxProf Class to read and
		//calculate standard deviations.
		iperts.auxs.profsigs(h, phi, thet, sth, sph, sdh, suh, svh, &st2, &sp2, &sd2,
			&su2, &sv2, &profwgt);

		sth = st2;
		sph = sp2;
		sdh = sd2;
		suh = su2;
		svh = sv2;

	}
	//Apply rpscale and ruscale to standard deviations
	sph = sph*iperts.inits1.rpscale;
	sdh = sdh*iperts.inits1.rpscale;
	sth = sth*iperts.inits1.rpscale;
	suh = suh*iperts.inits1.ruscale;
	svh = svh*iperts.inits1.ruscale;

	//Check for buell constraint violations
	if (sdh >= rlim*(sph + sdh)){
		sdh = rlim*(sph + sth);
	}
	else if (sth >= rlim*(sph + sdh)){
		sth = rlim*(sph + sdh);
	}
	else if (sph >= rlim*(sdh + sth)){
		sph = rlim*(sdh + sth);
	}

	//Compute large scale and small scale sigmas from total sigmas
	sphl = sqrt(plph)*sph;
	sphs = sqrt(1.0 - plph)*sph;
	sdhl = sqrt(dlph)*sdh;
	sdhs = sqrt(1.0 - dlph)*sdh;
	sthl = sqrt(tlph)*sth;
	sths = sqrt(1.0 - tlph)*sth;
	suhl = sqrt(ulph)*suh;
	suhs = sqrt(1.0 - ulph)*suh;
	svhl = sqrt(vlph)*svh;
	svhs = sqrt(1.0 - vlph)*svh;

	//Compute pressure-density correlation from the Buell relation
	rpd = (sph*sph + sdh*sdh - sth*sth) / (2.0*sph*sdh);

	//Get upper limit on rdt (which should be negative)
	if (h < 10.0){
		rdtlim = -0.9 + 0.05*h;
	}
	else if (h < 15.0){
		rdtlim = -0.4 - 0.08*(h - 10.0);
	}
	else if (h < 30.0){
		rdtlim = -0.8 + 0.04*(h - 15.0);
	}
	else {
		rdtlim = -0.2;
	}

	//Get upper limit on rpd from rdt limit
	rpdlim = (sdh + rdtlim*sth) / sph;
	if (abs(rpdlim) > 0.99) rpdlim = copysign(0.99, rpdlim);
	if (rpd > rpdlim) rpd = rpdlim;
	//Assume rpdl=rpds for perturbation p-d correlation
	rpdl = rpd*sph*sdh / (sphl*sdhl + sphs*sdhs);
	if (abs(rpdl) > 0.99)rpdl = copysign(0.99, rpdl);

	//Evaluate the perturbation model coefficients (a,b..aj,ak) for the small scale model
	coeff(sd1s, sdhs, sp1s, sphs, st1s, sths, su1s, suhs, sv1s, svhs, uds1, uds2, vds1, vds2, rpdl, rps, rds, rvs, &as, &bs, &cs, &ds, &es, &fs, &gs, &hs, &ais, &ajs, &aks);

	z2 = rvec[0];
	zd = iperts.ppnd(z2, &l);
	z2 = rvec[1];
	zt = iperts.ppnd(z2, &l);

	//Compute new small-scale perturbations for density (drhs), temperature (trhs) and pressure (prhs)
	//Set densfact to 1 if patchiness is off
	if (iperts.inits1.patchy == 0.0) densfact = 1.0;

	z2 = (as - 1.0)*rd1s + bs*zd*densfact;

	//Assure change in density during severe perturbation is not large compared to normal perturbation standard deviation
	if ((iperts.inits1.patchy != 0.0) & (abs(z2) > 3.0*sdhs)) z2 = copysign(3.0*sdhs, z2);
	drhs = rd1s + z2;

	z2 = (cs - 1.0)*rp1s + ds*drhs + es*zt*densfact;

	//Assure change in pressure during severe perturbations is not large compared to normal perturbation standard deviation
	if ((iperts.inits1.patchy != 0.0) & (abs(z2) > 3.0*sphs)) z2 = copysign(3.0*sphs, z2);
	prhs = rp1s + z2;
	spsr = abs(prhs / sphs);
	if (spsr > 4.0) prhs = copysign(abs(prhs) - 0.1*sphs, prhs);
	trhs = prhs - drhs;
	//Adjust temperature perturbation assuming rpdl p-d correlation
	stsC = sqrt(abs(sphs*sphs + sdhs*sdhs - 2.0*rpdl*sphs*sdhs));
	trhs = trhs*sths / stsC;

	z2 = rvec[2];
	zd = iperts.ppnd(z2, &l);
	z2 = rvec[3];
	zt = iperts.ppnd(z2, &l);
	//Compute new small-scale wind perturbations (urhs and vrhs)
	//Apply scalfact for normal or severe perturbed conditions
	scalfact = densfact;
	//Assure wind perturbation standard deviations are not large compared with severe turbulence sigmas from Figure 2, NASA TM 4168 (1990)
	if (iperts.inits1.patchy != 0.0){
		scalfact = min((4.0 + 0.7*h) / suhs, scalfact);
		scalfact = min((4.0 + 0.7*h) / svhs, scalfact);
		scalfact = max(densfact, scalfact);
	}
	z2 = (fs - 1.0)*ru1s + gs*drhs + hs*zd*scalfact;
	//Assure change in E-W wind during severe perturbations is not large compared to normal perturbation standard deviation
	if ((iperts.inits1.patchy != 0.0) & (abs(z2) > 3.0*suhs)) z2 = copysign(3.0*suhs, z2);
	urhs = ru1s + z2;

	z2 = (ais - 1.0)*rv1s + ajs*drhs + aks*zt*scalfact;
	//Assure change in N-S wind during severe perturbations is not large compared to normal perturbation standard deviation
	if ((iperts.inits1.patchy != 0.0) & (abs(z2) > 3.0*svhs)) z2 = copysign(3.0*svhs, z2);
	vrhs = rv1s + z2;

	a = 4.0*suhs*scalfact;
	b = 6.0*suhs;
	c = 200.0 + h;
	d = 400.0;

	e = min(a, b);
	f = min(c, d);

	umax = min(e, f);

	a = 4.0*svhs*scalfact;
	b = 6.0*svhs;
	c = 200.0 + h;
	d = 400.0;

	e = min(a, b);
	f = min(c, d);

	vmax = min(e, f);
	if (abs(urhs) > umax) urhs = copysign(abs(urhs) - 0.1*umax, urhs);
	if (abs(vrhs) > vmax) vrhs = copysign(abs(vrhs) - 0.1*vmax, vrhs);

	//Period and wavelengt parameters for large-scale perturbations
	wper = 2.0 + 60.0*iperts.waverand;
	wperv = 2.0 + 40.0*pow(iperts.waverand, 2);
	dphid = iperts.ppnd(iperts.waverand, &l);
	if (abs(dphid) > 3.0)dphid = copysign(3.0, dphid);
	av = 25.0 + 3.0*dphid;
	hwn = double(int(round(4.0 + 0.833*dphid)));

	//Compute vertical scale vll
	vll = av + bv*sqrt(abs(pow(h, 3)));
	vll = min(200.0, vll);
	vlls = av + bv*sqrt(abs(pow(iperts.hsrf1, 3)));
	//Horizontal wavenumbers for waver perturbations
	dxdlat = hwn*pi180;
	dxdlon = hwn*pi180;

	//Compute new large-scale perturbations for density (drhl), temperature (drhl) and pressure (prhl)
	alpha = dxdlon*thet + dxdlat*phi;
	alphav = alpha + twopi*elt / (wperv*86400.0);
	alpha = alpha + twopi*elt / (wper*86400.0);
	phidenz = iperts.phidens + twopi*h / vll;
	dphiu = acos(rpdl);
	drhl = dsqrt2*sdhl*cos(alpha + phidenz)*iperts.ampfact;
	prhl = dsqrt2*sphl*cos(alpha + phidenz + dphiu)*iperts.ampfact;
	trhl = prhl - drhl;

	//Adjust temperature perturbation assuming rpdl p-d correlation
	stsC = sqrt(abs(sphl*sphl + sdhl*sdhl - 2 * rpdl*sphl*sdhl));
	trhl = trhl*sthl / stsC;
	//Get phase angles for large-scale u and v components
	dphiu = acos(udl2);
	dphiv = uvt2 / (0.02 + 0.98*pow((suhl / suh), 2));
	if (abs(dphiv) > 1.0)dphiv = copysign(1.0, dphiv);
	dphiv = dphiu + copysign(acos(dphiv), iperts.waverand - 0.5);

	//Compute new large-scale wind perturbations (urhl and vrhl)
	urhl = dsqrt2*suhl*cos(alpha + phidenz + dphiu)*iperts.ampfact;
	vrhl = dsqrt2*svhl*cos(alphav + phidenz + dphiv)*iperts.ampfact;

	a = 4.0*suhl*scalfact;
	b = 6.0*suhl;
	c = 200.0 + h;
	d = 400.0;

	e = min(a, b);
	f = min(c, d);

	umax = min(e, f);

	a = 4.0*svhl*scalfact;
	b = 6.0*svhl;
	c = 200.0 + h;
	d = 400.0;

	e = min(a, b);
	f = min(c, d);

	vmax = min(e, f);
	if (abs(urhl) > umax) urhl = copysign(umax, urhl);
	if (abs(vrhl) > vmax) vrhl = copysign(vmax, vrhl);

	//Insure large-scale perturbations approach zero at the poles
	if (abs(phi) >= 85.0){
		polefac = cos(pi*(abs(phi) - 85.0) / 10.0);
		drhl = polefac*drhl;
		prhl = polefac*prhl;
		trhl = polefac*trhl;
		urhl = polefac*urhl;
		vrhl = polefac*vrhl;
	}

	//Calculate vertical wind perturbation
	if (h <= 10.0){
		//Get phase angles for large-scale U and V components at surface
		iperts.intr25(iperts.inits1.udl, iperts.inits1.uvt, iperts.hsrf1, phi, &iperts.udlsrf, &uvtsrf);
		dphiu = acos(iperts.udlsrf);
		iperts.intr25(iperts.inits1.ulp, iperts.inits1.vlp, iperts.hsrf1, phi, &ulphsrf, &vlphsrf);
		dphiv = iperts.uvtsrf1 / (0.02 + 0.98*vlphsrf);
		if (abs(dphiv) > 1.0)dphiv = copysign(1.0, dphiv);;
		dphiv = dphiu + copysign(acos(dphiv), iperts.waverand - 0.5);

	}

	//Calculate sigma-w from boundary layer (BL) model
	z2 = rvec[4];
	zd = iperts.ppnd(z2, &l);
	iperts.getsigw(h, phi, thet, alpha, alphav, vlls, dphiu, dphiv, elt);
	if (iperts.sw1 <= 0.0)iperts.sw1 = slim;
	if (iperts.swh <= 0.0)iperts.swh = slim;

	//Compute coefficients needed for vertical wind perturbations
	aw = rvs*iperts.swh / iperts.sw1;
	bw = iperts.swh*sqrt(1.0 - rvs*rvs);

	//Compute new vertical wind (small-scale) perturbation (wrh)
	z2 = (aw - 1.0)*iperts.rw1 + bw*zd*densfact;

	if (abs(z2) > 3.0*iperts.swh) z2 = copysign(3.0*iperts.swh, z2);
	wrh = iperts.rw1 + z2;

	//Total (small + large-scale) perturbations
	drh = drhs + drhl;
	trh = trhs + trhl;

	//Use 2nd order law for total and large-scale pressure perturbations
	prh = drh + trh + drh*trh;
	prhl = prh - prhs;
	//Avoid perturbation less than -90% of mean value
	prh = max(-0.9, prh);
	drh = max(-0.9, drh);
	trh = max(-0.9, trh);
	//Total (small + large-scale) wind perturbations
	urh = urhs + urhl;
	vrh = vrhs + vrhl;

	a = 4.0*suh*scalfact;
	b = 6.0*suh;
	c = 200.0 + h;
	d = 400.0;

	e = min(a, b);
	f = min(a, b);

	umax = min(e, f);

	a = 4.0*svh*scalfact;
	b = 6.0*svh;
	c = 200.0 + h;
	d = 400.0;

	e = min(a, b);
	f = min(c, d);

	vmax = min(e, f);

	if (abs(urh) > umax){
		urh = copysign(abs(urh) - 0.1*umax, urh);
		urhs = urh - urhl;
	}
	if (abs(vrh) > vmax){
		vrh = copysign(abs(vrh) - 0.1*vmax, vrh);
		vrhs = vrh - vrhl;
	}


	//Reset density-velocity corrlaltions
	uds1 = uds2;
	vds1 = vds2;

	//Set previous random perturbation in p,d,t to current perturbations,
	//for next cycle
	rp1s = prhs;
	rd1s = drhs;
	rt1s = trhs;
	rp1l = prhl;
	rd1l = drhl;
	rt1l = trhl;

	//Set previous wind perturbation values to current values, for next
	//cycle.
	ru1s = urhs;
	rv1s = vrhs;
	ru1l = urhl;
	rv1l = vrhl;
	iperts.rw1 = wrh;


	//Reset standard deviations
	sp1s = sphs;
	sd1s = sdhs;
	st1s = sths;
	su1s = suhs;
	sv1s = svhs;
	sp1l = sphl;
	sd1l = sdhl;
	st1l = sthl;
	su1l = suhl;
	sv1l = svhl;
	iperts.sw1 = iperts.swh;
	iperts.sd1 = sdh;
	iperts.sp1 = sph;
	iperts.st1 = sth;
	iperts.su1 = suh;
	iperts.sv1 = svh;


	if (densfact > 1.0){
		isev = 1;
	}
	else {
		isev = 0;
	}

RS:

	if (iupdate < 0){
		sph = iperts.sp1;
		sdh = iperts.sd1;
		sth = iperts.st1;
		suh = iperts.su1;
		svh = iperts.sv1;
	}

	//Set previous scale values and magnitudes to current values
	iperts.xl1 = xlh;
	iperts.sxl1 = sxlh;
	iperts.xbar1 = xbarh;
	iperts.zl1 = zlh;
	iperts.szl1 = szlh;
	iperts.zbar1 = zbarh;

	return;


}

void Pert::fair(int i, double hg, double pg, double dg, double tg, double ug, double vg, double wg,
	double hj, double pj, double dj, double tj, double uj, double vj, double wj, double h,
	double *p, double *d, double *t, double *u, double *v, double *w, double *czi)
{

	//fair member function from Pert class
	//Fairs at height h, between pg,dg,tg,ug,vg,wg and pj,dj,tj,uj,vj,wj.  Fairing is cosine**2,
	//such that all "g" values obtain at height hg and all "j" values obtain at height hj.  The
	//faired results are p,d,t,u,v,w.  Input i = 1 means density (d) is faired and pressure (p)
	//is computed from a faired value of the gas-law constant.  If i not equal to 1, both p and
	//d are faired directly.

	double r, rg, rj, szi, pi = 4 * atan(1.0), piby2 = pi / 2.0;

	//Fairs between values
	*czi = pow((cos(piby2*(h - hg) / (hj - hg))), 2);

	//Complement of fairing coefficient
	szi = 1.0 - *czi;

	//Faired temperature
	*t = tg**czi + tj*szi;

	if (i == 1){
		//Faired density
		*d = exp(log(dg)**czi + log(dj)*szi);

		//Faired gas constant and pressure
		rg = pg / (dg*tg);
		rj = pj / (dj*tj);
		r = *czi*rg + szi*rj;
		*p = r**d**t;
	}

	else {
		*d = dg**czi + dj*szi;
		*p = pg**czi + pj*szi;
	}

	//Faired wind components
	*u = ug**czi + uj*szi;
	*v = vg**czi + vj*szi;
	*w = wg**czi + wj*szi;

	return;

}

void Pert::coeff(double sd1, double sd2, double sp1, double sp2, double st1, double st2, double su1, double su2, double sv1, double sv2, double ud1, double ud2, double vd1, double vd2, double rpd, double rp, double rd, double rv, double *a, double *b, double *c, double *d, double *e, double *f, double *g, double *h, double *ai, double *aj, double *ak)
{

	//coeff member function from Pert class
	//Computes perturbation model coefficients (a,b...aj,ak) at previous 
	//positions (1) and current position (2) from the  

	double slim = 1.0e-5;
	double rlim = 0.99999;

	//Default values avoid division by zero
	if (sd1 <= 0.0)sd1 = slim;
	if (st1 <= 0.0)st1 = slim;
	if (sd2 <= 0.0)sd2 = slim;
	if (st2 <= 0.0)st2 = slim;
	if (rp <= 0.0)rp = slim;
	if (rd <= 0.0)rd = slim;
	if (rv <= 0.0)rv = slim;
	if (abs(ud1) <= 0.0) ud1 = slim;
	if (abs(vd1) <= 0.0) vd1 = slim;
	if (abs(su1) <= 0.0) su1 = slim;
	if (abs(sv1) <= 0.0) sv1 = slim;
	if (abs(ud1) >= rlim)ud1 = rlim*ud1 / abs(ud1);
	if (abs(vd1) >= rlim)vd1 = rlim*vd1 / abs(vd1);
	if (abs(ud2) <= 0.0)ud2 = slim;
	if (abs(vd2) <= 0.0)vd2 = slim;
	if (abs(su2) <= 0.0)su2 = slim;
	if (abs(sv2) <= 0.0)sv2 = slim;
	if (abs(ud2) >= rlim)ud2 = rlim*ud2 / abs(ud2);
	if (abs(vd2) >= rlim)vd2 = rlim*vd2 / abs(vd2);
	if (abs(rp) >= rlim)rp = rlim*rp / abs(rp);
	if (abs(rd) >= rlim)rd = rlim*rd / abs(rd);
	if (abs(rv) >= rlim)rv = rlim*rv / abs(rv);

	//Compute perturbation model coeffcients
	*a = rd*sd2 / sd1;
	*b = sd2*sqrt(1 - rd*rd);
	if (abs(rpd) <= 0.0) rpd = slim;
	if (abs(rpd) >= rlim) rpd = rlim*rpd / abs(rpd);
	*c = (sp2 / sp1) * (rp - (rd*rpd*rpd)) / (1.0 - (rpd*rpd*rd*rd));
	*d = rpd*(sp2 - (*c*rd*sp1)) / sd2;
	*e = sp2*sp2 - *c**c*sp1*sp1 - *d**d*sd2*sd2 - 2 * *c**d*rd*rpd*sp1*sd2;
	*e = max(0.0, *e);
	*e = sqrt(*e);
	*f = (su2 / su1)*(rv - rd*ud2*ud1) / (1.0 - rd*rd*ud1*ud1);
	*g = (ud2*su2 - *f*rd*ud1*su1) / sd2;
	*h = su2*su2 - *f**f*su1*su1 - *g**g*sd2*sd2 - 2.0**f**g*rd*ud1*sd2*su1;
	*h = max(0.0, *h);
	*h = sqrt(*h);
	*ai = (sv2 / sv1)*(rv - rd*vd2*vd1) / (1.0 - rd*rd*vd1*vd1);
	*aj = (vd2*sv2 - *ai*rd*vd1*sv1) / sd2;
	*ak = sv2*sv2 - *ai**ai*sv1*sv1 - *aj**aj*sd2*sd2 - 2.0**ai**aj*rd*vd1*sd2*sv1;
	*ak = max(0.0, *ak);
	*ak = sqrt(*ak);


}


double Pert::radll(double phi1, double thet1, double phi2, double thet2)
{

	//radll member function from Pert class
	//Returns great-circle distance (degrees) between two input lat-lon positions 
	//(in degrees)

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

void Pert::rig(double ch, double phir, double *g0, double *g, double *re, double *r0, double *ri)
{

	//rig member function from Pert class
	//Computes surface gravity g0 (m/s^2) and effective Earth radius re (km) from 
	//input geocentric latitude phir (radians).  Also computes gravity g (m/s^2)
	//and total radius ri (km) at input height ch (km).

	double g0a = 9.80616, g0b = 0.0026373, g0c = 0.0000059, rea = 3.085462e-3, reb = 2.27e-6, rec = 2.0e-9;
	double eps, c2phi, cphi2, c4phi;
	const double a = 6378.137, b = 6356.752314;

	//Parameters for computing effective Earth radius

	//Set up cosines:  cphi2 = [cos(phir)]**2, ...
	cphi2 = pow(cos(phir), 2);
	c2phi = 2.0*cphi2 - 1.0;
	c4phi = 8.0*cphi2*(cphi2 - 1.0) + 1.0;


	eps = 1.0 - pow(b / a, 2);

	//Compute Earth radius, r0
	*r0 = b / sqrt(1.0 - eps*cphi2);

	//Compute surface gravity
	*g0 = g0a*(1.0 - g0b*c2phi + g0c*pow(c2phi, 2));

	//Compute effective Earth radius, re
	*re = 2.0**g0 / (rea + reb*c2phi - rec*c4phi);

	//Compute g at height ch
	*g = *g0 / pow((1.0 + ch / (*re)), 2);

	//Compute radius at height ch
	*ri = *r0 + ch;

}

