//P. White
//Earth-GRAM 2016, Atmod Class
//The atmod class computes mean and standard deviation pressure, density and 
//winds for output.  As well as perturbed values and atmospheric constituents.


#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <string>
#include <cmath>
#include <stdio.h>
#include "Atmod1.h"
#include "iomanip"
#include <cstdlib>

using namespace std;


Atm1::Atm1()
{
	initializeMemberVariables();

}

Atm1::~Atm1()
{

}

void Atm1::initializeMemberVariables()
{

	//Initialize the member variables for the atmod class.

	h = 0.0;
	hj1 = 90.0;
	hj2 = 120.0;
	phi = 0.0;
	thet = 0.0;
	dz = 0.0;
	dphi = 0.0;
	dthet = 0.0;
	delt = 0.0;
	elt = 0.0;
	pi = 4 * atan(1.0);
	pi180 = pi / 180;
	yr = 2010;
	mn = 0;
	piby2 = pi / 2;
	pmean = 0.0;
	dmean = 0.0;
	tmean = 0.0;
	umean = 0.0;
	vmean = 0.0;
	wmean = 0.0;
	pmean1 = 0.0;
	dmean1 = 0.0;
	tmean1 = 0.0;
	umean1 = 0.0;
	vmean1 = 0.0;
	wmean1 = 0.0;
	pmean2 = 0.0;
	dmean2 = 0.0;
	tmean2 = 0.0;
	umean2 = 0.0;
	vmean2 = 0.0;
	wmean2 = 0.0;
	ppert = 0.0;
	dpert = 0.0;
	tpert = 0.0;
	upert = 0.0;
	vpert = 0.0;
	wpert = 0.0;
	pmpert = 0.0;
	dmpert = 0.0;
	tmpert = 0.0;
	umpert = 0.0;
	vmpert = 0.0;
	wmpert = 0.0;
	nr1 = 0;
	pstd = 0.0;
	dstd = 0.0;
	tstd = 0.0;
	ustd = 0.0;
	vstd = 0.0;
	wstd = 0.0;
	prramean = 0.0;
	drramean = 0.0;
	trramean = 0.0;
	urramean = 0.0;
	vrramean = 0.0;
	prrastd = 0.0;
	drrastd = 0.0;
	trrastd = 0.0;
	urrastd = 0.0;
	vrrastd = 0.0;
	spg = 0.0, sdg = 0.0, stg = 0.0, sug = 0.0, svg = 0.0;
	nm = 0;
	pnmean = 0.0, dnmean = 0.0, tnmean = 0.0, unmean = 0.0, vnmean = 0.0,
		wnmean = 0.0, Umean = 0.0, Vmean = 0.0;
	vpn = 0.0, rhn = 0.0, tdmean, stdd = 0.0, svpn = 0.0, srhn = 0.0;
	arnd = 0.0, hend = 0.0, hnd = 0.0, nnd = 0.0, o2nd = 0.0, n2nd = 0.0, wtmol = 0.0,
		waterg = 0.0, ond = 0.0, pstd1 = 0.0, dstd1 = 0.0, tstd1 = 0.0, ustd1 = 0.0, vstd1 = 0.0,
		wstd1 = 0.0, airmw;
	ppert1 = 0.0, dpert1 = 0.0, tpert1 = 0.0, upert1 = 0.0, vpert1 = 0.0, wpert1 = 0.0;
	wtair = 0.0, wth2o = 0.0, wto3 = 0.0, wtn2o = 0.0, wtco = 0.0, wtch4 = 0.0, wtco2 = 0.0,
		wtn2 = 0.0, wto2 = 0.0, wto = 0.0, wtar = 0.0, wthe = 0.0, wth = 0.0, wtn = 0.0;
	hgt1 = 0.0, time1 = 0.0, lat1 = 0.0, lon1 = 0.0, iupdate1 = 0;
	dtz1 = 0.0, dtzm = 0.0, dmdz1 = 0.0, dmdzm = 0.0;
}


void Atm1::initdata()
{

	//The member function initdata, of the Atmod Class, that initializes the data
	//needed to compute the atmospheric variables for output.

	string namefile;

	//Output statement displayed when program is executed.
	cout << "Earth-GRAM 2016" << '\n';

	//Call the namelist function from the Init Class for user-defined inputs.
	msisa.maps.perts.iperts.inits1.namelist();

	//Variable namefile for saving namefile directory and file name 
	namefile = msisa.maps.perts.iperts.inits1.home +
		msisa.maps.perts.iperts.inits1.namelst;

	//Set the variables for number seed (nr1) and month (mn) from the namelist file.
	nr1 = msisa.maps.perts.iperts.inits1.nr1;
	mn = msisa.maps.perts.iperts.inits1.mn;

	//Output statement when GRAM initialization begins.
	cout << "GRAM Initialization" << '\n';

	//Call the init1 member function from the Init class to initialize the 
	//atmosdat data into arrays.
	msisa.maps.perts.iperts.inits1.init1(mn);

	//Initialize NCEP data
	msisa.maps.perts.iperts.ncps.namelist(namefile);
	msisa.maps.perts.iperts.ncps.NCEPread();

	//Initialize Auxiliary Profile data
	if (msisa.maps.perts.iperts.inits1.iaux > 0){
		msisa.maps.perts.iperts.auxs.rdprof(msisa.maps.perts.iperts.inits1.sitenear,
			msisa.maps.perts.iperts.inits1.sitelim, msisa.maps.perts.iperts.inits1.home,
			msisa.maps.perts.iperts.inits1.profile);
	}

	//Initialize RRA data
	msisa.maps.perts.iperts.rras.init(namefile);

	//Initialize Thermosphere models
	msisa.maps.perts.iperts.mets.namelist(namefile);
	msisa.jb2008.namelist(namefile);
	msisa.namelist(namefile);

	//File open for the "special.txt" and "species.txt" output files.
	string nppath = msisa.maps.perts.iperts.inits1.home + msisa.maps.perts.iperts.inits1.nprpath;
	string conpath = msisa.maps.perts.iperts.inits1.home + msisa.maps.perts.iperts.inits1.conpath;

	output.open(nppath.c_str());
	species.open(conpath.c_str());

	//Headers for the "species.txt" and "special.txt" output files.
	species << "   Height   GcLat    Long.     Conc.     #Dens.   |    Conc.     #Dens.     Species" << '\n';
	species << "    (km)    (deg)    (deg)    (ppmv)    (#/m^3)   |    (ppmv)   (#/m^3)" << '\n';
	output << "      Time    Hgtkm GeocenLat  Lon(East)   DensMean   PresMean    Tmean "
		"  EWmean"
		"   NSmean   DensPert   PresPert    Tpert   EWpert  NSpert SDden% SDprs% "
		"SDtemK SDuwnd SDvwnd SDwwnd Wpert SpdAvg  SpdStd SOSmean SOSpert Sev" << '\n';

}

void Atm1::initpert()
{

	//The member function initpert from the Atmod class for computing
	//initial perturbations.

	//Call to member function initpert of InitP Class to calculate initial perturbations.
	msisa.maps.perts.iperts.initpert(h, phi, thet, elt, nr1);

	//Atmod member variables for initial perturbations and standard deviations for
	//pressure, density, temperature, east-west wind, north-south wind, and vertical
	//wind.
	ppert1 = msisa.maps.perts.iperts.prh1;
	dpert1 = msisa.maps.perts.iperts.drh1;
	tpert1 = msisa.maps.perts.iperts.trh1;
	upert1 = msisa.maps.perts.iperts.urh1;
	vpert1 = msisa.maps.perts.iperts.vrh1;
	wpert1 = msisa.maps.perts.iperts.rw1;
	pstd1 = msisa.maps.perts.iperts.sp1;
	dstd1 = msisa.maps.perts.iperts.sd1;
	tstd1 = msisa.maps.perts.iperts.st1;
	ustd1 = msisa.maps.perts.iperts.su1;
	vstd1 = msisa.maps.perts.iperts.sv1;
	wstd1 = msisa.maps.perts.iperts.swh;

}

void Atm1::map()
{

	//The member function map from the Atmod class used to calculate the MAP data.

	double pmm, dmm, tmm, umm, vmm, wmm, dtz;
	double g0, g, re, r0, ri1;
	double pi180 = pi / 180;
	double phir = phi * pi180;

	//The member function rig from the Pert Class used to compute surface gravity
	//and effective Earth radius.
	msisa.maps.perts.rig(h, phir, &g0, &g, &re, &r0, &ri1);

	//The member function mapmod from the Map class used to compute mean pressure,
	//density, temperture, east-west wind, north-south wind, and vertical wind
	//from the middle atmosphere program.
	msisa.maps.mapmod(h, phi, thet, ri1, g, &pmm, &dmm, &tmm, &umm, &vmm, &wmm, &dtz);

	pmean1 = pmm;
	dmean1 = dmm;
	tmean1 = tmm;
	umean1 = umm;
	vmean1 = vmm;
	wmean1 = wmm;
	dtzm = dtz;

}


void Atm1::pert()
{

	//The member function pert from the Atmod class used to calculate perturbations.
	//The member function pert1 from the Pert class used to calculate perturbations
	//for pressure, density, temperature, east-west wind, north-south wind, and
	//vertical wind.
	msisa.maps.perts.pert1(h, phi, thet, iupdate1, nm);

	ppert = msisa.maps.perts.prh;
	dpert = msisa.maps.perts.drh;
	tpert = msisa.maps.perts.trh;
	upert = msisa.maps.perts.urh;
	vpert = msisa.maps.perts.vrh;
	wpert = msisa.maps.perts.wrh;
	pstd = msisa.maps.perts.sph;
	dstd = msisa.maps.perts.sdh;
	tstd = msisa.maps.perts.sth;
	ustd = msisa.maps.perts.suh;
	vstd = msisa.maps.perts.svh;
	wstd = msisa.maps.perts.iperts.swh;


}


void Atm1::met()
{

	//The member function met from the Atmod class used to calculate atmospheric
	//variables with the MET model.

	double pmm, dmm, tmm, umm, vmm, wmm, arnd1, hend1, hnd1, o2nd1, n2nd1, ond1, wtmol1;

	//The member function jacmod from the MET class used to calculate mean pressure,
	//density, temperature, east-west wind, north-south wind, and vertical wind.  
	//Atmospheric constituents ar, he, h, o2, n2, and o are also calculated.
	msisa.maps.perts.iperts.mets.jacmod(h, phi, thet, elt, &pmm, &dmm, &tmm, &umm, &vmm, &wmm, &arnd1, &hend1, &hnd1, &o2nd1, &n2nd1, &ond1, &wtmol1);

	pmean2 = pmm;
	dmean2 = dmm;
	tmean2 = tmm;
	umean2 = umm;
	vmean2 = vmm;
	wmean2 = wmm;

	arnd = arnd1;
	hend = hend1;
	hnd = hnd1;
	o2nd = o2nd1;
	n2nd = n2nd1;
	ond = ond1;
	wtmol = wtmol1;


}

void Atm1::jb2008()
{

	//The member function jb2008 from the Atmod class used to calculate atmospheric
	//variables with the JB2008 model.

	double pmm, dmm, tmm, umm, vmm, wmm, n2nd1, o2nd1, ond1, arnd1, hend1,
		hnd1, wtmol1, dmdz2, nnd1;
	double pi = 3.1415926535897931, g0, g, re, r0, ri1;
	double pi180 = pi / 180;
	double phir = phi * pi180;


	//The member function rig from the Pert Class used to compute surface gravity
	//and effective Earth radius.
	msisa.maps.perts.rig(h, phir, &g0, &g, &re, &r0, &ri1);

	//The member function JB08mod from the JB2008 class used to calculate mean pressure,
	//density, temperature, east-west wind, north-south wind, and vertical wind.  
	//Atmospheric constituents ar, he, h, o2, n2, and o are also calculated.
	msisa.jb2008.JB08mod(h, phi, thet, elt, ri1, g, &pmm, &dmm, &tmm, &umm,
		&vmm, &wmm, &n2nd1, &o2nd1, &ond1, &arnd1, &hend1, &hnd1, &wtmol1,
		&dmdz2, &nnd1);

	pmean2 = pmm;
	dmean2 = dmm;
	tmean2 = tmm;
	umean2 = umm;
	vmean2 = vmm;
	wmean2 = wmm;
	arnd = arnd1;
	hend = hend1;
	hnd = hnd1;
	o2nd = o2nd1;
	n2nd = n2nd1;
	ond = ond1;
	nnd = nnd1;
	wtmol = wtmol1;
	dmdz1 = dmdz2;

}


void Atm1::msis()
{

	//The member function msis from the Atmod class used to calculate atmospheric
	//variables with the MSIS model.

	double pmm, dmm, tmm, umm, vmm, wmm, n2nd1, o2nd1, ond1, arnd1, hend1,
		hnd1, wtmol1, dmdz2, nnd1;
	double pi = 3.1415926535897931, g0, g, re, r0, ri1;
	double pi180 = pi / 180;
	double phir = phi * pi180;

	//The member function rig from the Pert Class used to compute surface gravity
	//and effective Earth radius.
	msisa.maps.perts.rig(h, phir, &g0, &g, &re, &r0, &ri1);

	//The member function msismod from the MSIS class used to calculate mean pressure,
	//density, temperature, east-west wind, north-south wind, and vertical wind.  
	//Atmospheric constituents ar, he, h, o2, n2, and o are also calculated.
	msisa.msismod(h, phi, thet, elt, ri1, g, &pmm, &dmm, &tmm, &umm, &vmm,
		&wmm, &n2nd1, &o2nd1, &ond1, &arnd1, &hend1, &hnd1, &wtmol1, &dmdz2,
		&nnd1);

	pmean2 = pmm;
	dmean2 = dmm;
	tmean2 = tmm;
	umean2 = umm;
	vmean2 = vmm;
	wmean2 = wmm;
	arnd = arnd1;
	hend = hend1;
	hnd = hnd1;
	o2nd = o2nd1;
	n2nd = n2nd1;
	nnd = nnd1;
	ond = ond1;
	wtmol = wtmol1;
	dmdz1 = dmdz2;

}

void Atm1::ncep()
{

	//The member function ncep from the Atmod class used to calculate atmospheric
	//variables with NCEP data.

	double pz1, rhoz1, tz1, uz1, vz1, wz1, tdz1, spz1, srhoz1, stz1, suz1,
		svz1, stdz1, rhn1, srhn1, vpn1, svpn1, spdavz1, spdsdz1, uvcorrz1;

	//The member function ncedmd from the NCEP class used to calculate mean and standard
	//deviation for pressure, density, temperature, east-west wind, north-south wind, 
	//vertical wind, dewpoint, vapor pressure, relative humidity and wind speed.  
	//The u-v correlation is also provided from NCEP data.
	msisa.maps.perts.iperts.ncps.ncepmd(h, phi, thet, &pz1, &rhoz1, &tz1, &uz1, &vz1, &wz1, &tdz1, &spz1,
		&srhoz1, &stz1, &suz1, &svz1, &stdz1, &rhn1, &srhn1, &vpn1, &svpn1, &spdavz1,
		&spdsdz1, &uvcorrz1);

	pnmean = pz1;
	dnmean = rhoz1;
	tnmean = tz1;
	unmean = uz1;
	vnmean = vz1;
	wnmean = wz1;
	spg = spz1;
	sdg = srhoz1;
	sug = suz1;
	svg = svz1;
	spdavg = spdavz1;
	spdsd = spdsdz1;
	vpn = vpn1;
	rhn = rhn1;
	tdmean = tdz1;
	stdd = stdz1;
	stg = stz1;
	svpn = svpn1;
	srhn = srhn1;
}


void Atm1::speconc()
{
	//The speconc member function from the Atmod class used the calculate species
	//concentration for Earth-GRAM.  
	double arppm = 9.34e+03, ckf = 273.15, heppm = 5.2,
		onec = 100.0, onemeg = 1.0e+6, pgw = 30000.0, r0 = 287.055,
		zero = 0.0, rlim = 0.99, tdgh, ws, w, dewpt,
		corr, seratio, stratio, std;

	wtair = 28.9649, wth2o = 18.01528, wto3 = 47.9982, wtn2o = 44.0129, wtco = 28.01, wtch4 = 16.043, wtco2 = 44.01, wtn2 = 28.0134,
		wto2 = 31.9988, wto = 15.9994, wtar = 39.948, wthe = 4.0026, wth = 1.00797,
		wtn = 14.0067;
	avn = 6.0221417e+26;
	//Evaluate concentrations
	//If h <= 130 km, GRAM calls the concvals member function from the MAP class
	//to evaluate concentrations.
	if (h <= 130.0) msisa.maps.concvals(h, phi, yr, pmean, waterg, &ond);

	//Heights < 120 then NCEP/LaRC/MAP/AFGL values or mixed values
	//(Thermosphere-NCEP/LaRC/MAP/AFGL
	airmw = wtmol;

	if (h < hj2){
		if (pmean >= pgw){
			//Mean H2O parameters from NCEP mean dewpoint, vapor pressure,
			//RH, and temperature if pressure > 300 mb
			eoft = vpn;
			rhp = rhn;
			airmw = wtair;
			eps = wth2o / wtair;
			rhov = eps*eoft / (r0*tmean);
			tdgh = tdmean;
		}
		else {
			//Mean H2O parameters from LaRC/MAP/AFGL concentrations (ppm)
			//if pressure < 300 mb
			//Molecular weight of air
			if (h <= hj1){
				airmw = wtair;
			}
			else if (h >= hj2){
				airmw = wtmol;
			}
			else {
				//Intermediate value if ch = 90-120 km
				airmw = r0*dmean*tmean*wtair / pmean;
			}
			//Mixing ratio, RH, dewpoint, vapor pressure, vapor density
			eps = wth2o / airmw;
			ws = mixrat(tmean, pmean);
			w = msisa.maps.perts.iperts.inits1.water / onemeg;
			rhp = 0.0;
			if (ws > 0.0) rhp = (w / ws)*(ws + eps) / (w + eps);
			dewpt = tdbuck(tmean - ckf, rhp, pmean);
			eoft = w*pmean / (w + eps);
			rhp = rhp*onec;
			tdmean = dewpt + ckf;
			rhov = eps*eoft / (r0*tmean);

			if (pmean > 10000.0){
				seoft = stdd*msisa.maps.dedt(tdmean, pmean);
				corr = 0.8;
				seratio = seoft / eoft;
				stratio = tstd*msisa.maps.dedt(tmean, pmean) / (onec*eoft);
				srhp = rhp*sqrt(pow(seratio, 2) + pow(seratio, 2) - 2 * corr*
					seratio*stratio);
				srhov = rhov*sqrt(pow(seratio, 2) + tstd);
			}
		}

		//Conversion factor from ppm to number density
		ppmtond = avn*dmean / (onemeg*airmw);

		//Water vapor concentrations
		ppmh2o = msisa.maps.perts.iperts.inits1.water;
		h2ond = msisa.maps.perts.iperts.inits1.water*ppmtond;

		//For standard deviation of water vapor use either:
		//(a) NCEP data (or RH extrapolation), surface to 100 mb
		//(b) MAP vol 31 sigma/mean for 100-0.01 mb
		//(c) A sigma/mean ratio of 0.36 (average of values in Table 1 of 
		// Harries, Rev. Geophys. Space Phys., vol 14(4), p 565, 1976) for
		// heights above 0.01 mb level

		if (pmean <= 10000.0){
			if (pmean < 1.0){
				//case with sigma/mean = 0.36 case
				seoft = 0.36*eoft;
			}
			else
			{
				seoft = (msisa.maps.perts.iperts.inits1.sigwater / msisa.maps.perts.iperts.inits1.water)*eoft;

			}

			//MAP sigma case
			if (eoft <= 0.0){
				seoft = 0.0;
				stdd = 0.0;
				srhp = 0.0;
				srhov = 0.0;
			}
			else {
				stdd = seoft / msisa.maps.dedt(tdmean, pmean);
				seratio = seoft / eoft;
				stratio = tstd*msisa.maps.dedt(tmean, pmean) / (onec*eoft);
				corr = 0.8;
				srhp = rhp*sqrt(pow(seratio, 2.0) + pow(stratio, 2.0) - 2.0*
					corr*seratio*stratio);
				srhov = rhov*sqrt(pow(seratio, 2.0) + pow((stdd / tmean), 2.0));
			}
		}
		else {
			//NCEP of RH-extrapolated case
			std = stdd;
			tstd = stg;
			if (pmean >= pgw){
				seoft = svpn;
				seratio = seoft / eoft;
				srhp = srhn;
			}
			else {
				seoft = std*msisa.maps.dedt(tdmean, pmean);
				seratio = seoft / eoft;
				stratio = tstd*msisa.maps.dedt(tmean, pmean) / (onec*eoft);
				corr = 0.8;
				srhp = rhov*sqrt(pow(seratio, 2.0) + pow(stratio, 2.0) -
					2.0 * corr*seratio*stratio);
			}
			srhov = rhov*sqrt(pow(seratio, 2.0) + pow((std / tstd), 2.0));
		}

		//Convert oxygen number density to ppm
		ppmo = ond / ppmtond;

		//Zero if ch<90 km; fixed He and A; AFGL O2 and N2, ch and N=0
		if (h <= hj1){
			hnd = zero;
			ppmh = zero;
			ppmhe = heppm;
			ppmar = arppm;
			hend = heppm*ppmtond;
			arnd = arppm*ppmtond;
			ppmo2 = msisa.maps.perts.iperts.inits1.moloxy;
			ppmn2 = msisa.maps.perts.iperts.inits1.nitrogen;
			nnd = zero;
			ppmn = zero;
		}
		else {
			//Thermosphere Ar, He if h = 90-120 km
			ppmar = arnd / ppmtond;
			ppmhe = hend / ppmtond;
			ppmh = hnd / ppmtond;
			ppmo2 = o2nd / ppmtond;
			ppmn2 = n2nd / ppmtond;
			ppmn = nnd / ppmtond;
			//Faired values for O2 and N2 if h = 90-120 km
			ppmo2 = ((h - hj1)*ppmo2 + (hj2 - h)*msisa.maps.perts.iperts.inits1.moloxy) / (hj2 - hj1);
			ppmn2 = ((h - hj1)*ppmn2 + (hj2 - h)*msisa.maps.perts.iperts.inits1.nitrogen) / (hj2 - hj1);
			o2nd = ppmo2*ppmtond;
			n2nd = ppmn2*ppmtond;
		}


		//AFGL concentrations for other species at all heights < 120
		ppmo3 = msisa.maps.perts.iperts.inits1.ozone;
		ppmn2o = msisa.maps.perts.iperts.inits1.nitrous;
		ppmco = msisa.maps.perts.iperts.inits1.carbmon;
		ppmch4 = msisa.maps.perts.iperts.inits1.methane;
		ppmco2 = msisa.maps.perts.iperts.inits1.carbdiox;
		o3nd = msisa.maps.perts.iperts.inits1.ozone*ppmtond;
		n2ond = msisa.maps.perts.iperts.inits1.nitrous*ppmtond;
		cond = msisa.maps.perts.iperts.inits1.carbmon*ppmtond;
		ch4nd = msisa.maps.perts.iperts.inits1.methane*ppmtond;
		co2nd = msisa.maps.perts.iperts.inits1.carbdiox*ppmtond;
	}
	else {
		//Heights > 120 km, concentrations (=0.0) in pure thermosphere range
		eoft = zero;
		rhov = zero;
		tdmean = zero;
		rhp = zero;
		srhp = zero;
		seoft = zero;
		srhov = zero;
		stdd = zero;
		ppmh2o = zero;
		ppmo3 = zero;
		ppmn2o = zero;
		ppmco = zero;
		ppmch4 = zero;
		ppmco2 = zero;
		h2ond = zero;
		o3nd = zero;
		n2ond = zero;
		cond = zero;
		ch4nd = zero;
		co2nd = zero;

		//Conversion factor from ppm to number density
		ppmtond = avn*dmean / (onemeg*wtmol);
		ppmn2 = n2nd / ppmtond;
		ppmo2 = o2nd / ppmtond;
		ppmo = ond / ppmtond;
		ppmar = arnd / ppmtond;
		ppmhe = hend / ppmtond;
		ppmh = hnd / ppmtond;
		ppmn = nnd / ppmtond;
	}





}


double Atm1::tdbuck(double t, double rh, double p)
{
	//tdbuck is a member function from the Atmod class.
	//Dewpoint temperature (C) and temperature T (C), relative humidity
	//rh (0-1), and pressure p (Pa), from the inverse of the Buck 4 formula;
	//See Table 2 of Elliott and Gaffen, Bull. Amer. Meterol. Soc., 72(10),
	//1507, 1991.

	double e0 = 6.1121, a = 1.0 / 227.3, b = 18.729, c = 257.87, onec = 100.0,
		ckf = 273.15, aa = 1.0007, bb = 3.46e-08;

	double alph, w, tdbuck_out;

	tdbuck_out = 0.0;
	w = wexler(t + ckf, p);
	if ((rh > 0.0) & (w > 0.0)){
		alph = log(rh*w / (onec*e0*(aa + bb*p)));
		tdbuck_out = 2.0*alph*c / ((b - alph) + sqrt((pow((b - alph), 2) - 4.0
			*alph*a*c)));
	}

	return tdbuck_out;



}

double Atm1::wexler(double t, double p)
{
	//wexler is a member function for the atmod class.
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

double Atm1::mixrat(double tx, double p)
{
	//mixrat is a member function for the atmod class
	//Water vapor volume mixing ratio (m**3 vapor/m**3 dry air)
	//Tx - dewpoint temperature in degrees Kelvin
	//p - pressure (N/m**2)

	double eoft, w, mixrat_out;
	double one = 1.0, two = 2.0;

	//Get vapor pressure in (N/m**2)
	eoft = wexler(tx, p);

	//Avoid cases that would give mixing ratio > 1
	if (eoft < p / two){
		w = eoft / (p - eoft);
	}
	else{
		w = one;
	}

	mixrat_out = w;

	return mixrat_out;

}

void Atm1::fair(int i, double hg, double pg, double dg, double tg, double ug, double vg, double wg,
	double hj, double pj, double dj, double tj, double uj, double vj, double wj, double h,
	double *p, double *d, double *t, double *u, double *v, double *w, double *czi)
{
	//fair is a member function for the atmod class
	//Fairs at height h, between pg,dg,tg,ug,vg,wg and pj,dj,tj,uj,vj,wj.  Fairing is cosine**2,
	//such that all "g" values obtain at height hg and all "j" values obtain at height hj.  The
	//faired results are p,d,t,u,v,w.  Input i = 1 means density (d) is faired and pressure (p)
	//is computed from a faired value of the gas-law constant.  If i not equal to 1, both p and
	//d are faired directly.

	double r, rg, rj, szi;

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

void Atm1::stdatm(double z, double *t, double *p, double *d)
{

	//Member function stdatm from the Atmod Class.
	//US-76 Standard Atmosphere values of the temperature (K), pressure (N/m^2),
	//and density (kg/m^3) at height z (km).  Uses vertical interpolation between
	//values in data table s (zs = height, tms = temperature, wms = molecular 
	//weight, ps = pressure).

	// DATA
	const double zs[49] =
	{
		0.0, 11.019, 20.063, 32.162, 47.35, 51.413, 71.802, 86.0, 91.0, 94.0, 97.0, 100.0, 103.0,
		106.0, 108.0, 110.0, 112.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0,
		155.0, 160.0, 165.0, 170.0, 180.0, 190.0, 210.0, 230.0, 265.0, 300.0, 350.0, 400.0,
		450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0
	};

	const double tms[49] =
	{
		288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.95, 186.87, 187.74, 190.40,
		195.08, 202.23, 212.89, 223.29, 240.0, 264.0, 300.00, 360.00, 417.23, 469.27, 516.59,
		559.63, 598.78, 634.39, 666.80, 696.29, 723.13, 747.57, 790.07, 825.31, 878.84, 915.78,
		955.20, 976.01, 990.06, 995.83, 998.22, 999.24, 999.67, 999.85, 999.93, 999.97, 999.99,
		999.99, 1000.0, 1000.0, 1000.0, 1000.0
	};

	const double wms[49] =
	{
		28.9644, 28.9644, 28.9644, 28.9644, 28.9644, 28.9644, 28.9644, 28.9522, 28.889, 28.783,
		28.620, 28.395, 28.104, 27.765, 27.521, 27.268, 27.020, 26.680, 26.205, 25.803, 25.436,
		25.087, 24.749, 24.422, 24.103, 23.792, 23.488, 23.192, 22.902, 22.342, 21.809, 20.825,
		19.952, 18.688, 17.726, 16.735, 15.984, 15.247, 14.330, 13.092, 11.505, 9.718, 7.998,
		6.579, 5.543, 4.849, 4.404, 4.122, 3.940
	};

	const double ps[49] =
	{
		1013.25, 226.32, 54.7487, 8.68014, 1.10905, 0.66938, 0.039564, 3.7338e-3, 1.5381e-3,
		9.0560e-4, 5.3571e-4, 3.2011e-4, 1.9742e-4, 1.2454e-4, 9.3188e-5, 7.1042e-5, 5.5547e-5,
		4.0096e-5, 2.5382e-5, 1.7354e-5, 1.25054e-5, 9.3568e-6, 7.2028e-6, 5.6691e-6, 4.5422e-6,
		3.6930e-6, 3.0395e-6, 2.5278e-6, 2.1210e-6, 1.5271e-6, 1.1266e-6, 6.4756e-7, 3.9276e-7,
		1.7874e-7, 8.7704e-8, 3.4498e-8, 1.4518e-8, 6.4468e-9, 3.0236e-9, 1.5137e-9,
		8.2130e-10, 4.8865e-10, 3.1908e-10, 2.2599e-10, 1.7036e-10, 1.3415e-10, 1.0873e-10,
		8.9816e-11, 7.5138e-11
	};

	// Local Variables
	int i, j;
	double alp0, alp1, alp2, alpa, alpb, g, gm, go, hm, ht, ro, rs, wm, wma, wmb,
		wmo, xi, z0, z1, z2, zl, zlm, zu;

	// Set values to 0 if z < 0 (or z > 1000 km)
	if (z < 0.0) {
		*t = 0.0;
		*p = 0.0;
		*d = 0.0;
	}
	else {
		ro = 6356.766;
		go = 9.80665;
		wmo = 28.9644;
		rs = 8314.472;

		// Do height interpolation for temperature t (K)
		if (z < 86.0) {

			for (i = 0; i < 7; ++i) {

				if ((zs[i] <= z) && (z < zs[i + 1])) {
					break;
				}
			}

			zl = ro*zs[i] / (ro + zs[i]);
			zu = ro*zs[i + 1] / (ro + zs[i + 1]);
			zlm = zl*1000.0;
			wm = wmo;
			ht = (ro*z) / (ro + z);
			hm = ht*1000.0;

			// Compute temperature gradient (deg/km)
			g = (tms[i + 1] - tms[i]) / (zu - zl);
			gm = g*0.001;

			if (fabs(g) <= 0.0001) {
				*p = ps[i] * exp(-(go*wmo*(hm - zlm)) / (rs*tms[i]))*100.0;
			}
			else {
				*p = ps[i] * (pow(tms[i] / (tms[i] + g*(ht - zl)), (go*wmo) / (rs*gm)))*100.0;
			}

			*t = tms[i] + g*(ht - zl);

			goto L25;

		}

		for (i = 7; i < 48; ++i) {

			if (zs[i] <= z && z < zs[i + 1]) {
				goto L8;
			}
		}

		i = 47;
		if (z > 1000.0) {
			*t = 0.0;
			*p = 0.0;
			*d = 0.0;

			return;
		}
	}

L8:
	if (i == 7) {
		*t = tms[8];
	}
	else if (i >= 15 && i < 18) {
		*t = 240.0 + 12.0*(z - 110.0);
	}
	else if (i < 15) {
		*t = 263.1905 - 76.3232*sqrt(1.0 - pow((z - 91.0) / 19.9429, 2.0));
	}
	else {
		xi = (z - 120.0)*(ro + 120.0) / (ro + z);
		*t = (1000.0 - 640.0*exp(-0.01875*xi));
	}

	j = i;

	// Do height interpolation for molecular weight, wm
	if (i == 47) {
		j = i - 1;
	}

	z0 = zs[j];
	z1 = zs[j + 1];
	z2 = zs[j + 2];

	wma = wms[j] * (z - z1)*(z - z2) / ((z0 - z1)*(z0 - z2))
		+ wms[j + 1] * (z - z0)*(z - z2) / ((z1 - z0)*(z1 - z2))
		+ wms[j + 2] * (z - z0)*(z - z1) / ((z2 - z0)*(z2 - z1));

	alp0 = log(ps[j]);
	alp1 = log(ps[j + 1]);
	alp2 = log(ps[j + 2]);
	alpa = alp0*(z - z1)*(z - z2) / ((z0 - z1)*(z0 - z2))
		+ alp1*(z - z0)*(z - z2) / ((z1 - z0)*(z1 - z2))
		+ alp2*(z - z0)*(z - z1) / ((z2 - z0)*(z2 - z1));
	alpb = alpa;

	wmb = wma;

	if (i != 7 && i != 47) {
		j = j - 1;

		z0 = zs[j];
		z1 = zs[j + 1];
		z2 = zs[j + 2];

		alp0 = log(ps[j]);
		alp1 = log(ps[j + 1]);
		alp2 = log(ps[j + 2]);

		alpb = alp0*(z - z1)*(z - z2) / ((z0 - z1)*(z0 - z2))
			+ alp1*(z - z0)*(z - z2) / ((z1 - z0)*(z1 - z2))
			+ alp2*(z - z0)*(z - z1) / ((z2 - z0)*(z2 - z1));

		wmb = wms[j] * (z - z1)*(z - z2) / ((z0 - z1)*(z0 - z2))
			+ wms[j + 1] * (z - z0)*(z - z2) / ((z1 - z0)*(z1 - z2))
			+ wms[j + 2] * (z - z0)*(z - z1) / ((z2 - z0)*(z2 - z1));
	}

	// Convert pressure in mb to N/m^2
	*p = 100.0*exp((alpa + alpb) / 2.0);
	wm = (wma + wmb) / 2.0;

L25:
	// Compute density from perfect gas law and molecular weight
	*d = (wm*(*p)) / (rs*(*t));

	return;
}



void Atm1::traj(double hgt, double lat, double lon, double time, int iupdate, int initonce, double *dmout,
	double *pmout, double *tmout, double *umout, double *vmout, double *wmout, double *dpout, double *ppout,
	double *tpout, double *upout, double *vpout, double *wpout, double *dsout, double *psout,
	double *tsout, double *usout, double *vsout, double *wsout, double *dsmall, double *psmall,
	double *tsmall, double *usmall, double *vsmall, double *wsmall, double *sosmean, double *sospert)
{

	//traj member function for the atmod class
	//Main member function to calculate mean and standard deviation of pressure,
	//density, temperature, and winds.  As well as perturbed values and atmospheric
	//constituents.  This member function allows the capability to plug-in to model.

	using namespace std;


	double pgh, tgh, dgh, ugh, vgh, hg1 = 0.0, hg2 = 0.0;
	double g0 = 0.0, g = 0.0, re = 0.0, r0 = 0.0, ri1 = 0.0, phir = 0.0, ri = 0.0;
	double czi, pmm, dmm, tmm, umm, vmm, wmm, water, dest,
		rhrra, shrra, sigest, srhf, corr, srhrra, tdsrf,
		stdsrf, totppm, ppmscale, totnd, dphi0 = 0.0, FS, spdgh, sdsph, su, sv,
		sp = 0.0, st = 0.0, sd = 0.0, sw = 0.0, swg = 0.0, spj, stj, sdj, suj, svj, swj;
	double r0a = 287.055, onemeg = 1.0e+06, gasconst, cpmn, csp0, csp, mwnd, t1, p1, d1;
	double pghp, dghp, tghp, php, dhp, thp, onec = 100.0, prhp, drhp, trhp, prsp, drsp,
		trsp, prlp, drlp, trlp, spsp, sdsp, stsp, splp, sdlp, stlp, sphp, sdhp, sthp;
	double gdlat, req = 6378.137, rpo = 6356.752314, gdhgt, ffac;
	string rra_id;
	string nprp = msisa.maps.perts.iperts.inits1.home + msisa.maps.perts.iperts.inits1.nprpath;
	string conp = msisa.maps.perts.iperts.inits1.home + msisa.maps.perts.iperts.inits1.conpath;
	double rrawt = 0.0, grmwt = 0.0, rrawts, grmwts, hgtp, hgtd;
	int f1 = 0, i, ihr = msisa.maps.perts.iperts.inits1.ihro;
	phir = phi * pi180;
	iupdate1 = iupdate;

	if (initonce > 0){

		h = hgt;
		phi = lat;
		thet = lon;
		elt = time;

		//If the monte carlo parameter, mc == 0, set to 1 to run through at least 1 
		//trajectory
		if (msisa.maps.perts.iperts.inits1.mc == 0){
			msisa.maps.perts.iperts.inits1.mc = 1;
		}

		//Loop through number of monte carlo runs.
		for (int j = 0; j < msisa.maps.perts.iperts.inits1.mc; j++){

			i = 0;
			nm = 0;


			//Interpret input height as radius if h > 6000
			if (h > 6000.0){
				//Get local Earth radius r0
				phir = phi*pi180;
				msisa.maps.perts.rig(0.0, phir, &g0, &g, &re, &r0, &ri);
				h = h - r0;
			}

			//Terminate if height < -30 m below sea level
			if (h < (-0.031)){
				cout << "Height below -30m sea level, h = " << h << '\n';
				system("pause");
				exit(1);

			}

			//Output seed value if mc > 1
			if (msisa.maps.perts.iperts.inits1.mc > 1){
				cout << "Seed value for Monte Carlo:  " << nr1 << '\n';
			}



			spdavg = 0.0;

			//Adjust phi, thet if trajectory cross the pole.
			if (abs(phi) > 90.0){
				phi = copysign(180.0 - abs(phi), phi);
				thet = thet + 180.0;
				if (thet >= 180.0) thet = thet - 360.0;
			}
			if (thet < (-180.0)) thet = thet + 360.0;

			//Calculate the initial perturbations at initial state.
			initpert();

			//Call pert member function to calculate xlh and zlh for next perturbed value. 
			pert();
			//Increment first number seed (nr1) for Monte Carlo runs
			nr1 = nr1 + 1;

			//Call member function initsigs from Pert class to save the initial standard
			//deviations needed for calculating perturbation.
			msisa.maps.perts.initsigs();


			//If h < 40.0 call the gethgs member function from the NCEP class to get
			//the height at the 10 mb (hg2) and 20 mb (hg1) level
			if (h < 40.0){
				msisa.maps.perts.iperts.ncps.gethgs(phi, thet, &hg1, &hg2);
			}


			//If height <= height of 10 mb level use NCEP data 
			if (h <= hg2){

				//Compute NCEP data values by calling ncep member function
				ncep();
				dtz1 = msisa.maps.perts.iperts.ncps.dtz;
				//Get mean and sigma water vapor volume mixing ratio for NCEP
				waterg = onemeg*vpn / (pnmean - vpn);
				msisa.maps.perts.iperts.inits1.sigwater = waterg*svpn / vpn;

				//Use only NCEP values if below height hg1
				if (h < hg1){
					pmean = pnmean;
					dmean = dnmean;
					tmean = tnmean;
					umean = unmean;
					vmean = vnmean;
					wmean = wnmean;
				}
				else {
					//Fair between NCEP and middle-atmosphere values if height
					//between hg1 and hg2

					//Call map member function for calculating values from MAP data
					map();

					//Call intruv member function from InitP class to calculate standard
					//deviation of winds from MAP data.  Used to calculate spdavg and 
					//spdsd.
					msisa.maps.perts.iperts.intruv(msisa.maps.perts.iperts.inits1.ur,
						msisa.maps.perts.iperts.inits1.vr, h, phi, &su, &sv);
					su = sqrt(abs(su));
					sv = sqrt(abs(sv));

					//call fair member function to fair between NCEP and MAP
					fair(1, hg1, pnmean, dnmean, tnmean, unmean, vnmean, wnmean, hg2,
						pmean1, dmean1, tmean1, umean1, vmean1, wmean1, h, &pmm,
						&dmm, &tmm, &umm, &vmm, &wmm, &czi);
					fair(2, hg1, spg, sdg, stg, sug, svg, swg, hg2, sp, sd, st, su, sv,
						sw, h, &spj, &sdj, &stj, &suj, &svj, &swj, &czi);
					pmean = pmm;
					dmean = dmm;
					tmean = tmm;
					umean = umm;
					vmean = vmm;
					wmean = wmm;

					if ((tmean1 - tnmean) != 0.0){
						ffac = (tmean - tnmean) / (tmean1 - tnmean);
						dtz1 = msisa.maps.perts.iperts.ncps.dtz + ffac*(dtzm -
							msisa.maps.perts.iperts.ncps.dtz);
					}

					//Apply ruscale to standard deviation of u and v wind components.
					suj = msisa.maps.perts.iperts.inits1.ruscale*suj;
					svj = msisa.maps.perts.iperts.inits1.ruscale*svj;

					//Calculate a faired value for wind speed and wind speed standard deviation.
					FS = pow(suj, 2.0) + pow(svj, 2.0);
					spdavg = czi*spdavg + (1.0 - czi)*sqrt(pow(umean, 2.0) + pow(vmean, 2.0) + 0.605*FS);
					spdsd = czi*spdsd + (1.0 - czi)*sqrt(0.395*FS);
				}

			}
			//Use MET, MSIS, or JB2008 thermosphere model if height above hj1, 90 km.
			else if (h >= hj1){

				//Call thermosphere model based on user-selected input 'itherm'.
				if (msisa.maps.perts.iperts.inits1.itherm == 1){
					met();
					dtz1 = msisa.maps.perts.iperts.mets.dtz;
					dmdz1 = msisa.maps.perts.iperts.mets.dmdz;
				}
				else if (msisa.maps.perts.iperts.inits1.itherm == 2){
					msis();
					dtz1 = msisa.dtz;
				}
				else if (msisa.maps.perts.iperts.inits1.itherm == 3){
					jb2008();
					dtz1 = msisa.jb2008.dtz;
				}


				//Use only thermosphere model data if height above hj2, 120 km.
				if (h >= hj2){
					pmean = pmean2;
					dmean = dmean2;
					tmean = tmean2;
					umean = umean2;
					vmean = vmean2;
					wmean = wmean2;
					Umean = umean2;
					Vmean = vmean2;
				}
				else {
					//Fair between thermosphere and middle-atmosphere values
					//if height between hj1 and hj2.
					map();
					fair(1, hj1, pmean1, dmean1, tmean1, umean1, vmean1, wmean1, hj2,
						pmean2, dmean2, tmean2, umean2, vmean2, wmean2, h, &pmm, &dmm,
						&tmm, &umm, &vmm, &wmm, &czi);

					pmean = pmm;
					dmean = dmm;
					tmean = tmm;
					umean = umm;
					vmean = vmm;
					wmean = wmm;
					Umean = umm;
					Vmean = vmm;

					if ((tmean1 - tmean2) != 0.0){
						ffac = (tmean - tmean2) / (tmean1 - tmean2);
						dtz1 = dtz1 + ffac*(dtzm - dtz1);
					}

				}
			}
			else {
				//Use only middle-atmosphere values if height between hg2 and 
				//hj1.
				map();
				pmean = pmean1;
				dmean = dmean1;
				tmean = tmean1;
				umean = umean1;
				vmean = vmean1;
				wmean = wmean1;
				Umean = umean;
				Vmean = vmean;
				dtz1 = dtzm;
			}

			//Call speconc member function to calculate species concentration data.
			speconc();

			//Use RRA data if height <= 70.0 and user-selected input variable iurra > 0.
			int iurra = msisa.maps.perts.iperts.inits1.iurra;
			if ((h <= 70.0) & (iurra > 0)){

				phir = phi*pi180;
				msisa.maps.perts.rig(h, phir, &g0, &g, &re, &r0, &ri1);

				pgh = pmean;
				dgh = dmean;
				tgh = tmean;
				ugh = umean;
				vgh = vmean;
				//Call rramean memeber function from RRA class to calculate RRA mean data.
				msisa.maps.perts.iperts.rras.rramean(h, phi, thet, mn, r0, &pgh,
					&dgh, &tgh, &ugh, &vgh, rra_id, &rrawt);

				if (msisa.maps.perts.iperts.rras.z1[msisa.maps.perts.iperts.rras.num1 - 1]){
					umean = ugh;
					vmean = vgh;
				}

				if (msisa.maps.perts.iperts.rras.z2[msisa.maps.perts.iperts.rras.num2 - 1]){
					pmean = pgh;
					dmean = dgh;
					tmean = tgh;
				}

				//Revise GRAM vapor pressure, dew point, RH, mixing ratio
				//and standard deviations by weighted average with RRA data
				int num3a = msisa.maps.perts.iperts.rras.num3;

				if (h <= msisa.maps.perts.iperts.rras.z3[num3a - 1]){
					grmwt = 1.0 - msisa.maps.perts.iperts.rras.rrawt3;
					eoft = rrawt*msisa.maps.perts.iperts.rras.vprra + grmwt*eoft;
					seoft = rrawt*msisa.maps.perts.iperts.rras.svprra + grmwt*seoft;
					tdmean = rrawt*msisa.maps.perts.iperts.rras.tdrra + grmwt*tdmean;
					stdd = rrawt*msisa.maps.perts.iperts.rras.stdrra + grmwt*stdd;
					water = 1.0e+06*eoft / (pmean - eoft);

					dest = 0.5*msisa.maps.d2edt2(msisa.maps.perts.iperts.rras.trra,
						msisa.maps.perts.iperts.rras.prra)*msisa.maps.perts.iperts.rras.strra;
					rhrra = 100.0*msisa.maps.perts.iperts.rras.vprra / (msisa.maps.wexler(msisa.maps.perts.iperts.rras.trra,
						msisa.maps.perts.iperts.rras.prra) + dest);
					rhrra = max(3.0, rhrra);
					rhrra = min(100.0, rhrra);

					rhp = rrawt*rhrra + grmwt*rhp;
					rhov = eps*eoft / (r0a*tmean);

					if (eoft <= 0.0){
						msisa.maps.perts.iperts.inits1.sigwater = 0.0;
						shrra = 0.0;
						srhov = 0.0;
					}
					else{
						msisa.maps.perts.iperts.inits1.sigwater = water*seoft / eoft;
						sigest = msisa.maps.dedt(tmean, pmean)*tstd1*tmean / eoft;
						srhov = rhov*sqrt(pow((seoft / eoft), 2.0) +
							pow((stdd / tdmean), 2.0));
						srhf = 1.43 + 0.315*log10(msisa.maps.perts.iperts.rras.prra / 1.0e+05);
						corr = 0.85 - 0.7*pow(cos(pi180*phi), 6.0);
						srhrra = rhrra*sqrt(pow((seoft / eoft), 2.0) + pow(sigest, 2.0) -
							2.0*corr*sigest*seoft / eoft) / srhf;
						srhrra = min(50.0, srhrra);
						srhrra = min(1.05*rhrra, srhrra);
						srhrra = max(0.05*rhrra, srhrra);
						srhp = rrawt*srhrra + grmwt*srhp;
					}


					ppmh2o = water;
					h2ond = water*ppmtond;

				}
				//If height <= 27.0, calcuate a weighted RRA surface dewpoint mean and
				//standard deviation.
				if (h <= 27.0){
					tdsrf = grmwt*tdmean + rrawt*msisa.maps.perts.iperts.rras.tdsrf;
					stdsrf = grmwt*stdd + rrawt*msisa.maps.perts.iperts.rras.stdsrf;
				}

			}


			//Rescale non-water-vapor concentrations to correct for the effects of 
			//interpolation, fairing, and water vapor.
			totppm = ppmo3 + ppmn2o + ppmco + ppmch4 + ppmco2 + ppmn2 + ppmo2 +
				ppmo + ppmar + ppmhe + ppmh + ppmn;

			ppmscale = (1.0e+06 - ppmh2o) / totppm;

			ppmo3 = ppmo3*ppmscale;
			ppmn2o = ppmn2o*ppmscale;
			ppmco = ppmco*ppmscale;
			ppmch4 = ppmch4*ppmscale;
			ppmco2 = ppmco2*ppmscale;
			ppmn2 = ppmn2*ppmscale;
			ppmo2 = ppmo2*ppmscale;
			ppmo = ppmo*ppmscale;
			ppmar = ppmar*ppmscale;
			ppmhe = ppmhe*ppmscale;
			ppmh = ppmh*ppmscale;
			ppmn = ppmn*ppmscale;
			ppmtond = avn*pmean / (8314.472e+6*tmean);
			h2ond = ppmh2o*ppmtond;
			o3nd = ppmo3*ppmtond;
			n2ond = ppmn2o*ppmtond;
			cond = ppmco*ppmtond;
			ch4nd = ppmch4*ppmtond;
			co2nd = ppmco2*ppmtond;
			n2nd = ppmn2*ppmtond;
			o2nd = ppmo2*ppmtond;
			ond = ppmo*ppmtond;
			arnd = ppmar*ppmtond;
			hend = ppmhe*ppmtond;
			hnd = ppmh*ppmtond;
			nnd = ppmn*ppmtond;
			totnd = h2ond + o3nd + n2ond + cond + ch4nd + co2nd + n2nd + o2nd + ond
				+ arnd + hend + hnd + nnd;
			mwnd = (h2ond*wth2o + o3nd*wto3 + n2ond*wtn2o + cond*wtco + ch4nd*wtch4 +
				co2nd*wtco2 + n2nd*wtn2 + o2nd*wto2 + ond*wto + arnd*wtar + hend*wthe +
				hnd*wth + nnd*wtn) / totnd;


			//Variables needed for using auxiliary profile data.
			double profout, tout, pout, dout, uout, vout,
				sitenear = msisa.maps.perts.iperts.inits1.sitenear,
				sitelim = msisa.maps.perts.iperts.inits1.sitelim;
			string home1 = msisa.maps.perts.iperts.inits1.home;
			string profile1 = msisa.maps.perts.iperts.inits1.profile;

			//If user-selected variable iaux > 0 use available auxiliary profile data.
			if (msisa.maps.perts.iperts.inits1.iaux > 0){

				//call rdprof and profterp member functions from the AuxProf Class to 
				//read profile data calculate mean profile data.
				msisa.maps.perts.iperts.auxs.rdprof(sitenear, sitelim, home1, profile1);
				msisa.maps.perts.iperts.auxs.profterp(h, phi, thet, tmean, pmean, dmean, umean,
					vmean, &tout, &pout, &dout, &uout, &vout, &profout);

				pmean = pout;
				dmean = dout;
				tmean = tout;
				vmean = vout;
				umean = uout;

			}

			//Calculate total perturbed values added to mean data.
			pmpert = pmean * (1 + ppert1);
			dmpert = dmean * (1 + dpert1);
			tmpert = tmean * (1 + tpert1);
			umpert = umean + upert1;
			vmpert = vmean + vpert1;
			wmpert = wmean + wpert1;


			//Calculate average wind speed and wind speed standard deviation.
			//Use wind speed data if it is available.
			if (spdavg > 0){
				spdgh = spdavg;
				sdsph = spdsd;
			}
			else{
				//Calculate wind speed if data is not available.
				msisa.maps.perts.iperts.intruv(msisa.maps.perts.iperts.inits1.ur,
					msisa.maps.perts.iperts.inits1.vr, h, phi, &su, &sv);
				su = sqrt(abs(su))*msisa.maps.perts.iperts.inits1.ruscale;
				sv = sqrt(abs(sv))*msisa.maps.perts.iperts.inits1.ruscale;
				FS = pow(su, 2.0) + pow(sv, 2.0);
				spdgh = sqrt(pow(Umean, 2.0) + pow(Vmean, 2.0) + 0.605*FS);
				sdsph = sqrt(0.3950*FS);

			}
			//Weight RRA wind speed data when it is available.
			if (msisa.maps.perts.iperts.rras.sitewgt > 0.0 && h <= msisa.maps.perts.iperts.rras.z1[msisa.maps.perts.iperts.rras.num1 - 1]){
				rrawts = msisa.maps.perts.iperts.rras.rrawt1;
				grmwts = 1.0 - rrawts;
				spdgh = rrawts*msisa.maps.perts.iperts.rras.avspdrra + grmwts*spdgh;
				sdsph = rrawts*msisa.maps.perts.iperts.rras.sdspdrra + grmwts*sdsph;
			}


			if ((msisa.maps.perts.iperts.inits1.iaux) && (profout > 0.0)){
				FS = pow(ustd1, 2.0) + pow(vstd1, 2.0);
				spdgh = (sqrt(pow(umean, 2.0) + pow(vmean, 2.0) + 0.605*FS)*profout) + spdavg*(1 - profout);
				sdsph = (sqrt(0.3950*FS)*profout) + spdsd*(1 - profout);
			}

			//Gas constant and ratio of specific heats
			gasconst = pmean / (dmean*tmean);
			cpmn = (1.4 / 0.4)*gasconst;
			//Sound speed (m/s) from mean and perturbed temperature
			csp0 = sqrt(1.4*gasconst*tmean);
			csp = sqrt(1.4*gasconst*tmpert);

			//Pressure scale height (m), density scale height (m)
			msisa.maps.perts.rig(h, phir, &g0, &g, &re, &r0, &ri);
			hgtp = pmean / (dmean*g);
			if (h <= hg2) dtz1 = dtz1 / 1000.0;
			if (h <= hj1) {
				hgtd = hgtp / (1.0 + hgtp*dtz1 / tmean);
			}
			else {
				hgtd = hgtp / (1.0 + hgtp*dtz1 / tmean - hgtp*dmdz1 / wtmol);
			}

			//Output species concentration data to "species.txt"
			species << fixed << setw(9) << setprecision(3) << h << setw(8) << setprecision(3) <<
				phi << setw(9) << setprecision(3) << thet << scientific << setw(12) << setprecision(3) <<
				ppmh2o << setw(11) << setprecision(3) << h2ond << " |" << setw(11) << setprecision(3) <<
				ppmo3 << setw(11) << setprecision(3) << o3nd << "  H2O |  O3" <<
				'\n';
			species << scientific << setw(38) << setprecision(3) << ppmn2o << setw(11) <<
				setprecision(3) << n2ond << " |" << setw(11) << setprecision(3) << ppmco <<
				setw(11) << setprecision(3) << cond << "  N2O |  CO" <<
				'\n';
			species << scientific << setw(38) << setprecision(3) << ppmch4 << setw(11) <<
				setprecision(3) << ch4nd << " |" << setw(11) << setprecision(3) << ppmco2 <<
				setw(11) << setprecision(3) << co2nd << "  CH4 | CO2" <<
				'\n';
			species << scientific << setw(38) << setprecision(3) << ppmn2 << setw(11) <<
				setprecision(3) << n2nd << " |" << setw(11) << setprecision(3) << ppmo2 <<
				setw(11) << setprecision(3) << o2nd << "   N2 |  O2" <<
				'\n';
			species << scientific << setw(38) << setprecision(3) << ppmo << setw(11) <<
				setprecision(3) << ond << " |" << setw(11) << setprecision(3) << ppmar <<
				setw(11) << setprecision(3) << arnd << "    O |  Ar" <<
				'\n';
			species << scientific << setw(38) << setprecision(3) << ppmhe << setw(11) <<
				setprecision(3) << hend << " |" << setw(11) << setprecision(3) << ppmh <<
				setw(11) << setprecision(3) << hnd << "   He |   H" <<
				'\n';
			species << scientific << setw(38) << setprecision(3) << ppmn << setw(11) <<
				setprecision(3) << nnd << " | " << "MW=" << fixed << setw(6) << setprecision(3) << mwnd <<
				scientific << setw(12) << setprecision(3) << totnd << "    N | Tot" <<
				'\n';

			species << " -------- ------- --------  ---------- ----------   ---------- ----------  ---- ---- " << '\n';

			//Output data to "Special.txt".
			output << fixed << setw(10) << setprecision(2) << elt << setw(9) << setprecision(3) <<
				h << setw(10) << setprecision(5) << phi << setw(11) << setprecision(5) << thet <<
				scientific << setw(12) << setprecision(4) << dmean << setw(12) << setprecision(4) <<
				pmean << fixed << setw(8) << setprecision(2) << tmean << setw(8) << setprecision(2) << umean <<
				setw(8) << setprecision(2) << vmean << scientific << setw(12) << setprecision(4) << dmpert <<
				setw(12) << setprecision(4) << pmpert << fixed << setw(8) << setprecision(2) << tmpert <<
				setw(8) << setprecision(2) << umpert << setw(8) << setprecision(2) << vmpert << setw(7) <<
				setprecision(2) << dstd1 * 100 << setw(7) << setprecision(2) << pstd1 * 100 << setw(7) << setprecision(2) <<
				tmean*tstd1 << setw(7) << setprecision(2) << ustd1 << setw(7) << setprecision(2) << vstd1 <<
				setw(7) << setprecision(2) << wstd1 << setw(7) << setprecision(2) << wmpert << setw(7) << setprecision(2) <<
				spdgh << setw(7) << setprecision(2) << sdsph << setw(8) << setprecision(2) <<
				csp0 << setw(8) << setprecision(2) << csp << setw(4) << msisa.maps.perts.isev <<
				'\n';

			//Process data for Output data to Output.txt

			//stdatm member function from Atmod Class to calculate standard atmosphere
			stdatm(h, &t1, &p1, &d1);

			//Calculate percent deviation of GRAM mean from stadard atmosphere
			if (p1*d1*t1 > 0){
				pghp = onec*(pmean - p1) / p1;
				dghp = onec*(dmean - d1) / d1;
				tghp = onec*(tmean - t1) / t1;
				php = onec*(pmpert - p1) / p1;
				dhp = onec*(dmpert - d1) / d1;
				thp = onec*(tmpert - t1) / t1;
			}
			else {
				pghp = 0.0;
				dghp = 0.0;
				tghp = 0.0;
				php = 0.0;
				dhp = 0.0;
				thp = 0.0;
			}




			//Converts random p,d,t and standard deviations to percent
			prhp = onec*msisa.maps.perts.iperts.prh1;
			drhp = onec*msisa.maps.perts.iperts.drh1;
			trhp = onec*msisa.maps.perts.iperts.trh1;
			prsp = onec*msisa.maps.perts.iperts.rp1s;
			drsp = onec*msisa.maps.perts.iperts.rd1s;
			trsp = onec*msisa.maps.perts.iperts.rt1s;
			prlp = onec*msisa.maps.perts.iperts.rp1l;
			drlp = onec*msisa.maps.perts.iperts.rd1l;
			trlp = onec*msisa.maps.perts.iperts.rt1l;
			sphp = onec*msisa.maps.perts.iperts.sp1;
			sdhp = onec*msisa.maps.perts.iperts.sd1;
			sthp = onec*msisa.maps.perts.iperts.st1;
			spsp = onec*msisa.maps.perts.iperts.sp1s;
			sdsp = onec*msisa.maps.perts.iperts.sd1s;
			stsp = onec*msisa.maps.perts.iperts.st1s;
			splp = onec*msisa.maps.perts.iperts.sp1l;
			sdlp = onec*msisa.maps.perts.iperts.sd1l;
			stlp = onec*msisa.maps.perts.iperts.st1l;

			phir = phi *pi180;

			msisa.maps.perts.rig(h, phir, &g0, &g, &re, &r0, &ri);

			msisa.maps.perts.iperts.rras.geocenttogeodet(ri*cos(phir), ri*sin(phir),
				&gdlat, &gdhgt, req, rpo);



			msisa.maps.perts.iperts.inits1.out << fixed << setw(9) << setprecision(3) << h
				<< setw(8) << phi << setw(9) << thet << scientific << setw(11) << pmean <<
				setw(12) << dmean << fixed << setw(7) << setprecision(1) << tmean <<
				setw(7) << umean << setw(7) << vmean << setw(7) << setprecision(3) <<
				wmean << " Mean" << '\n' << setw(9) << ri << setw(8) << gdlat << setw(17) <<
				setprecision(2) << pghp << "%" << setw(11) << dghp << "%" << setw(8) <<
				tghp << "%" << setw(26) << "M-76" << '\n' << setw(10) << setprecision(1) << elt <<
				setw(7) << rra_id << setw(9) << setprecision(4) << rrawt << setw(8) << setprecision(2) <<
				prsp << "%" << setw(11) << drsp << "%" << setw(8) << trsp << "%" << setw(7) <<
				setprecision(1) << msisa.maps.perts.iperts.ru1s << setw(7) << msisa.maps.perts.iperts.rv1s <<
				setw(12) << "ranS" << '\n' << setw(34) << setprecision(2) << spsp << "%" <<
				setw(11) << sdsp << "%" << setw(8) << stsp << "%" << setw(7) << setprecision(1) << msisa.maps.perts.iperts.su1s <<
				setw(7) << msisa.maps.perts.sv1s << setw(12) << "sigS" << '\n' <<
				" Wind and SoS" << setw(21) << setprecision(2) << prlp << "%" << setw(11) <<
				drlp << "%" << setw(8) << trlp << "%" << setw(7) << setprecision(1) << msisa.maps.perts.iperts.ru1l <<
				setw(7) << msisa.maps.perts.iperts.rv1l << setw(12) << "ranL" << '\n' << " ------------" <<
				setw(21) << setprecision(2) << splp << "%" << setw(11) << sdlp << "%" <<
				setw(8) << stlp << "%" << setw(7) << setprecision(1) << msisa.maps.perts.iperts.su1l <<
				setw(7) << msisa.maps.perts.iperts.sv1l << setw(12) << "sigL" << '\n' << " Ruv =  " <<
				setw(5) << setprecision(3) << msisa.maps.perts.iperts.uvth << setw(21) <<
				setprecision(2) << prhp << "%" << setw(11) << drhp << "%" << setw(8) <<
				trhp << "%" << setw(7) << setprecision(1) << msisa.maps.perts.iperts.urh1 <<
				setw(7) << msisa.maps.perts.iperts.vrh1 << setw(7) << setprecision(2) << msisa.maps.perts.iperts.rw1 << " ranT" << '\n' << " SpdAv= "
				<< setw(6) << setprecision(1) << spdgh << setw(20) << setprecision(2) << sphp << "%" << setw(11) <<
				sdhp << "%" << setw(8) << sthp << "%" << setw(7) << setprecision(1) << msisa.maps.perts.iperts.su1 <<
				setw(7) << msisa.maps.perts.iperts.sv1 << setw(7) << setprecision(2) << msisa.maps.perts.iperts.sw1 << " sigL" << '\n' <<
				" SpdSd=" << setw(6) << setprecision(1) << sdsph << setw(24) << scientific << setprecision(3) << pmpert << setw(12) << dmpert <<
				setw(7) << fixed << setprecision(1) << tmpert << setw(7) << umpert << setw(7) << vmpert << setw(7) << setprecision(2)
				<< wmpert << " Tot." << '\n' << " SoSav=" << setw(6) << setprecision(1) << csp0 << setw(21) << setprecision(2) <<
				php << "%" << setw(11) << dhp << "%" << setw(8) << thp << "%" << setw(26) << "T-76" <<
				'\n' << " SoSpt=" << setw(6) << setprecision(1) << csp << setw(24) << setprecision(3) <<
				scientific << eoft << setw(12) << rhov << fixed << setprecision(1) << setw(7) << tdmean <<
				setw(20) << rhp << "%" << " H2O" << '\n' << setw(37) << scientific << setprecision(3) <<
				seoft << setw(12) << srhov << fixed << setprecision(1) << setw(7) << stdd << setw(20) <<
				srhp << "%" << " sigH" <<
				'\n';


			*dmout = dmean;
			*pmout = pmean;
			*tmout = tmean;
			*umout = umean;
			*vmout = vmean;
			*wmout = wmean;
			*dpout = dmpert;
			*ppout = pmpert;
			*tpout = tmpert;
			*upout = umpert;
			*vpout = vmpert;
			*wpout = wmpert;
			*psout = pstd1 * 100;
			*dsout = dstd1 * 100;
			*tsout = tstd1*tmean;
			*usout = ustd1;
			*vsout = vstd1;
			*wsout = wstd1;
			*psmall = msisa.maps.perts.iperts.rp1s*100;
			*dsmall = msisa.maps.perts.iperts.rd1s*100;
			*tsmall = msisa.maps.perts.iperts.rt1s*100;
			*usmall = msisa.maps.perts.iperts.ru1s;
			*vsmall = msisa.maps.perts.iperts.rv1s;
			*wsmall = msisa.maps.perts.iperts.rw1;
			*sosmean = csp0;
			*sospert = csp;
		}


	}



	else if (initonce == 0) {

		msisa.maps.perts.iperts.inits1.dhgt = hgt - hgt1;
		msisa.maps.perts.iperts.inits1.delt = time - time1;
		msisa.maps.perts.iperts.inits1.dphi = lat - lat1;
		msisa.maps.perts.iperts.inits1.dthet = lon - lon1;
		h = hgt;
		phi = lat;
		thet = lon;
		elt = time;
		//dphi = msisa.maps.perts.iperts.inits1.dphi;
		//dthet = msisa.maps.perts.iperts.inits1.dthet;
		
		//Interpret input height as radius if h > 6000
		if (h > 6000.0){
			//Get local Earth radius r0
			phir = phi*pi180;
			msisa.maps.perts.rig(0.0, phir, &g0, &g, &re, &r0, &ri);
			h = h - r0;
		}

		//Terminate if height < -30 m below sea level
		if (h < (-0.031)){
			cout << "Height below -30 m sea level, h = " << h << '\n';
			exit(1);
		}

		//Set nm = 1, to calculate perturbations normally.
		nm = 1;
		spdavg = 0.0;
		
		//Adjust phi, thet, and dphi if trajectory crossed the pole.
		if (abs(phi) > 90.0){
			phi = copysign(180.0 - abs(phi), phi);
			dphi = -dphi;
			thet = thet + 180.0;
		
		}
		
		if (thet < (-180.0)) thet = thet + 360.0;
		if (thet >= 180.0) thet = thet - 360.0;
		
		//If h < 40.0 call the gethgs member function from the NCEP class to get
		//the height at the 10 mb (hg2) and 20 mb (hg1) level
		if (h < 40.0){
			msisa.maps.perts.iperts.ncps.gethgs(phi, thet, &hg1, &hg2);
		}

		//If height <= height of 10 mb level use NCEP data
		if (h <= hg2){

			//Compute NCEP data values by calling ncep member function
			ncep();
			dtz1 = msisa.maps.perts.iperts.ncps.dtz;
			//Get mean and sigma water vapor volume mixing ratio for NCEP
			waterg = onemeg*vpn / (pnmean - vpn);
			msisa.maps.perts.iperts.inits1.sigwater = waterg*svpn / vpn;

			//Use only NCEP values if below height hg1
			if (h < hg1){
				pmean = pnmean;
				dmean = dnmean;
				tmean = tnmean;
				umean = unmean;
				vmean = vnmean;
				wmean = wnmean;

			}
			else {
				//Fair between NCEP and middle-atmosphere values if height
				//between hg1 and hg2

				//Call map member function for calculating values from MAP data.
				map();

				//Call intruv member function from InitP class to calculate standard
				//deviation of winds from MAP data.  Used to calculate spdavg and 
				//spdsd.
				msisa.maps.perts.iperts.intruv(msisa.maps.perts.iperts.inits1.ur,
					msisa.maps.perts.iperts.inits1.vr, h, phi, &su, &sv);
				su = sqrt((abs(su)));
				sv = sqrt(abs(sv));

				//call fair member function to fair between NCEP and MAP
				fair(1, hg1, pnmean, dnmean, tnmean, unmean, vnmean, wnmean, hg2,
					pmean1, dmean1, tmean1, umean1, vmean1, wmean1, h, &pmm,
					&dmm, &tmm, &umm, &vmm, &wmm, &czi);
				fair(2, hg1, spg, sdg, stg, sug, svg, swg, hg2, sp, sd, st, su, sv,
					sw, h, &spj, &sdj, &stj, &suj, &svj, &swj, &czi);

				pmean = pmm;
				dmean = dmm;
				tmean = tmm;
				umean = umm;
				vmean = vmm;
				wmean = wmm;

				if ((tmean1 - tnmean) != 0.0){
					ffac = (tmean - tnmean) / (tmean1 - tnmean);
					dtz1 = msisa.maps.perts.iperts.ncps.dtz + ffac*(dtzm -
						msisa.maps.perts.iperts.ncps.dtz);
				}

				//Apply ruscale to standard deviation of u and v wind components.
				suj = msisa.maps.perts.iperts.inits1.ruscale*suj;
				svj = msisa.maps.perts.iperts.inits1.ruscale*svj;

				//Calculate a faired value for wind speed and wind speed standard deviation.
				FS = pow(suj, 2.0) + pow(svj, 2.0);
				spdavg = czi*spdavg + (1.0 - czi)*sqrt(pow(umean, 2.0) + pow(vmean, 2.0) + 0.605*FS);
				spdsd = czi*spdsd + (1.0 - czi)*sqrt(0.395*FS);

			}

		}
		//Use MET, MSIS, or JB2008 thermosphere model if height above hj1, 90 km.
		else if (h >= hj1){

			//Call thermosphere model based on user-selected input 'itherm'.
			if (msisa.maps.perts.iperts.inits1.itherm == 1){
				met();
				dtz1 = msisa.maps.perts.iperts.mets.dtz;
				dmdz1 = msisa.maps.perts.iperts.mets.dmdz;
			}
			else if (msisa.maps.perts.iperts.inits1.itherm == 2){
				msis();
				dtz1 = msisa.dtz;
			}
			else if (msisa.maps.perts.iperts.inits1.itherm == 3){
				jb2008();
				dtz1 = msisa.jb2008.dtz;
			}

			//Use only thermosphere model data if height above hj2, 120 km.
			if (h >= hj2){
				pmean = pmean2;
				dmean = dmean2;
				tmean = tmean2;
				umean = umean2;
				vmean = vmean2;
				wmean = wmean2;
				Umean = umean2;
				Vmean = vmean2;
			}
			else {
				//Fair between thermosphere and middle-atmosphere values
				//if height between hj1 and hj2.
				map();
				fair(1, hj1, pmean1, dmean1, tmean1, umean1, vmean1, wmean1, hj2,
					pmean2, dmean2, tmean2, umean2, vmean2, wmean2, h, &pmm, &dmm,
					&tmm, &umm, &vmm, &wmm, &czi);


				pmean = pmm;
				dmean = dmm;
				tmean = tmm;
				umean = umm;
				vmean = vmm;
				wmean = wmm;
				Umean = umm;
				Vmean = vmm;

				if ((tmean1 - tmean2) != 0.0){
					ffac = (tmean - tmean2) / (tmean1 - tmean2);
					dtz1 = dtz1 + ffac*(dtzm - dtz1);
				}
			}
		}
		else {
			//Use only middle-atmosphere values if height between hg2 and 
			//hj1.
			map();
			pmean = pmean1;
			dmean = dmean1;
			tmean = tmean1;
			umean = umean1;
			vmean = vmean1;
			wmean = wmean1;
			Umean = umean1;
			Vmean = vmean1;
			dtz1 = dtzm;
		}


		//Call speconc member function to calculate species concentration data.
		speconc();

		pert();

		//Use RRA data if height <= 70.0 and user-selected input variable iurra > 0.
		int iurra = msisa.maps.perts.iperts.inits1.iurra;
		if ((h <= 70.0) & (iurra > 0)){

			phir = phi*pi180;
			msisa.maps.perts.rig(h, phir, &g0, &g, &re, &r0, &ri1);

			pgh = pmean;
			dgh = dmean;
			tgh = tmean;
			ugh = umean;
			vgh = vmean;

			//Call rramean memeber function from RRA class to calculate RRA mean data.
			msisa.maps.perts.iperts.rras.rramean(h, phi, thet, mn, r0, &pgh,
				&dgh, &tgh, &ugh, &vgh, rra_id, &rrawt);

			if (msisa.maps.perts.iperts.rras.z1[msisa.maps.perts.iperts.rras.num1 - 1]){
				umean = ugh;
				vmean = vgh;
			}

			if (msisa.maps.perts.iperts.rras.z2[msisa.maps.perts.iperts.rras.num2 - 1]){
				pmean = pgh;
				dmean = dgh;
				tmean = tgh;
			}


			//Revise GRAM vapor pressure, dew point, RH, mixing ratio
			//and standard deviations by weighted average with RRA data
			int num3b = msisa.maps.perts.iperts.rras.num3;
			if (h <= msisa.maps.perts.iperts.rras.z3[num3b - 1]){
				grmwt = 1.0 - msisa.maps.perts.iperts.rras.rrawt3;

				eoft = rrawt*msisa.maps.perts.iperts.rras.vprra + grmwt*eoft;
				seoft = rrawt*msisa.maps.perts.iperts.rras.svprra + grmwt*seoft;
				tdmean = rrawt*msisa.maps.perts.iperts.rras.tdrra + grmwt*tdmean;
				stdd = rrawt*msisa.maps.perts.iperts.rras.stdrra + grmwt*stdd;
				water = 1.0e+06*eoft / (pmean - eoft);

				dest = 0.5*msisa.maps.d2edt2(msisa.maps.perts.iperts.rras.trra,
					msisa.maps.perts.iperts.rras.prra)*msisa.maps.perts.iperts.rras.strra;
				rhrra = 100.0*msisa.maps.perts.iperts.rras.vprra / (msisa.maps.wexler(msisa.maps.perts.iperts.rras.trra,
					msisa.maps.perts.iperts.rras.prra) + dest);
				rhrra = max(3.0, rhrra);
				rhrra = min(100.0, rhrra);
				rhp = rrawt*rhrra + grmwt*rhp;
				rhov = eps*eoft / (r0a*tmean);
				if (eoft <= 0.0){
					msisa.maps.perts.iperts.inits1.sigwater = 0.0;
					shrra = 0.0;
					srhov = 0.0;
				}
				else{
					msisa.maps.perts.iperts.inits1.sigwater = water*seoft / eoft;
					sigest = msisa.maps.dedt(tmean, pmean)*msisa.maps.perts.sth*tmean / eoft;

					srhov = rhov*sqrt(pow((seoft / eoft), 2.0) +
						pow((stdd / tdmean), 2.0));
					srhf = 1.43 + 0.315*log10(msisa.maps.perts.iperts.rras.prra / 1.0e+05);
					corr = 0.85 - 0.7*pow(cos(pi180*phi), 6.0);
					srhrra = rhrra*sqrt(pow((seoft / eoft), 2.0) + pow(sigest, 2.0) -
						2.0*corr*sigest*seoft / eoft) / srhf;
					srhrra = min(50.0, srhrra);
					srhrra = min(1.05*rhrra, srhrra);
					srhrra = max(0.05*rhrra, srhrra);
					srhp = rrawt*srhrra + grmwt*srhp;
				}

				ppmh2o = water;
				h2ond = water*ppmtond;


			}
			//If height <= 27.0, calcuate a weighted RRA surface dewpoint mean and
			//standard deviation.
			if (h <= 27.0){
				tdsrf = grmwt*tdmean + rrawt*msisa.maps.perts.iperts.rras.tdsrf;
				stdsrf = grmwt*stdd + rrawt*msisa.maps.perts.iperts.rras.stdsrf;
			}

		}




		//Rescale non-water-vapor concentrations to correct for the effects of 
		//interpolation, fairing, and water vapor
		totppm = ppmo3 + ppmn2o + ppmco + ppmch4 + ppmco2 + ppmn2 + ppmo2 +
			ppmo + ppmar + ppmhe + ppmh + ppmn;

		ppmscale = (1.0e+06 - ppmh2o) / totppm;

		ppmo3 = ppmo3*ppmscale;
		ppmn2o = ppmn2o*ppmscale;
		ppmco = ppmco*ppmscale;
		ppmch4 = ppmch4*ppmscale;
		ppmco2 = ppmco2*ppmscale;
		ppmn2 = ppmn2*ppmscale;
		ppmo2 = ppmo2*ppmscale;
		ppmo = ppmo*ppmscale;
		ppmar = ppmar*ppmscale;
		ppmhe = ppmhe*ppmscale;
		ppmh = ppmh*ppmscale;
		ppmn = ppmn*ppmscale;
		ppmtond = avn*pmean / (8314.472e+6*tmean);

		h2ond = ppmh2o*ppmtond;
		o3nd = ppmo3*ppmtond;
		n2ond = ppmn2o*ppmtond;
		cond = ppmco*ppmtond;
		ch4nd = ppmch4*ppmtond;
		co2nd = ppmco2*ppmtond;
		n2nd = ppmn2*ppmtond;
		o2nd = ppmo2*ppmtond;
		ond = ppmo*ppmtond;
		arnd = ppmar*ppmtond;
		hend = ppmhe*ppmtond;
		hnd = ppmh*ppmtond;
		nnd = ppmn*ppmtond;

		totnd = h2ond + o3nd + n2ond + cond + ch4nd + co2nd + n2nd + o2nd + ond
			+ arnd + hend + hnd + nnd;
		mwnd = (h2ond*wth2o + o3nd*wto3 + n2ond*wtn2o + cond*wtco + ch4nd*wtch4 +
			co2nd*wtco2 + n2nd*wtn2 + o2nd*wto2 + ond*wto + arnd*wtar + hend*wthe +
			hnd*wth + nnd*wtn) / totnd;

		//If user-selected variable iaux > 0 use available auxiliary profile data.

		double profout, tout, pout, dout, uout, vout,
			sitenear = msisa.maps.perts.iperts.inits1.sitenear,
			sitelim = msisa.maps.perts.iperts.inits1.sitelim;
		string home1 = msisa.maps.perts.iperts.inits1.home;
		string profile1 = msisa.maps.perts.iperts.inits1.profile;

		if (msisa.maps.perts.iperts.inits1.iaux > 0){
			msisa.maps.perts.iperts.auxs.rdprof(sitenear, sitelim, home1, profile1);
			msisa.maps.perts.iperts.auxs.profterp(h, phi, thet, tmean, pmean, dmean, umean,
				vmean, &tout, &pout, &dout, &uout, &vout, &profout);

			pmean = pout;
			dmean = dout;
			tmean = tout;
			vmean = vout;
			umean = uout;

		}

		//Call the member function pert to calculate perturbation at current 
		//state.


		//Calculate total perturbed values added to mean data.
		pmpert = pmean * (1 + ppert);
		dmpert = dmean * (1 + dpert);
		tmpert = tmean * (1 + tpert);
		umpert = umean + upert;
		vmpert = vmean + vpert;
		wmpert = wmean + wpert;

		//Calculate average wind speed and wind speed standard deviation.
		//Use wind speed data if it is available.
		if (spdavg > 0){
			spdgh = spdavg;
			sdsph = spdsd;
		}
		else{
			//Calculate wind speed if data is not available.
			msisa.maps.perts.iperts.intruv(msisa.maps.perts.iperts.inits1.ur,
				msisa.maps.perts.iperts.inits1.vr, h, phi, &su, &sv);
			su = sqrt(abs(su))*msisa.maps.perts.iperts.inits1.ruscale;
			sv = sqrt(abs(sv))*msisa.maps.perts.iperts.inits1.ruscale;
			FS = pow(su, 2.0) + pow(sv, 2.0);
			spdgh = sqrt(pow(Umean, 2.0) + pow(Vmean, 2.0) + 0.605*FS);
			sdsph = sqrt(0.3950*FS);
		}
		//Weight RRA wind speed data when it is available.
		if (msisa.maps.perts.iperts.rras.sitewgt > 0.0 && h <= msisa.maps.perts.iperts.rras.z1[msisa.maps.perts.iperts.rras.num1 - 1]){
			rrawts = msisa.maps.perts.iperts.rras.rrawt1;
			grmwts = 1.0 - rrawts;
			spdgh = rrawts*msisa.maps.perts.iperts.rras.avspdrra + grmwts*spdgh;
			sdsph = rrawts*msisa.maps.perts.iperts.rras.sdspdrra + grmwts*sdsph;

		}

		if ((msisa.maps.perts.iperts.inits1.iaux) && (profout > 0.0)){
			FS = pow(msisa.maps.perts.suh, 2.0) + pow(msisa.maps.perts.svh, 2.0);
			spdgh = (sqrt(pow(umean, 2.0) + pow(vmean, 2.0) + 0.605*FS)*profout) + spdavg*(1 - profout);
			sdsph = (sqrt(0.3950*FS)*profout) + spdsd*(1 - profout);
		}

		//Gas constant and ratio of specific heats
		gasconst = pmean / (dmean*tmean);
		cpmn = (1.4 / 0.4)*gasconst;
		//Sound speed (m/s) from mean and perturbed temperature
		csp0 = sqrt(1.4*gasconst*tmean);
		csp = sqrt(1.4*gasconst*tmpert);

		//Pressure scale height (m), density scale height (m)
		msisa.maps.perts.rig(h, phir, &g0, &g, &re, &r0, &ri);
		hgtp = pmean / (dmean*g);
		if (h <= hg2) dtz1 = dtz1 / 1000.0;
		if (h <= hj1) {
			hgtd = hgtp / (1.0 + hgtp*dtz1 / tmean);
		}
		else {
			hgtd = hgtp / (1.0 + hgtp*dtz1 / tmean - hgtp*dmdz1 / wtmol);
		}

		//Output species concentration data to "species.txt"
		species << fixed << setw(9) << setprecision(3) << h << setw(8) << setprecision(3) <<
			phi << setw(9) << setprecision(3) << thet << scientific << setw(12) << setprecision(3) <<
			ppmh2o << setw(11) << setprecision(3) << h2ond << " |" << setw(11) << setprecision(3) <<
			ppmo3 << setw(11) << setprecision(3) << o3nd << "  H2O |  O3" <<
			'\n';
		species << scientific << setw(38) << setprecision(3) << ppmn2o << setw(11) <<
			setprecision(3) << n2ond << " |" << setw(11) << setprecision(3) << ppmco <<
			setw(11) << setprecision(3) << cond << "  N2O |  CO" <<
			'\n';
		species << scientific << setw(38) << setprecision(3) << ppmch4 << setw(11) <<
			setprecision(3) << ch4nd << " |" << setw(11) << setprecision(3) << ppmco2 <<
			setw(11) << setprecision(3) << co2nd << "  CH4 | CO2" <<
			'\n';
		species << scientific << setw(38) << setprecision(3) << ppmn2 << setw(11) <<
			setprecision(3) << n2nd << " |" << setw(11) << setprecision(3) << ppmo2 <<
			setw(11) << setprecision(3) << o2nd << "   N2 |  O2" <<
			'\n';
		species << scientific << setw(38) << setprecision(3) << ppmo << setw(11) <<
			setprecision(3) << ond << " |" << setw(11) << setprecision(3) << ppmar <<
			setw(11) << setprecision(3) << arnd << "    O |  Ar" <<
			'\n';
		species << scientific << setw(38) << setprecision(3) << ppmhe << setw(11) <<
			setprecision(3) << hend << " |" << setw(11) << setprecision(3) << ppmh <<
			setw(11) << setprecision(3) << hnd << "   He |   H" <<
			'\n';
		species << scientific << setw(38) << setprecision(3) << ppmn << setw(11) <<
			setprecision(3) << nnd << " | " << "MW=" << fixed << setw(6) << setprecision(3) << mwnd <<
			scientific << setw(12) << setprecision(3) << totnd << "    N | Tot" <<
			'\n';

		species << " -------- ------- --------  ---------- ----------   ---------- ----------  ---- ---- " << '\n';

		//Output data to "Special.txt".
		output << fixed << setw(10) << setprecision(2) << elt << setw(9) << setprecision(3) <<
			h << setw(10) << setprecision(5) << phi << setw(11) << setprecision(5) << thet <<
			scientific << setw(12) << setprecision(4) << dmean << setw(12) << setprecision(4) <<
			pmean << fixed << setw(8) << setprecision(2) << tmean << setw(8) << setprecision(2) << umean <<
			setw(8) << setprecision(2) << vmean << scientific << setw(12) << setprecision(4) << dmpert <<
			setw(12) << setprecision(4) << pmpert << fixed << setw(8) << setprecision(2) << tmpert <<
			setw(8) << setprecision(2) << umpert << setw(8) << setprecision(2) << vmpert << setw(7) <<
			setprecision(2) << msisa.maps.perts.sdh * 100 << setw(7) << setprecision(2) << msisa.maps.perts.sph * 100 << setw(7) << setprecision(2) <<
			tmean*msisa.maps.perts.sth << setw(7) << setprecision(2) << msisa.maps.perts.suh << setw(7) << setprecision(2) << msisa.maps.perts.svh <<
			setw(7) << setprecision(2) << msisa.maps.perts.iperts.swh << setw(7) << setprecision(2) << wmpert << setw(7) << setprecision(2) <<
			spdgh << setw(7) << setprecision(2) << sdsph << setw(8) << setprecision(2) <<
			csp0 << setw(8) << setprecision(2) << csp << setw(4) << msisa.maps.perts.isev <<
			'\n';

		//Process data for Output data to Output.txt

		//stdatm member function from Atmod Class to calculate standard atmosphere
		stdatm(h, &t1, &p1, &d1);
		//Calculate percent deviation of GRAM mean from stadard atmosphere
		if (p1*d1*t1 > 0){
			pghp = onec*(pmean - p1) / p1;
			dghp = onec*(dmean - d1) / d1;
			tghp = onec*(tmean - t1) / t1;
			php = onec*(pmpert - p1) / p1;
			dhp = onec*(dmpert - d1) / d1;
			thp = onec*(tmpert - t1) / t1;
		}
		else {
			pghp = 0.0;
			dghp = 0.0;
			tghp = 0.0;
			php = 0.0;
			dhp = 0.0;
			thp = 0.0;
		}

		//Converts random p,d,t and standard deviations to percent
		prhp = onec*msisa.maps.perts.prh;
		drhp = onec*msisa.maps.perts.drh;
		trhp = onec*msisa.maps.perts.trh;
		prsp = onec*msisa.maps.perts.prhs;
		drsp = onec*msisa.maps.perts.drhs;
		trsp = onec*msisa.maps.perts.trhs;
		prlp = onec*msisa.maps.perts.prhl;
		drlp = onec*msisa.maps.perts.drhl;
		trlp = onec*msisa.maps.perts.trhl;
		sphp = onec*msisa.maps.perts.sph;
		sdhp = onec*msisa.maps.perts.sdh;
		sthp = onec*msisa.maps.perts.sth;
		splp = onec*msisa.maps.perts.sphl;
		sdlp = onec*msisa.maps.perts.sdhl;
		stlp = onec*msisa.maps.perts.sthl;
		spsp = onec*msisa.maps.perts.sphs;
		sdsp = onec*msisa.maps.perts.sdhs;
		stsp = onec*msisa.maps.perts.sths;

		phir = phi *pi180;

		msisa.maps.perts.rig(h, phir, &g0, &g, &re, &r0, &ri);

		msisa.maps.perts.iperts.rras.geocenttogeodet(ri*cos(phir), ri*sin(phir),
			&gdlat, &gdhgt, req, rpo);

		msisa.maps.perts.iperts.inits1.out << "-------- ------- -------- --------- ---------- ------ ------ ------ ------ ----" << '\n';

		msisa.maps.perts.iperts.inits1.out << fixed << setw(9) << setprecision(3) << h
			<< setw(8) << phi << setw(9) << thet << scientific << setw(11) << pmean <<
			setw(12) << dmean << fixed << setw(7) << setprecision(1) << tmean <<
			setw(7) << umean << setw(7) << vmean << setw(7) << setprecision(3) <<
			wmean << " Mean" << '\n' << setw(9) << ri << setw(8) << gdlat << setw(17) <<
			setprecision(2) << pghp << "%" << setw(11) << dghp << "%" << setw(8) <<
			tghp << "%" << setw(26) << "M-76" << '\n' << setw(10) << setprecision(1) << elt <<
			setw(7) << rra_id << setw(9) << setprecision(4) << rrawt << setw(8) << setprecision(2) <<
			prsp << "%" << setw(11) << drsp << "%" << setw(8) << trsp << "%" << setw(7) <<
			setprecision(1) << msisa.maps.perts.urhs << setw(7) << msisa.maps.perts.vrhs <<
			setw(12) << "ranS" << '\n' << setw(34) << setprecision(2) << spsp << "%" <<
			setw(11) << sdsp << "%" << setw(8) << stsp << "%" << setw(7) << setprecision(1) << msisa.maps.perts.suhs <<
			setw(7) << msisa.maps.perts.svhs << setw(12) << "sigS" << '\n' <<
			" Wind and SoS" << setw(21) << setprecision(2) << prlp << "%" << setw(11) <<
			drlp << "%" << setw(8) << trlp << "%" << setw(7) << setprecision(1) << msisa.maps.perts.urhl <<
			setw(7) << msisa.maps.perts.vrhl << setw(12) << "ranL" << '\n' << " ------------" <<
			setw(21) << setprecision(2) << splp << "%" << setw(11) << sdlp << "%" <<
			setw(8) << stlp << "%" << setw(7) << setprecision(1) << msisa.maps.perts.suhl <<
			setw(7) << msisa.maps.perts.svhl << setw(12) << "sigL" << '\n' << " Ruv =  " <<
			setw(5) << setprecision(3) << msisa.maps.perts.uvt2 << setw(21) <<
			setprecision(2) << prhp << "%" << setw(11) << drhp << "%" << setw(8) <<
			trhp << "%" << setw(7) << setprecision(1) << msisa.maps.perts.urh <<
			setw(7) << msisa.maps.perts.vrh << setw(7) << setprecision(2) << msisa.maps.perts.wrh << " ranT" << '\n' << " SpdAv= "
			<< setw(6) << setprecision(1) << spdgh << setw(20) << setprecision(2) << sphp << "%" << setw(11) <<
			sdhp << "%" << setw(8) << sthp << "%" << setw(7) << setprecision(1) << msisa.maps.perts.suh <<
			setw(7) << msisa.maps.perts.svh << setw(7) << setprecision(2) << msisa.maps.perts.iperts.sw1 << " sigT" << '\n' <<
			" SpdSd=" << setw(6) << setprecision(1) << sdsph << setw(24) << scientific << setprecision(3) << pmpert << setw(12) << dmpert <<
			setw(7) << fixed << setprecision(1) << tmpert << setw(7) << umpert << setw(7) << vmpert << setw(7) << setprecision(2)
			<< wmpert << " Tot." << '\n' << " SoSav=" << setw(6) << setprecision(1) << csp0 << setw(21) << setprecision(2) <<
			php << "%" << setw(11) << dhp << "%" << setw(8) << thp << "%" << setw(26) << "T-76" <<
			'\n' << " SoSpt=" << setw(6) << setprecision(1) << csp << setw(24) << setprecision(3) <<
			scientific << eoft << setw(12) << rhov << fixed << setprecision(1) << setw(7) << tdmean <<
			setw(20) << rhp << "%" << " H2O" << '\n' << setw(37) << scientific << setprecision(3) <<
			seoft << setw(12) << srhov << fixed << setprecision(1) << setw(7) << stdd << setw(20) <<
			srhp << "%" << " sigH" <<
			'\n';

		*dmout = dmean;
		*pmout = pmean;
		*tmout = tmean;
		*umout = umean;
		*vmout = vmean;
		*wmout = wmean;
		*dpout = dmpert;
		*ppout = pmpert;
		*tpout = tmpert;
		*upout = umpert;
		*vpout = vmpert;
		*wpout = wmpert;
		*psout = msisa.maps.perts.sph * 100;
		*dsout = msisa.maps.perts.sdh * 100;
		*tsout = msisa.maps.perts.sth*tmean;
		*usout = msisa.maps.perts.suh;
		*vsout = msisa.maps.perts.svh;
		*wsout = msisa.maps.perts.iperts.swh;
		*psmall = msisa.maps.perts.prhs * 100;
		*dsmall = msisa.maps.perts.drhs * 100;
		*tsmall = msisa.maps.perts.trhs * 100;
		*usmall = msisa.maps.perts.urhs;
		*vsmall = msisa.maps.perts.vrhs;
		*wsmall = wpert;
		*sosmean = csp0;
		*sospert = csp;


	}

	hgt1 = hgt;
	lat1 = lat;
	lon1 = lon;
	time1 = time;

	

	

	
	
}

/*int main(){

#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <string>
#include <cmath>
#include <stdio.h>
#include "Atmod1.h"
#include "iomanip"
#include <cstdlib>
	
	Atm1 atms1;
	
	double pm1, dm1, tm1, um1, vm1, wm1, pp1, dp1, tp1, up1, vp1, wp1, ps1, ds1, ts1,
		us1, vs1, ws1, psmall, dsmall, tsmall, usmall, vsmall, wsmall, sos, sosp;
	
	atms1.initdata();

	double h = 140.0, phi = 0.45, thet = -164.53, time = 0.0, time1 = 0.0;
	int iupdate = 1, initonce = 1;
	
	atms1.traj(h, phi, thet, time, iupdate, initonce, &dm1, &pm1, &tm1, &um1, &vm1, &wm1,
		&dp1, &pp1, &tp1, &up1, &vp1, &wp1, &ds1, &ps1, &ts1, &us1, &vs1, &ws1, &dsmall,
		&psmall, &tsmall, &usmall, &vsmall, &wsmall, &sos, &sosp);

	for (int i = 0; i < 70; i++){
		time = time + 60.0;
		initonce = 0;
		phi = phi + 0.4;
		thet = thet + 1.2;
		h = h - 2.0;

		atms1.traj(h, phi, thet, time, iupdate, initonce, &dm1, &pm1, &tm1, &um1, &vm1, &wm1,
			&dp1, &pp1, &tp1, &up1, &vp1, &wp1, &ds1, &ps1, &ts1, &us1, &vs1, &ws1, &dsmall,
			&psmall, &tsmall, &usmall, &vsmall, &wsmall, &sos, &sosp);


	}

	//Close output files
	atms1.output.close();
	atms1.species.close();
	atms1.msisa.maps.perts.iperts.inits1.ibl.close();
	atms1.msisa.maps.perts.iperts.inits1.out.close();

	system("pause");

}*/