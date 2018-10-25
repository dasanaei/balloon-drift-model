//RRA class for reading and processing Range Reference Atmosphere (RRA) 
//database sites.
//P. White

#include <iostream>
#include <fstream>
#include <sstream>
#include "RRA.h"
#include <iomanip>
#include <cstdlib>
#include <cmath>

using namespace std;


RRA::RRA()
{
	initializeMemberVariables();
}


RRA::~RRA()
{

	delete[] u;
	delete[] v;
	delete[] su;
	delete[] sv;
	delete[] ruvt;
	delete[] avspd;
	delete[] sdspd;
	delete[] p;
	delete[] sp;
	delete[] t;
	delete[] st;
	delete[] d;
	delete[] sd;
	delete[] vp;
	delete[] svp;
	delete[] td;
	delete[] std;

}

void RRA::initializeMemberVariables(){


	//Initialize member variables for RRA class

	for (int i = 0; i < 100; i++){
		rralat[i] = 0.0;
		rralon[i] = 0.0;
		rrasfc[i] = 0.0;
		zmax[i] = 0.0;

	}

	u = new float[300];
	v = new float[300];
	su = new float[300];
	sv = new float[300];
	ruvt = new float[300];
	avspd = new float[300];
	sdspd = new float[300];
	p = new float[300];
	sp = new float[300];
	t = new float[300];
	st = new float[300];
	d = new float[300];
	sd = new float[300];
	vp = new float[300];
	svp = new float[300];
	td = new float[300];
	std = new float[300];


	for (int i = 0; i < 300; i++){
		z1[i] = 0.0;
		u[i] = 0.0;
		su[i] = 0.0;
		v[i] = 0.0;
		sv[i] = 0.0;
		z2[i] = 0.0;
		p[i] = 0.0;
		sp[i] = 0.0;
		t[i] = 0.0;
		st[i] = 0.0;
		d[i] = 0.0;
		sd[i] = 0.0;
		z3[i] = 0.0;
		vp[i] = 0.0;
		svp[i] = 0.0;
		td[i] = 0.0;
		std[i] = 0.0;
	}

	for (int i = 0; i < 99; i++){
		rrasite[i] = "0";
		nsrf1[i] = 0;
		nsrf2[i] = 0;
		nsrf3[i] = 0;
	}

	pi = 4.0*atan(1.0);
	pi180 = pi / 180.0;
	nrra = 0;
	sitelim = 0.0;
	sitenear = 0.0;
	isitered = 0;
	rrascale = 1.0;
	rrsrf = 0.0, usrf = 0.0, vsrf = 0.0, susrf = 0.0, svsrf = 0.0, shsrf = 0.0,
		uvtsrf = 0.0, spdavsrf = 0.0, spdsdsrf = 0.0;
	psrf = 0.0, spsrf = 0.0, sdsrf = 0.0, dsrf = 0.0, tsrf = 0.0,
		strrasrf = 0.0, stsrf = 0.0, uvrra = 0.0;
	rrawt3 = 0.0, rrawt1 = 0.0, rrawt2 = 0.0;
	num1 = 0, num2 = 0, num3 = 0, sitewgt = 0.0;
}

void RRA::init(string namef)
{
	//init member function from RRA class
	//Read RRA site file for initializing RRA data

	string rracode, rrayear, wmo;
	char rraname[50];
	double phidet, phicent, gclat;
	const double fac = atan(1.0) / 45;
	double g0, re, r0, a, b, g, ri, rtot, xy, z;

	namef1 = namef;
	namelist();

	int i;
	ifstream rrasites;

	//Open file for Range Reference Atmosphere (RRA) site data
	lenpath = rrapath.find(" ");

	rrasites.open((rrapath.substr(0, lenpath) + "rrasites.txt").
		c_str());

	//File open error check
	if (rrasites.is_open()){
	}
	else {
		cout << "File Open Error!  " << rrapath << '\n';
		system("pause");
		exit(1);
	}

	//Read RRA site data, load into arrays
	for (i = 0; i < 11; ++i) {
		rrasites >> rracode;

	}

	nrra = 1;

	//Read file to end of file
	while (!rrasites.eof()){

		rrasites >> rracode >> rrayear >> phidet >> gclat >> rralon[nrra - 1] >> rrasfc[nrra - 1] >> zmax[nrra - 1] >> wmo;

		rrasites.get(rraname, 50);

		if (rrayear.substr(2, 3) != rrayr1) {
			continue;
		}

		if (rrasites.eof()){
			break;
		}

		//Get local Earth radius
		rig(0.0, phidet*fac, &g0, &g, &re, &r0, &ri, &a, &b);

		//Compute and store geocentric latitude and rrasite
		geodettogeocent(phidet, 0.0, &phicent, &rtot, &xy, &z, a, b);
		rralat[nrra - 1] = phicent;
		rrasite[nrra - 1] = rracode + rrayear.substr(2, 3);
		nrra++;

	}

	//Close RRA site data file to re-use unit for RRA data files
	rrasites.close();
	isitered = 0;
	
	return;
}

void RRA::namelist()
{

	//namelist member function from RRA class
	//Open and read namelist input file for setting model parameters.

	ifstream namelist;
	string dummy, NCEPpath, NCEPmn, atmpath, trapath, prtpath, nprpath,
		conpath, rndpath, profile, NCEPpath1;
	double h1, phi1, thet1, s10, s10b, xm10, xm10b, y10,
		y10b, dstdtc, dphi, dthet, dhgt, delt, rpscale, ruscale, rwscale,
		rdinit, rtinit, ruinit, rvinit, rwinit, patchy,
		z0in, f10, f10b, ap, seco;
	int iopt, NCEPyr, NCEPhr, nr1,
		iurra, initpert, itherm, ibltest, nmax, ida, iyr,
		mino, ihro, mn, iaux, mc;


	string mn1 = "12";
	NCEPmn = "Nb9008" + mn1 + ".bin";



	namelist.open(namef1.c_str());

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



	if (iyrrra == 1){
		rrayr1 = "83";
	}
	else if (iyrrra == 2){
		rrayr1 = "06";
	}
	else {
		rrayr1 = "13";
	}



}

void RRA::readrra1(string &rra_id, double xlat, double xlon, int m, int isite)
{

	//readrra1 member function from RRA class
	//Reads Range Reference Atmosphere (RRA) wind components (u and v) and 
	//standard deviations (su and sv) for site at geocentric latitude = xlat
	//longitude = xlon from RRA data file T1yyrra.txt.  Month is m.

	string lat, lon, ns, ew, header;
	stringstream convert_lat, convert_lon;
	double ztop, zlat, zlon, zlati = 0.0, zloni = 0.0, zprev, zcurr, uwn, usd, ruv, vwn,
		vsd, ws, sws, skew, xnobs;
	int maxnum, count = 0, head_count = 0, pop_header = 1, m1;

	//Highest RRA data allowed is 70 km; max number = 300
	ztop = zmax[isite - 1];
	maxnum = 300;
	
	cout << "Reading RRA data" << '\n';

	ifstream rratable1;

	//Initialize number of RRA heights
	num1 = 0;

	//Open file for reading RRA data
	rratable1.open((rrapath.substr(0, lenpath) +
		"T1" + rra_id + ".txt").c_str());

	//File open error check
	if (rratable1.is_open()){
	}
	else{
		cout << "File Open Error!  " << rrapath << '\n';
		system("pause");
		exit(1);
	}

	//Read lat-lon on data file
	rratable1.seekg(17);
	rratable1 >> lat >> lon;

	convert_lat << lat.substr(0, lat.length() - 1) << " " << lat.substr(lat.length() - 1, lat.length())
		<< std::endl;
	convert_lon << lon.substr(0, lon.length() - 1) << " " << lon.substr(lon.length() - 1, lon.length())
		<< std::endl;

	convert_lat >> zlat >> ns;
	convert_lon >> zlon >> ew;

	if (ew.length() > 1) {
		ew.resize(1);
	}

	//Make north latitudes positive and south latitudes negative or 
	//write error message and stop otherwise
	if (ns == "N") {
		zlati = zlat;
	}
	else if (ns == "S") {
		zlati = -zlat;
	}
	else {
		cout << " Bad lat code " << xlat << "  " << zlati << "  " << zlat << "  "
			<< ns << '\n';
		cout << " RRA 1 Bad lat code " << '\n';
		system("pause");
		exit(1);
	}

	// Make east longitudes positive and west longitudes negative or 
	//write error message and stop otherwise
	if (ew == "E") {
		zloni = zlon;
	}
	else if (ew == "W") {
		zloni = -zlon;
	}
	else {
		cout << " Bad lon code " << xlon << "  " << zloni << "  " << zlon << "  " << ew << '\n';
		cout << "RRA 1 Bad lon code" << '\n';
		system("pause");
		exit(1);
	}

	//Write error message and stop if RRA latitude from file does not agree with 
	//expected latitude.
	if (abs(xlat - zlati) > 0.005){
		cout << "Bad latitude on input " << xlat << "  " << zlati << "  " << zlat << "  " << ns << '\n';
		cout << " RRA 1 Bad lat input" << '\n';
		system("pause");
		exit(1);
	}

	//Write error message and stop if RRA longitude from file does not agree with
	//expected longitude
	if (abs(xlon - zloni) > 0.005){
		cout << "Bad longitude on input " << xlon << "  " << zloni << "  " << zlon << "  " << ew << '\n';
		cout << " RRA 1 Bad lon input" << '\n';
		system("pause");
		exit(1);
	}

	// Initialize month read from file and previous height read (first data in file are annual mean
	// values, month = 0)
	if (rra_id.substr(3, 4) == "83") {
		m1 = 0;
	}
	else {
		m1 = 1;
	}

	zprev = -99.0;

	nsrf1[isite - 1] = 0;

	// Read a data line from the RRA file; ignore any lines with character data (or other condition
	// that causes a read error
	while (m1 <= 12){

		++count;
		if (count > 5000) {
			break;
		}

		if (pop_header == 1) {
			rratable1 >> header;

			if ((header == "M/S") || (header == "m/s")) {
				++head_count;
			}

			if (head_count == 6) {
				head_count = 0;
				pop_header = 0;
			}

			continue;
		}

		if (pop_header == 0) {
			rratable1 >> zcurr >> uwn >> usd >> ruv >> vwn >> vsd >> ws >>
				sws >> skew >> xnobs;
			
		}
		rratable1 >> header;
		
		if ((header == "T") || (header == "TABLE")) {
			pop_header = 1;
		}
		else if (header == "1") {

			rratable1 >> header;


			if (header == "T") {
				pop_header = 1;
			}
			else {
				rratable1.seekg(-8, std::ifstream::cur);
			}
		}
		else {
			rratable1.seekg(-8, std::ifstream::cur);
		}

		//Ignore any heights above maximum allowable RRA height
		if (zcurr > ztop){
			continue;
		}

		//Ignore any data lines with too few observed data points
		if (xnobs <= 10.0){
			continue;
		}

		//Ignore invalid data points
		if ((fabs(uwn) <= 0.0) && (fabs(vwn) <= 0.0)){
			continue;
		}

		//Set minimum allowed standard deviations
		if (usd < 0.01) {
			usd = 0.01;
		}

		if (vsd < 0.01){
			vsd = 0.01;
		}

		//If current height < previous height, a new month is being read
		if (zcurr < zprev){

			//Increment the month counter
			m1++;

			//Bypass loading RRA data arrays if month > 12 (or eof)
			if (m1 > 12){
				break;
			}
		}

		//Re-set previous height read to current height
		zprev = zcurr;

		//Load data into RRA arrays if month = month to be read (m)
		if (m1 == m) {

			//Increment array counter
			++num1;

			//Stop if too many data in RRA file
			if (num1 > maxnum){
				//hsSendMsg(TS_HS_FATAL, "ENV", "Too many data in RRA Table 1 file.");
			}

			//Load RRA array values
			z1[num1-1] = zcurr;
			u[num1-1] = float(uwn);
			su[num1-1] = float(usd);
			v[num1-1] = float(vwn);
			sv[num1-1] = float(vsd);
			ruvt[num1-1] = float(ruv);
			avspd[num1-1] = float(ws);
			sdspd[num1-1] = float(sws);
			
			

			if (fabs(1000.0*zcurr - rrasfc[isite - 1]) < 1.0){
				nsrf1[isite - 1] = num1;
			}
		}

		//Re-cycle to read next data line
		//Terminate if eof or > 12 months read
	}

	//Close the RRA input file (to re-use the RRA unit number)
	rratable1.close();

	return;

}

void RRA::readrra2(string &rra_id, double xlat, double xlon, int m, int isite)
{

	//readrra2 member function from RRA class
	//Reads Range Reference Atmosphere (RRA) pressure (p), temperature (t), density (d)
	//and standard deviations (sp, st, sd) for site at latitude = xlat, longitude = xlon
	//from RRA data file T2yyrra.txt.  Month is m.

	string lat, lon, ns, ew, header;
	stringstream convert_lat, convert_lon;
	double ztop, zlat, zlon, zlati = 0.0, zloni = 0.0, zprev, zcurr, pr, spr, skp,
		tmp, stmp, skt, dns, sdns, skd, xnp, xnt, xnd;
	int maxnum, month, count = 0, head_count = 0, pop_header = 1;


	ifstream rratable2;

	//Highest RRA data allowed is 70 km; max number = 300
	ztop = zmax[isite - 1];
	maxnum = 300;

	//Initalize RRA number of heights
	num2 = 0;

	//Open file for reading RRA data
	rratable2.open((rrapath.substr(0, lenpath) + "T2" + rra_id +
		".txt").c_str());

	//File open error check
	if (rratable2.is_open()){
	}
	else{
		cout << "File Open Error!  " << rrapath << '\n';
		system("pause");
		exit(1);
	}

	//Read RRA lat-lon on data file
	rratable2.seekg(17);
	rratable2 >> lat >> lon;

	convert_lat << lat.substr(0, lat.length() - 1) << " " << lat.substr(lat.length() - 1, lat.length())
		<< std::endl;
	convert_lon << lon.substr(0, lon.length() - 1) << " " << lon.substr(lon.length() - 1, lon.length())
		<< std::endl;

	convert_lat >> zlat >> ns;
	convert_lon >> zlon >> ew;

	if (ew.length() > 1) {
		ew.resize(1);
	}

	// Make north latitudes positive and south latitudes negative or 
	//write error message and stop otherwise
	if (ns == "N") {
		zlati = zlat;
	}
	else if (ns == "S") {
		zlati = -zlat;
	}
	else {
		cout << " Bad lat code " << xlat << "  " << zlati << "  " << zlat << "  " << ns << '\n';
		cout << " RRA 2 Bad lat code" << '\n';
		system("pause");
		exit(1);
	}

	// Make east longitudes positive and west longitudes negative or 
	//write error message and stop otherwise
	if (ew == "E") {
		zloni = zlon;
	}
	else if (ew == "W") {
		zloni = -zlon;
	}
	else {

		cout << " Bad lon code " << xlon << "  " << zloni << "  " << zlon << "  " << ew << '\n';
		cout << " RRA 2 Bad lon code " << '\n';
		system("pause");
		exit(1);
	}


	//Write error message and stop if RRA latitude from file does not agree 
	//with expected latitude
	if (fabs(xlat - zlati) > 0.005) {
		cout << "Bad latitude on input " << xlat << "  " << zlati << "  " << zlat << "  " << ns << '\n';
		cout << " RRA 2 Bad lat input" << '\n';
		system("pause");
		exit(1);
	}


	// Write error message and stop if RRA longitude from file does not agree with expected longitude
	if (fabs(xlon - zloni) > 0.005) {
		cout << "Bad longitude on input " << xlon << "  " << zloni << "  " << zlon << "  " << ew << '\n';
		cout << "RRA 2 Bad lon input" << '\n';
		system("pause");
		exit(1);
	}


	// Initialize month read from file and previous height read 
	//(first data in file are annual mean values, month = 0)

	if (rra_id.substr(3, 4) == "83") {
		month = 0;
	}
	else {
		month = 1;
	}

	zprev = -99.0;


	nsrf2[isite - 1] = 0;

	// Read a data line from the RRA file; ignore any lines with character data (or other condition
	// that causes a read error
	while (!rratable2.eof()) {

		++count;
		if (count > 5000) {
			break;
		}


		if (pop_header == 1) {
			rratable2 >> header;

			if (header == "G/M3") {
				++head_count;
			}

			if (head_count == 2) {
				head_count = 0;
				pop_header = 0;
			}

			continue;
		}


		if (pop_header == 0) {
			rratable2 >> zcurr >> pr >> spr >> skp >> tmp >> stmp >> skt >> dns
				>> sdns >> skd >> xnp >> xnt >> xnd;
		}

		rratable2 >> header;

		if ((header == "T") || (header == "TABLE")) {
			pop_header = 1;
		}
		else if (header == "1") {

			rratable2 >> header;

			if (header == "T") {
				pop_header = 1;
			}
			else {
				rratable2.seekg(-8, std::ifstream::cur);
			}
		}
		else {
			rratable2.seekg(-8, std::ifstream::cur);
		}

		// Ignore any heights above maximum allowable RRA height
		if (zcurr > ztop) {
			continue;
		}


		//Ignore any data lines with too few observed data points
		if ((xnp <= 10.0) || (xnt <= 10.0) || (xnd <= 10.0)){
			continue;
		}

		//Ignore invalid data points
		if ((tmp > 999.0) || (stmp > 99.0)){
			continue;
		}

		if ((pr <= 0.0) || (dns <= 0.0) || (tmp <= 0.0)){
			continue;
		}

		//Set minimum allowed standard deviations
		if (spr < 0.001){
			spr = 0.001;
		}

		if (stmp < 0.01){
			stmp = 0.01;
		}

		if (sdns < 0.001){
			sdns = 0.001;
		}


		//If current height < previous height, a new month is being read
		if (zcurr < zprev){
			//Increment the month counter
			++month;

			//Bypass loading RRA data arrays if month > 12 (or eof)

			if (month > 12){
				break;
			}
		}

		//Re-set previous height read to current height
		zprev = zcurr;

		//Load data into RRA arrays if month = month to be read (m)

		if (month == m){

			//Increment array counter
			++num2;

			//Insure pressure standard deviations are not too large
			if (spr > 0.3*pr) {
				cout << " Large RRA pressure sigma (reset) " << zcurr << "  "
					<< pr << "  " << spr << '\n';
				spr = 0.3*pr;
			}

			//Insure temperature standard deviations are not too large
			if (stmp > 0.3*tmp){
				cout << " Large RRA temp. sigma (reset) " << zcurr << "  "
					<< tmp << "  " << stmp << '\n';
				stmp = 0.3*tmp;
			}


			//Insure density standard deviations are not too large
			if (sdns > 0.3*dns){

				cout << " Large RRA density sigma (reset) " << zcurr << "  " << dns << "  " << sdns << '\n';

				sdns = 0.3*dns;

			}

			// Stop if too many data in RRA file
			if (num2 > maxnum) {
				cout << " Too many data in RRA Table 2 file." << '\n';
				system("pause");
				exit(1);
			}


			//Load RRA array values
			z2[num2 - 1] = zcurr;
			p[num2 - 1] = float(pr*100.0);
			sp[num2 - 1] = float(spr*100.0 / p[num2 - 1]);
			t[num2 - 1] = float(tmp);
			st[num2 - 1] = float(stmp / t[num2 - 1]);
			d[num2 - 1] = float(dns / 1000.0);
			sd[num2 - 1] = float(sdns / (1000.0*d[num2 - 1]));

			if (fabs(1000.0*zcurr - rrasfc[isite - 1]) < 1.0) {
				nsrf2[isite - 1] = num2;
			}

		}

		//Re-cycle to read next data line
		//Terminate if eof or > 12 months read

	}

	//Close the RRA input file (to reuse RRA unit number
	rratable2.close();

	return;

}


void RRA::readrra3(string &rra_id, double xlat, double xlon, int m, int isite)
{
	//readrra3 member function from RRA class
	//Read Range Reference Atmosphere (RRA) vapor pressure (vp), dewpoint 
	//temperature (td) and standard deviations (svp, std) for a site at latitude
	//= xlat, longitude = xlon from RRA data file T3yyrra.txt.  Month is m.

	string lat, lon, ns, ew, header;
	stringstream convert_lat, convert_lon;
	double ztop, zlat, zlon, zlati = 0.0, zloni = 0.0, zprev, zcurr, vpr, svpr, skvp, tvir,
		stvir, sktv, tdp, stdp, sktd, dum, xntp, xntt;
	int maxnum, month, count = 0, head_count = 0, pop_header = 1;

	ifstream rratable3;

	//Highest RRA data allowed is 70 km; max number = 300
	ztop = 30.0;
	maxnum = 300;

	//Initialize RRA number of heights
	num3 = 0;

	//Open file for reading RRA data
	rratable3.open((rrapath.substr(0, lenpath) + "T3" + rra_id + ".txt").c_str());

	//File open error check
	if (rratable3.is_open()){
	}
	else{
		cout << "File Open Error!  " << rrapath << '\n';
		system("pause");
		exit(1);
	}

	//Read RRA lat-lon on data file
	rratable3.seekg(17);
	rratable3 >> lat >> lon;

	convert_lat << lat.substr(0, lat.length() - 1) << " " << lat.substr(lat.length() - 1, lat.length())
		<< std::endl;
	convert_lon << lon.substr(0, lon.length() - 1) << " " << lon.substr(lon.length() - 1, lon.length())
		<< std::endl;

	convert_lat >> zlat >> ns;
	convert_lon >> zlon >> ew;

	if (ew.length() > 1) {
		ew.resize(1);
	}

	//Make north latitudes positive and south latitudes negative or write error
	//message to stop otherwise

	if (ns == "N"){
		zlati = zlat;
	}
	else if (ns == "S"){
		zlati = -zlat;
	}
	else {
		cout << " Bad lat code " << xlat << "  " << zlati << "  " << zlat << "  " << ns << '\n';
		cout << " RRA 3 Bad lat code" << '\n';
		system("pause");
		exit(1);
	}


	//Make east longitudes positive and west longitudes negative or write error
	//message and stop otherwise

	if (ew == "E"){
		zloni = zlon;
	}
	else if (ew == "W"){
		zloni = -zlon;
	}
	else {
		cout << " Bad lon code " << xlon << "  " << zloni << "  " << zlon << "  " << ew << '\n';
		cout << " RRA 3 Bad lon code " << '\n';
		system("pause");
		exit(1);
	}

	// Write error message and stop if RRA latitude from file does not agree with expected latitude
	if (fabs(xlat - zlati) > 0.005) {
		cout << "Bad latitude on input " << xlat << "  " << zlati << "  " << zlat << "  " << ns << '\n';
		cout << " RRA 3 Bad lat input" << '\n';
		system("pause");
		exit(1);
	}

	// Write error message and stop if RRA longitude from file does not agree with expected
	// longitude
	if (fabs(xlon - zloni) > 0.005) {
		cout << "Bad longitude on input " << xlon << "  " << zloni << "  " << zlon << "  " << ew << '\n';
		cout << "RRA 2 Bad lon input" << '\n';
		system("pause");
		exit(1);
	}


	//Initialize month read from file and previous height read (first data in 
	//file are annual mean values, month = 0)

	if (rra_id.substr(3, 4) == "83") {
		month = 0;
	}
	else {
		month = 1;
	}


	zprev = -99.0;


	nsrf3[isite - 1] = 0;



	//Read a data line from the RRA file; ignore any lines with character data
	//(or other condition that causes read error)

	while (!rratable3.eof()) {

		++count;
		if (count > 5000)	{
			break;
		}


		if (pop_header == 1){
			rratable3 >> header;

			if (header == "KELVIN"){
				++head_count;
			}

			if ((header == "DEG") || (header == "Deg")){
				rratable3 >> header;
				++head_count;
			}

			if (head_count == 4){
				head_count = 0;
				pop_header = 0;
			}

			continue;
		}


		if (pop_header == 0) {
			if ((rra_id == "fad83") || (rra_id == "nel83") || (rra_id == "shm83") ||
				(rra_id == "thu83") || (rra_id == "wak83")) {

				rratable3 >> zcurr >> vpr >> svpr >> skvp >> tvir >> stvir >>
					sktv >> tdp >> stdp >> sktd >> xntp >> xntt >> dum;
			}
			else{
				rratable3 >> zcurr >> vpr >> svpr >> skvp >> tvir >> stvir >>
					sktv >> tdp >> stdp >> sktd >> xntp >> xntt;
			}

			rratable3 >> header;

			if ((header == "T") || (header == "TABLE")){

				pop_header = 1;
			}
			else if (header == "1"){
				rratable3 >> header;

				if (header == "T"){
					pop_header = 1;
				}
				else {
					rratable3.seekg(-8, std::ifstream::cur);

				}
			}
			else {
				rratable3.seekg(-8, std::ifstream::cur);
			}

			//Ignore any heights above maximum allowable RRA height
			if (zcurr > ztop){
				continue;
			}

			//Ignore any data lines with invalid data points
			if ((tdp > 999.9) || (stdp > 99.9)){
				continue;
			}

			if ((vpr > 99.9) || (svpr > 99.9)){
				continue;
			}

			if ((tdp <= 0.0) || (vpr <= 0.0)){
				continue;
			}

			if ((xntp <= 10.0) || (xntt <= 10.0)){
				continue;
			}

			//Set minimum allowed standard deviations
			if (svpr < 0.001){
				svpr = 0.001;
			}

			if (stdp < 0.01){
				stdp = 0.01;
			}

			//If current height < previous height, a new month is being read
			if (zcurr < zprev){

				//Increment month counter
				++month;

				//Bypass loading RRA data arrays if month > 12 (or eof)
				if (month > 12){
					break;
				}

			}

			//Reset previous height read to current height
			zprev = zcurr;

			//Load data into RRA arrays if month = month to be read (m)
			if (month == m){

				//Increment array counter
				++num3;

				//Insure vapor pressure standard deviations are not too large
				if (svpr > 1.08*vpr){
					cout << "Large RRA vapor pressure sigma (reset)" << "  " <<
						zcurr << "  " << vpr << "  " << svpr << '\n';
					svpr = 1.08*vpr;
				}

				//Insure dewpoint temperature standard deviations are not too large
				if (stdp > 0.3*tdp){
					cout << "Large RRA dewpoint sigma (reset)" << "  " <<
						zcurr << "  " << tdp << "  " << stdp << '\n';
					stdp = 0.3*tdp;
				}

				//Stop if too many data in RRA file
				if (num3 > maxnum){
					cout << "Too many data in RRA Table 3 file." << '\n';
					system("pause");
					exit(1);
				}

				//Load RRA array values
				z3[num3 - 1] = zcurr;
				vp[num3 - 1] = float(vpr*100.0);
				svp[num3 - 1] = float(svpr*100.0);
				td[num3 - 1] = float(tdp);
				std[num3 - 1] = float(stdp);

				if (fabs(1000.0*zcurr - rrasfc[isite - 1]) < 1.0) {
					nsrf3[isite - 1] = num3;
				}
			}
			//Recycle to read next data line
			//Terminate if eof or >12 month read
		}

	}

	//Close the RRA input file (to reuse the RRA unit number)
	rratable3.close();

	return;
}


double RRA::heightwt(double h, double *z, int num)
{
	//heightwt member function from RRA class
	//Weighting factor for smooth altitude transition from Range Reference 
	//Atmosphere (RRA) data to GRAM data 

	double heightwt_out;

	//Height weight factor = 1 up to RRA height next to top, 0 if above top
	//RRA height, linear variation between top and next to top height

	if (h <= z[num - 1]){
		heightwt_out = 1.0;
	}
	else if (h >= z[num]){
		heightwt_out = 0.0;
	}
	else {
		heightwt_out = (z[num] - h) / (z[num] - z[num - 1]);
	}

	return heightwt_out;
}


int RRA::numfind(double h, double *z, int num)
{
	//numfind member function from RRA class
	//Find array index number (1-num) for Range Reference Atmosphere (RRA) 
	//height array (z) for input altitude (h)

	int i, numfind_out;

	numfind_out = 1;

	//Finds array index i for which z[i] < h < z[i+1]
	for (i = 1; i < num; i++){
		if ((h > z[i - 1]) && (h <= z[i])){
			numfind_out = i;
		}
	}

	return numfind_out;
}


double RRA::weightfact(double delta, double sitelim, double sitenear)
{
	//weightfact member function from RRA class
	//Horizontal weight factor for Range Reference Atmosphere (RRA) data.  
	//Inputs are delta, the latitude-longitude distance from location to RRA
	//site, sitelim, the maximum lat-lon distance allowed for adjustment to 
	//take place, and sitenear.  The lat-lon limit within which RRA data is 
	//used with full weight of 1.

	double pi, weightfact_out;


	pi = 4 * atan(1.0);

	if (delta < sitenear){
		//RRA weight factor = 1 if location is near RRA site
		weightfact_out = 1.0;
	}
	else if (delta > sitelim){
		//RRA weight factor = 0 if location outside sitelim
		weightfact_out = 0.0;
	}
	else {
		//RRA weight factor decrease to 0 (as cosine squared) between 'sitenear'
		//limit and distance = sitelim
		weightfact_out = pow(cos(pi*(delta - sitenear) / (2.0*(sitelim - sitenear))), 2);
	}

	return weightfact_out;
}


void RRA::rig(double ch, double phir, double *g0, double *g, double *re, double *r0, double *ri, double *a, double *b)
{

	//rig member function from RRA class
	//Computes surface gravity g0 (m/s^2) and effective Earth radius re (km) 
	//from input geocentric latitude phir (radians).  Also computes gravity 
	//g (m/s^2) and total radius ri (km) at input height ch (km).

	using namespace std;

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
}


double RRA::radll(double phi1, double thet1, double phi2, double thet2)
{
	//radll member function from radll class
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


void RRA::geodettogeocent(double fidet, double h, double *ficent, double *rtot,
	double *xy, double *z, double a, double b)
{
	//geodettogeocent member function from RRA class
	//Program to transform geodetic latitude, height to geocentric latitude,
	//radius.  Method from page K12 of Astronomical Almanac.  Input:  fidet,
	//height = geodetic latitude (deg), height (km);  a, b = equatorial and 
	//polar planetary radii (km).  Output:  ficent = geocentric lat (deg);
	//rtot = geocentric total radius (km);  xy, z = equatorial, polar cartesian
	//components (km)

	double omf, sphi, cphi, c;

	//Special case for the poles
	if (fabs(fidet) == 90.0){
		*ficent = copysign(90.0, fidet);
		*rtot = b + h;
		*xy = 0.0;
		*z = copysign(*rtot, fidet);
	}
	else {
		// 1 - flattening
		omf = b / a;

		//Sin and Cos of geodetic latitude
		sphi = sin(pi180*fidet);
		cphi = cos(pi180*fidet);

		//Computational factor c for cartesian coordinates
		c = 1.0 / sqrt(pow(cphi, 2) + pow(omf*sphi, 2));

		//Polar and equatorial cartesian coordinates
		*z = (a*c*pow(omf, 2) + h)*sphi;
		*xy = (a*c + h)*cphi;

		//Total geocentric radius
		*rtot = sqrt(pow(*xy, 2) + pow((*z), 2));

		//Geocentric latitude, deg
		*ficent = atan((*z) / (*xy)) / pi180;
	}

}

void RRA::rrasigs(double h, double phi, double thet, int mn, double hsrf1, double usrf1,
	double vsrf1, double tsrf1, double uvtsrf1, double susrf1, double svsrf1, double stsrf1,
	double spdavsrf1, double spdsdsrf1, double r0, double uvth, double *sph, double *sdh, double *sth, double *suh, double *svh,
	string &rra_id, double *rrawt)
{
	//rrasigs member function from RRA class
	//Based on a given height (h), geocentric latitude (phig), and longitude 
	//(thet) determine if GRAM values need to be adjusted by Range Reference 
	//Atmosphere (RRA) values (given location must be within latitude-longitude 
	//"radius" < sitelim for adjustment to be used).  If more than one RRA 
	//site is within adjustment range, the RRA site with minimum distance from 
	//the given location is used.  Computes modified atmospheric values, based 
	//on weighted average of original GRAM value and RRA value.

	int isite, i, n1, n2;
	double delmin, delsite, xlat, g0, g, re, rzero, rlocal, rpole,
		requa, gdlat, gdhgt, xlon, urra, vrra, surra, svrra, hgtwt, grmwt,
		sprra, sdrra,
		uvtrra, avspdrra, sdspdrra;

	//Find RRA site to read (if any), by looking for RRA site with minimum
	//deviation from current lat, lon.  

	//Set initial values of isite, rra_id, sitewgt, and delmin
	uvrra = uvth;
	isite = 0;
	rra_id = "  GRM";
	sitewgt = 0.0;
	delmin = 999.0;
	hgtwt = 0.0;

	//Step through the RRA sites
	for (i = 0; i < nrra; i++){

		//Compute great circle distance of given lat-lon from RRA site
		delsite = radll(phi, thet, rralat[i], rralon[i]);

		//Find minimum radius that is within lat-lon limit (if any)
		if ((delsite <= sitelim) && (delsite < delmin)){
			delmin = delsite;
			isite = i + 1;

			//RRA site name year & code
			rra_id = rrasite[i];

			//RRA site latitude and longitude
			xlat = rralat[i];

			//Get local Earth radius
			rig(0.0, xlat*pi180, &g0, &g, &re, &rzero, &rlocal, &rpole, &requa);

			//Get geodetic latitude
			geocenttogeodet(rlocal*cos(pi180*xlat), rlocal*sin(pi180*xlat),
				&gdlat, &gdhgt, rpole, requa);
			xlat = gdlat;
			xlon = rralon[i];

			//Horizontal weighting factor for RRA site
			sitewgt = weightfact(delsite, sitelim, sitenear);

		}

	}

	//Bypass RRA modification if all RRA sites are outside the lat-lon limit
	//(sitelim) range
	if (isite == 0){
		*rrawt = 0.0;
		return;
	}

	//Read RRA data unless this is the RRA site that was last read in
	if (isite != isitered){

		//Read Table 1 RRA data (wind components and standard deviations)
		readrra1(rra_id, xlat, xlon, mn, isite);

		//Read Table 2 RRA data (pressure, density, temperature, and their
		//standard deviations)
		readrra2(rra_id, xlat, xlon, mn, isite);

		//Read Table 3 RRA data (vapor pressure and dewpoint temperature and 
		//their standard deviations)
		readrra3(rra_id, xlat, xlon, mn, isite);

		//Store RRA site number of site just read in
		isitered = isite;

	}

	//Modify GRAM wind values if in height range
	if (h <= z1[num1 - 1]){

		//Find height index values for vertical interpolation
		n1 = numfind(h, z1, num1);
		n2 = n1 + 1;

		//Do vertical interpolation to get wind components
		interz(u[n1 - 1], v[n1 - 1], ruvt[n1 - 1], z1[n1 - 1], u[n2 - 1], v[n2 - 1], ruvt[n2 - 1],
			z1[n2 - 1], &urra, &vrra, &uvtrra, h);

		//Do vertical interpolation to get speed statistics
		interw(avspd[n1 - 1], sdspd[n1 - 1], z1[n1 - 1], avspd[n2 - 1], sdspd[n2 - 1],
			z1[n2 - 1], &avspdrra, &sdspdrra, h);

		if (h <= z1[0]){
			urra = u[0];
			vrra = v[0];
			uvrra = ruvt[0];
			avspdrra = avspd[0];
			sdspdrra = sdspd[0];
		}

		else {
			uvrra = uvtrra;
		}

		//Do vertical interpolation to get wind standard deviations
		interw(su[n1 - 1], sv[n1 - 1], z1[n1 - 1], su[n2 - 1], sv[n2 - 1],
			z1[n2 - 1], &surra, &svrra, h);

		if (h <= z1[0]){
			surra = su[0];
			svrra = sv[0];
		}

		else {
			surra = surra*rrascale;
			svrra = svrra*rrascale;
		}

		//Compute height weighting factor for RRA values
		hgtwt = heightwt(h, z1, num1 - 1);

		//Total RRA weighting factor = products of height weight and horizontal
		//weight factors
		*rrawt = sitewgt*hgtwt;
		rrawt1 = *rrawt;

		//GRAM weighting factor = 1 - RRA weight factor
		grmwt = 1.0 - *rrawt;

		//Revise GRAM winds and standard deviations by weighted average with 
		//RRA data
		*suh = *rrawt*surra + grmwt**suh;
		*svh = *rrawt*svrra + grmwt**svh;
		uvrra = *rrawt*uvrra + grmwt*uvth;
	}

	else {
		rrawt1 = 0.0;
	}

	if (h <= 27.0){
		*rrawt = sitewgt;
		grmwt = 1.0 - *rrawt;
		usrf = grmwt*usrf1 + *rrawt*u[0];
		vsrf = grmwt*vsrf1 + *rrawt*v[0];
		susrf = grmwt*susrf1 + *rrawt*su[0];
		svsrf = grmwt*svsrf1 + *rrawt*sv[0];
		rrsrf = grmwt*hsrf1 + *rrawt*z1[0];
		shsrf = grmwt*shsrf;
		uvtsrf = grmwt*uvtsrf1 + *rrawt*ruvt[0];
		spdavsrf = grmwt*spdavsrf1 + *rrawt*avspd[0];
		spdsdsrf = grmwt*spdsdsrf1 + *rrawt*sdspd[0];
	}


	//Modify GRAM pressure, density and temperature values with RRA values if in
	//height range
	if (h <= z2[num2 - 1]){

		//Find height index values for vertical interpolation
		n1 = numfind(h, z2, num2);
		n2 = n1 + 1;

		//Do vertical interpolation to get pressure, density, and temperature
		inter2(p[n1 - 1], d[n1 - 1], t[n1 - 1], z2[n1 - 1], p[n2 - 1], d[n2 - 1],
			t[n2 - 1], z2[n2 - 1], &prra, &drra, &trra, h);

		interz(sp[n1 - 1], sd[n1 - 1], st[n1 - 1], z2[n1 - 1], sp[n2 - 1],
			sd[n2 - 1], st[n2 - 1], z2[n2 - 1], &sprra, &sdrra, &strra, h);

		if (h <= z2[1]){
			sprra = sp[1] * rrascale;
			sdrra = sd[1] * rrascale;
			strra = st[1] * rrascale;
		}

		else {
			sprra = sprra*rrascale;
			sdrra = sdrra*rrascale;
			strra = strra*rrascale;
		}

		if ((rra_id == "dug83") || (rra_id == "wsm83") || (rra_id == "nel83")){
			if (h <= z2[2]){
				sprra = sp[2] * rrascale;
				sdrra = sd[2] * rrascale;
				strra = st[2] * rrascale;
			}
		}

		//Compute height weighting factor for RRA values
		hgtwt = heightwt(h, z2, num2 - 1);

		//Total RRA weighting factor = product of height weight and horizontal
		//weight factors
		*rrawt = sitewgt*hgtwt;
		rrawt2 = *rrawt;
		//GRAM weight factor = 1 - RRA weight factor
		grmwt = 1.0 - *rrawt;

		//Revise GRAM pressure, density, temperature and standard deviations
		//by weighted average with RRA data
		*sph = *rrawt*sprra + grmwt**sph;
		*sdh = *rrawt*sdrra + grmwt**sdh;
		*sth = *rrawt*strra + grmwt**sth;

		//Compute height weighting factor for surface RRA values
		if (h <= 27.0){
			*rrawt = sitewgt;
			grmwt = 1.0 - *rrawt;
			if ((psrf > 0.0) & (dsrf > 0.0)){
				spsrf = grmwt*spsrf / psrf + *rrawt*sp[nsrf2[isite - 1]];
				sdsrf = grmwt*sdsrf / dsrf + *rrawt*sd[nsrf2[isite - 1]];
				psrf = exp(grmwt*log(psrf)) + *rrawt*log(p[nsrf2[isite - 1]]);
				dsrf = exp(grmwt*log(dsrf)) + *rrawt*log(d[nsrf2[isite - 1]]);
				spsrf = spsrf*psrf;
				sdsrf = sdsrf*dsrf;
			}
			else{
				psrf = 0.0;
				dsrf = 0.0;
				spsrf = 0.0;
				sdsrf = 0.0;
			}
			if (rra_id == "dug83" || rra_id == "wsm83" || rra_id == "nel83"){
				tsrf = grmwt*tsrf1 + *rrawt*t[nsrf2[isite] + 2];
				strrasrf = st[nsrf2[isite] + 1] * t[nsrf2[isite] + 2];
				stsrf = grmwt*stsrf1 + *rrawt*strrasrf;
			}
			else{
				tsrf = grmwt*tsrf1 + *rrawt*t[nsrf2[isite] + 1];
				strrasrf = st[nsrf2[isite] + 1] * t[nsrf2[isite] + 1];
				stsrf = grmwt*stsrf1 + *rrawt*strrasrf;
			}
		}

	}
	else {
		rrawt2 = 0.0;
	}

	//Set rra_id to " GRM", weight to 0 if no RRA adjustment done
	*rrawt = sitewgt*hgtwt;

	if ((h >= z1[num1 - 1]) || (h >= z2[num2 - 1])){
		rra_id = " GRM";
		*rrawt = 0.0;
	}
	return;
}

void RRA::rramean(double h, double phi, double thet, int mn, double r0,
	double *pgh, double *dgh, double *tgh, double *ugh, double *vgh,
	string &rra_id, double *rrawt)
{
	//rramean member function from RRA class
	//Based on a given height (h), geocentric latitude (phig), and longitude (thet)
	//determine if GRAM values need to be adjusted by Range Reference Atmosphere
	//(RRA) values (given location must be within latitude-longitude "radius"
	//< sitelim for adjustment to be used).  If more than one RRA site is within
	//adjustment range, the RRA site with minimum distance from the given 
	//location is used.  Computes modified atmospheric values, based on weighted
	//average of original GRAM value and RRA value.

	int isite, i, n1, n2;
	double delmin, delsite, xlat, g0, g, re, rzero, rlocal, rpole,
		requa, gdlat, gdhgt, xlon, urra, vrra, hgtwt, grmwt, hgtwt3;


	//Find RRA site to read (if any), by looking for RRA site with minimum
	//deviation from current lat, lon.  

	//Set initial values of isite, rra_id, sitewgt, and delmin
	isite = 0;
	rra_id = "  GRM";
	sitewgt = 0.0;
	delmin = 999.0;
	hgtwt = 0.0;

	//Step through the RRA sites
	for (i = 0; i < nrra; i++){

		//Compute great circle distance of given lat-lon from RRA site
		delsite = radll(phi, thet, rralat[i], rralon[i]);

		//Find minimum radius that is within lat-lon limit (if any)
		if ((delsite <= sitelim) && (delsite < delmin)){
			delmin = delsite;
			isite = i + 1;

			//RRA site name year & code
			rra_id = rrasite[i];

			//RRA site latitude and longitude
			xlat = rralat[i];

			//Get local Earth radius
			rig(0.0, xlat*pi180, &g0, &g, &re, &rzero, &rlocal, &rpole, &requa);

			//Get geodetic latitude
			geocenttogeodet(rlocal*cos(pi180*xlat), rlocal*sin(pi180*xlat), &gdlat,
				&gdhgt, rpole, requa);

			xlat = gdlat;
			xlon = rralon[i];

			//Horizontal weighting factor for RRA site
			sitewgt = weightfact(delsite, sitelim, sitenear);

		}

	}

	//Bypass RRA modification if all RRA sites are outside the lat-lon limit
	//(sitelim) range
	if (isite == 0){
		*rrawt = 0.0;
		return;
	}

	//Read RRA data unless this is the RRA site that was last read in
	if (isite != isitered){

		//Read Table 1 RRA data (wind components and standard deviations)
		readrra1(rra_id, xlat, xlon, mn, isite);

		//Read Table 2 RRA data (pressure, density, temperature, and their
		//standard deviations)
		readrra2(rra_id, xlat, xlon, mn, isite);

		//Read Table 3 RRA data (vapor pressure and dewpoint temperature and 
		//their standard deviations)
		readrra3(rra_id, xlat, xlon, mn, isite);

		//Store RRA site number of site just read in
		isitered = isite;

	}


	//Modify GRAM wind values if in height range
	if (h <= z1[num1 - 1]){

		//Find height index values for vertical interpolation
		n1 = numfind(h, z1, num1);
		n2 = n1 + 1;
		
		//Do vertical interpolation to get wind components
		interw(u[n1 - 1], v[n1 - 1], z1[n1 - 1], u[n2 - 1], v[n2 - 1], z1[n2 - 1],
			&urra, &vrra, h);


		//Do vertical interpolation to get wind statistics
		interw(avspd[n1 - 1], sdspd[n1 - 1], z1[n1 - 1], avspd[n2 - 1], sdspd[n2 - 1],
			z1[n2 - 1], &avspdrra, &sdspdrra, h);

		if (h <= z1[0]){
			urra = u[0];
			vrra = v[0];
			avspdrra = avspd[0];
			
		}
		
		if (h <= z1[0]){
			sdspdrra = sdspd[0];
		}

		else {
			sdspdrra = sdspdrra;
		}

		//Compute height weighting factor for RRA values
		hgtwt = heightwt(h, z1, num1 - 1);

		//Total RRA weighting factor = products of height weight and horizontal
		//weight factors
		*rrawt = sitewgt*hgtwt;
		rrawt1 = *rrawt;


		//GRAM weighting factor = 1 - RRA weight factor
		grmwt = 1.0 - *rrawt;

		//Revise GRAM winds and standard deviations by weighted average with 
		//RRA data
		*ugh = *rrawt*urra + grmwt**ugh;
		*vgh = *rrawt*vrra + grmwt**vgh;

	}

	else {
		rrawt1 = 0.0;
	}

	//Modify GRAM pressure, density and temperature values with RRA values if in
	//height range

	if (h <= z2[num2 - 1]){

		//Find height index values for vertical interpolation
		n1 = numfind(h, z2, num2);
		n2 = n1 + 1;

		//Do vertical interpolation to get pressure, density, and temperature
		inter2(p[n1 - 1], d[n1 - 1], t[n1 - 1], z2[n1 - 1], p[n2 - 1], d[n2 - 1],
			t[n2 - 1], z2[n2 - 1], &prra, &drra, &trra, h);


		//Compute height weighting factor for RRA values
		hgtwt = heightwt(h, z2, num2 - 1);

		//Total RRA weighting factor = product of height weight and horizontal
		//weight factors
		*rrawt = sitewgt*hgtwt;
		rrawt2 = *rrawt;
		//GRAM weight factor = 1 - RRA weight factor
		grmwt = 1.0 - *rrawt;

		//Revise GRAM pressure, density, temperature and standard deviations
		//by weighted average with RRA data
		*pgh = *rrawt*prra + grmwt**pgh;
		*dgh = *rrawt*drra + grmwt**dgh;
		*tgh = *rrawt*trra + grmwt**tgh;

	}

	else {
		rrawt2 = 0.0;
	}
	//Modify GRAM vapor pressure and dewpoint temperature values with RRA values
	//if in height range

	if (h <= z3[num3 - 1]){

		//Find height index values for vertical interpolation
		n1 = numfind(h, z3, num3);
		n2 = n1 + 1;

		//Do vertical interpolation to get vapor pressure and dewpoint
		interw(vp[n1 - 1], td[n1 - 1], z3[n1 - 1], vp[n2 - 1], td[n2 - 1],
			z3[n2 - 1], &vprra, &tdrra, h);

		//Do vertical interpolations to get standard deviations
		interw(svp[n1 - 1], std[n1 - 1], z3[n1 - 1], svp[n2 - 1], std[n2 - 1],
			z3[n2 - 1], &svprra, &stdrra, h);


		if (h <= z3[nsrf3[isite] + 1]){
			svprra = svp[nsrf3[isite] + 1];
			stdrra = std[nsrf3[isite] + 1];
		}

		//Compute height weighting factor for RRA values
		hgtwt3 = heightwt(h, z3, num3 - 1);
		rrawt3 = sitewgt*hgtwt3;

		//Total RRA weighting factor = product of height weight and horizontal
		//weight factors

		//GRAM weight factor = 1 - RRA weight factor
		grmwt = 1.0 - rrawt3;

	}

	else {
		rrawt3 = 0.0;
	}

	if (h <= 27.0){
		tdsrf = td[nsrf3[isite] + 1];
		stdsrf = std[nsrf3[isite] + 1];
	}

	//Set rra_id to " GRM", weight to 0 if no RRA adjustment done
	*rrawt = sitewgt*hgtwt;

	if ((h >= z1[num1 - 1]) || (h >= z2[num2 - 1])){
		rra_id = " GRM";
		*rrawt = 0.0;
	}

	return;

}

void RRA::inter2(double p1, double d1, double t1, double z1, double p2, double d2,
	double t2, double z2, double *p, double *d, double *t, double z)
{
	//inter2 member function from RRA class
	//Interpolates between p1,d1,t1 at height z1 to p2,d2,t2 at height z2 to 
	//output values of p,d,t and height z.  Checks t1,d1,t2,d2 product=0, for 
	//gas constant interpolation.
	double a, tz, b, pb, pz, r1, r2, r;

	//Set p=d=t=0 if some input values are negative or zero
	if (p1*d1*t1*p2*d2*t2 <= 0.0 || p2*p1 <= 0.0){
		*p = 0.0;
		*d = 0.0;
		*t = 0.0;
		return;
	}
	//Sets p,d,t = p1,d1,t1 if z1 = z2
	if (fabs(z1 - z2) <= 0.001){
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
	if (fabs(t2 - t1) <= 0.001){
		b = (z2 - z1) / log(p1 / p2);
		*p = p1*exp((z1 - z) / b);
		*d = *p / (r*tz);
		*t = tz;
	}
	//Power law interpolation on pressure
	else {
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


void RRA::interz(double p1, double d1, double t1, double z1, double p2, double d2,
	double t2, double z2, double *p, double *d, double *t, double z)
{
	//interz member function from RRA class
	//Linear interpolation between p1,d1,t1 at height z1 and p2,d2,t2 at 
	//height z1 to output values of p,d,t at height z.

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


void RRA::interw(double u1, double v1, double z1, double u2, double v2, double z2, double *u, double *v, double z)
{
	//interw member function from RRA class
	//Linear interpolation between u1,v1 at z1 and u2,v2 at height z1.  Output
	//is u,v at height z.

	double a;

	if (fabs(z1 - z2) <= 0.001) {
		*u = u1;
		*v = v1;
	}
	else{
		a = (z - z1) / (z2 - z1);
		*u = u1 + (u2 - u1)*a;
		*v = v1 + (v2 - v1)*a;
	}
}


void RRA::geocenttogeodet(double r, double zin, double *fi, double *h, double a,
	double b)
{
	//geocenttogeodet member function from RRA class
	//Program to transform Cartesian to geodetic coordinates.  Code adapted from 
	//Polish version of K. M.  Borkowski, Bull. Geod. vol 63, pp. 50-56
	//Input:  r, zin = equitorial and polar Cartesian components (km); a, b = 
	//equitorial and polar planetary radii (km)
	//Output:  fi, h, geodetic coord's (latitude (deg), height (km))


	double z, rabs, zlim, e, f, p, q, d, s, v, g, t, yfi, xfi, srfi;

	z = fabs(zin);
	rabs = fabs(r);
	zlim = 1.0e-6;

	if (rabs < zlim*z){
		*fi = 90.0 - atan(rabs / z) / pi180;
		*h = z - b;
	}


	else {

		//Analytical solution for non-polar case
		//See also page K12 of Astronomical Almanac for iterative solution
		e = ((z + b)*b / a - a) / rabs;
		f = ((z - b)*b / a + a) / rabs;
		p = (e*f + 1.0)*4.0 / 3.0;
		q = (e*e - f*f) * 2;
		d = p*p*p + q*q;

		if (d >= 0.0) {
			s = sqrt(d) + q;
			s = copysign(exp(log(fabs(s)) / 3.0), s);
			v = p / s - s;
			v = -(q + q + v*v*v) / (3.0*p);
		}
		else {
			v = 2.0*sqrt(-p)*cos(acos(q / p / sqrt(-p)) / 3.0);
		}

		g = 0.5*(e + sqrt(e*e + v));
		t = sqrt(g*g + (f - v*g) / (g + g - e)) - g;
		yfi = (1.0 - t*t)*a;
		xfi = 2.0*b*t;

		if (t <= zlim){
			*fi = 90.0 - atan(rabs / z) / pi180;
			*h = z - b;
		}
		else {
			srfi = sqrt(xfi*xfi + yfi*yfi);
			*h = ((rabs - a*t)*xfi + (z - b)*yfi / srfi);
			*fi = atan(yfi / xfi) / pi180;
		}

	}

	if (zin < 0.0){
		*fi = -*fi;
	}

	return;

}
