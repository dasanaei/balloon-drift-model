//Init Class for initializing data for the model prior to runs.
//P. White
#include "Init.h"
#include <iostream>
#include <cmath>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <iomanip>


using namespace std;

Init::Init()
{
	initializeMemberVariables();
}


Init::~Init()
{



}


void Init::initializeMemberVariables()
{

	//Initialize Init class member variables

	psp = new float**[15];
	dsp = new float**[15];
	tsp = new float**[15];
	usp = new float**[15];
	vsp = new float**[15];
	for (int i = 0; i < 15; ++i){
		psp[i] = new float*[19];
		dsp[i] = new float*[19];
		tsp[i] = new float*[19];
		usp[i] = new float*[19];
		vsp[i] = new float*[19];
		for (int j = 0; j < 19; ++j){
			psp[i][j] = new float[18];
			dsp[i][j] = new float[18];
			tsp[i][j] = new float[18];
			usp[i][j] = new float[18];
			vsp[i][j] = new float[18];
		}
	}

	for (int i = 0; i < 15; ++i) {
		for (int j = 0; j < 19; ++j) {
			for (int k = 0; k < 18; ++k) {
				psp[i][j][k] = 0.0;     // Pa    Stationary perturbation data for pressure
				dsp[i][j][k] = 0.0;     // kg/M3 Stationary perturbation data for density
				tsp[i][j][k] = 0.0;     // K     Stationary perturbation data for temperature
				usp[i][j][k] = 0.0;     // M/s   Stationary perturbation data for Eastward wind
				vsp[i][j][k] = 0.0;     // M/s   Stationary perturbation data for Northward wind
			}
		}
	}


	pg = new float*[21];
	dg = new float*[21];
	tg = new float*[21];
	ug = new float*[21];
	for (int i = 0; i < 21; ++i){
		pg[i] = new float[19];
		dg[i] = new float[19];
		tg[i] = new float[19];
		ug[i] = new float[19];
	}


	for (int i = 0; i < 21; ++i) {
		for (int j = 0; j < 19; ++j) {
			pg[i][j] = 0.0;             // Pa    Zonal pressure data
			dg[i][j] = 0.0;             // kg/M3 Zonal density data
			tg[i][j] = 0.0;             // K     Zonal mean temperature data
			ug[i][j] = 0.0;             // M/s   Zonal average wind data
		}
	}

	pr = new float*[29];
	dr = new float*[29];
	tr = new float*[29];
	ur = new float*[29];
	vr = new float*[29];
	for (int i = 0; i < 29; ++i){
		pr[i] = new float[19];
		dr[i] = new float[19];
		tr[i] = new float[19];
		ur[i] = new float[19];
		vr[i] = new float[19];
	}

	for (int i = 0; i < 29; i++){
		for (int j = 0; j < 19; j++){
			pr[i][j] = 0.0;
			dr[i][j] = 0.0;
			tr[i][j] = 0.0;
			ur[i][j] = 0.0;
			vr[i][j] = 0.0;
		}
	}

	dlp = new float*[25];
	plp = new float*[25];
	tlp = new float*[25];
	ulp = new float*[25];
	vlp = new float*[25];
	uds = new float*[25];
	vds = new float*[25];
	udl = new float*[25];
	uvt = new float*[25];
	for (int i = 0; i < 25; ++i){
		dlp[i] = new float[10];
		plp[i] = new float[10];
		tlp[i] = new float[10];
		ulp[i] = new float[10];
		vlp[i] = new float[10];
		uds[i] = new float[10];
		vds[i] = new float[10];
		udl[i] = new float[10];
		uvt[i] = new float[10];
	}

	for (int i = 0; i < 25; i++){
		for (int j = 0; j < 10; j++){
			dlp[i][j] = 0.0;
			plp[i][j] = 0.0;
			tlp[i][j] = 0.0;
			ulp[i][j] = 0.0;
			vlp[i][j] = 0.0;
			uds[i][j] = 0.0;
			vds[i][j] = 0.0;
			udl[i][j] = 0.0;
			uvt[i][j] = 0.0;
		}
	}

	for (int i = 0; i < 29; i++){
		xlbar[i] = 0.0;
		zlbar[i] = 0.0;
		xsigl[i] = 0.0;
		zsigl[i] = 0.0;
		xlmin[i] = 0.0;
		zlmin[i] = 0.0;
		xscale[i] = 0.0;
		zscale[i] = 0.0;
		wr[i] = 0.0;
	}

	h2ol = new float*[10];
	for (int i = 0; i < 10; ++i){
		h2ol[i] = new float[35];
	}

	for (int i = 0; i < 10; ++i) {
		for (int j = 0; j < 35; ++j) {
			h2ol[i][j] = 0.0;       // ppmv Atmospheric constituent of water
			if (i == 0) {
				hgtl[j] = 0.0;
			}
		}
	}

	for (int i = 0; i < 50; i++){
		hgta[i] = 0.0;
		co2[i] = 330.0;
		o2[i] = 209470.0;
		n2[i] = 780830.0;
	}

	co2[41] = 328.0, co2[42] = 320.0, co2[43] = 310.0, co2[44] = 270.0,
		co2[45] = 195.0, co2[46] = 110.0, co2[47] = 60.0, co2[48] = 40.0,
		co2[49] = 35.0;

	o2[42] = 200000.0, o2[43] = 190000.0, o2[44] = 180000.0,
		o2[45] = 160000.0, o2[46] = 140000.0, o2[47] = 120000.0,
		o2[48] = 94000.0, o2[49] = 72500.0;


	n2[43] = 780000.0, n2[44] = 779000.0, n2[45] = 777000.0,
		n2[46] = 774000.0, n2[47] = 770000.0, n2[48] = 765000.0,
		n2[49] = 760000.0;




	h2oa = new float*[8];
	o3 = new float*[8];
	n2o = new float*[8];
	co = new float*[8];
	ch4 = new float*[8];
	for (int i = 0; i < 8; ++i){
		h2oa[i] = new float[50];
		o3[i] = new float[50];
		n2o[i] = new float[50];
		co[i] = new float[50];
		ch4[i] = new float[50];
	}


	for (int i = 0; i < 8; i++){
		for (int j = 0; j < 50; j++){
			h2oa[i][j] = 0.0;
			o3[i][j] = 0.0;
			n2o[i][j] = 0.0;
			co[i][j] = 0.0;
			ch4[i][j] = 0.0;
		}
	}

	o3map = new float*[19];
	for (int i = 0; i < 19; ++i){
		o3map[i] = new float[24];
	}


	h2omap = new float*[8];
	sh2omap = new float*[8];
	for (int i = 0; i < 8; ++i){
		h2omap[i] = new float[19];
		sh2omap[i] = new float[19];
	}

	n2omap = new float*[19];
	ch4map = new float*[19];
	for (int i = 0; i < 19; ++i){
		n2omap[i] = new float[17];
		ch4map[i] = new float[17];
	}

	oxymap = new float*[19];
	for (int i = 0; i < 19; ++i){
		oxymap[i] = new float[19];
	}


	for (int i = 0; i < 19; i++){
		for (int j = 0; j < 24; j++){
			if (i == 0){
				po3map[j] = 0.0;
			}
			o3map[i][j] = 0.0;
		}
	}

	for (int i = 0; i < 8; i++){
		for (int j = 0; j < 19; j++){
			if (i == 0){
				ph2omap[j] = 0.0;
			}
			h2omap[i][j] = 0.0;
			sh2omap[i][j] = 0.0;
		}
	}

	for (int i = 0; i < 19; i++){
		for (int j = 0; j < 17; j++){
			if (i == 0){
				pmap31[j] = 0.0;
			}
			n2omap[i][j] = 0.0;
			ch4map[i][j] = 0.0;
		}
	}

	for (int i = 0; i < 19; i++){
		for (int j = 0; j < 19; j++){
			oxymap[i][j] = 0.0;
			if (i == 0){
				mapzox[j] = 0;
			}
		}
	}



	water = 0.0, ozone = 0.0, nitrous = 0.0, carbmon = 0.0, methane = 0.0,
		carbdiox = 0.0, moloxy = 0.0, nitrogen = 0.0, sigwater = 0.0,
		wtmol = 0.0, ond = 0.0;

	ztopo = new int*[360];
	for (int i = 0; i < 360; ++i){
		ztopo[i] = new int[180];
	}

	landcd = new int*[360];
	for (int i = 0; i < 360; ++i){
		landcd[i] = new int[180];
	}


	for (int i = 0; i < 360; i++){
		for (int j = 0; j < 180; j++){
			ztopo[i][j] = 0;
			landcd[i][j] = 0;
		}
	}

	f = 0;
	dz1 = 0.0;
	hgt1 = 0.0, dphi1 = 0.0, dthet1 = 0.0, lat1 = 0.0, lon1 = 0.0;
	time = 0.0, time1 = 0.0, delt1 = 0.0, dtime1 = 0.0;

	return;

}

void Init::rtran1(ifstream &atmosdat, string *i1, int *i2, int *i3, int *ndata, int *i5)
{

	//rtran1 member function from Init class
	//Member function reads the "atmosdat" data file

	int i;

	atmosdat >> *i1 >> *i2 >> *i3;

	for (i = 0; i < 19; i++){
		atmosdat >> ndata[i];
	}

	atmosdat >> *i5;

}

void Init::rtran2(ifstream &atmosdat, string *i1, int *i2, int *i3, int *i4, int *ndata1){

	//rtran2 member function from Init class
	//Member function reads the "atmosdat" data file

	int i;
	atmosdat >> *i1 >> *i2 >> *i3 >> *i4;
	for (i = 0; i < 18; i++){
		atmosdat >> ndata1[i];
	}

}

void Init::rtran3(ifstream &atmosdat, string *i1, int *i2, int *i3, int *ndata2)
{

	//rtran3 member function from Init class
	//Member function reads the "atmosdat" data file

	int i;

	atmosdat >> *i1 >> *i2 >> *i3;
	for (i = 0; i < 15; i++){
		atmosdat >> ndata2[i];
	}

}

void Init::rtran4(ifstream &atmosdat, string *i1, int *i2, int *i3, int *ndata3)
{

	//rtran4 member function from Init class
	//Member function reads the "atmosdat" data file


	int i;

	atmosdat >> *i1 >> *i2 >> *i3;
	for (i = 0; i < 10; i++){
		atmosdat >> ndata3[i];
	}

}

void Init::init1(int mn)
{

	//init1 member function from Init class
	//Loads in the data for the specific month form the "atmosdat" data file

	using namespace std;
	int n, j, i, ish, ip[5], id[5], it[5], k;
	int i2, i3, i4, ihra, ndata1[18], ndata2[15], ndata3[10];
	int ndata[19], i5;
	string i1;
	float tenx;
	ifstream atmosdat;
	int lenpath1;

	float rpscale = 1.0;
	string rscode;
	double zin;

	lenpath1 = atmpath.find(" ");

	//Open the atmosdat data file
	atmosdat.open((atmpath.substr(0, lenpath1) + "atmosdat_E10.txt").c_str());
	
	if (atmosdat.is_open()){
	}
	else {
		cout << "File Open Error!  " << atmpath << '\n';
		system("pause");
		exit(1);
	}

	for (n = 0; n < 252; n++){

		rtran1(atmosdat, &i1, &i2, &i3, ndata, &i5);
		//Read zonal pressure data
		if (i2 == mn){
			tenx = powf(10.0, i5);
			i3 = (i3 - 15) / 5;

			for (i = 0; i < 19; i++){
				pg[i3 - 1][i] = ndata[i] * tenx;
			}
		}

	}

	for (n = 0; n < 252; n++){

		rtran1(atmosdat, &i1, &i2, &i3, ndata, &i5);
		//Reads zonal density data
		if (i2 == mn){
			tenx = powf(10.0, i5);
			i3 = (i3 - 15) / 5;

			for (i = 0; i < 19; i++){
				dg[i3 - 1][i] = ndata[i] * tenx;
			}
		}

	}
	
	for (n = 0; n < 252; n++){

		rtran1(atmosdat, &i1, &i2, &i3, ndata, &i5);
		//Reads zonal mean temperature data
		if (i2 == mn){
			tenx = powf(10.0, i5);
			i3 = (i3 - 15) / 5;

			for (i = 0; i < 19; i++){
				tg[i3 - 1][i] = ndata[i] * tenx;
			}
		}

	}

	for (n = 0; n < 252; n++){

		rtran1(atmosdat, &i1, &i2, &i3, ndata, &i5);
		//Reads zonal avg. zonal wind data
		if (i2 == mn){
			tenx = powf(10.0, i5);
			i3 = (i3 - 15) / 5;
			for (i = 0; i < 19; i++){
				ug[i3 - 1][i] = ndata[i] * tenx;
			}
		}

	}

	for (n = 0; n < 3060; n++){

		rtran2(atmosdat, &i1, &i2, &i3, &i4, ndata1);
		//Reads stationary perturbation data for pressure (to be stored in psp array)
		if (i2 == mn){
			ish = (i3 - 15) / 5;
			k = (i4 + 100) / 10;
			for (i = 0; i < 18; i++){
				psp[ish - 1][k - 1][i] = float(ndata1[i] / 1000.0);
			}
		}
	}

	for (n = 0; n < 3060; n++){

		rtran2(atmosdat, &i1, &i2, &i3, &i4, ndata1);
		//Reads stationary perturbation data for density (to be stored in dsp array)
		if (i2 == mn){
			ish = (i3 - 15) / 5;
			k = (i4 + 100) / 10;
			for (i = 0; i < 18; i++){
				dsp[ish - 1][k - 1][i] = float(ndata1[i] / 1000.0);
			}
		}
	}

	for (n = 0; n < 3060; n++){

		rtran2(atmosdat, &i1, &i2, &i3, &i4, ndata1);
		//Reads stationary perturbation data for temperature (to be stored in tsp array)
		if (i2 == mn){
			ish = (i3 - 15) / 5;
			k = (i4 + 100) / 10;
			for (i = 0; i < 18; i++){
				tsp[ish - 1][k - 1][i] = float(ndata1[i] / 1000.0);
			}
		}
	}

	for (n = 0; n < 3060; n++){

		rtran2(atmosdat, &i1, &i2, &i3, &i4, ndata1);
		//Reads stationary perturbation data for zonal wind (to be stored in usp array)
		if (i2 == mn){
			ish = (i3 - 15) / 5;
			k = (i4 + 100) / 10;
			for (i = 0; i < 18; i++){
				usp[ish - 1][k - 1][i] = float(ndata1[i] / 10.0);
			}
		}
	}

	for (n = 0; n < 3060; n++){

		rtran2(atmosdat, &i1, &i2, &i3, &i4, ndata1);
		//Reads stationary perturbation data for meridional wind (to be stored in vsp array)
		if (i2 == mn){
			ish = (i3 - 15) / 5;
			k = (i4 + 100) / 10;
			for (i = 0; i < 18; i++){
				vsp[ish - 1][k - 1][i] = float(ndata1[i] / 10.0);
			}
		}
	}


	for (i = 0; i < 300; i++){
		//Reads random perturbations in pressure
		rtran1(atmosdat, &i1, &i2, &i3, ndata, &i5);
		if (i2 == mn){
			ihra = (i3 + 5) / 5;
			if (i3 > 120){
				ihra = 19 + i3 / 20;
			}
			for (k = 0; k < 19; k++){
				pr[ihra - 1][k] = powf(ndata[k] * rpscale / 1000, 2);
			}

		}
	}

	for (i = 0; i < 300; i++){
		//Reads random perturbations in density
		rtran1(atmosdat, &i1, &i2, &i3, ndata, &i5);
		if (i2 == mn){
			ihra = (i3 + 5) / 5;
			if (i3 > 120){
				ihra = 19 + i3 / 20;
			}
			for (k = 0; k < 19; k++){
				dr[ihra - 1][k] = powf(ndata[k] * rpscale / 1000, 2);
			}

		}
	}

	for (i = 0; i < 300; i++){
		//Reads random perturbations in temperature
		rtran1(atmosdat, &i1, &i2, &i3, ndata, &i5);
		if (i2 == mn){
			ihra = (i3 + 5) / 5;
			if (i3 > 120){
				ihra = 19 + i3 / 20;
			}
			for (k = 0; k < 19; k++){
				tr[ihra - 1][k] = powf(ndata[k] * rpscale / 1000, 2);
			}

		}
	}

	for (i = 0; i < 348; i++){
		//Reads random perturbations in zonal wind
		rtran1(atmosdat, &i1, &i2, &i3, ndata, &i5);
		if (i2 == mn){
			ihra = (i3 + 5) / 5;
			if (i3 > 120){
				ihra = 19 + i3 / 20;
			}
			for (k = 0; k < 19; k++){
				ur[ihra - 1][k] = powf(ndata[k] * rpscale / 10, 2);
			}

		}
	}

	for (i = 0; i < 348; i++){
		//Reads random perturbations in meridional wind
		rtran1(atmosdat, &i1, &i2, &i3, ndata, &i5);
		if (i2 == mn){
			ihra = (i3 + 5) / 5;
			if (i3 > 120){
				ihra = 19 + i3 / 20;
			}
			for (k = 0; k < 19; k++){
				vr[ihra - 1][k] = powf(ndata[k] * rpscale / 10, 2);
			}

		}
	}

	for (i = 0; i < 25; i++){
		rtran3(atmosdat, &i1, &i2, &i3, ndata2);

		//split into 3 smaller arrays
		for (k = 0; k < 5; k++){
			ip[k] = ndata2[k];
			id[k] = ndata2[5 + k];
			it[k] = ndata2[10 + k];
		}

		for (j = 0; j < 5; j++){
			plp[i][j + 5] = float(ip[j] / 1000.0);
			plp[i][6 - (j + 2)] = float(ip[j] / 1000.0);
			dlp[i][j + 5] = float(id[j] / 1000.0);
			dlp[i][6 - (j + 2)] = float(id[j] / 1000.0);
			tlp[i][j + 5] = float(it[j] / 1000.0);
			tlp[i][6 - (j + 2)] = float(it[j] / 1000.0);
		}

	}

	for (i = 0; i < 25; i++){
		rtran4(atmosdat, &i1, &i2, &i3, ndata3);

		//split into 2 small arrays
		for (k = 0; k < 5; k++){
			ip[k] = ndata3[k];
			id[k] = ndata3[k + 5];
		}

		//assume large-scale fraction for u equal that for v
		for (j = 0; j < 5; j++){
			ulp[i][j + 5] = float(id[j] / 1000.0);
			ulp[i][6 - (j + 2)] = float(id[j] / 1000.0);
			vlp[i][j + 5] = float(id[j] / 1000.0);
			vlp[i][6 - (j + 2)] = float(id[j] / 1000.0);
		}
	}

	for (i = 0; i < 25; i++){
		rtran4(atmosdat, &i1, &i2, &i3, ndata3);
		//split into 2 small arrays
		for (k = 0; k < 5; k++){
			ip[k] = ndata3[k];
			id[k] = ndata3[k + 5];
		}


		for (j = 0; j < 5; j++){
			uds[i][j + 5] = float(ip[j] / 1000.0);
			uds[i][6 - (j + 2)] = float(ip[j] / 1000.0);
			vds[i][j + 5] = float(id[j] / 1000.0);
			vds[i][6 - (j + 2)] = float(id[j] / 1000.0);
		}
	}


	for (i = 0; i < 25; i++){
		rtran4(atmosdat, &i1, &i2, &i3, ndata3);
		//split into 2 small arrays
		for (k = 0; k < 5; k++){
			ip[k] = ndata3[k];
			id[k] = ndata3[k + 5];
		}


		for (j = 0; j < 5; j++){
			udl[i][j + 5] = float(ip[j] / 1000.0);
			udl[i][6 - (j + 2)] = float(ip[j] / 1000.0);
			uvt[i][j + 5] = float(id[j] / 1000.0);
			uvt[i][6 - (j + 2)] = float(id[j] / 1000.0);
		}
	}


	for (i = 0; i < 29; ++i) {

		atmosdat >> rscode >> zin
			>> xlbar[i] >> xsigl[i] >> xlmin[i] >> xscale[i]
			>> zlbar[i] >> zsigl[i] >> zlmin[i] >> zscale[i]
			>> wr[i];

		wr[i] = rwscale*wr[i];

	}


	//Initialize concentrations data
	concinit(atmosdat, mn);
	mapinit(atmosdat, mn);

	//Store topography data in array
	topo();

	atmosdat.close();


}

void Init::concinit(ifstream &atmosdat, int imon)
{

	//concinit member function from Init class
	//Initializes concentrations data by store atmosdat data into arrays

	double const alph = 8.0, bet = 5.0, bet2 = bet + 2.0, betoh = bet + 0.5, alph3q = alph + 0.75,
		alphoh = alph + 1.5, qalph = 0.25*alph + 0.5, tqalph = 0.75*alph + 0.5, bet1 = bet + 1.0,
		alph1 = alph + 1.0, alph1h = alph1 / 2.0;



	int j, l, m, i, imonp, mp, k;
	double xl[10], xa[5];
	string lcode, afcode;
	string lseas[4] = { "LDJF", "LMAM", "LJJA", "LSON" };
	string acode[4] = { "AFMS", "AFMW", "AFSS", "AFSW" };

	double facl[4][12] = {
		{ bet2, betoh, 0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.5, betoh },
		{ -0.5, 0.5, betoh, bet2, betoh, 0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, -0.5, 0.5, betoh, bet2, betoh, 0.5, -0.5, 0.0, 0.0 },
		{ -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.5, betoh, bet2, betoh, 0.5 }
	};

	double faca[2][12] = {

		{ -0.5, 0.25, qalph, alph1h, tqalph, alph3q, alphoh, alph3q, tqalph, alph1h,
		qalph, 0.25 },
		{ alphoh, alph3q, tqalph, alph1h, qalph, 0.25, -0.5, 0.25, qalph, alph1h,
		tqalph, alph3q }

	};

	imonp = imon + 5;

	if (imonp > 11){
		imonp = imonp - 11;
	}


	//Loop through the 4 season of LaRC H2O data
	for (i = 0; i < 4; i++){
		//Compute monthly weighting factors for LaRC H2O data
		for (j = 0; j < 12; j++){
			facl[i][j] = float(facl[i][j] / bet1);
		}

		//Read the LaRC H2O data (lat = -70 to 70)
		for (k = 0; k < 35; k++){
			atmosdat >> lcode >> hgtl[k];
			for (j = 1; j < 9; j++){
				atmosdat >> xl[j];
			}

			//Compute lat = -90 and lat = +90
			xl[0] = exp((4.0*log(xl[1]) - log(xl[2])) / 3.0);
			xl[9] = exp((4.0*log(xl[8]) - log(xl[7])) / 3.0);

			for (l = 0; l < 10; l++){
				h2ol[l][k] = h2ol[l][k] + float(facl[i][imon - 1] * xl[l]);
			}
		}
	}


	//Compute monthly weighting factors for AFGL data
	for (i = 0; i < 2; i++){
		for (j = 0; j < 12; j++){
			faca[i][j] = faca[i][j] / alph1;
		}
	}

	//Read the tropical (lat = 15) AFGL concentration data
	for (k = 0; k < 50; k++){
		atmosdat >> afcode >> hgta[k] >> xa[0] >> xa[1] >> xa[2] >> xa[3] >> xa[4];

		//Load data as both +15 and -15 latitude
		for (l = 3; l < 5; l++){
			h2oa[l][k] = float(xa[0]);
			o3[l][k] = float(xa[1]);
			n2o[l][k] = float(xa[2]);
			co[l][k] = float(xa[3]);
			ch4[l][k] = float(xa[4]);
		}
	}

	//Loop through the 2 latitudes of AFGL concentration data
	for (m = 5; m < 7; m++){
		mp = 7 - m;

		//Loop through the 2 seasons of AFGL concentration data
		for (i = 0; i < 2; i++){
			//Read the AFGL data (lat 45 & 60, summer/winter)
			for (k = 0; k < 50; k++){
				atmosdat >> afcode >> hgta[k] >> xa[0] >> xa[1] >> xa[2] >> xa[3] >> xa[4];

				//Sum the weighted AFGL values for the month given Northern hemisphere data
				h2oa[m][k] = h2oa[m][k] + float(xa[0] * faca[i][imon - 1]);
				o3[m][k] = o3[m][k] + float(xa[1] * faca[i][imon - 1]);
				n2o[m][k] = n2o[m][k] + float(xa[2] * faca[i][imon - 1]);
				co[m][k] = co[m][k] + float(xa[3] * faca[i][imon - 1]);
				ch4[m][k] = ch4[m][k] + float(xa[4] * faca[i][imon - 1]);

				//Southern hemisphere data (6 months displaced)
				h2oa[mp][k] = h2oa[mp][k] + float(xa[0] * faca[i][imonp]);
				o3[mp][k] = o3[mp][k] + float(xa[1] * faca[i][imonp]);
				n2o[mp][k] = n2o[mp][k] + float(xa[2] * faca[i][imonp]);
				co[mp][k] = co[mp][k] + float(xa[3] * faca[i][imonp]);
				ch4[mp][k] = ch4[mp][k] + float(xa[4] * faca[i][imonp]);

			}
		}
	}

	//Compute lat +90 and lat -90 AFGL concentration values 
	for (k = 0; k < 50; k++){
		h2oa[7][k] = float(exp((9.0*log(h2oa[6][k]) - 4.0*log(h2oa[5][k])) / 5.0));
		o3[7][k] = float(exp((9.0*log(o3[6][k]) - 4.0*log(o3[5][k])) / 5.0));
		n2o[7][k] = float(exp((9.0*log(n2o[6][k]) - 4.0*log(n2o[5][k])) / 5.0));
		co[7][k] = float(exp((9.0*log(co[6][k]) - 4.0*log(co[5][k])) / 5.0));
		ch4[7][k] = float(exp((9.0*log(ch4[6][k]) - 4.0*log(ch4[5][k])) / 5.0));

		h2oa[0][k] = float(exp((9.0*log(h2oa[1][k]) - 4.0*log(h2oa[2][k])) / 5.0));
		o3[0][k] = float(exp((9.0*log(o3[1][k]) - 4.0*log(o3[2][k])) / 5.0));
		n2o[0][k] = float(exp((9.0*log(n2o[1][k]) - 4.0*log(n2o[2][k])) / 5.0));
		co[0][k] = float(exp((9.0*log(co[1][k]) - 4.0*log(co[2][k])) / 5.0));
		ch4[0][k] = float(exp((9.0*log(ch4[1][k]) - 4.0*log(ch4[2][k])) / 5.0));
	}




}


void Init::mapinit(ifstream &atmosdat, int imon)
{

	//mapinit member functio from Init class
	//Open and read the map species data from the atmosdat file and load into arrays

	//Local Variables
	int i, j, k, nox[17], nmap[15], m, iexp;

	double xwv[5], swv[5];

	const double ten = 10.0;

	string ocode, code;

	//Read the map ozone data at latitudes -80 to +80
	for (i = 0; i < 12; i++){
		for (j = 0; j < 24; j++){

			atmosdat >> ocode >> m >> po3map[j] >> nox[0] >> nox[1] >> nox[2] >> nox[3] >>
				nox[4] >> nox[5] >> nox[6] >> nox[7] >> nox[8] >> nox[9] >> nox[10] >>
				nox[11] >> nox[12] >> nox[13] >> nox[14] >> nox[15] >> nox[16] >> iexp;

			if (i + 1 == imon){
				for (k = 1; k < 18; k++){
					o3map[k][j] = float(nox[k - 1] * pow(ten, iexp));
				}

				//Fill in lat = -90 and +90
				o3map[0][j] = float(exp((4.0*log(o3map[1][j]) - log(o3map[2][j])) / 3.0));
				o3map[18][j] = float(exp((4.0*log(o3map[17][j]) - log(o3map[16][j])) / 3.0));
			}
		}
	}

	//Read the map h2o data at latitudes -60, -45, +/-15, +45, +60 for pressure levels 1.5 to 100 mb

	for (i = 0; i < 12; i++){
		for (j = 8; j < 19; j++){
			atmosdat >> code >> m >> ph2omap[j] >> xwv[0] >> xwv[1] >> xwv[2] >> xwv[3] >>
				xwv[4] >> swv[0] >> swv[1] >> swv[2] >> swv[3] >> swv[4];

			if (i + 1 == imon){
				//Fill in the h2o and sigma h2o arrays
				for (k = 0; k < 2; k++){
					h2omap[1 + k][j] = float(xwv[k]);
					sh2omap[1 + k][j] = float(swv[k]);
					h2omap[5 + k][j] = float(xwv[3 + k]);
					sh2omap[5 + k][j] = float(swv[3 + k]);
					h2omap[3 + k][j] = float(xwv[2]);
					sh2omap[3 + k][j] = float(swv[2]);
				}
			}

		}
	}


	//Read the annual map h2o data
	for (j = 0; j < 8; j++){
		atmosdat >> code >> m >> ph2omap[j] >> xwv[0] >> xwv[1] >> xwv[2] >> xwv[3] >>
			xwv[4] >> swv[0] >> swv[1] >> swv[2] >> swv[3] >> swv[4];

		//Fill in the h2o and sigma h2o arrays
		for (k = 0; k < 2; k++){
			h2omap[1 + k][j] = float(xwv[k]);
			sh2omap[1 + k][j] = float(swv[k]);
			h2omap[5 + k][j] = float(xwv[3 + k]);
			sh2omap[5 + k][j] = float(swv[3 + k]);
			h2omap[3 + k][j] = float(xwv[2]);
			sh2omap[3 + k][j] = float(swv[2]);
		}
	}

	//Fill in the lat = -90 and +90 values
	for (j = 0; j < 19; j++){
		h2omap[0][j] = float(exp((9.0*log(h2omap[1][j]) - 4.0*log(h2omap[2][j])) / 5.0));
		sh2omap[0][j] = sh2omap[1][j] * h2omap[0][j] / h2omap[1][j];
		h2omap[7][j] = float(exp((9.0*log(h2omap[6][j]) - 4.0*log(h2omap[5][j])) / 5.0));
		sh2omap[7][j] = sh2omap[6][j] * h2omap[7][j] / h2omap[6][j];
	}



	//Read the map n2o data at latitudes -70 to +70
	for (i = 0; i < 12; i++){
		for (j = 0; j < 17; j++){
			atmosdat >> code >> m >> pmap31[j] >> nmap[0] >> nmap[1] >> nmap[2] >>
				nmap[3] >> nmap[4] >> nmap[5] >> nmap[6] >> nmap[7] >> nmap[8] >>
				nmap[9] >> nmap[10] >> nmap[11] >> nmap[12] >> nmap[13] >>
				nmap[14] >> iexp;

			if (i + 1 == imon){
				for (k = 2; k < 17; k++){
					n2omap[k][j] = float((nmap[k - 2])*pow(ten, iexp));
				}
			}

			//Fill in lat = -90, -80, +80, +90
			n2omap[0][j] = float(exp((9.0*log(n2omap[2][j]) - 4.0*log(n2omap[3][j])) / 5.0));
			n2omap[1][j] = float(exp((8.0*log(n2omap[2][j]) - 3.0*log(n2omap[3][j])) / 5.0));
			n2omap[18][j] = float(exp((9.0*log(n2omap[16][j]) - 4.0*log(n2omap[15][j])) / 5.0));
			n2omap[17][j] = float(exp((8.0*log(n2omap[16][j]) - 3.0*log(n2omap[15][j])) / 5.0));
		}
	}

	//Read the map ch4 data latitudes -70 to +70
	for (i = 0; i < 12; i++){
		for (j = 0; j < 17; j++){
			atmosdat >> code >> m >> pmap31[j] >> nmap[0] >> nmap[1] >> nmap[2] >> nmap[3] >>
				nmap[4] >> nmap[5] >> nmap[6] >> nmap[7] >> nmap[8] >> nmap[9] >> nmap[10] >>
				nmap[11] >> nmap[12] >> nmap[13] >> nmap[14] >> iexp;

			if (i + 1 == imon){
				for (k = 2; k < 17; k++){
					ch4map[k][j] = float((nmap[k - 2]) * pow(ten, iexp));
				}
			}

			//Fill in lat = -90, -80, +80, and +90
			ch4map[0][j] = float(exp((9.0*log(ch4map[2][j]) - 4.0*log(ch4map[3][j])) / 5.0));
			ch4map[1][j] = float(exp((8.0*log(ch4map[2][j]) - 3.0*log(ch4map[3][j])) / 5.0));
			ch4map[18][j] = float(exp((9.0*log(ch4map[16][j]) - 4.0*log(ch4map[15][j])) / 5.0));
			ch4map[17][j] = float(exp((8.0*log(ch4map[16][j]) - 3.0*log(ch4map[15][j])) / 5.0));
		}
	}

	//Read the map ox data at latitudes -80 to +80
	for (i = 0; i < 12; i++){
		for (j = 0; j < 19; j++){
			atmosdat >> ocode >> m >> mapzox[j] >> nox[0] >> nox[1] >> nox[2] >> nox[3] >>
				nox[4] >> nox[5] >> nox[6] >> nox[7] >> nox[8] >> nox[9] >> nox[10] >>
				nox[11] >> nox[12] >> nox[13] >> nox[14] >> nox[15] >> nox[16] >> iexp;

			if (i + 1 == imon){
				for (k = 0; k < 18; k++){
					oxymap[k][j] = float(nox[k - 1] * pow(ten, iexp));
				}
			}

			//Fill in lat = -90 and +90
			oxymap[0][j] = float(exp((4.0*log(oxymap[1][j]) - log(oxymap[2][j])) / 3.0));
			oxymap[18][j] = float(exp((4.0*log(oxymap[17][j]) - log(oxymap[16][j])) / 3.0));
		}
	}

	//Convert pressure levels to log (base e) N/m^2 scale
	for (j = 0; j < 19; j++){
		ph2omap[j] = float(log(100.0*ph2omap[j]));
	}

	for (j = 0; j < 24; j++){
		po3map[j] = float(log(100.0*po3map[j]));
	}

	for (j = 0; j < 17; j++){
		pmap31[j] = float(log(100.0*pmap31[j]));
	}





}


void Init::topo(){


	//topo member function from init class
	//Reads topography and land type data and stores in arrays

	ifstream tops;

	int j, i, ltopom;
	double xlo, xla;
	int lenpath2;
	lenpath2 = atmpath.find(" ");

	//Open topography file
	tops.open((atmpath.substr(0, lenpath2) + "topo.txt").c_str());

	//File open error check
	if (tops.is_open()){
	}
	else {
		cout << "File Open Error" << atmpath << '\n';
	}

	//Read topography and land code data and store in arrays
	for (j = 179; j >= 0; j--){
		for (i = 0; i < 360; i++){
			tops >> xlo >> xla >> landcd[i][j] >> ltopom;
			ztopo[i][j] = ltopom;

		}
	}



	tops.close();


}


void Init::traj()
{

	//traj member function from Init class
	//Reads trajectory dat

	string code;


	//read time, height, lat, lon
	traj1 >> time >> hgt >> lat >> lon;
	
	//Calculate deltas
	dz1 = (hgt - hgt1);
	dphi1 = (lat - lat1);
	dthet1 = (lon - lon1);
	delt1 = (time - time1);

	
	hgt1 = hgt;
	lon1 = lon;
	lat1 = lat;
	time1 = time;
	if (traj1.eof()){
		f = 0;
		traj1.close();
	}
	else {
		f = 1;
	}

}


void Init::trajopen()
{

	//trjopen member function from Init class
	//Open trajectory file

	string trjp = home + trapath;

	traj1.open(trjp.c_str());
	if (traj1.is_open()){
	}
	else {
		cout << "File Open Error!  " << home + trapath << '\n';
		system("pause");
		exit(1);
	}
	out << "Positions generated from trajectory file:" << '\n';
	out << trapath << '\n';

}


void Init::outopen()
{

	//outopen member function from Init class
	//Opens output.txt file and write info to file

	string outp = home + prtpath;

	out.open(outp.c_str());

	out << " **** Earth Global Reference Atmospheric Model - 2016 (Earth-GRAM 2016) ****" << '\n';
	out << " MM/DD/YYYY = " << mn << "/" << ida << "/" << iyr << "HH:MM:SS(UTC) = " << ihro << ":" << ":" << mino << ":" << seco << '\n';
	out << "F10.7 = " << f10 << " Mean F10.7 = " << f10b << " ap Index = " << ap << '\n';
	out << "Range Reference Atmosphere (RRA) data availave for " << Ryear << '\n';
	out << "NCEP Global Climatology Data:  POR = " << NCEPyr << " NCEPhr = " << NCEPhr << '\n';
	out << "Path = " << NCEPpath << '\n';
	out << "Thermospheric conditions from " << therm << " model" << '\n';
	out << "1st Random No. = " << nr1 << " Random Scale Factors = " << rpscale << " " << ruscale << " " << rwscale << '\n';
	out << "Patchy Turbulence Option = " << patch << '\n';
	out << '\n';
	out << " Mean-76 and Total-76 are percent deviations from 1976 US Standard Atmosphere." << '\n';
	out << " Other deviations in percent are with respect to mean values.  RH is relative" << '\n';
	out << " humidity in percent.  Zeroes for H2O indicate no estimate available." << '\n';
	out << " E-W wind positive toward East; N-S wind positive toward North." << '\n';
	out << '\n';



}


void Init::namelist()
{

	//namelist member function from Init class
	//Open and read namelist input file for setting model parameters.

	double iday1[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
	ifstream namelist;
	string dummy, NCEPmn;
	string mn1 = "12";

	

	NCEPmn = "Nb9008" + mn1 + ".bin";

	cout << "Enter directory:" << '\n';
	getline(cin, home);

	cout << "Enter filename:" << '\n';
	getline(cin, namelst);

	string namel = home + namelst;
	namelist.open(namel.c_str());


	string iblo = home + "BLTest.txt";

	//File open error check
	if (namelist.is_open()){
	}
	else {
		cout << "File Open Error!  " << namel << '\n';
		system("pause");
		exit(1);
	}

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

	namelist >> dummy >> dummy >> iyrra;

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

	if (ibltest > 1){
		ibl.open(iblo.c_str());
		ibl << "   Day    LST    Hgt    Lat      Lon   hsrf spdsrf"
			" LC     z0   Elmn   El   Elmd     sha blfct  nri   S"
			"       ool  ustar    BVfsq    hN   hbl   chb  sigrat"
			"  swb   swh spdavsrf sdpsdsrf   tsrf  stsrf" << '\n';
	}

	//Check bounds on z0in value
	if (z0in > 0.0){
		if (z0in < 1.0e-05) z0in = 1.0e-05;
		if (z0in > 3.0) z0in = 3.0;
	}

	//Check sitelim and sitenear
	if (sitelim <= 0.0 || sitenear >= sitelim){
		cout << "Bad RRA site limits (sitenear must be < sitelim)  " << sitenear <<
			"  " << sitelim << '\n';
		system("pause");
		exit(1);
	}

	//Terminate if month out of range
	if (mn < 1 || mn > 12){
		cout << "Month out of range.  " << mn << '\n';
		system("pause");
		exit(1);
	}


	//Terminate if day out of range
	if ((ida > iday1[mn - 1]) || (ida < 0)){
		cout << "Bad day input.  " << ida << '\n';
		system("pause");
		exit(1);
	}

	//Terminate initial hour if out of range
	if ((ihro > 23) || (ihro < 0)){
		cout << "Bad hour input.  " << ihro << '\n';
		system("pause");
		exit(1);
	}

	//Terminate if initial minute out of range
	if ((mino > 60) || (mino < 0)){
		cout << "Bad minute input.  " << mino << '\n';
		system("pause");
		exit(1);
	}

	//Terminate if initial second out of range
	if ((seco > 60.0) || (seco < 0.0)){
		cout << "Bad second input.  " << seco << '\n';
		system("pause");
		exit(1);
	}

	//Terminate if out of range NCEPhr
	if (NCEPhr < 0 || NCEPhr > 5){
		cout << "Bad NCEPhr value." << '\n';
		system("pause");
		exit(1);
	}
	//Terminate if auxiliary profile and RRA data turned on
	if (iaux > 0 && iurra > 0){
		cout << "Can not use RRA data and auxiliary profile data simultaneously."
			<< '\n';
		system("pause");
		exit(1);
	}


	if (itherm == 1) therm = "MET";
	if (itherm == 2) therm = "MSIS";
	if (itherm == 3) therm = "JB2008";

	if (patchy > 0) patch = "On";
	if (patchy == 0) patch = "Off";

	if (iyrra == 1) Ryear = "1983";
	if (iyrra == 2) Ryear = "2006";
	if (iyrra == 3) Ryear = "2013";

	outopen();

}

