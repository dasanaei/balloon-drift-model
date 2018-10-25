//Class for reading and processing NCEP data, lower atmosphere data
//P. White

#include "NCEP.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

NCEPmods::NCEPmods()
{

	initializeMemberVariables();

}

NCEPmods::~NCEPmods()
{

	for (int i = 0; i < 19; i++){

		for (int j = 0; j < 73; j++){
			delete[] temp[i][j];
			delete[] dens[i][j];
			delete[] dewp[i][j];
			delete[] uwnd[i][j];
			delete[] vwnd[i][j];
			delete[] geop[i][j];
			delete[] stmp[i][j];
			delete[] sden[i][j];
			delete[] sdewp[i][j];
			delete[] suwd[i][j];
			delete[] svwd[i][j];
			delete[] sprs[i][j];
			delete[] rhum[i][j];
			delete[] srhum[i][j];
			delete[] vprs[i][j];
			delete[] svprs[i][j];
			delete[] spdav[i][j];
			delete[] spdsd[i][j];
			delete[] uvcor[i][j];
		}
		delete[] temp[i];
		delete[] dens[i];
		delete[] dewp[i];
		delete[] uwnd[i];
		delete[] vwnd[i];
		delete[] geop[i];
		delete[] stmp[i];
		delete[] sden[i];
		delete[] sdewp[i];
		delete[] suwd[i];
		delete[] svwd[i];
		delete[] sprs[i];
		delete[] rhum[i];
		delete[] srhum[i];
		delete[] vprs[i];
		delete[] svprs[i];
		delete[] spdav[i];
		delete[] spdsd[i];
		delete[] uvcor[i];
	}

	for (int j = 0; j < 73; j++){
		delete[] slpn[j];
		delete[] sfcp[j];
		delete[] sslp[j];
		delete[] ssfcp[j];
	}

	delete[] temp;
	delete[] dens;
	delete[] dewp;
	delete[] uwnd;
	delete[] vwnd;
	delete[] geop;
	delete[] slpn;
	delete[] sfcp;
	delete[] stmp;
	delete[] sden;
	delete[] sdewp;
	delete[] suwd;
	delete[] svwd;
	delete[] sprs;
	delete[] sslp;
	delete[] ssfcp;
	delete[] rhum;
	delete[] srhum;
	delete[] vprs;
	delete[] svprs;
	delete[] spdav;
	delete[] spdsd;
	delete[] uvcor;

}

void NCEPmods::initializeMemberVariables()
{


	//Initialize member variables for NCEPmods class

	int i, j, k;


	temp = new float**[19];
	for (i = 0; i < 19; ++i){
		temp[i] = new float*[73];
		for (j = 0; j < 73; ++j){
			temp[i][j] = new float[145];
		}
	}



	dens = new float**[19];
	for (i = 0; i < 19; ++i){
		dens[i] = new float*[73];
		for (j = 0; j < 73; ++j){
			dens[i][j] = new float[145];
		}
	}

	dewp = new float**[19];
	for (i = 0; i < 19; ++i){
		dewp[i] = new float*[73];
		for (j = 0; j < 73; ++j){
			dewp[i][j] = new float[145];
		}
	}

	uwnd = new float**[19];
	for (i = 0; i < 19; ++i){
		uwnd[i] = new float*[73];
		for (j = 0; j < 73; ++j){
			uwnd[i][j] = new float[145];
		}
	}

	vwnd = new float**[19];
	for (i = 0; i < 19; ++i){
		vwnd[i] = new float*[73];
		for (j = 0; j < 73; ++j){
			vwnd[i][j] = new float[145];
		}
	}

	geop = new float**[19];
	for (i = 0; i < 19; ++i){
		geop[i] = new float*[73];
		for (j = 0; j < 73; ++j){
			geop[i][j] = new float[145];
		}
	}

	slpn = new float*[73];
	for (i = 0; i < 73; ++i){
		slpn[i] = new float[145];
	}

	sfcp = new float*[73];
	for (i = 0; i < 73; ++i){
		sfcp[i] = new float[145];
	}

	stmp = new float**[19];
	for (i = 0; i < 19; ++i){
		stmp[i] = new float*[73];
		for (j = 0; j < 73; ++j){
			stmp[i][j] = new float[145];
		}
	}

	sden = new float**[19];
	for (i = 0; i < 19; ++i){
		sden[i] = new float*[73];
		for (j = 0; j < 73; ++j){
			sden[i][j] = new float[145];
		}
	}

	sdewp = new float**[19];
	for (i = 0; i < 19; ++i){
		sdewp[i] = new float*[73];
		for (j = 0; j < 73; ++j){
			sdewp[i][j] = new float[145];
		}
	}

	suwd = new float**[19];
	for (i = 0; i < 19; ++i){
		suwd[i] = new float*[73];
		for (j = 0; j < 73; ++j){
			suwd[i][j] = new float[145];
		}
	}

	svwd = new float**[19];
	for (i = 0; i < 19; ++i){
		svwd[i] = new float*[73];
		for (j = 0; j < 73; ++j){
			svwd[i][j] = new float[145];
		}
	}

	sprs = new float**[19];
	for (i = 0; i < 19; ++i){
		sprs[i] = new float*[73];
		for (j = 0; j < 73; ++j){
			sprs[i][j] = new float[145];
		}
	}

	sslp = new float*[73];
	for (i = 0; i < 73; ++i){
		sslp[i] = new float[145];
	}

	ssfcp = new float*[73];
	for (i = 0; i < 73; ++i){
		ssfcp[i] = new float[145];
	}

	rhum = new float**[19];
	for (i = 0; i < 19; ++i){
		rhum[i] = new float*[73];
		for (j = 0; j < 73; ++j){
			rhum[i][j] = new float[145];
		}
	}

	srhum = new float**[19];
	for (i = 0; i < 19; ++i){
		srhum[i] = new float*[73];
		for (j = 0; j < 73; ++j){
			srhum[i][j] = new float[145];
		}
	}

	vprs = new float**[19];
	for (i = 0; i < 19; ++i){
		vprs[i] = new float*[73];
		for (j = 0; j < 73; ++j){
			vprs[i][j] = new float[145];
		}
	}

	svprs = new float**[19];
	for (i = 0; i < 19; ++i){
		svprs[i] = new float*[73];
		for (j = 0; j < 73; ++j){
			svprs[i][j] = new float[145];
		}
	}

	spdav = new float**[19];
	for (i = 0; i < 19; ++i){
		spdav[i] = new float*[73];
		for (j = 0; j < 73; ++j){
			spdav[i][j] = new float[145];
		}
	}

	spdsd = new float**[19];
	for (i = 0; i < 19; ++i){
		spdsd[i] = new float*[73];
		for (j = 0; j < 73; ++j){
			spdsd[i][j] = new float[145];
		}
	}

	uvcor = new float**[19];
	for (i = 0; i < 19; ++i){
		uvcor[i] = new float*[73];
		for (j = 0; j < 73; ++j){
			uvcor[i][j] = new float[145];
		}
	}

	for (i = 0; i < 19; ++i){
		for (j = 0; j < 73; ++j){
			for (k = 0; k < 145; ++k){
				temp[i][j][k] = 0.0;
				dens[i][j][k] = 0.0;
				dewp[i][j][k] = 0.0;
				uwnd[i][j][k] = 0.0;
				vwnd[i][j][k] = 0.0;
				geop[i][j][k] = 0.0;
				stmp[i][j][k] = 0.0;
				sden[i][j][k] = 0.0;
				sdewp[i][j][k] = 0.0;
				suwd[i][j][k] = 0.0;
				svwd[i][j][k] = 0.0;
				sprs[i][j][k] = 0.0;
				rhum[i][j][k] = 0.0;
				srhum[i][j][k] = 0.0;
				vprs[i][j][k] = 0.0;
				svprs[i][j][k] = 0.0;
				spdav[i][j][k] = 0.0;
				spdsd[i][j][k] = 0.0;
				uvcor[i][j][k] = 0.0;
			}
		}
	}



	for (j = 0; j < 73; ++j){
		for (k = 0; k < 145; ++k){
			slpn[j][k] = 0.0;
			sfcp[j][k] = 0.0;
			sslp[j][k] = 0.0;
			ssfcp[j][k] = 0.0;
		}
	}


	t = new float*[4];
	p = new float*[4];
	h = new float*[4];
	td = new float*[4];
	u = new float*[4];
	v = new float*[4];
	w = new float*[4];
	st = new float*[4];
	std = new float*[4];
	su = new float*[4];
	sv = new float*[4];
	sp = new float*[4];
	rho = new float*[4];
	srho = new float*[4];
	vp = new float*[4];
	svp = new float*[4];
	spdavl = new float*[4];
	spdsdl = new float*[4];
	uvcorrl = new float*[4];
	rh = new float*[4];
	srh = new float*[4];
	for (i = 0; i < 4; ++i){
		t[i] = new float[19];
		p[i] = new float[19];
		h[i] = new float[19];
		td[i] = new float[19];
		u[i] = new float[19];
		v[i] = new float[19];
		w[i] = new float[19];
		st[i] = new float[19];
		std[i] = new float[19];
		su[i] = new float[19];
		sv[i] = new float[19];
		sp[i] = new float[19];
		rho[i] = new float[19];
		srho[i] = new float[19];
		vp[i] = new float[19];
		svp[i] = new float[19];
		spdavl[i] = new float[19];
		spdsdl[i] = new float[19];
		uvcorrl[i] = new float[19];
		rh[i] = new float[19];
		srh[i] = new float[19];
	}

	for (i = 0; i < 4; i++){
		for (j = 0; j < 19; j++){
			t[i][j] = 0.0;
			p[i][j] = 0.0;
			h[i][j] = 0.0;
			td[i][j] = 0.0;
			u[i][j] = 0.0;
			v[i][j] = 0.0;
			w[i][j] = 0.0;
			st[i][j] = 0.0;
			std[i][j] = 0.0;
			rho[i][j] = 0.0;
			srho[i][j] = 0.0;
			su[i][j] = 0.0;
			sv[i][j] = 0.0;
			sp[i][j] = 0.0;
			vp[i][j] = 0.0;
			svp[i][j] = 0.0;
			spdavl[i][j] = 0.0;
			spdsdl[i][j] = 0.0;
			uvcorrl[i][j] = 0.0;
			rh[i][j] = 0.0;
			srh[i][j] = 0.0;

		}
	}


	pi = 3.1415926535897931;
	pi180 = pi / 180, dtz = 0.0;

}


void NCEPmods::NCEPread()
{

	//NCEPread member function from NCEPmods class
	//Reads NCEP data from binary files.  Temperature (temp), density (dens),
	//dewpoint (dewp),  east-west wind (uwnd), north-south wind (vwnd),
	//geopotential height (geop), sea-level pressure (slpn), surface-level
	//pressure (sfcp), standard deviation of temperature (stmp), standard
	//deviation of density (sden), standard deviation of dewpoint (sdewp), 
	//standard deviation of east-west wind (suwd), standard deviation of 
	//north-south wind (svwd), standard deviation of pressure (sprs), standard
	//deviation of sea-level pressure (sslp), standard deviation of surface 
	//pressure, relative humidity (rhum), standard deviation of relative humidity
	//(srhum), vapor pressure (vprs), standard deviation of vapor pressure (svprs)
	//wind speed (spdav), standard deviation of wind speed (spdsd), u-v correlation
	//(uvcor) read into dynamic arrays using character pointer.  


	float monk, monk1, monk2, monk3, monk4, monk5;

	int i, j, k, m, h;

	ifstream file;

	cout << "Reading NCEP data" << '\n';

	file.open(NCEPpath1.c_str(), ios::in | ios::binary | ios::ate);

	if (file.is_open())
	{
		file.seekg(0, ios::beg);

		//Read first set of data for NCEPhr = 1, 00 UT
		//Read file buffer
		file.read((char *)(&monk), 4);


		//Read temperature  
		for (i = 0; i < 19; i++){
			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&temp[i][j][k]), 4);
				}
			}
		}
		
		for (m = 0; m < 2; m++){
			file.read((char *)(&monk1), 4);

		}

		//Read density 
		for (i = 0; i < 19; i++){
			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&dens[i][j][k]), 4);
				}
			}
		}

		for (m = 0; m < 2; m++){
			file.read((char *)(&monk2), 4);
		}

		//Read dewpoint 
		for (i = 0; i < 19; i++){
			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&dewp[i][j][k]), 4);
				}
			}
		}


		for (m = 0; m < 2; m++){
			file.read((char *)(&monk3), 4);
		}

		//Read east-west wind 
		for (i = 0; i < 19; i++){
			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&uwnd[i][j][k]), 4);
				}
			}
		}

		for (m = 0; m < 2; m++){
			file.read((char *)(&monk4), 4);
		}

		//Read north-south wind 
		for (i = 0; i < 19; i++){
			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&vwnd[i][j][k]), 4);
				}
			}
		}

		for (m = 0; m < 2; m++){
			file.read((char *)(&monk5), 4);
		}

		//Read geopotential height 
		for (i = 0; i < 19; i++){
			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&geop[i][j][k]), 4);
				}
			}
		}


		for (m = 0; m < 2; m++){
			file.read((char *)(&monk5), 4);
		}

		//Read sea-level pressure
		for (j = 0; j < 73; j++){
			for (k = 0; k < 145; k++){
				file.read((char *)(&slpn[j][k]), 4);
			}
		}

		for (m = 0; m < 2; m++){
			file.read((char *)(&monk5), 4);
		}

		//Read surface pressure
		for (j = 0; j < 73; j++){
			for (k = 0; k < 145; k++){
				file.read((char *)(&sfcp[j][k]), 4);
			}
		}

		for (m = 0; m < 2; m++){
			file.read((char *)(&monk5), 4);
		}

		//Read standard deviation of temperature
		for (i = 0; i < 19; i++){
			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&stmp[i][j][k]), 4);
				}
			}
		}

		for (m = 0; m < 2; m++){
			file.read((char *)(&monk5), 4);
		}

		//Read standard deviation of density
		for (i = 0; i < 19; i++){
			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&sden[i][j][k]), 4);
				}
			}
		}

		for (m = 0; m < 2; m++){
			file.read((char *)(&monk5), 4);
		}

		//Read standard deviation of dewpoint
		for (i = 0; i < 19; i++){
			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&sdewp[i][j][k]), 4);
				}
			}
		}

		for (m = 0; m < 2; m++){
			file.read((char *)(&monk5), 4);
		}

		//Read standard deviation of east-west wind
		for (i = 0; i < 19; i++){
			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&suwd[i][j][k]), 4);
				}
			}
		}

		for (m = 0; m < 2; m++){
			file.read((char *)(&monk5), 4);
		}

		//Read standard deviation of north-south wind
		for (i = 0; i < 19; i++){
			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&svwd[i][j][k]), 4);
				}
			}
		}

		for (m = 0; m < 2; m++){
			file.read((char *)(&monk5), 4);
		}

		//Read standard deviation of pressure
		for (i = 0; i < 19; i++){
			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&sprs[i][j][k]), 4);
				}
			}
		}

		for (m = 0; m < 2; m++){
			file.read((char *)(&monk5), 4);
		}

		//Read standard deviation of sea-level pressure
		for (j = 0; j < 73; j++){
			for (k = 0; k < 145; k++){
				file.read((char *)(&sslp[j][k]), 4);
			}
		}

		for (m = 0; m < 2; m++){
			file.read((char *)(&monk5), 4);
		}

		//Read standard deviation of surface pressure
		for (j = 0; j < 73; j++){
			for (k = 0; k < 145; k++){
				file.read((char *)(&ssfcp[j][k]), 4);
			}
		}

		for (m = 0; m < 2; m++){
			file.read((char *)(&monk5), 4);
		}

		//Read relative humidity
		for (i = 0; i < 19; i++){
			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&rhum[i][j][k]), 4);
				}
			}
		}

		for (m = 0; m < 2; m++){
			file.read((char *)(&monk5), 4);
		}

		//Read standard deviation of relative humidity
		for (i = 0; i < 19; i++){
			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&srhum[i][j][k]), 4);
				}
			}
		}


		for (m = 0; m < 2; m++){
			file.read((char *)(&monk5), 4);
		}

		//Read vapor pressure
		for (i = 0; i < 19; i++){
			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&vprs[i][j][k]), 4);
				}
			}
		}

		for (m = 0; m < 2; m++){
			file.read((char *)(&monk5), 4);
		}

		//Read standard deviation of vapor pressure
		for (i = 0; i < 19; i++){
			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&svprs[i][j][k]), 4);
				}
			}
		}

		for (m = 0; m < 2; m++){
			file.read((char *)(&monk5), 4);
		}

		//Read wind speed
		for (i = 0; i < 19; i++){
			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&spdav[i][j][k]), 4);
				}
			}
		}

		for (m = 0; m < 2; m++){
			file.read((char *)(&monk5), 4);
		}

		//Read standard deviation of wind speed
		for (i = 0; i < 19; i++){
			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&spdsd[i][j][k]), 4);
				}
			}
		}

		for (m = 0; m < 2; m++){
			file.read((char *)(&monk5), 4);
		}


		//Read u-v correlation
		for (i = 0; i < 19; i++){
			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&uvcor[i][j][k]), 4);
				}
			}
		}

		//If NCEPhr = 1, exit reading NCEP data
		if (NCEPhr == 1){
			return;
		}

		//Loop to desired NCEPhr from input
		for (h = 1; h < NCEPhr; h++){

			for (m = 0; m < 2; m++){
				file.read((char *)(&monk5), 4);
			}

			for (i = 0; i < 19; i++){
				for (j = 0; j < 73; j++){
					for (k = 0; k < 145; k++){
						file.read((char *)(&temp[i][j][k]), 4);
					}
				}
			}


			for (m = 0; m < 2; m++){
				file.read((char *)(&monk1), 4);
			}


			for (i = 0; i < 19; i++){
				for (j = 0; j < 73; j++){
					for (k = 0; k < 145; k++){
						file.read((char *)(&dens[i][j][k]), 4);
					}
				}
			}

			for (m = 0; m < 2; m++){
				file.read((char *)(&monk2), 4);
			}

			for (i = 0; i < 19; i++){
				for (j = 0; j < 73; j++){
					for (k = 0; k < 145; k++){
						file.read((char *)(&dewp[i][j][k]), 4);
					}
				}
			}


			for (m = 0; m < 2; m++){
				file.read((char *)(&monk3), 4);
			}

			for (i = 0; i < 19; i++){
				for (j = 0; j < 73; j++){
					for (k = 0; k < 145; k++){
						file.read((char *)(&uwnd[i][j][k]), 4);
					}
				}
			}

			for (m = 0; m < 2; m++){
				file.read((char *)(&monk4), 4);
			}

			for (i = 0; i < 19; i++){
				for (j = 0; j < 73; j++){
					for (k = 0; k < 145; k++){
						file.read((char *)(&vwnd[i][j][k]), 4);
					}
				}
			}

			for (m = 0; m < 2; m++){
				file.read((char *)(&monk5), 4);
			}

			for (i = 0; i < 19; i++){
				for (j = 0; j < 73; j++){
					for (k = 0; k < 145; k++){
						file.read((char *)(&geop[i][j][k]), 4);
					}
				}
			}


			for (m = 0; m < 2; m++){
				file.read((char *)(&monk5), 4);
			}

			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&slpn[j][k]), 4);
				}
			}

			for (m = 0; m < 2; m++){
				file.read((char *)(&monk5), 4);
			}

			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&sfcp[j][k]), 4);
				}
			}

			for (m = 0; m < 2; m++){
				file.read((char *)(&monk5), 4);
			}


			for (i = 0; i < 19; i++){
				for (j = 0; j < 73; j++){
					for (k = 0; k < 145; k++){
						file.read((char *)(&stmp[i][j][k]), 4);
					}
				}
			}

			for (m = 0; m < 2; m++){
				file.read((char *)(&monk5), 4);
			}

			for (i = 0; i < 19; i++){
				for (j = 0; j < 73; j++){
					for (k = 0; k < 145; k++){
						file.read((char *)(&sden[i][j][k]), 4);
					}
				}
			}

			for (m = 0; m < 2; m++){
				file.read((char *)(&monk5), 4);
			}

			for (i = 0; i < 19; i++){
				for (j = 0; j < 73; j++){
					for (k = 0; k < 145; k++){
						file.read((char *)(&sdewp[i][j][k]), 4);
					}
				}
			}

			for (m = 0; m < 2; m++){
				file.read((char *)(&monk5), 4);
			}

			for (i = 0; i < 19; i++){
				for (j = 0; j < 73; j++){
					for (k = 0; k < 145; k++){
						file.read((char *)(&suwd[i][j][k]), 4);
					}
				}
			}

			for (m = 0; m < 2; m++){
				file.read((char *)(&monk5), 4);
			}

			for (i = 0; i < 19; i++){
				for (j = 0; j < 73; j++){
					for (k = 0; k < 145; k++){
						file.read((char *)(&svwd[i][j][k]), 4);
					}
				}
			}

			for (m = 0; m < 2; m++){
				file.read((char *)(&monk5), 4);
			}

			for (i = 0; i < 19; i++){
				for (j = 0; j < 73; j++){
					for (k = 0; k < 145; k++){
						file.read((char *)(&sprs[i][j][k]), 4);
					}
				}
			}

			for (m = 0; m < 2; m++){
				file.read((char *)(&monk5), 4);
			}

			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&sslp[j][k]), 4);
				}
			}

			for (m = 0; m < 2; m++){
				file.read((char *)(&monk5), 4);
			}


			for (j = 0; j < 73; j++){
				for (k = 0; k < 145; k++){
					file.read((char *)(&ssfcp[j][k]), 4);
				}
			}

			for (m = 0; m < 2; m++){
				file.read((char *)(&monk5), 4);
			}

			for (i = 0; i < 19; i++){
				for (j = 0; j < 73; j++){
					for (k = 0; k < 145; k++){
						file.read((char *)(&rhum[i][j][k]), 4);
					}
				}
			}

			for (m = 0; m < 2; m++){
				file.read((char *)(&monk5), 4);
			}

			for (i = 0; i < 19; i++){
				for (j = 0; j < 73; j++){
					for (k = 0; k < 145; k++){
						file.read((char *)(&srhum[i][j][k]), 4);
					}
				}
			}

			for (m = 0; m < 2; m++){
				file.read((char *)(&monk5), 4);
			}

			for (i = 0; i < 19; i++){
				for (j = 0; j < 73; j++){
					for (k = 0; k < 145; k++){
						file.read((char *)(&vprs[i][j][k]), 4);
					}
				}
			}

			for (m = 0; m < 2; m++){
				file.read((char *)(&monk5), 4);
			}

			for (i = 0; i < 19; i++){
				for (j = 0; j < 73; j++){
					for (k = 0; k < 145; k++){
						file.read((char *)(&svprs[i][j][k]), 4);
					}
				}
			}

			for (m = 0; m < 2; m++){
				file.read((char *)(&monk5), 4);
			}

			for (i = 0; i < 19; i++){
				for (j = 0; j < 73; j++){
					for (k = 0; k < 145; k++){
						file.read((char *)(&spdav[i][j][k]), 4);
					}
				}
			}

			for (m = 0; m < 2; m++){
				file.read((char *)(&monk5), 4);
			}

			for (i = 0; i < 19; i++){
				for (j = 0; j < 73; j++){
					for (k = 0; k < 145; k++){
						file.read((char *)(&spdsd[i][j][k]), 4);
					}
				}
			}

			for (m = 0; m < 2; m++){
				file.read((char *)(&monk5), 4);
			}

			for (i = 0; i < 19; i++){
				for (j = 0; j < 73; j++){
					for (k = 0; k < 145; k++){
						file.read((char *)(&uvcor[i][j][k]), 4);
					}
				}
			}


		}



		//Close NCEP data file
		file.close();



	}
	else {
		cout << "File Open Error!  " << NCEPpath1 << '\n';
		system("pause");
		exit(1);
	}

}


void NCEPmods::ncepterp(double alpha, double beta, double value[], double *val)
{

	//Member function for lat-lon interpolation of 2.5 degree NCEP data 
	//from grid "square" of values, with interpolation coefficients alpha
	//& beta

	double alphap, betap;

	//Interpolation coefficients interpolation across lat-lon "sqaure"
	alphap = 1.0 - alpha;
	betap = 1.0 - beta;

	//val[0] = value at lat, lon; val[1] = value at lon+2.5, lat;
	//val[2] = value at lon, lat+2.5; val[3] = value at lon+2.5, lat+2.5

	//2-D interpolation across the lat-lon "square"
	*val = alphap*betap*value[0] + alpha*betap*value[1] +
		alphap*beta*value[2] + alpha*beta*val[3];

	return;
}


double NCEPmods::gascon(double td, double p)
{
	//gascon member function from NCEPmods class
	//gas constant R for dewpont temperature Td (deg C) at pressure p (mb)

	const double  r0 = 287.055, omeps = 0.37803, ckf = 273.15,
		onec = 100.0;
	double e, gascon_out;

	e = wexler(td + ckf, onec*p) / onec;

	gascon_out = r0 / (1.0 - omeps*e / p);

	return gascon_out;



}


double NCEPmods::wexler(double t, double p)
{

	//wexler member function from NCEPmods class
	//Wexler formulation for saturation vapor pressure (Pa) as a function
	//of temperature T (kelvin) and pressure p (Pa), as given by Flatau
	//et al., J. Appl. Meteorol., 31(12), 1507, Dec., 1992, with 
	//correction factor (f3) as given by Buck, J. Appl. Meteorol., 
	//20(12), 1527, Dec. 1981.
	//wexler(T dry bulb) gives saturation vapor pressure.
	//wexler(T dew point) gives actual vapor pressure.
	//Relative humidity (0-1) is wexler(Td)/wexler(T)

	const double g0 = -0.29912729e4, g1 = -0.60170128e4,
		g2 = 18.87643854e0, g3 = -0.028354721e0, g4 = 0.17838301e-4,
		g5 = -0.84150417e-9, g6 = 0.44412543e-12, g7 = 2.858487e0,
		a = 1.0007e0, b = 3.46e-08;
	double wexler_out;

	if (t <= 75.0){
		wexler_out = 1.0e-23;
	}
	else {
		wexler_out = exp((g0 + (g1 + (g2 + g7*log(t) + (g3 + (g4 + (g5 +
			g6*t)*t)*t)*t)*t)*t) / pow(t, 2))*(a + b*p);
	}

	return wexler_out;

}


double NCEPmods::ztoh(double g0, double r, double z)
{
	//member function ztoh from NCEPmods class
	//Converts input geometric height z into output geopotential
	//height H.  Parameters are:  g0 - local surface gravity (m/s**2),
	//r - local effective Earth radius, gref - reference value of surface
	//gravity (9.80665 m/s**2)
	//Note: r, z, and H must all be in the same units (e.g. all in meters
	//or all in km

	const double gref = 9.80665;
	double ztoh_out;

	ztoh_out = r*z*g0 / (gref*(r + z));

	return ztoh_out;

}


void NCEPmods::gethgs(double phi, double thet, double *hg1, double *hg2)
{

	//member function gethgs from NCEPmods class
	//Compute geometric height of lowest 10 mb level and highest 20 mb
	//level from four 2.5 by 2.5 degree NCEP profiles surrounding given 
	//lat-lon

	const double grid = 0.4, gref = 9.80665;

	int i1, j1, id, jd, i, j;
	double elon, phir, g0, gg0, re, r0, ri0, a, b, h10, h20, z10, z20;

	//Find i, j indexes for next lower lat-lon grid point
	j1 = int((phi + 90.0)*grid);
	if (j1 == 73) j1 = 72;
	elon = thet;
	if (elon < 0.0) elon = elon + 360.0;
	i1 = int(elon*grid);
	if (i1 > 144) i1 = i1 - 144;

	//Convert latitude to radians and get local gravity and effective
	//radius
	phir = phi*pi180;
	rig(0.0, phir, &g0, &gg0, &re, &r0, &ri0, &a, &b);
	*hg1 = 0.0;
	*hg2 = 99.9;

	//Get geometric heights at corners of NCEP lat-lon "sqaure"
	for (jd = 0; jd < 2; jd++){
		j = j1 + jd;
		for (id = 0; id < 2; id++){
			i = i1 + id;
			//Get geopotential heights of 10mb and 20mb levels
			h10 = double(geop[18][j][i]);
			h20 = double(geop[17][j][i]);
			//Get geometric heights at 10mb and 20mb levels
			z10 = (h10*gref*re / (re*g0 - h10*gref));
			z20 = (h20*gref*re / (re*g0 - h20*gref));
			//Find hg2 = min z10 and hg1 = max z20
			if (z10 < *hg2) *hg2 = z10;
			if (z20 > *hg1) *hg1 = z20;
		}
	}

	return;
}

void NCEPmods::rig(double ch, double phir, double *g0, double *g,
	double *re, double *r0, double *ri, double *a, double *b)
{

	//rig member function from NCEPmods class
	//Computes surface gravity g0 (m/s^2) and effective Earth radius re (km) from 
	//input geocentric latitude phir (radians).  Also computes gravity g (m/s^2) and 
	//total radius ri (km) at input height ch (km).

	const double g0a = 9.80616, g0b = 0.0026373, g0c = 0.0000059,
		rea = 3.085462e-3, reb = 2.27e-6, rec = 2.0e-9;
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

void NCEPmods::sortlevel(int lsort[][19], double zsort[][19])
{

	//sortlevel member function from NCEPmods class
	//Sorts level numbers, to put surface and sea level in right order.

	int i, lvl, j, n;
	double qsort[4][19];
	double zsort0, zsort1;

	//Step through "corners" of lat-lon "square"
	for (n = 0; n < 4; n++){
		//Initialize level numbers
		for (i = 0; i < 19; i++){
			lsort[n][i] = i;
		}


		//Store surface (level 1) and sea level (level 0) altitudes
		zsort1 = zsort[n][1];
		zsort0 = zsort[n][0];

		//Sort to put surface level number into right order
		for (i = 2; i < 19; i++){
			if (zsort1 > zsort[n][i]){
				lsort[n][i - 1] = lsort[n][i];
				lsort[n][i] = 1;
			}
			else{
				break;
			}
		}

		//Store altitudes with surface altitude in right order
		for (j = 0; j < 19; j++){
			lvl = lsort[n][j];
			qsort[n][j] = zsort[n][lvl];
		}

		//Sort to put sea level number into right order
		for (i = 1; i < 19; i++){
			if (zsort0 > qsort[n][i]){
				lsort[n][i - 1] = lsort[n][i];
				lsort[n][i] = 0;
			}
			else{
				break;
			}
		}
	}
	return;

}

void NCEPmods::ncepzterp(double g0, double re, double z, double *pzo,
	double *tzo, double *tdzo, double *rhozo, double *uzo,
	double *vzo, double *wzo, double *spzo, double *stzo,
	double *stdzo, double *srhozo, double *suzo, double *svzo,
	double *dtdzo, double *rhno, double *srhno, double *vpno,
	double *svpno, double *spdavzo, double *spdsdzo,
	double *uvcorrzo, double alpha, double beta)
{

	//ncepzterp member function from NCEPmods class
	//Interpolates from array of values at pressure levels 0-18 to any
	//height z, 0-27 km.

	int k1, k2, i, n;
	int lsort[4][19];
	double pz[4], tz[4], tdz[4], rhoz[4], uz[4], vz[4], wz[4], spz[4],
		stz[4], stdz[4], srhoz[4], suz[4], svz[4], dtdz[4], rhn[4],
		srhn[4], vpn[4], svpn[4], spdavz[4], spdsdz[4], uvcorrz[4];
	double hden, hz, dh, rz, r1, r2, hvp, tzo1 = 0.0;
	double zsort[4][19];
	double pzo1 = 0.0, tdzo1 = 0.0, rhozo1 = 0.0, uzo1 = 0.0, vzo1 = 0.0, spzo1 = 0.0,
		stzo1 = 0.0, srhozo1 = 0.0, suzo1 = 0.0, svzo1 = 0.0, dtdzo1 = 0.0,
		rhno1 = 0.0, srhno1 = 0.0, vpno1 = 0.0, svpno1 = 0.0, spdavzo1 = 0.0,
		spdsdzo1 = 0.0, uvcorrzo1 = 0.0, stdzo1 = 0.0;
	double alphap, betap;


	//Covert input geometric height to geopotential height, for interpolation
	hz = ztoh(g0, re, 1000.0*z) / 1000.0;

	//Step through 4 "corners" of lat-lon "square"
	for (n = 0; n < 4; n++){
		//Find index values k1 and k2 for interpolation
		zsort[n][0] = 0.0;
		zsort[n][1] = h[n][1];
		for (i = 2; i < 19; i++){
			zsort[n][i] = h[n][i];
		}
		//Sort levels to put surface and sea level in right order
		sortlevel(lsort, zsort);
		if (hz < 0.0){
			k1 = lsort[n][0];
			k2 = lsort[n][1];
		}
		else {
			for (i = 0; i < 18; i++){
				k1 = lsort[n][i];
				k2 = lsort[n][i + 1];
				if ((hz < h[n][k1]) || (hz >= h[n][k2])) continue;
				break;
			}
		}

		//Density scale height, km
		if (abs(h[n][k1] - h[n][k2]) == 0.0){
			dtdz[n] = 0.0;
			dh = 0.0;
		}
		else {
			dtdz[n] = (t[n][k2] - t[n][k1]) / (h[n][k2] - h[n][k1]);
			dh = (hz - h[n][k1]) / (h[n][k2] - h[n][k1]);
		}
		if (abs(rho[n][k1] - rho[n][k2]) == 0.0 || abs(h[n][k1] - h[n][k2]) == 0.0){
			hden = 1.0;
			rhoz[n] = rho[n][k1];
		}
		else {
			hden = (h[n][k2] - h[n][k1]) / log(rho[n][k1] / rho[n][k2]);
			//Logarithmic interpolation for density
			rhoz[n] = rho[n][k1] * exp((h[n][k1] - hz) / hden);
		}

		//Linear interpolation on temperature, dewpoint, and gas constant
		tz[n] = t[n][k1] + (t[n][k2] - t[n][k1])*dh;
		tdz[n] = td[n][k1] + (td[n][k2] - td[n][k1])*dh;
		r1 = p[n][k1] / (rho[n][k1] * t[n][k1]);
		r2 = p[n][k2] / (rho[n][k2] * t[n][k2]);
		rz = r1 + (r2 - r1)*dh;
		//Pressure (mb) from gas law
		pz[n] = rhoz[n] * rz*tz[n];
		//Interpolate winds and standard deviations linearly
		uz[n] = u[n][k1] + (u[n][k2] - u[n][k1])*dh;
		vz[n] = v[n][k1] + (v[n][k2] - v[n][k1])*dh;
		wz[n] = w[n][k1] + (w[n][k2] - w[n][k1])*dh;
		spz[n] = sp[n][k1] + (sp[n][k2] - sp[n][k1])*dh;
		stz[n] = st[n][k1] + (st[n][k2] - st[n][k1])*dh;
		stdz[n] = std[n][k1] + (std[n][k2] - std[n][k1])*dh;
		srhoz[n] = srho[n][k1] + (srho[n][k2] - srho[n][k1])*dh;
		suz[n] = su[n][k1] + (su[n][k2] - su[n][k1])*dh;
		svz[n] = sv[n][k1] + (sv[n][k2] - sv[n][k1])*dh;
		rhn[n] = rh[n][k1] + (rh[n][k2] - rh[n][k1])*dh;
		srhn[n] = srh[n][k1] + (srh[n][k2] - srh[n][k1])*dh;
		if (vp[n][k1] - vp[n][k2] == 0.0){
			hvp = 1.0e10;
		}
		else {
			hvp = (h[n][k2] - h[n][k1]) / log(vp[n][k1] / vp[n][k2]);
		}
		vpn[n] = vp[n][k1] * exp((h[n][k1] - hz) / hvp);
		svpn[n] = svp[n][k1] + (svp[n][k2] - svp[n][k1])*dh;
		spdavz[n] = spdavl[n][k1] + (spdavl[n][k2] - spdavl[n][k1])*dh;
		spdsdz[n] = spdsdl[n][k1] + (spdsdl[n][k2] - spdsdl[n][k1])*dh;
		uvcorrz[n] = uvcorrl[n][k1] + (uvcorrl[n][k2] - uvcorrl[n][k1])*dh;

	}

	//Lat-lon interpolation to get output values at input height
	alphap = 1.0 - alpha;
	betap = 1.0 - beta;

	*pzo = alphap*betap*pz[0] + alpha*betap*pz[1] +
		alphap*beta*pz[2] + alpha*beta*pz[3];

	*tzo = alphap*betap*tz[0] + alpha*betap*tz[1] +
		alphap*beta*tz[2] + alpha*beta*tz[3];

	*tdzo = alphap*betap*tdz[0] + alpha*betap*tdz[1] +
		alphap*beta*tdz[2] + alpha*beta*tdz[3];

	*rhozo = alphap*betap*rhoz[0] + alpha*betap*rhoz[1] +
		alphap*beta*rhoz[2] + alpha*beta*rhoz[3];

	*uzo = alphap*betap*uz[0] + alpha*betap*uz[1] +
		alphap*beta*uz[2] + alpha*beta*uz[3];

	*vzo = alphap*betap*vz[0] + alpha*betap*vz[1] +
		alphap*beta*vz[2] + alpha*beta*vz[3];

	*wzo = alphap*betap*wz[0] + alpha*betap*wz[1] +
		alphap*beta*wz[2] + alpha*beta*wz[3];

	*spzo = alphap*betap*spz[0] + alpha*betap*spz[1] + alphap*beta*spz[2]
		+ alpha*beta*spz[3];

	*stzo = alphap*betap*stz[0] + alpha*betap*stz[1] + alphap*beta*stz[2]
		+ alpha*beta*stz[3];

	*stdzo = alphap*betap*stdz[0] + alpha*betap*stdz[1] +
		alphap*beta*stdz[2] + alpha*beta*stdz[3];

	*srhozo = alphap*betap*srhoz[0] + alpha*betap*srhoz[1] +
		alphap*beta*srhoz[2] + alpha*beta*srhoz[3];

	*suzo = alphap*betap*suz[0] + alpha*betap*suz[1] +
		alphap*beta*suz[2] + alpha*beta*suz[3];

	*svzo = alphap*betap*svz[0] + alpha*betap*svz[1] +
		alphap*beta*svz[2] + alpha*beta*svz[3];

	*dtdzo = alphap*betap*dtdz[0] + alpha*betap*dtdz[1] +
		alphap*beta*dtdz[2] + alpha*beta*dtdz[3];

	*rhno = alphap*betap*rhn[0] + alpha*betap*rhn[1] +
		alphap*beta*rhn[2] + alpha*beta*rhn[3];

	*srhno = alphap*betap*srhn[0] + alpha*betap*srhn[1] +
		alphap*beta*srhn[2] + alpha*beta*srhn[3];

	*vpno = alphap*betap*vpn[0] + alpha*betap*vpn[1] +
		alphap*beta*vpn[2] + alpha*beta*vpn[3];

	*svpno = alphap*betap*svpn[0] + alpha*betap*svpn[1] +
		alphap*beta*svpn[2] + alpha*beta*svpn[3];

	*spdavzo = alphap*betap*spdavz[0] + alpha*betap*spdavz[1] +
		alphap*beta*spdavz[2] + alpha*beta*spdavz[3];

	*spdsdzo = alphap*betap*spdsdz[0] + alpha*betap*spdsdz[1] +
		alphap*beta*spdsdz[2] + alpha*beta*spdsdz[3];

	*uvcorrzo = alphap*betap*uvcorrz[0] + alpha*betap*uvcorrz[1] +
		alphap*beta*uvcorrz[2] + alpha*beta*uvcorrz[3];

	return;

}

void NCEPmods::ncepmd(double z, double phi, double thet, double *pz, double *rhoz,
	double *tz, double *uz, double *vz, double *wz, double *tdz,
	double *spz, double *srhoz, double *stz, double *suz, double *svz,
	double *stdz, double *rhn, double *srhn, double *vpn, double *svpn,
	double *spdavz, double *spdsdz, double *uvcorrz)
{

	//ncepmd member function from NCEPmods class
	//Compute mean pressure, density, temperature, wind components, 
	//dewpoint temperature and their standard deviations, for the NCEP data

	const double grid = 0.4, onek = 1000.0, gref = 9.80665;

	int i1, j1, k1, k, id, jd, n, i, j;
	double phir, a, b, elon, alpha,
		beta, g0, gg0, re, r0, ri0, dy, dx, cp, dtd;
	double ri = 6395.7026853182988, ri1, re1, g;
	double dtdx[19], dtdy[19], dtdz1[19], pn[19];
	double pncep[19] = { 9999.0, 9999.0, 1000.0, 925.0, 850.0, 700.0,
		600.0, 500.0, 400.0, 300.0, 250.0, 200.0, 150.0, 100.0, 70.0, 50.0,
		30.0, 20.0, 10.0 };
	double pz1, tz1, tdz1, rhoz1, uz1, vz1, wz1, spz1, stz1, stdz1,
		srhoz1, suz1, svz1, vpn1, svpn1, spdavz1, spdsdz1,
		uvcorrz1, rhn1, srhn1;
	for (i = 0; i < 19; i++){
		pn[i] = 100.0*pncep[i];
		dtdz1[i] = 0.0;
	}

	//Terminate if phi or thet too large
	if ((abs(phi) <= 90.0) & (abs(thet) <= 360.0)){
		//Find i, j indexes for next lower lat-lon grid point
		j1 = int((phi + 90.0)*grid);
		if (j1 == 73) j1 = 72;
		elon = thet;
		if (elon < 0.0) elon = elon + 360.0;
		i1 = int(elon*grid);
		if (i1 > 144) i1 = i1 - 144;
		alpha = elon*grid - int(elon*grid);
		beta = (phi + 90.0 - (j1 - 1.0) / grid)*grid - 1.0;

		//Get distances dx, dy across lat-lon "square"
		phir = phi*pi180;
		rig(0.0, phir, &g0, &gg0, &re, &r0, &ri0, &a, &b);

		//Convert re to meters
		re = onek*re;

		//Compute vertical winds frirom Montgomery stream function
		rig(z, phir, &g0, &g, &re1, &r0, &ri1, &a, &b);
		dy = onek*pi180*ri1 / grid;
		dx = 0.01*dy;

		if (abs(phi) < 89.0) dx = dy*cos(pi180*phi);

		//Get horizontal temperature gradients
		for (k = 2; k < 19; k++){
			dtdx[k] = double(temp[k][j1][i1 + 1] - temp[k][j1][i1]) / dx;
			dtdy[k] = double(temp[k][j1 + 1][i1] - temp[k][j1][i1]) / dy;
		}

		//Get profiles at "corners" of lat-lon "square" for interpolation
		n = -1;
		for (jd = 0; jd < 2; jd++){
			j = j1 + jd;
			for (id = 0; id < 2; id++){
				n = n + 1;
				i = i1 + id;
				//Get pressure, sigma-pressure, and vertical wind at sea level
				//and surface
				p[n][0] = (slpn[j][i]);
				p[n][1] = (sfcp[j][i]);
				sp[n][0] = (sslp[j][i]);
				sp[n][1] = (ssfcp[j][i]);
				w[n][0] = 0.0;
				w[n][1] = 0.0;
				//Gert "corner" profiles from full NCEP arrays
				for (k = 0; k < 19; k++){
					if (k > 1){
						p[n][k] = float(pn[k]);
						sp[n][k] = (sprs[k][j][i]);
					}
					t[n][k] = (temp[k][j][i]);
					td[n][k] = (dewp[k][j][i]);
					rho[n][k] = (dens[k][j][i]);
					u[n][k] = (uwnd[k][j][i]);
					v[n][k] = (vwnd[k][j][i]);
					h[n][k] = (geop[k][j][i]);
					st[n][k] = (stmp[k][j][i]);
					std[n][k] = (sdewp[k][j][i]);
					srho[n][k] = (sden[k][j][i]);
					su[n][k] = (suwd[k][j][i]);
					sv[n][k] = (svwd[k][j][i]);
					rh[n][k] = (rhum[k][j][i]);
					srh[n][k] = (srhum[k][j][i]);
					vp[n][k] = (vprs[k][j][i]);
					svp[n][k] = (svprs[k][j][i]);
					spdavl[n][k] = (spdav[k][j][i]);
					spdsdl[n][k] = (spdsd[k][j][i]);
					uvcorrl[n][k] = (uvcor[k][j][i]);
				}


				//Get vertical wind at "corners" of lat-lon "square"
				for (k = 2; k < 19; k++){
					k1 = k;
					if (k == 18) k1 = k - 1;
					dtdz1[k] = (t[n][k1 + 1] - t[n][k1]) / (h[n][k1 + 1] - h[n][k1]);
					cp = 7.0*gascon(td[n][k], pn[k]) / 2.0;
					dtdz1[k] = dtdz1[k] + g0 / cp;
					dtdz1[k] = max(1.0e-4, dtdz1[k]);
					w[n][k] = float(-(u[n][k] * dtdx[k] + v[n][k] * dtdy[k]) / dtdz1[k]);
				}
			}
		}


		//Do vertical and lat-lon interpolation for values at current 
		//height 
		ncepzterp(g0, re, z, &pz1, &tz1, &tdz1, &rhoz1, &uz1, &vz1, &wz1,
			&spz1, &stz1, &stdz1, &srhoz1, &suz1, &svz1, &dtd, &rhn1,
			&srhn1, &vpn1, &svpn1, &spdavz1, &spdsdz1, &uvcorrz1, alpha, beta);


		*pz = pz1;
		*tz = tz1;
		*tdz = tdz1;
		*rhoz = rhoz1;
		*uz = uz1;
		*vz = vz1;
		*wz = wz1;
		*spz = spz1;
		*stz = stz1;
		*stdz = stdz1;
		*srhoz = srhoz1;
		*suz = suz1;
		*svz = svz1;
		*rhn = rhn1;
		*srhn = srhn1;
		*vpn = vpn1;
		*svpn = svpn1;
		*spdavz = spdavz1;
		*spdsdz = spdsdz1;
		*uvcorrz = uvcorrz1;
		dtz = dtd;

		//Convert standard deviations to variances relative to mean
		*spz = pow((*spz / *pz), 2.0);
		*srhoz = pow((*srhoz / *rhoz), 2.0);
		*stz = pow((*stz / *tz), 2.0);


	}
}


void NCEPmods::namelist(string namef)
{

	//namelist member function from NCEPmods class
	//Open and read namelist input file for setting model parameters.

	ifstream namelist;
	string dummy, NCEPpath, NCEPmn, atmpath, trapath, prtpath, nprpath,
		conpath, rndpath, rrapath, rralist, profile, mn1, NCEPyr1;
	double h1, phi1, thet1, f10, f10b, ap, s10, s10b, xm10, y10,
		y10b, dstdtc, seco, dphi, dthet, dhgt, delt, xm10b;
	int mn, ida, iyr, ihro, mino, nmax, iopt, NCEPyr, iaux;





	//Open and read namelist file
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

	namelist.close();


	if (NCEPhr == 0){
		NCEPhr = 1 + int(round(ihro / 6.0));
		if (NCEPhr > 4) NCEPhr = 1;
	}

	ostringstream convert;
	ostringstream convert1;
	convert << mn;
	convert1 << NCEPyr;
	mn1 = convert.str();
	NCEPyr1 = convert1.str();


	if (mn >= 1 && mn <= 9){
		NCEPmn = "Nb" + NCEPyr1 + "0" + mn1 + ".bin";
	}
	else{
		NCEPmn = "Nb" + NCEPyr1 + mn1 + ".bin";
	}
	NCEPpath1 = NCEPpath + NCEPmn;


}



