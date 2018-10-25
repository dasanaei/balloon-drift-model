            NATIONAL CENTERS FOR ENVIRONMENTAL PREDICTION (NCEP) 
           REANALYSIS PROJECT CLIMATOLOGY DATA AND HOW TO USE IT

         Adapted from NCEP/NCAR Reanalysis Data Archive Description
                      Version Revised: 26 September 2003 

________________________________________________________________________________

Overview -

The NCEP/NCAR Reanalysis Project is a joint project between the National 
Centers for Environmental Prediction (NCEP, formerly "NMC") and the National 
Center for Atmospheric Research (NCAR). The goal of this joint effort is to 
produce new atmospheric analyses using historical data (1948 onwards) and as 
well to produce analyses of the current atmospheric state (Climate Data 
Assimilation System, CDAS).  Until recently, the meteorological community has 
had to use analyses that supported the real-time weather forecasting. These 
analyses are very inhomogeneous in time as there have been big improvements in 
the data assimilation systems. The quality and utility of the re-analyses 
should be superior to NCEP's original analyses because:

* a state-of-the-art data assimilation is used 
* more observations are used 
* quality control has been improved 
* the model/data assimilation procedure remains unchanged during the project 
* many more fields are being saved 
* global (some older analyses were hemispheric) 
* better vertical resolution (stratosphere) 

More information about the reanalysis project and data are available from
several sources, including (Kalnay et al., 1996) and at:

http://www.cdc.noaa.gov/cdc/data.nmc.reanalysis.html

Data used in Earth-GRAM 2010 were downloaded from the NOAA Air Resources Lab
(ARL) archives at

http://www.arl.noaa.gov/archives.php
________________________________________________________________________________

Availability -

A subset of the full NCEP data is available from ARL in a format suitable for  
many applications, by selecting "Reanalysis" in the meteorological data set 
selection pull-down menu.

	Data Period Availability: 1948 ==> 2008 

The directory contains data files with the following syntax:

R{S|P}{YEAR}{MONTH}.{gbl|tbd}

Where R indicates "Reanalysis", S or P indicates that the data are on Sigma or 
Pressure surfaces, YEAR is a four digit year, and MONTH is a two digit month.
The "RP" (pressure level data) is what is used in Earth-GRAM 2010.  The file 
suffix identifies the projection as either the 2.5 degree global latitude
-longitude projection (gbl, as used in Earth-GRAM 2010), or a regional conformal 
map projection.

The sigma level data were obtained from NCEP's internal spectral coefficient
archive, and are not used by Earth-GRAM 2010. The pressure level data were 
obtained from the NOAA-CIRES Climate Diagnostics Center, Boulder, Colorado, USA.

________________________________________________________________________________

Additional Data Set Details

  Pressure Level Data

	2.5 degree latitude-longitude global grid
	144x73 points from 90N-90S, 0E-357.5E 
	1/1/1948 - present with output every 6 hours 
	Levels (hPa): 1000,925,850,700,600,500,400,300,250,200,
                      150,100,70,50,30,20,10 
	Surface or near the surface (.995 sigma level) winds and temperature 
	Precipitation 

	Model Type:            LAT-LON
	Vert Coord:            2
	Numb X pt:           144
	Numb Y pt:            73
	Numb Levels:          18
	Sfc Variables:         5 PRSS T02M U10M V10M TPP6
	Upper Levels:          6 HGTS TEMP UWND VWND WWND RELH


Example Header

92 1 1 0 0 099INDX   0  .0000000E+00  .0000000E+00CDC1  0 0  90.00 357.50   2.50   2.50
.00    .00    .00   1.00   1.00 -90.00    .00    .00144 73 18 2 996.00000 5PRSS199 T02M232
U10M 77 V10M176 TPP6195 1000.0 6HGTS 60 TEMP183 UWND212 VWND 49 WWND 42 RELH149 925.00 6HGTS
16 TEMP110 UWND215 VWND 88 WWND131 RELH  3 850.00 6HGTS247 TEMP185 UWND234 VWND100 WWND182 
RELH 52 700.00 6HGTS161 TEMP 20 UWND123 VWND 72 WWND187 RELH 96 600.00 6HGTS152 TEMP 60 
UWND211 VWND148 WWND137 RELH 38 500.00 6HGTS136 TEMP193 UWND144 VWND 99 WWND130 RELH161
400.00 6HGTS107 TEMP212 UWND159 VWND117 WWND136 RELH139 300.00 6HGTS 80 TEMP127 UWND131
VWND  6 WWND198 RELH255 250.00 5HGTS 25 TEMP 66 UWND101 VWND247 WWND 57 200.00 5HGTS188
TEMP153 UWND169 VWND 92 WWND 23 150.00 5HGTS 86 TEMP  3 UWND133 VWND 98 WWND149 100.00
5HGTS201 TEMP 71 UWND135 VWND136 WWND 12 70.000 4HGTS 46 TEMP 82 UWND153 VWND208 50.000 
4HGTS122 TEMP  4 UWND230 VWND162 30.000 4HGTS 20 TEMP132 UWND117 VWND211 20.000 4HGTS 37 
TEMP119 UWND212 VWND165 10.000 4HGTS108 TEMP157 UWND101 VWND159

Data Packing Format

NCEP typically saves their model output in GRIB format. However, at ARL the 
data are reprocessed and stored in a 1-byte packing algorithm. This 1-byte 
packing is a bit more compact than GRIB and can be directly used on a variety 
of computing platforms with direct access I/O. 

The data array is packed and stored into one byte characters. To preserve 
as much data precision as possible the difference between the values at grid 
points is saved and packed rather than the actual values. The grid is then 
reconstructed by adding the differences between grid values starting with 
the first value, which is stored in unpacked ASCII form in the header record. 
To illustrate the process, assume that a grid of real data, R, of dimensions 
i,j is given by the below example. 


  1,j       2,j        ....    i-1,j      i,j
  1,j-1     2,j-1      ....    i-1,j-1    i,j-1
  ....      ....       ....    ....       ....
  1,2       2,2        ....    i-1,2      i,2
  1,1       2,1        ....    i-1,1      i,1


The packed value, P, is then given by 

     Pi,j = (Ri,j  - Ri-1,j)* (2**(7-N)),


where the scaling exponent 


     N = ln dRmax / ln 2 .


The value of dRmax is the maximum difference between any two adjacent grid 
points for the entire array. It is computed from the differences along each 
i index holding j constant. The difference at index (1,j) is computed from 
index (1,j-1), and at 1,1 the difference is always zero. The packed values 
are one byte unsigned integers, where values from 0 to 126 represent -127 
to -1, 127 represents zero, and values of 128 to 254 represent 1 to 127. 
Each record length is then equal in bytes to the number of array elements 
plus 50 bytes for the header label information. The 50 byte label field 
precedes each packed data field and contains the following ASCII data: 


Field          Format    Description                             
Year           I2        Greenwich date for which data valid
Month          I2                       "
Day            I2                       "
Hour           I2                       "
Forecast*      I2        Hours forecast, zero for analysis
Level          I2        Level from the surface up
Grid           I2        Grid identification
Variable       A4        Variable label
Exponent       I4        Scaling exponent needed for unpacking
Precision      E14.7     Precision of unpacked data
Value 1,1      E14.7     Unpacked data value at grid point 1,1
*Forecast hour is -1 for missing data.  

Example: 2 1 1 0 0 099PRSS   8  .1007874E+01  .6867000E+03

ARL provides a Fortran-90 program that can be used to unpack and read the 
first few elements of the data array for each record of an ARL packed 
meteorological file. This file, called CHK_DATA.F, was adapted for reading 
and processing the ARL data files for use in Earth-GRAM 2010.


Reference:
Kalnay et al.,The NCEP/NCAR 40-year reanalysis project, Bull. Amer. Meteor. 
Soc., 77, 437-470, 1996. 

----------------------------------------------------------------------------

                         NCEP Data in Earth-GRAM 2010
                             
NCEP climatology is provided in ASCII format in folder NCEPdata\FixedASCII.
Files for each month are named Nfy1y2mm.txt, where the period of record (POR)
covers years y1 through y2 (e.g. 9008 is for POR 1990 through 2008), and 
month is mm.  Surface and near-surface NCEP data in these files have been
adjusted according to procedures described in file NCEPmods.pdf, in the
DOCUMENTATION folder.  ASCII data in these files have been prepared in a compact 
format, with decimals removed, and some spaces removed.  An example file 
header and first few data lines is given by:


h L J  I   N   Hav  Hsd Uav Us Vav Vs Tav Ts Dav Dsn RHaRHs TdaTds eTa eTsn Sav Ss Ruv
1 1 1  1 589  6828   56  12 27 -20 302467 3696401321 4572742355107357627202  42 21  12
1 2 1  1 589    16  597   8 23 -30 262656 381310 190 7451692615 51270910241  42 20 -89
1 3 1  1 589  6111  573   7 23 -29 262617 371231 180 7451692576 501989 7611  41 19 -91
1 4 1  1 589 12621  562   7 22 -29 262576 371149 160 7451692538 491439 5571  41 19 -80
1 5 1  1 589 27172  591  14 25 -32 282483 3398181311 7451692448 46645624912  46 22-137
1 6 1  1 589 38380  651  22 42 -29 442455 278512 951 4562732344 96295619212  62 33-230
1 7 1  1 589 51290  745  10 48 -19 512381 287314 861 4392272282 761417 8292  65 35-240
1 8 1  1 589 66535  869   6 58   3 632286 256097 661 3352162164 84418830893  75 42-265
1 9 1  1 589 85419  967   1 65  25 692220 214708 441  941321996 67533373864  84 50-314
110 1  1 589 97290  952  -3 61  33 642232 303903 531  861182003 62535966994  81 50-350
111 1  1 589111984  932  -7 53  38 572264 333078 451  791042025 62721790274  73 45-326
112 1  1 589131199  962  -8 44  37 482291 262281 261  72 902042 59 87310453  65 37-244
113 1  1 589158571 1021  -8 34  32 432320 171502 111  64 762061 54105511593  55 32-182
114 1  1 589182933 1012 -11 31  29 392343 151041  71  57 622074 49119711983  51 30-180
115 1  1 589206107  978 -11 28  28 332359 177383 542  50 472082 44126411293  45 27-128
116 1  1 589241563  948 -11 25  20 292384 214383 392  43 332096 38142310613  38 23 -38
117 1  1 589270016  990 -13 24  13 262409 222892 262  37 202107 301557 8273  35 20  -6
118 1  1 589319695 1189  -9 24  21 242500 211392 122  31 202166 183065 7283  36 19 -12


where 
 h   = UT hour of day: 1=00 UT, 2=06UT, 3=12UT, 4=18UT.  Statistics for all
       hours of day combined (h=5) are added in the conversion from ASCII
       format to binary format, by utility program NCEPbinF.f90.
 L   = Pressure levels 1-18.  In mb, these are surface,1000,925,850,700,600,
       500,400,300,250,200,150,100,70,50,30,20, and 10. Data for sea-level
       are computed in the conversion from ASCII format to binary format,
       by utility program NCEPbinF.f90.
 J   = Latitude index (1-73)
 I   = East longitude index (1-144).
 N   = Number of data in the period of record (e.g. if there are no missing
       data, POR 1990-2008 has 589 values for a 31-day month)
 Hav = Monthly average geopotential altitude for given UT hour (m times 10).
       For level 1, Hav is average surface pressure (mb times 10)
 Hsd = Standard deviation for geopotential altitude (m times 10).  For level
       1, Hsd is standard deviation of surface pressure (mb times 10)
 Uav = Monthly average Eastward wind for given UT hour (m/s times 10)
 Us  = Standard deviation for Eastward wind (m/s times 10)
 Vav = Monthly average Northward wind for given UT hour (m/s times 10)
 Vs  = Standard deviation for Northward wind (m/s times 10)
 Tav = Monthly average temperature for given UT hour (K times 10)
 Ts  = Standard deviation for temperature (K times 10)
 Dav = Monthly average density for given UT hour
 Ds  = Standard deviation for density
 n   = Density exponent.  Densities in kg/m**3 are given by Dav*10**(-(3+n))
       and Ds*10**(-(3+n))
 RHa = Monthly average relative humidity for given UT hour (percent times 10)
 RHs = Standard deviation for relative humidity (percent times 10)
 Tda = Monthly average dewpoint temperature for given UT hour (K times 10)
 Tds = Standard deviation for dewpoint temperature (K times 10)
 eTa = Monthly average vapor pressure for given UT hour
 eTs = Standard deviation for vapor pressure
 n   = Vapor pressure exponent.  Vapor pressures in N/m**2 are given by
       eTa*10**(-n) and eTs*10**(-n)
 Sav = Monthly average wind speed (m/s times 10)
 Ss  = Monthly standard deviation of wind speed (m/s times 10)
 Ruv = Cross-correlation between Eastward and Northward wind components (times 1000)
 
NCEP ASCII data format is (I1,2I2,I3,I4,I6,I5,4(I4,I3),I1,2(I4,I3),2I4,I1,I4,I3,I4).

Binary format NCEP data are provided in PC form. Files for each month are named 
Nby1y2mm.bin, where the period of record is for years y1 through y2 (e.g. 9008 is 
for POR 1990 through 2008), and month is mm.  Utility program NCEPbinF is used to 
generate similar binary format files on any non-PC platform.  See README8.txt for 
additional details.  In addition to building binary format versions of the NCEP data,
utility program NCEPbinF can also be used at any time to output easily-readable 
vertical profile data and horizontal grids ("maps") of NCEP data values.  The 
NCEPbinF utility program MUST be used to generate binary version NCEP data for use 
by Earth-GRAM 2010 on non-PC computer platforms (See README0.txt and README8.txt).