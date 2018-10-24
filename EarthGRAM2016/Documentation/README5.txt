                OPTIONAL RANGE REFERENCE ATMOSPHERE (RRA) DATA


    An option exists to use data (in the form of vertical profiles) from a set 
of Range Reference Atmospheres (RRA), as an alternate to the usual GRAM 
climatology, at a set of RRA site locations.  With this feature it is possible, 
for example, to simulate a flight profile that takes off from the location of 
one RRA site (e.g. Edwards AFB, using the Edwards RRA atmospheric data), to 
smoothly transition into an atmosphere characterized by the GRAM climatology, 
then smoothly transition into an atmosphere characterized by a different RRA 
site (e.g. White Sands, NM), to be used as the landing site in the simulation.
    RRA data includes information on both monthly means and standard deviations
of the various parameters at the RRA site (see description below).  Under
the RRA option, when a given trajectory point is sufficiently close to an 
RRA site (latitude-longitude radius from site less than "sitenear", see 
below), then the mean RRA data replace the mean values of the conventional 
GRAM climatology, and the RRA standard deviations replace the conventional 
GRAM standard deviations in the perturbation model computations.  In Earth-GRAM 
2007 a feature was introduced to replace GRAM surface data with surface data  
from the appropriate RRA site, when the RRA option is used.
    A total of 21 sites are available for 2006 RRA data, 17 sites for the 1983 RRA data, and 4
sites for the 2013 RRA data.  After discovery of an error in standard deviations of pressure and 
density in RRA data from the Roosevelt Roads Puerto Rico site (code rrd), data from this RRA 
site are no longer included.  The 2006 RRA data period-of-record varies from site to site, 
but is generally 1990 to 2002 for the 2006 RRA data.  Exceptions are El Paso (1990-95), Great 
Falls (1990-94), Taguac (1990-99), China Lake (1948-2000), and White Sands (1949-1993).  El Paso 
RRA data (0-30 km) are augmented with White Sands rocket-sonde data from 30-70km.  
White Sands RRA file includes only 0-30 km data.  The user can also prepare (in the 
appropriate format, as described below) data for any other site desired, for use 
in the RRA mode.  Another mode to use alternatives for the conventional GRAM 
climatology is the "auxiliary profile" option (discussed in README7.txt).
Auxiliary profiles also allow substitution of alternate data for both means and
standard deviations.  The format for preparing auxiliary profiles is much simpler than 
the required formatting for additional RRA data files (see file README7.txt).
   In addition to the RRA data files, a file called rrasites.txt is provided. 
This file gives file-code identifier, geodetic and geocentric latitude, longitude, 
surface altitude, maximum altitude, and WMO site number, and site name, for each of 
the available 2006 RRA sites.  Alternate lists of RRA files may be prepared by the user,
if new sites or a different set of sites are desired.  The RRA site list file to be
used by Earth-GRAM 2010 for a specific program run is specified by input parameter rralist. 


               List of RRA site data provided, file "rrasites.txt":

Code Year  GdLat  GcLat Lon(E+) Hgt(m) Zmax WMO #  Site Name
 asc 1983  -7.93  -7.88  -14.42    20.  70. 619020 Ascension Island, Atlantic                 
 bar 1983  22.03  21.90 -159.78     5.  70. 911620 Barking Sands, Hawaii                      
 cap 1983  28.47  28.31  -80.55     3.  70. 747940 Cape Canaveral, Florida                    
 dug 1983  40.77  40.58 -111.97  1288.  70. 725720 Dugway Proving Ground (Salt Lake City), UT 
 eaf 1983  34.92  34.74 -117.90   705.  70. 723810 Edwards Air Force Base, California         
 egl 1983  30.48  30.31  -86.52    20.  30. 722210 Eglin AFB, Florida                         
 kmr 1983   8.73   8.67  167.75     2.  70. 913660 Kwajalein Missile Range, Pacific           
 ptu 1983  34.12  33.94 -119.12     4.  70. 723910 Point Mugu Naval Air Weapons Center, CA    
 tag 1983  13.55  13.46  144.85   111.  30. 912170 Taguac, Guam                               
 vaf 1983  34.75  34.57 -120.57   100.  70. 723930 Vandenberg AFB, California                 
 wal 1983  37.85  37.66  -75.48     3.  70. 724020 Wallops Island, Virginia                   
 wsm 1983  32.38  32.21 -106.48  1246.  70. 722696 White Sands, New Mexico                    
 fad 1983  64.82  64.67 -147.87   135.  30. 702610 Fairbanks, Alaska                          
 nel 1983  36.62  36.44 -116.02  1007.  30. 723870 Nellis AFB, Nevada                         
 shm 1983  52.72  52.53  174.12    39.  70. 704140 Shemya, Alaska                             
 thu 1983  76.52  76.43  -68.50    59.  70. 042020 Thule, Greenland                           
 wak 1983  19.28  19.16  166.65     5.  30. 912450 Wake Island, Pacific            
 anf 2006  47.62  47.43  -52.73   140.  30. 718010 Argentia, Newfoundland (St. Johns Airport) 
 asc 2006  -7.93  -7.88  -14.42    79.  70. 619020 Ascension Island, Atlantic                 
 bar 2006  21.98  21.85 -159.34    31.  30. 911650 Barking Sands, Hawaii (Lihue)              
 cap 2006  28.47  28.31  -80.55     3.  70. 747940 Cape Canaveral, Florida                    
 chl 2006  35.68  35.50 -117.68   665.  30. 746120 China Lake Naval Air Weapons Center, CA    
 dug 2006  40.77  40.58 -111.97  1288.  30. 725720 Dugway Proving Ground (Salt Lake City), UT 
 eaf 2006  34.92  34.74 -117.90   724.  30. 723810 Edwards Air Force Base, California         
 egl 2006  30.48  30.31  -86.52    20.  30. 722210 Eglin AFB, Florida                         
 elp 2006  31.81  31.64 -106.38  1199.  70. 722700 El Paso, Texas                             
 fad 2006  64.80  64.65 -147.88   135.  30. 702610 Fairbanks, Alaska                          
 fha 2006  32.12  31.95 -110.93   787.  30. 722740 Ft. Huachuca Elec Prvng Grnd (Tucson), AZ  
 gtf 2006  47.47  47.28 -111.38  1118.  30. 727750 Great Falls, MT                            
 kmr 2006   8.73   8.67  167.75     2.  70. 913660 Kwajalein Missile Range, Pacific           
 ncf 2006  43.87  43.68    4.40    62.  30. 076450 Nimes-Courbessac, France (STS TAL Site)    
 nel 2006  36.62  36.44 -116.02  1007.  30. 723870 Nellis AFB, Nevada (Mercury)               
 ptu 2006  34.12  33.94 -119.12     2.  70. 723910 Point Mugu Naval Air Weapons Center, CA    
 tag 2006  13.55  13.46  144.85    78.  30. 912170 Taguac, Guam (Anderson AFB)                
 vaf 2006  34.75  34.57 -120.57   121.  30. 723930 Vandenberg AFB, California                 
 wal 2006  37.85  37.66  -75.48    13.  30. 724020 Wallops Island, Virginia (NASA)            
 wsm 2006  32.38  32.21 -106.48  1207.  30. 722690 White Sands Missile Range, New Mexico       
 ysd 2006  32.87  32.69 -117.14   134.  30. 722930 Yuma Proving Ground, AZ (San Diego, CA)
 cap 2013  28.47  28.31  -80.57     4.  30. 747940 Cape Canaveral, Florida                    
 eaf 2013  34.91  34.74 -117.88   723.  30. 723810 Edwards Air Force Base, California 
 vaf 2013  34.73  34.53 -120.58   100.  30. 722930 Vandenberg AFB, California
 wsm 2013  32.94  32.88 -106.42  1251.  30. 722690 White Sands Missle Range, New Mexico    
           
NOTE: Latitudes in the RRA data files are geodetic.  In the rras2006.txt file, latitudes 
GdLat are geodetic and latitudes GcLat are geocentric.  Earth-GRAM input latitudes are 
geocentric.  Internally, Earth-GRAM works with geocentric RRA site latitudes if the 
RRA data option is used. Earth-GRAM output file lists both geodetic and geocentric 
RRA site latitudes.

NOTE: For Kwajalein Missile Range (kmr), data are available above 30 km for the
annual mean case only.  None of the individual months have sufficient data (number
must be greater than 10) to be used, so all months have missing data codes (e.g. 99.99)
between 30 and 70 km.

    Use of the RRA data option in Earth-GRAM 2010 is controlled by several input 
parameters (see file README2.txt): (1) rrapath gives the name of the directory
containing the RRA data; (2) rralist gives the name of the file containing RRA
sites to be used; (3) iurra is the unit number to be used by the program for reading 
the RRA data;  and (4) sitelim and (4) sitenear control the size of the latitude-
longitude region in the vicinity of any RRA site for which data from that RRA site is 
to be used (the RRA zone of influence).  For any location having a radial distance (in 
latitude-longitude terms) of less than the value given by sitenear, the RRA data (with 
a full weight of 1) is used.  For any location outside a latitude-longitude radius 
given by sitelim, the GRAM climatology data is used (i.e. a weight of 0 for the RRA 
data).  Between radial distances of sitenear and sitelim, a weighted average of RRA 
and GRAM climatology data is used, insuring a smooth transition from RRA data to GRAM 
data.
   Nominal (default) values are sitenear = 0.5 degrees and sitelim = 2.5 degrees.
For these values, then any latitude-longitude within a radius of 0.5 degrees from
any of the RRA sites, will use data from that RRA site.  Any location beyond a 
radius of 2.5 degrees would use the GRAM data. Between 0.5 and 2.5 degrees 
radius, a weighted average of RRA data and GRAM data would be used, with the RRA 
data weight smoothly changing from 1 at a radius of 0.5 degrees to a weight of 0 
at a radius of 2.5 degrees.
   Depending on the value of sitelim, and the proximity of the various RRA sites 
used, it may be possible that a given trajectory location is in the vicinity of 
more than one RRA site (e.g. for locations near Point Mugu, Edwards AFB, and 
Vandenberg AFB).  If a given trajectory location could be influenced by more than
one RRA site, only data from the NEAREST (highest weight) site is used.  NOTE 
that if, in such a situation, the user desires to ALWAYS use a specific RRA site 
(e.g. Edwards) and NEVER use a nearby RRA site (e.g. Point Mugu), then the name 
and information for the undesired nearby RRA site should be removed from the 
file list, and the modified site list used, as specified by input parameter 
rralist.
   RRA data apply from 0 to (at most) 70 km (above mean sea level).  There is 
also a smooth interpolation process used to transition from the RRA data to GRAM
data as the top of the RRA data is approached.  This transition takes place between
the next-to-highest RRA altitude (at which RRA weight = 1) and the highest RRA 
altitude (at which the RRA weight becomes 0).
   RRA date files for a given site consist of three data files, T1sssyy.txt,
T2sssyy.txt, and T3sssyy.txt, where sss is the three-character site code from
the list of sites given above, and yy=06 for the 2006 RRA data.  The RRA data 
files are in the format given in the series of Range Reference Atmosphere reports 
(e.g. Document 361-83, Cape Canaveral Range Reference Atmosphere 0-70 km Altitude, 
February, 1983, Meteorology Group, Range Commanders Council). Files Txsssyy.txt 
correspond to Table x (x = 1-3) in the RRA reports.
    Each Txsssyy.txt file contains an annual average data set, as well as 12 
monthly data sets.  Only the RRA data for the desired month are used by GRAM, 
and the annual average data are always ignored on read-in. For 2006 RRA data, the 
annual average data must follow the 12 monthly data sets.  Annual data preceded
the 12 monthly data sets in the 1983 RRA data files (no longer used in Earth-GRAM
2010).  
    Table-1 data contain wind statistical parameters: height, mean Eastward wind, 
standard deviation of Eastward wind, correlation between Eastward and Northward wind, 
mean Northward wind, standard deviation of Northward wind, mean wind speed, standard 
deviation of wind speed, skewness of wind speed(*), and number of observations(*). 
Asterisks denote parameters that are not used or output by GRAM (although any RRA
data having number of data less than 10 is ignored by the RRA read routines).
   Table-2 data contain thermodynamic statistical parameters: height, mean 
pressure, standard deviation of pressure, skewness of pressure(*), mean tempera-
ture, standard deviation of temperature, skewness of temperature(*), mean 
density, standard deviation of density, skewness of density(*), number of 
pressure observations(*), number of temperature observations(*), and number of
density observations(*).   Again, any data having fewer than 10 observations are
ignored by the RRA read routine.
   Table-3 data contain moisture related statistical parameters: height, mean 
vapor pressure, standard deviation of vapor pressure, skewness of vapor 
pressure(*), mean virtual temperature(*), standard deviation of virtual tempera-
ture(*), skewness of virtual temperature(*), mean dewpoint temperature, standard 
deviation of dewpoint temperature, skewness of dewpoint temperature(*), number 
of observations of vapor pressure and dewpoint temperature(*), and number of 
observations of virtual temperature(*). RRA Table-4 data are not used in GRAM.
   User-provided "RRA" data can also be used if the following conditions are 
adhered to.  Each new "RRA" site must be entered into the "rralist" file (maximum 
total number of sites allowed is 99). The site code can be any 3-character code 
not already being used.  Heights for any new RRA data must be in the range from 0 
to altitude Zmax, as given in the "rralist" file. Heights should be given in 
ascending order in the Txsssyy.txt files, with 300 or fewer heights entered (height 
increments can be any value, and fixed height increments do not have to be used).  
The first data line of each Txsssyy.txt may have descriptive information (such as 
site name). However, the first data line of each file MUST contain the site 
latitude and longitude. Latitude is given as xx.xxN or xx.xxS; longitude is given 
as xxx.xxE or xxx.xxW, in format (17X,F5.2,A1,2X,F6.2,A1).  Latitude and longitude 
values from the first data line are compared with the latitude and longitude in the
"rralist" file, to insure that the appropriate site data are being used. In the
"rralist" file, north latitudes are positive (and south latitudes are negative), 
while east longitudes are positive (and west longitudes are negative).  Each file 
T1sssyy.txt or T2sssyy.txt must have FIVE "header" lines preceding each monthly
set of data values.  Each file T3sssyy.txt must have SIX "header" limes preceding 
each monthly set of data values. Monthly data values for all 12 months must be 
provided.  The data lines may be in free-field (list directed) format, but MUST 
contain a numerical value for each of the parameters expected, depending on the 
Table type.  Parameters not used by GRAM (those indicated by asterisks, above) may 
be input as zero values (except for number of observations, which can be any number 
greater than 10). Missing values (i.e. those that will be ignored) may be indicated 
by using 99.99 or 999.99.
