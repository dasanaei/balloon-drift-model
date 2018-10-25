        DISCUSSION OF SEVERAL FEATURES AND OPTIONS IN EARTH-GRAM 2010
               

DIFFERENCES BETWEEN PC-VERSION AND UNIX-VERSION CODE

With the introduction of NCEP climatology data to represent the lower
altitude region (0-27 km), Earth-GRAM 2010 now has no differences between 
UNIX-version and PC-version code.  For discussion of the NCEP data, see
file README4.txt.


REVISED RANGE REFERENCE ATMOSPHERE DATA

In 2006, the Air Force Combat Climatology Center (AFCCC) developed a set of 
revised Range Reference Atmosphere (RRA) data. Earth-GRAM 2010 has the option 
of using the 2006 revised RRA data as a replacement for conventional GRAM (NCEP/
MAP etc.) climatology.  Details of the RRA data are provided in file README5.txt.
After discovery of an error in standard deviations of pressure and density in RRA 
data from the Roosevelt Roads Puerto Rico site (code rrd), data from this RRA site 
are no longer included.  The option to use was 1983 RRA database was added in 
Earth-GRAM2010V3.0. 
                  
                  
OPTIONAL AUXILIARY PROFILE INPUT   

In addition to RRA options, an "auxiliary profile" feature has been implemented.
This allows the user to input a data profile of pressure, density, temperature, 
and/or winds versus altitude, with the auxiliary profile values used in place
of conventional climatology (NCEP/MAP/etc.) values. This option is controlled 
by setting parameters profile, sitenear, and sitelim in the NAMELIST input
file.  Parameter profile gives the file name containing the profile data
values.  Parameter sitenear is the latitude-longitude radius (in degrees) 
within which weight for the auxiliary profile is 1.0.  Parameter sitelim
is the latitude-longitude radius (in degrees) beyond which weight for the
auxiliary profile is 0.0.  A weighting factor for the profile data (profwgt),
having values between 0 and 1, is applied between radii sitelim and sitenear.
Mean conditions are as given in the profile file if the desired point is within
a lat-lon radius of sitenear from the profile lat-lon at the given altitude; 
mean conditions are as given by the original NCEP/MAP/etc. data if the desired 
point is beyond a lat-lon radius of sitelim from the lat-lon of the profile 
at the given altitude.  If sitenear = 0, then profile data are NOT used.
The profile weight factor (profwgt) for the auxiliary profile also varies 
between 0 at the first profile altitude level and 1 at the second profile 
altitude level (and between 1 at the next-to-last profile altitude level and 
0 at the last profile altitude level).  First and second profile points (and 
next-to-last and last profile points) should therefore be selected widely 
enough apart in altitude that a smooth transition can occur as profwgt changes 
from 0 to 1 near these profile end points.  Characteristics and format of the 
auxiliary profile input file are given in file README2.txt. Example auxiliary 
profile files, generated from annual average Range Reference Atmosphere statistics
(files RRAsssAnn.txt for RRA site sss, in folder AuxProfiles) are also provided.  
NOTE: Auxiliary Profile option and RRA data option cannot both be invoked 
simultaneously.  Height values greater than 6000 km in the profile input are 
interpreted as radius values.

                  
THERMOSPHERE MODELS

Earth-GRAM 2007 included several updates to the Marshall Engineering Thermosphere
(MET) model, now denoted MET-07. MET updates included:
  (1) corrections for inconsistency between constituent number density and mass 
      density [described by Justus, C. G., Aleta Duvall, and Vernon W. Keller, 
      "Earth Global Reference Atmospheric Model (GRAM-99) and Trace Constituents", 
      Paper C4.1-0002-04, Presented at 35th COSPAR Scientific Assembly Paris, 
      France July 18-25, 2004].
  (2) Representation of gravity above an oblate spheroid Earth shape, rather than 
      using a spherical Earth approximation.
  (3) Treatment of day-of-year as a continuous variable in the semi-annual term, 
      rather than as an integer day.
  (4) Treatment of year as either 365 or 366 days in length (as appropriate), 
      rather than all years having length = 365.2422 days.
  (5) Allows continuous variation of time input, rather than limiting time 
      increments to integer minutes.
For Earth-GRAM 2010, the only changes to MET were to convert to Fortran 90, and make
the code double precision.

As an alternate to MET, the option is now available for Earth-GRAM to use either 
the Naval Research Labs Mass Spectrometer, Incoherent Scatter Radar Extended 
Model for the thermosphere (NRL MSIS E-00), or the Jacchia-Bowman 2008 thermosphere
model (JB2008).  If the MSIS or JB2008 model option is used, the NRL Harmonic Wind 
Model (HWM 93) is used to compute thermospheric winds. No further updates to these
thermosphere models were made for Earth-GRAM 2010.  For further discussion, see 
files README6.txt, MSISfix.txt, and HWMfix.txt.


COORDINATE SYSTEMS AND EARTH REFERENCE ELLIPSOID 
                 
Equatorial and polar Earth radii for the "sea-level" reference ellipsoid have 
been updated to World Geodetic System (WGS 84) values.  Previous (GRAM-99) 
radius values were from IAU 76.  WGS 84 values are used by the GPS navigation 
system.  These are also equivalent (to 10 significant figures) to the Geodetic 
Reference System (GRS 80) values.  Other recent values that could be used include 
the International Earth Rotation & Reference System (IERS 1989) values.  Earth 
radius values are set (by values of parameters a = equatorial radius, km, and b = 
polar radius, km) in subroutine RIG_E10 in the models_E10.f file.  In Earth-GRAM 
2010, input values of altitude greater than 6000 km are treated as geocentric 
radius values, rather than heights.  Both radius and height are now given on the 
output file.  Although all input latitudes are geocentric, GRAM now gives both 
geocentric and geodetic values on the output file. A subroutine (radll_E10) 
has been added, which computes horizontal distance from great-circle distance 
between two input latitude-longitude positions.  This subroutine is used to
calculate lat-lon "radius" of current position from Range Reference Atmosphere
site locations, and to compute horizontal step size in the perturbation model.
                  

EXAMPLE USES OF GRAM IN TRAJECTORY CODE

Several features of the GRAM perturbation models are illustrated by sample 
"trajectory" codes provided with Earth-GRAM 2010.  These include:
(1) Feature to update atmospheric mean values without updating perturbation 
    values (trajopts, Option 3).
(2) A multiple-trajectory driver routine (multtraj), that allows multiple 
    trajectories and perturbations to be simulated in one GRAM run.
(3) A multiple-profile driver routine (corrtraj), that allows multiple 
    profiles and perturbations to be simulated in one stand-alone GRAM run, with 
    both small-scale and large-scale time correlations maintained between the 
    profiles.
(4) A multiple-profile driver routine (multbody), that allows multiple 
    profiles and perturbations to be simulated with GRAM used in a trajectory 
    code, with small-scale space and time correlations maintained between the 
    separate bodies.
    
These sample trajectory codes are provided to illustrate how to use GRAM as a 
subroutine in a user-provided trajectory program, and how to implement these 
features of the perturbation model. These sample trajectory programs include:
   gramtraj_E10.f90  A subroutine for use in user-provided trajectory programs
   trajdemo_E10.f90  A main driver (replaces gram_E10.f) for use in a simple 
                     trajectory demonstration program
   trajopts_E10.f90  A main driver (replaces gram_E10.f) for use in an example 
                     trajectory program that demonstrates several special options
                     of the perturbation model.
   multtraj_E10.f90  A main driver (replaces gram_E10.f) program that illustrates 
                     how to call Earth-GRAM 2010 to evaluate multiple trajectories 
                     in one program run, with independent (un-correlated) small-
                     scale correlations between trajectories. 
   trajcalc_E10.f90  Subroutines used in these trajectory demonstration programs,
                     illustrating how the user's trajectory program should
                     calculate and provide position and velocity information
                     to GRAM.
   corrtraj_E10.f90  A main driver (replaces gram_E10.f) program illustrating 
                     how to call Earth-GRAM 2010 to evaluate multiple profiles
                     in one program run, with small-scale correlations preserved
                     between the profiles. 
   multbody_E10.f90  Main driver illustrating how to call GRAM to evaluate
                     perturbations along trajectories of multiple bodies which
                     start out together and at various times separate from 
                     each other. The perturbations for each body start out the 
                     same, and then remain correlated over both position and 
                     time along the body trajectories, and across the spatial 
                     separation between the bodies.
                   
Program corrtraj is used to simulate a sequence of profiles over a period of 
(typically) a few hours to a few days.  Time variation of wave phases of the 
large-scale perturbation model insures that each profile in the time sequence 
has appropriate large-scale correlation with the previous profile in the sequence.  
Specific correlations between small-scale perturbations in the profiles are
maintained by perturbation correlation routines built into the corrtraj
driver code.  Program corrtraj is meant to be run interactively, with input
values requested for: (1) number of profiles to generate, (2) time step (hours)
between successive profiles, (3) initial position (height, geocentric latitude, 
longitude) for each profile, and (4) steps in time, height, geocentric latitude 
and longitude along each individual profile.  The corrtraj program can be used 
to simulate a sequence of observed profiles during day-of-launch operations, 
when conditions can be observed only up to some time before launch, and actual 
conditions encountered at launch time are correlated with, but not identical to, 
the last set of observed conditions.

The multbody program is similar to program corrtraj, with two notable differences:
  (1) corrtraj maintains time correlations between a sequence of simulated 
      profiles, while multbody maintains both time and space correlations
  (2) corrtraj is designed for stand-alone use only, while multbody is 
      designed to be dropped into a user's trajectory program.  Both corrtraj 
      and multbody make use of the "wrapper" routine gramtraj.  However, multbody 
      has moved much of the code from the main driver into subroutines, in order 
      to simplify use in a trajectory program.  A subroutine (crosscorr) has 
      been added to compute the spatial cross correlations. 

Applications for multbody include bodies which separate on entry (e.g. a heat 
shield that separates from the main spacecraft body) or which separate during 
launch (e.g. boosters from the final stage, or escape modules from the main 
spacecraft).  Simulations could also be done for separate parts of an intact
spacecraft that are assumed to experience different (but correlated) pertur-
bations (e.g. the nose, c.g. and tail of a long spacecraft).  The user can 
specify (with input parameter nmulti) up to 10 separate bodies to be tracked.  
Each body can have a different separation time and different separation rate
parameters (as specified in sample routine firstpos). Up until its separation 
time, each body is assumed to be at the same position (and to therefore 
experience the same perturbations) as the main body.

Both corrtraj and multtraj write a sequence of cross-correlated profiles to
files name Profnn.txt, where nn = 01, 02, ... to the number of correlated
profiles generated.  Headers for these files are:

    Header      Description
   =========    ==================================================    
    Hgtkm     = Height (km)
    GcLat     = Geocentric latitude (deg)
    ELon      = East longitude (deg)
    TpertK    = Perturbed (mean+perturbation) temperature (K)
    PpertNm2  = Perturbed (mean+perturbation) pressure (N/m**2)
    Dpertkgm3 = Perturbed (mean+perturbation) density (kg/m**3)
    UPert     = Perturbed (mean+perturbation) eastward wind (m/s)
    Vpert     = Perturbed (mean+perturbation) northward wind (m/s)
    TimeSec   = Time (sec) from initial point on first profile
    delT%     = Total deviation of temperature from mean (%)
    delD%     = Total deviation of density from mean (%)
    delUms    = Total deviation of eastward wind from mean (m/s)
    delVms    = Total deviation of northward wind from mean (m/s)
    Wpert     = Vertical wind (m/s)
    delTlrgK  = Large-scale temperature deviation (K)
    delUlrg   = Large-scale eastward wind deviation (m/s)
    delVlrg   = Large-scale northward wind deviation (m/s)
    delTsmlK  = Small-scale temperature deviation (K)
    delUsml   = Small-scale eastward wind deviation (m/s)
    delVsml   = Small-scale northward wind deviation (m/s)

NOTE: The "print format" and "special format" files from programs corrtraj and 
multtraj contain values for uncorrelated profiles (i.e. data before cross-
correlation effects between profiles are accounted for). Unless otherwise needed,
the print format and special format files can be suppressed by setting Namelist
paremeters iup and iopp to 0.

For further details on this and the other trajectory demonstration driver
programs, see comments within these code files.
                  
                 
NEW PARAMETERS IN NAMELIST FORMAT INPUT

Earth-GRAM 2010 has 10 new input parameters that are provided through the 
NAMELIST input file:

rralist:  File name for list of RRA sites (optional)
y10:      Solar X-Ray & Lya index scaled to F10 (for JB2008)
y10b:     Solar X-Ray & Lya 81-day avg. centered index (for JB2008)
dstdtc:   Temperature change computed from Dst index (for JB2008)
NCEPyr:   Period-of-record (POR=y1y2) for NCEP climatology. POR includes 
            years y1 through y2 (e.g. NCEPyr=9008 for POR = 1990 through 
            2008). NCEP monthly climatology is determined by input value 
            of month (mn) in initial time input
NCEPhr:   Code for UT hour of day for NCEP climatology: 1=00 UT, 2=06UT,
            3=12UT, 4=18UT, 5=all times of day combined, or 0 to use NCEP 
            time-of-day based on input UTC hour (ihro)
ruscale:  random perturbation scale for horizontal winds; nominal=1.0, 
            maximum=2.0, minimum=0.1
rwscale:  random perturbation scale for vertical winds; nominal=1.0, 
            maximum=2.0, minimum=0.1
z0in:     surface roughness (z0) for sigma-w model [ < 0 to use 1-by-1 deg 
           lat-lon surface data, from new file atmosdat_E10.txt;   = 0 for
           speed-dependent z0 over water; or enter a value between 1.0e-5 
           and 3 for user-specified z0 value ].   For more information, 
           see vertical wind model section, below.
ibltest:  unit number for boundary layer (B)L model output file (bltest.txt), 
            or 0 for no BL model output (see description below)
            
Two previous input parameters (iguayr and iyrrra) are no longer used, and 
should not be included in the namelist input file.

Further details are also described in file README2.txt.


NEW INFORMATION ON STANDARD FORMATTED OUTPUT

    A reference NAMELIST input file (NameRef.txt) and resulting standard
formatted output file (OutputRef.txt) "special-formatted" output file
(SpecialRef.txt), and species output file (SpeciesRef.txt) are supplied (see 
files README0.txt, README2.txt, and README3.txt). New information appearing on 
the standard formatted output file include:
(1) Ruv:    Cross-correlation between Eastward and Northward wind components
(2) SpdAv:  Monthly average wind speed (m/s)
(3) SpdSd:  Standard deviation of wind speed about monthly average (m/s)
(4) SoSav:  Average speed of sound (m/s)
(5) SoSpt:  Perturbed speed of sound (mean-plus-perturbation, m/s)


REVISED VERTICAL WIND MODEL (added for Earth-GRAM 2007 Ver 1.2, Updated for 2010)

This new model computes standard deviation for boundary-layer (BL) vertical wind 
(sigma-w) as a function of surface type (water or various land types), and 
"surface" horizontal wind (at 10 m height). Surface height (above Mean Sea Level, 
MSL) is interpolated from a 1-by-1 degree topographic data base in the new atmosdat_E10.txt 
file [or from Range Reference Atmosphere (RRA) surface altitude, if within the "zone of 
influence" of any RRA site].  BL depth (height of the top of the BL above the surface) is 
computed from a new time-of-day and stability-dependent model.  Land cover type is taken 
from a 1-by-1 degree resolution from data also in the new atmosdat_E10.txt file (DeFries, 
R. S. and J. R. G. Townshend, 1994, "NDVI-derived land cover classification at a global 
scale.", International Journal of Remote Sensing, 15:3567-3586; and Gates, W.L., and A.B. 
Nelson, 1975, "A New (Revised) Tabulation of the Scripps Topography on a 1deg Global Grid. 
Part I: Terrain Heights; Part II: Ocean Depths", Repts. 1276,1277, Rand Corp., Santa 
Monica, Calif.)

Surface Type codes are given in the following table.  Surface roughness length (z0) values 
assumed for each surface type were computed as the geometric mean value from a variety of 
sources.  Values for z0 are also given in the following table.  An option for user input
of surface roughness length is also provided.

     Code   Land Cover Class                                         z0 (m)
     -----  ---------------------------------------------------   ------------
     0      water                                                  u-dependent
     1      broadleaf evergreen forest                                0.6
     2      coniferous evergreen forest and woodland                  0.48
     3      high latitude deciduous forest and woodland               0.42
     4      tundra                                                    0.0056
     5      mixed coniferous forest and woodland                      0.45
     6      wooded grassland                                          0.12
     7      grassland                                                 0.046
     8      bare ground                                               0.015
     9      shrubs and bare ground                                    0.042
     10     cultivated crops                                          0.065
     11     broadleaf deciduous forest and woodland                   0.45
     12     data unavailable (re-assigned to Codes 4, 6, or 13,        
              as appropriate)                                         -.-
     13     ice                                                      3.2E-4

Wind dependence is on "surface hourly average" wind speed, computed from wind 
components given by monthly mean wind at 10 m plus GRAM large-scale perturbed wind 
components at the "surface". Since GRAM large-scale wind perturbations change from 
profile-to-profile in a Monte-Carlo run, then surface wind speed changes from 
profile-to-profile.

As given above, surface roughness (z0) is based on surface type array (from the 
atmosdat_E10.txt file), or for water uses  z0 = alpha x u*^2 / g, where alpha is a 
constant, u* is surface friction velocity, and g is acceleration of gravity.  
Friction velocity is computed from  

                 u*  = 0.4 x U10 / [Alog(10/z0) - psi(10/L) ] , 

where U10 is the "surface" wind speed (at 10 m height), L is the Monin-Obukhov scaling
length (stability-dependent), and the BL wind profile function psi is given by

       psi(10/L)  =  -50 / L                        if  1/ L > 0    (stable)
                       0.0                          if  1 / L = 0   (neutral)
                     1.0496 x ( -10 / L )**0.4591   if  1/L < 0     (unstable).
                 
The unstable formulation for psi is from Hsu et al. (1999), and is a simplification of
am often-used, but more complicated expression derived by Paulson (1970).

For water, the simultaneous dependence of z0 on u* and u* on z0 is solved by the 4-step 
iteration process.  File atmosdat_E10.txt also contains values of sigma-w at 5 km 
intervals (above MSL).  

Inverse of the Monin-Obukhov (1/L) is calculated by a multi-step process, based on
information in Tables 4-7 and 4-8 plus Figure 4-9 of Justus (1978):
(1) Compute a net radiation index (nri) that depends on solar elevation angle and time-of
    day (or night), where nri ranges from -3.5 (strong outgoing net radiation) to + 4.5 
    (strong incoming net radiation).  See Justus (1978), Table 4-7.
(2) Compute a wind-speed factor F(U10) from empirically-derived functions 

               F(U10) = (1 - U10 / 7.5)           if U10 < 6 m/s
                         0.2                      if U10 = 6 m/s
                         0.2 x Exp(12 - 2 U10 )   if U10 > 6 m/s

(3) Compute a stability category S as a function of F(U10) and nri by
   
               S  =  4.229 - nri x F
     
    which is an empirical fit to Justus (1978) Table 4-7.  Values of S are limited to
    0.5 on the low side (most unstable) and 7.5 on the high side (most stable).

(4) Compute inverse Monin-Obukhov length as a function of stability category S and
    surface roughness length z0, from
    
              1/L  =  (1/4) x ( -0.2161 + 0.0511 S ) Log10(10/z0)
              
     which is an empirical fit to Figure 4-9 of Justus (1978)

Standard deviation of vertical wind, sigma-w, is computed as a function of height
above the surface (z), and stability-dependent Monin-Obukhov length by relations:

         sigma-w  =  1.25 u* (1 + 0.2 z/L)           (stable; 1/L > 0)
                     1.25 u*                         (neutral; 1/L = 0)
                     1.25 u* ( 1 - 3 z/L)**(1/3)     (unstable; 1/L < 0)

where the stable relation is from Equation 1.33 of Kaimal and Finnigan (1994) and 
Pahlow et al. (2001), and the unstable relation is from Equation (2), Page 161, of
Panofsky and Dutton (1984).  For stable cases, sigma-w is limited to a value of
3.75 u*, while for unstable cases, it is limited by the magnitude of the convective
velocity w*, to a value of sigma-w < 0.62 w*, where w* is given by

                   w* = [ - h / (0.4 L) ]**(1/3)
                   
(Equation 4 of Panofsky, 1978), where h is the stability-dependent boundary layer 
depth, calculated from simplifications of methodologies given by Sugiyama and Nasstrom
(1999), Batchvarova and Gryning (1991), and Seibert (2000).

Subject to these limiting values, the above height-dependent equations are used to 
compute sigma-w from the surface to the top of the boundary layer.  Between the top of 
the BL and the next height for which sigma-w is available from the atmosdat file, 
linear interpolation is used to estimate sigma-w.

Calculation of the vertical wind perturbations is not changed, only the methodology
for computing standard deviation of vertical wind.  Values of sigma-w are constrained 
to be 0.1 m/s or greater, since the perturbation calculation methodology does not 
work properly if sigma-w is zero.

For more details on the new vertical wind model, including complete reference for
the sources noted above, see file VerticalWindModel.doc, on the GRAM distribution DVD.


VERTICAL WIND MODEL BOUNDARY LAYER OUTPUT FILE (BLTest.txt)

A new (optional) output file, called BLTest.txt, is provided, for examination of,
or graphing of parameters from the new vertical wind boundary layer (BL) model.  
Headers in BLTest.txt are as follows:

Header     Description
------     -----------------------------------------------------------------------
Day:       Day of month (including fractional part) 
LST:       Local solar time (hours)   
Hgt:       Height above mean sea level (MSL, km)    
Lat:       Geocentric latitude (degrees)      
Lon:       East longitude (degrees)  
hsrf:      Surface height (MSL, km) 
spdsrf:    Surface (10m) wind speed (from mean-plus-large-scale-perturbation, m/s)   
LC:        Land surface type code (see table above)   
z0:        Surface roughness length (m)   
Elmn:      Solar elevation at midnight (degrees)   
El:        Current solar elevation angle (degrees)   
Elmd:      Solar elevation at mid-day (degrees)     
sha:       Current solar hour angle (degrees) 
blfct:     Height factor for unstable BL during early daytime (unitless) 
nri:       Net radiation index (unitless), for stability calculation  
S:         Atmospheric stability category (unitless)       
ool:       Inverse of Monin-Obukhov scale length (1/meters)  
ustar:     Surface friction velocity (u*, m/s)    
BVfsq:     Square of Brunt-Vaisala frequency (sec**-2)   
hN:        Depth of neutral boundary layer (meters)   
hbl:       Current depth of boundary layer (meters)   
chb:       Current height above surface (m) or height to top of boundary layer (for 
             evaluating height variation of standard deviation of vertical wind)
sigrat:    Ratio of sigma-w to ustar at current height
swb:       Vertical wind standard deviation at top of boundary layer (m/s)
swh:       Current standard deviation of vertical winds (m/s)
spdavsrf:  Monthly average wind speed at surface (10m) (m/s)
spdsdsrf:  Standard deviation of surface wind speed (m/s)
tsrf:      Monthly average surface temperature (K)
stsrf:     Standard deviation of surface temperature (K)

For additional information, see vertical wind model description, above, or
file VerticalWindModel.doc, in the DOCUMENTATION folder.



      FEATURES AND OPTIONS INTRODUCED IN Earth-GRAM 2007 OR EARLIER


REVISED PERTURBATION MODEL

In order to produce more obvious "intermittency" in the perturbation
model, earlier GRAM models introduced revisions to details of calculations 
of the variable-scale random perturbation model.
(1) The time-series simulations of the variable length scale is used to
categorize the turbulence as "normal" or "perturbed".  The turbulence
(wind, density, temperature, etc.) is in disturbed conditions whenever
the length scale drops below a prescribed "minimum" value (described in the
GRAM-95 report, NASA TM 4715).  The probability of being in disturbed 
conditions is taken from statistics in NASA Tech. Memo. 4168 (1990), and 
varies from 1 to 2.5% near the surface to about 0.15% near 25 km altitude, 
to about 2% near 75 km, and back to about 1% above 120 km altitude.  The 
values for standard deviation of the length scales were modified in the 
"atmosdat" data file (described in NASA TM 4715) to get these appropriate 
probability values.
(2) The same total standard deviation values , sigma(total), are retained 
(for density, temperature, wind etc.), but they are now divided into two
values, a "small" standard deviation for normal conditions, sigma(small),
and a "large" standard deviation, sigma(large), for disturbed conditions.
The standard deviations follow the relation

 [sigma(total)}**2 = Pdist*[sigma(large)]**2 + (1-Pdist)*[sigma(small)]**2

where Pdist is the above-mentioned probability for disturbed conditions.
A parameterization (based on data in NASA TM 4168) specifies the ratio
of sigma(large) to sigma(total).  From these two conditions, values for
both sigma(large) and sigma(small) are computed.
(3)  Instead of allowing the perturbation length scale to vary contin-
uously, length scales are set to either a small value (the minimum length
scale, Lmin) if conditions are perturbed, or a large value (Lmax) if
conditions are normal.  The value of Lmin is taken directly from the Lmin
value in the "atmosdat" data file (see Section 4.3 of NASA TM 4715).  The
value of Lmax is determined from Lmin and Lbar (the average scale size
from the atmosdat data file, by requiring that

       Lbar  =  Pdist*Lmin  +  (1 - Pdist)*Lmax

where Pdist is the above-mentioned probability for disturbed conditions.
(4)  Perturbation magnitudes are larger during disturbed periods, because
of the larger standard deviation that applies. "Shears", changes in  wind
(or density or temperature etc.), are also larger during disturbed
conditions than in normal periods, both because the standard deviation is
larger and because the length scale is smaller (being the minimum value).
For a discussion of the effect of length scale on shears, see Section 2.6
of NASA TM 4715.


OPTIONAL SURFACE DATA OUTPUT

    As did GRAM-95, GRAM-99 provides a mechanism for the user to easily
modify the code to produce a user-defined "special format" output file. 
Features and output variable selections are discussed in Section 4.6 of 
NASA TM-4715. An Earth-GRAM 2010 feature is that data values at the terrain
surface can also easily be output.   In GRAM-99, surface data were
always provided from the GUACA data, even if the Range Reference Atmosphere
(RRA) option was being used.  For Earth-GRAM 2010, NCEP data provide the
surface values, with RRA data providing surface values if the RRA option 
is selected and the current location is within the "zone of influence" of
an RRA site.  Details of all available surface data output parameters are 
given in file README3.txt.


WAVE MODEL FOR LARGE-SCALE PERTURBATIONS

A two-scale perturbation model is used in GRAM.  Prior to GRAM-99, perturba-
tions at both large and small scales were both computed by a one-step 
Markov process (a 1st order autoregressive approach, equivalent to the
1st order autoregressive model of Hickey in NASA Contractor Report 4605). 
The combination of small-scale and large-scale one-step Markov processes 
in GRAM is somewhat more general than the 2nd order autoregressive model 
used by Hickey at high latitudes.  However, because the one-step Markov 
process for large-scales can yield somewhat inaccurate results if small 
steps are taken along a trajectory, a  wave model was developed to use
for the large scale perturbations in GRAM-99.  This model uses a cosine 
wave, with specified horizontal and vertical wavelengths. In Earth-GRAM
2007, a process was introduced for randomly varying horizontal and vertical 
wavelengths used in the large-scale perturbations.


OPTION FOR USER-SELECTED INITIAL PERTURBATIONS

Earlier versions of GRAM (through GRAM-90) required the user to specify
initial values for all of the random perturbation variables.  For typical
Monte-Carlo applications, most users found it difficult to choose different
starting values for each set of perturbed profiles to be generated.  Using
the same starting perturbation values for all profiles (e.g. 0.0) was
feasible, but not always the most appropriate thing to do.  GRAM-95 was
changed to have the program automatically select (with the appropriate
range of variability) a random starting value for each perturbation
parameter and each Monte-Carlo perturbation profile. GRAM-99 allowed, as a 
user-controlled option, the input of user-selected initial perturbation 
values.  This option is controlled by the input parameter initpert (with 0, 
the default value, meaning GRAM-selected random initial values; and 1 
triggering user-selected initial values).  An example application for user-
selected initial perturbations would be the following:  Suppose a measured 
profile (e.g. a day-of-launch atmospheric sounding) is to be used as an 
actual (perturbed) profile up to the highest measured altitude, and that an 
Earth-GRAM perturbed profile (or profiles) is desired for higher altitudes.  
The actual atmospheric values from the measured profile can be used to 
compute perturbed values with which to initialize Earth-GRAM (by methods 
discussed more fully below).  With these initial perturbation values, GRAM 
perturbed profile(s) for the higher altitude region will have complete 
continuity with the measured profile from the lower altitude region.  For
additional discussion, see comments in the trajopts_E10.f code file.


EXAMPLE CALCULATION OF USER-SELECTED INITIAL PERTURBATIONS

Suppose, in the example application of the previous section, that the
measured profile extends to an altitude of 25 km.  Further suppose that
the measured values at 25 km are: density = 0.044 kg/m**3 and eastward wind
= 14 m/s (along with other measured values).  Next do an Earth-GRAM 2010 run, 
starting at 25 km in this example, with the same latitude, longitude, and 
month as the measured profile.  Suppose this Earth-GRAM run gives a mean 
density of 0.040 kg/m**3 and mean eastward wind of 20 m/s.  The values of 
user-selected initial density and eastward wind perturbations to use in 
subsequent Earth-GRAM runs would be:

   initial density perturbation (%) =  100*(measured density - Earth-GRAM 
                 mean)/Earth-GRAM mean = 100*(0.044-0.040)/0.040 = 10 (%)

   initial eastward wind perturbation (m/s) = measured eastward wind -
                    Earth-GRAM mean eastward wind =  14 - 20  = -6 m/s.

Similar calculations would apply to the pressure, temperature, and other
wind components.  For subsequent Earth-GRAM runs (in this example case), 
start at 25 km, using initpert = 1, rdinit = 10.0, ruinit = -6.0 (and 
whatever values apply to the other initial perturbation components). These 
values will insure that the highest altitude (25 km) data point in the 
measured (perturbation) profile is consistent with the starting value (at 
25 km) of the perturbed profile from Earth-GRAM.  If multiple perturbed 
profiles are desired, for a Monte Carlo application, each of them will be 
initialized with the same values (rdinit = 10.0 and ruinit = -6.0, etc., 
in this case).


