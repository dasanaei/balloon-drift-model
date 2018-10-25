
                    DESCRIPTION OF EARTH-GRAM 2010 INPUT FILES

SAMPLE NAMELIST FORMAT INPUT FILE

    Following is a sample NAMELIST format input file (NameRef.txt) for use 
by Earth-GRAM 2010.  Note that some compilers use different formats for the 
beginning and ending lines of the file.  This sample file contains values 
required to produce the reference output data (file OutputRef.txt). These
input values are also the default values, set if none are provided.  Only 
values that differ from these default values actually need to be input in the 
NAMELIST file.
    Definitions and discussion of Earth-GRAM 2010 input parameters that are 
the same as GRAM-99 input parameters are given in Appendix D of NASA/TM-1999-
209630.

See also general introductory discussion about Range Reference Atmospheres 
in file README5.txt, about MET, MSIS, and JB2008  thermosphere models in 
README6.txt, about auxiliary profiles and about other new features and 
options in file README7.txt.

 $namein_e10
  atmpath = 'D:\GRAMs\EarthGRAM2010V2.0\PC_IOfiles\atmosdat_E10.txt'
  NCEPpath = 'D:\GRAMs\EarthGRAM2010V2.0\NCEPdata\BinData\'  
  trapath = 'null'
  prtpath = 'output.txt'
  nprpath = 'special.txt'
  conpath = 'species.txt'
  rndpath = 'null'
  rrapath = 'D:\GRAMs\EarthGRAM2010V2.0\RRAdata\'
  rralist = 'rras2006.txt'
  profile = 'null'
  h1 = 140.
  phi1 = 0.45
  thet1 = -164.53
  f10 = 230.
  f10b = 230.
  ap = 20.3
  s10     = 0. 
  s10b    = 0.
  xm10    = 0.
  xm10b   = 0.
  y10     = 0.
  y10b    = 0.
  dstdtc  = 0.
  mn = 1
  ida = 1
  iyr = 2010
  ihro = 0
  mino = 0
  seco = 0.0
  dphi = 0.4
  dthet = 1.2
  dhgt = -2.0
  nmax = 71
  delt = 60.0
  iopt = 0
  iopp = 17
  iu0 = 0
  iup = 6
  ius = 3
  iuc = 4
  iug = 22
  NCEPyr = 9008
  NCEPhr = 5
  iopr = 1
  nr1 = 1234
  iun = 0
  rpscale = 1.0
  ruscale = 1.0
  rwscale = 1.0
  iurra = 35
  iyrrra = 1
  sitelim = 2.5
  sitenear = 0.5
  initpert = 0
  rdinit = 0.
  rtinit = 0.
  ruinit = 0.
  rvinit = 0.
  rwinit = 0.
  patchy = 0.
  itherm = 1
  z0in   = -1.
  ibltest = 99
 $End
Parameter Descriptions:
atmpath  = path name for "atmosdat" atmospheric data file
NCEPpath = path name for NCEP data files
trapath  = path name for trajectory input file ('null' if none)
prtpath  = path name for standard formatted output file  ('null' if none)
nprpath  = path name for the "special" format output file ('null' if none)
conpath  = path name for species concentration output file ('null' if none)
rndpath  = path name for file containing more  random number seeds (optional)
rrapath  = DIRECTORY for Range Reference Atmosphere (RRA) data (optional)
rralist  = File name for list of RRA sites (optional)
profile  = path name for auxiliary profile data ('null' if none)
h1       = initial height (km).  Heights > 6000 km are interpreted as radius.
phi1     = initial geocentric latitude (degrees, N positive)
thet1    = initial longitude (degrees, East positive)
f10      = daily 10.7-cm flux
f10b     = mean 10.7-cm flux
ap       = geomagnetic index (Note: Valid ap must be used if JB2008 selected,
             for use in HWM wind model)
s10      = EUV index (26-34 nm) scaled to F10 units (0.0 -> s10=f10)
s10b     = EUV 81-day center-averaged index (0.0 -> s10b = f10b)
xm10     = MG2 index scaled to F10 units (0.0 -> xm10 = f10)
xm10b    = MG2 81-day center-averaged index (0.0 -> xm10b = f10b)
y10      = Solar X-Ray & Lya index scaled to F10 (0.0 -> y10=f10)
y10b     = Solar X-Ray & Lya 81-day avg. centered index (0.0 -> y10b=f10b)
dstdtc   = Temperature change computed from Dst index (for JB2008)
mn       = month (1-12)
ida      = day of month
iyr      = 4-digit year, or 2-digit year: >56=19xx <57=20XX
ihro     = initial UTC (Greenwich) time hour (0-23)
mino     = initial UTC (Greenwich) time minutes (0-59)
seco     = initial UTC (Greenwich) time seconds (0.0-60.0)
dphi     = geocentric latitude increment (degrees, Northward positive)
dthet    = longitude increment (degrees, Eastward positive)
dhgt     = height increment (km, upward positive). If radius input is used,
             dhgt is interpreted as a radius increment.
nmax     = maximum number of positions (including initial one; 0 means read
             trajectory input file)
delt     = time increment between positions (real seconds)
iopt     = trajectory option (0 = no trajectory data; otherwise unit number
             for trajectory input file)
iopp     = "special" output option (0 = no "special" output; otherwise unit
             number of "special" output file)
iu0      = unit number for screen output (normally 6 or 0)
iup      = unit number for standard formatted output file (0 for none)
ius      = unit number for atmosdat data
iuc      = unit number for concentrations output (0 for none)
iug      = unit for NCEP data files
NCEPyr   = y1y2  to use NCEP climatology for period-of-record (POR) from year 
             y1 through year y2 (e.g. NCEPyr=9008 for POR = 1990 through 
             2008). NCEP monthly climatology is determined by input value 
             of month (mn) in initial time input
NCEPhr   = Code for UT hour of day if NCEP climatology is used: 1=00 UT, 
             2=06UT, 3=12UT, 4=18UT, 5=all times of day combined, or 0 to 
             use NCEP time-of-day based on input UTC hour (ihro)
iopr     = random output option (1 = random output, 2 = none)
nr1      = first starting random number (1 to 9 * 10**8)
iun      = unit number for more starting random numbers (0 for none)
rpscale  = random perturbation scale for density, temperature and pressure; 
             nominal=1.0, max=2.0, min=0.1
ruscale  = random perturbation scale for horizontal winds; nominal=1.0, 
             maximum=2.0, minimum=0.1
rwscale  = random perturbation scale for vertical winds; nominal=1.0, 
             maximum=2.0, minimum=0.1
iurra    = unit number for Range Reference Atmosphere (RRA) data (0 if none)
iyrra    = 1 for 1983 RRAs, 2 for 2006 RRAs, 3 for 2013 RRAs
sitelim  = lat-lon radius (deg) from RRA site, outside which RRA data are
             NOT used.  Also used, with a similar meaning, for auxiliary
             profile input.  Note that RRA and auxiliary profile input
             cannot be used simultaneously.
sitenear = lat-lon radius (deg) from RRA site, inside which RRA data is
             used with full weight of 1 (smooth transition of weight factor
             from 1 to 0 between sitenear and sitelim). Also used, with a 
             similar meaning, for auxiliary profile input.
initpert = Use 1 for user-selected initial perturbations or 0 (default) for
             GRAM-derived, random initial perturbation values
rdinit   = initial density perturbation value (% of mean)
rtinit   = initial temperature perturbation value (% of mean). Note - initial
            pressure perturbation is computed from rdinit and rtinit.
ruinit   = initial eastward velocity perturbation (m/s)
rvinit   = initial northward velocity perturbation (m/s)
rwinit   = initial upward velocity perturbation (m/s)
patchy   = not equal 0 for patchiness; 0 to suppress patchiness in 
             perturbation model
itherm   = 1 for MET (Jacchia), 2 for MSIS, or 3 for JB2008 thermosphere   
z0in     = surface roughness (z0) for sigma-w model [ < 0 to use 1-by-1 deg 
            lat-lon surface data, from file atmosdat_E10.txt;   = 0 for
            speed-dependent z0 over water; or enter a value between 1.0e-5 
            and 3 for user-specified z0 value ].   For more information, 
            see file README7.txt.
ibltest  = unit number for BL model output file (bltest.txt), or 0 for no  
            BL model output


OPTIONAL TRAJECTORY INPUT FILES

An optional pre-computed input trajectory file of positions and times 
can be provided by setting the trajectory file pathname (trapath) in the 
NAMELIST input file to whatever file name is desired, setting the number 
of computed positions (nmax) to zero, and setting the trajectory file unit 
number (iopt) to any non-conflicting value.  The trajectory file contains:
one time and position per line, having values of:
(1) elapsed time (seconds), 
(2) height (km), 
(3) geocentric latitude (degrees, North positive), and
(4) longitude (degrees, East positive).  
Trajectory input heights greater than 6000 km are treated as geocentric 
radius values, rather than as altitudes above reference ellipsoid.  Example 
programs illustrating how to drop Earth-GRAM into trajectory code cannot 
be run with trajectory file input, since trajectory positions are generated 
"on the fly" in these programs.  Any additional REFERENCE information (e.g. 
orbit number, measured density, etc.) included at the end of each line, is 
ignored.  Trajectory positions in these files do not have to be at small 
time or space increments.  For example, a trajectory file may consist of 
successive periapsis times and positions for a simulated (or observed) 
aerobraking operation.  Trajectory files may also contain arrays of 
locations used for computing height-latitude cross sections or latitude-
longitude cross sections.  Such trajectory input files can be as built by 
program bldtraj.  For additional discussion of the bldtraj program, see 
file README8.txt.

A special trajectory file (RRAtraj.txt, in folder PC_IOfiles) is provided which 
contains height, latitude and longitude information for all of the Range 
Reference Atmosphere (RRA) sites.  See further description of trajectory files 
in README8.txt.  Using RRAtraj.txt as an input trajectory file provides an easy 
method for computing (in one Earth-GRAM run) output values for means and standard 
deviations from all of the RRA sites, and comparing these with values computed 
from a similar run using NCEP climatology (described in README4.txt).


OPTIONAL AUXILIARY PROFILE INPUT FILE

As an option, data read from an auxiliary profile may be used to replace
data from the conventional (NCEP/MAP/etc.) climatology. This option is 
controlled by setting parameters profile, sitenear, and sitelim in the 
NAMELIST input file.  Each line of the auxiliary profile input file consists 
of: (1) height, in km [height values greater than 6000 km are interpreted
as radius values, in km], (2) geocentric latitude, in degrees, (3) longitude, in 
degrees (East positive), (4) mean temperature, in K, (5) mean pressure, in N/m**2, 
(6) mean density, in kg/m**3, (7) mean Eastward wind, in m/s, (8) mean Northward 
wind, in m/s, (9) standard deviation (sigma) temperature, in K, (10) sigma 
pressure, in N/m**2, (11) sigma density, in kg/m**3, (12) sigma Eastward wind, 
in m/s, (13) sigma Northward wind, in m/s.  Heights are relative to the 
reference ellipsoid, except that values greater than 6000 km are interpreted 
as radius values, rather than altitudes. Latitudes are geocentric.  Regular 
climatological values (NCEP/MAP data etc.) are used if any of the three values
for temperature, pressure, or density data are input as zero in the auxiliary 
profile.  Regular climatological values of wind components are used if BOTH 
wind components are zero in the auxiliary profile file.  A sample auxiliary 
profile file (named KSCannRRA.txt) is provided.  See additional discussion in 
file README7.txt.

A special set of auxiliary profiles is provided (in folder AuxProfiles), 
constructed from values of annual means and standard deviations out of the Range 
Reference Atmosphere (RRA) data tables.  These files provide an easy method for 
users to simulate annual-average conditions at any of the RRA sites.
