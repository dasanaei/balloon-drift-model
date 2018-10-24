               UTILITY PROGRAMS FOR USE WITH EARTH-GRAM 2010
            
Below are described four utility programs, for use with Earth-GRAM 2010:
(1) NCEPbinF, (2) bldtraj, (3) timetraj, and (4) makerand. Program source 
code and PC-executables are provided, in the Utilities_Code folder.

                              Program NCEPbinF

This program reads the ASCII format NCEP files (adjusted according to procedures 
described in NCEPmods.pdf, in the DOCUMENTATION folder) and writes them in 
binary, for use by Earth GRAM-2010.  See file README4.txt for additional 
discussion.  NCEP binary data files are provided on the distribution DVD.
These should work as-is on a PC when using the PC executable programs provided. 
If necessary, program NCEPbinF is used to generate NCEP binary files from the 
NCEP ASCII files provided.  To generate NCEP binary files, compile
NCEPbinF.f90 with commands such as (on a Unix platform):

   f90 -static NCEPbinF.f90
   mv a.out NCEPbinF.x
   rm NCEPbinF.o
   
Note: the "-static" option may or may not be necessary, depending on whether
the array sizes in NCEPbinF exceed the default stack size on your computer.
If the static option is required, see instructions for your specific compiler
on how to set this option.  Run program NCEPbinF.x for each of the 12 months 
(input files Nfyyyymm.txt, where yyyy = period-of-record, e.g. 9008, and mm=
month), and put the resultant binary files (named Nbyyyymm.bin) in the folder 
whose path name is specified by input parameter NCEPpath in the NAMELIST 
input file.
              
                              Program bldtraj

This is a program to build pseudo-trajectory file for using in Earth-
GRAM 2010, to compute output for maps or cross-sections.

It is frequently desirable to produce Earth-GRAM output for graphing as a
map (i.e. lat-lon cross section at a given height) or other cross-section
(e.g. height-lat cross section at a given longitude).  Program bldtraj
generates a "trajectory" file (with input lines containing time, height,
geocentric latitude, and longitude) that can be used as Earth-GRAM input for
generation of such maps or cross-sections.  Program bldtraj is interactive
and prompts the user to input starting values, ending values, and step
sizes for height (z1,z2,dz), geocentric latitude (lat1,lat2,dlat), and 
longitude (lon1,lon2,dlon).  The program also prompts for a value of time 
increment which is applied between each "trajectory" step (the time increment 
may be 0, if all trajectory points are at the same time). Time values in the
trajectory file are time (seconds) from the start time specified by date
and time information provided in the Earth-GRAM NAMELIST-format input file.
NOTE: If heights > 6000 km are used, they are interpreted as planeto-centric
radius values (km). A new input option for bldtraj is for input height values
to be in either above ground level (AGL) or above mean sea level (MSL). If the
AGL option is selected, then surface heights can be selected to come only from 
the 1-by-1 degree topographic data base, or from the topo data base modified by
surface heights of any Range Reference Atmosphere sites that are within selected
"zone of influence" limits.  If the AGL option is selected, an input file named
"PathNames.txt" must be provided, containing: (1) the path to the atmosdat_E10.txt 
file, (2) the path to the folder containing the Range Reference Atmosphere (RRA) 
data, and (3) the file name (in the RRA folder) containing information of 
available RRA sites.


Example:
For a lat-lon map at height 10 km (above mean sea level), between geocentric 
latitudes -30 and 30 degrees (in steps of 5 degrees), and longitudes 0 to 180 
degrees (in steps of 20 degrees), enter 10 10 0 for z1, z2, dz; enter -30 30 5 
for lat1, lat2, dlat; and enter 0 180 20 for lon1, lon2, dlon.  All of these
input quantities are of type real, and can be entered to one or more significant 
digits beyond the decimal.  For use with the AGL input height option, a sample
"PathNames.txt" file is provided in the Utilities_Code folder.

              
                              Program timetraj  

This is a program to generate a pseudo-trajectory file consisting of any 
number of profiles separated by a user-selected time step.
                
This program generates a sequence of "trajectory" profiles, each separated by
a specified time step.  Note that it is not necessary to generate such a 
sequence of profiles to run the correlated profile program (corrtraj), since
corrtraj generates its own profile positions and time from user-provided
inputs.

The timetraj program runs interactively.  User inputs required are:

(a) Time step (seconds) between successive output "trajectory" profiles.
(b) Number of time steps in the sequence of profiles to be generated (one step
    generates two profiles; two steps generates three total profiles, etc.).
(c) Whether the input conditions are provided from a "pseudo-trajectory" file,
    (consisting of a sequence of time, height, geocentric latitude, and longitude 
    values), or are generated as a "linear" profile from user-specified steps in 
    height, geocentric latitude, and longitude.
    
If the program starts from an input "trajectory" file, further inputs required
are:

(a) Name of input trajectory file
(b) Name of output time-sequence "trajectory" file

If the program takes input from user-specified profile parameters, the 
additional inputs required are:

(a) Name of output time-sequence "trajectory" file
(b) Initial (starting) time for first profile in the sequence
(c) Starting and ending heights (km) and height step (km) for each profile
(d) Starting geocentric latitude and geocentric latitude increment for each 
     profile (degrees)
(e) Starting longitude and longitude increment for each profile (degrees)


                               Program makerand

This is a program to generate and write a file of random seed numbers for the 
GRAM perturbation model.  These files are used in Monte-Carlo simulations, with
file name specified as input parameter rndpath in the namelist input file.
                 
The makerand program runs interactively.  User inputs required are:

(a) Number of seed values to generate (number of seed values written to output 
    file is one less than this number because the first seed value is taken 
    from the namelist input file).
(b) First seed number to use (integer)
(c) Step size (integer) to use in generating seed values.  Final seed value
    [first seed + (step size) times (number of steps)] should not be larger
    than 900,000,000
(d) Name of output file for random seed values
