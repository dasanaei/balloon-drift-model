          Instructions for Compiling and Running Earth GRAM-2010


General Introduction

Earth-GRAM 2010 can be run in stand-alone mode (i.e. as an independent
program), or as a set of subroutines in a user-provided program, such
as a trajectory code.  The following material gives several examples
of how to compile and run Earth-GRAM 2010 in each of these modes.
Examples are given for both PC and UNIX machines.  If you wish to run
stand-alone GRAM on a PC, and do not need to make any code changes,
you can use PC executable files provided on the distribution DVD. It
is necessary to compile the Earth-GRAM 2010 source code -
(1) if you wish to run on a UNIX or other non-PC platform,
(2) if you make changes in the code (e.g. by selecting other/additional 
    parameters for output, as described in README3.txt), or
(3) if you want to incorporate Earth-GRAM 2010 as subroutines in your 
    own code, such as a trajectory program.
Unless otherwise noted below, it is assumed in the following material 
that all files from the distribution DVD have been copied to folders as 
named in the instructions for how to set-up Earth-GRAM 2010, given in 
file README0.txt.  Brief descriptions of all program and data files are 
also given in file README0.txt.


Command-Line Compiling on a PC

One example of how to compile Earth-GRAM 2010 on a PC is to use command
line mode.  Command line mode is initiated by clicking on "Start", then
"Run", then entering "CMD.EXE" (without the quotes) in the run window.
For the GNU PC Fortran complier (gfortran), stand-alone Earth-GRAM can  
be compiled (to an executable file named gram_E10.exe) by entering the 
following commands -

gfortran -c Cfiles_E10_C.F90 
gfortran -c Ifiles_E10_I.F90
gfortran -c GRAMSUBS_E10.F90
gfortran -c INITIAL_E10.F90
gfortran -c MET07PRG_E10.F90
gfortran -c MODELS_E10.F90
gfortran -c RANDOM_E10.F90
gfortran -c SPECONC_E10.F90
gfortran -c RRAmods_E10.F90
gfortran -c MSISsubs_E10.F90
gfortran -c HWMsubs_E10.F90
gfortran -c JB2008_E10.F90
gfortran -c NCEPsubs_E10.F90
gfortran GRAM_E10.F90 Cfiles_E10_C.o Ifiles_E10_I.o GRAMSUBS_E10.o INITIAL_E10.o 
  MET07PRG_E10.o MODELS_E10.o RANDOM_E10.o SPECONC_E10.o RRAmods_E10.o 
  MSISsubs_E10.o HWMsubs_E10.o JB2008_E10.o NCEPsubs_E10.o -o gram_E10.exe

Alternatively, these commands may be put into an executable batch (.bat) file, 
in which case it is advisable to include a pause statement at the end of the
batch file, so the temporary command window remains open after compilation has
completed, allowing the user to view any compilation error messages that
may have been produced.  Optionally, the last gfortran call may include any of 
several switches to get various run-time error messages, such as

gfortran -fbacktrace -std=f95 -Wextra -Wall -pedantic -fbounds-check ...

Other compilers will have different commands and compile switches for doing 
command-line compile operations.  Compilers also generally provide a "Programmer's
Workbench" or similar graphical-user-interface mode for doing program compile-
and-link operations.

Subroutines trajcalc_E10.f90 and gramtraj_E10.f90 and a dummy main driver
(trajdemo_E10.f90) are provided, to give an example of how to incorporate 
and use Earth-GRAM 2010 as subroutines in your own code.  To compile this example 
"trajectory" code (using gfortran), enter each of the following command lines 
(or incorporate them into an executable batch file) -

To compile trajdemo:

gfortran -c  Cfiles_E10_C.F90 
gfortran -c Ifiles_E10_I.F90
gfortran -c gramtraj_E10.F90
gfortran -c trajcalc_E10.F90
gfortran -c GRAMSUBS_E10.F90
gfortran -c INITIAL_E10.F90
gfortran -c MET07PRG_E10.F90
gfortran -c MODELS_E10.F90
gfortran -c RANDOM_E10.F90
gfortran -c SPECONC_E10.F90
gfortran -c RRAmods_E10.F90
gfortran -c MSISsubs_E10.F90
gfortran -c HWMsubs_E10.F90
gfortran -c JB2008_E10.F90
gfortran -c NCEPsubs_E10.F90
gfortran trajdemo_E10.F90 gramtraj_E10.o trajcalc_E10.o Cfiles_E10_C.o 
  Ifiles_E10_I.o GRAMSUBS_E10.o INITIAL_E10.o MET07PRG_E10.o MODELS_E10.o 
  RANDOM_E10.o SPECONC_E10.o RRAmods_E10.o MSISsubs_E10.o HWMsubs_E10.o 
  JB2008_E10.o NCEPsubs_E10.o -o trajdemo_E10.exe

NOTE - After successful compilation to executable programs, all object files and 
module files can be deleted, such as by the commands

erase *.o
erase *.mod

These commands may also be incorporated into the executable compile batch files.

Program trajdemo illustrates a basic trajectory program implementation.  
Other example trajectory driver programs, offering special options and
features, are available upon approved request.  These include: (1) Program 
trajopts, which shows how to implement various trajectory options, (2) Program 
multtraj, which illustrates implementation for calculating multiple trajectories 
in one program run (without cross-correlation between perturbations on the
separate trajectories), (3) Program corrtraj, which evaluates multiple profiles 
in one program run, with time cross-correlations preserved between the 
profiles, and (4) Program multbody (including associated subroutines), which 
provides as an example of how to use GRAM in a trajectory code to compute space 
and time-correlated perturbations influencing multiple spacecraft bodies.
Extensive comments in these example trajectory program driver routines provide 
more details.


Compiling on a UNIX Machine

The UNIX f90 command can be used to compile Earth-GRAM 2010, to
an executable named gram_E10.x, by entering the following -

f90 Cfiles_E10_C.f90 Ifiles_E10_I.f90 gram_E10.f90 gramsubs_E10.f90 /
   NCEPsubs_E10.f90 initial_E10.f90 MET07prg_E10.f90 models_E10.f90 /
   random_E10.f90 speconc_E10.f90 rramods_E10.f90 MSISsubs_E10.f90 /
   HWMsubs_E10.f90 JB2008_E10.f90
mv a.out gram_E10.x

Similarly, the f90 command may also be used to compile the trajdemo
program that provides an example of how to incorporate Earth-GRAM 2010
in user code -

f90 Cfiles_E10_C.f90 Ifiles_E10_I.f90 trajdemo_E10.f90 trajcalc_E10.f90 /
   gramtraj_E10.f90 gramsubs_E10.f90 NCEPsubs_E10.f90 initial_E10.f90 /
   MET07prg_E10.f90 models_E10.f90 random_E10.f90 speconc_E10.f90 /
   rramods_E10.f90 MSISsubs_E10.f90 HWMsubs_E10.f90 JB2008_E10.f90
mv a.out trajdemo_E10.x

NOTE - The forward slash character indicates a continuation line.  
NOTE - After successful compilation to executable programs, all object files 
and module files can be deleted, such as by the commands

rm -f *.o
rm -f *.mod

The above compilation commands and rm commands may also be incorporated
into an executable batch file.

Note: If the fullwarn compile option is used on the UNIX f90 compiler, several
CAUTION or NOTE messages may result.  Specifically, there may be several 
CAUTIONS about the VAST_KIND_PARAM module having "already been directly or
indirectly use associated into" the current scope. These CAUTION and NOTE 
messages can safely be ignored.  However, no Error or Warning messages should 
be produced by the compilation process.

Running in Command-Line Mode on a PC

If the setup procedure discussed in README0.txt is used, PC executable
files are placed in folder C:\EarthGRAM10\PC_Executables, and example 
input/output files are placed in folder C:\EarthGRAM10\PC_IOfiles.
Earth-GRAM 2010 can easily be run from any directory.  For example,
to run from a directory named C:\MyTest, do the following -

(1) Copy (and edit and/or rename, if desired) the file 
    C:\EarthGRAM10\PC_IOfiles\NameRef.txt to the C:\MyTest folder.
    Note that NAMELIST input file names are limited to 99  characters.
(2) Open a PC command-line window by Clicking Start and Run and
    entering the command CMD.EXE
(3) Change directory by entering the commands C: and CD \MyTest
(4) Execute Earth-GRAM 2010 by entering the command
     C:\EarthGRAM10\PC_Executables\gram_E10.exe
(5) Enter the name of the input file (from step 1) when prompted
(6) Output files will appear in the C:\MyTest folder

If file open errors are encountered, make sure that file pathnames
are properly set in the NAMELIST input file.  For files residing
in folders as described in README0.txt, set the following path names in 
the NAMELIST input file:
  atmpath = 'C:\EarthGRAM10\PC_IOfiles\atmosdat_E10.txt'
  NCEPpath = 'C:\EarthGRAM10\NCEPdata\'
    if the NCEP data are to be read from the C: drive, or
  NCEPpath = 'E:\NCEPdata\'
    if the NCEP data are to be read from the DVD on the E: drive.
    
For test or demonstration purposes, it is also possible to run Earth-GRAM
2010 directly from the distribution DVD.  To do this, from the C:\MyTest 
directory, as above, simply replace "C:\EarthGRAM10\" with "E:\" (assuming 
that E: is the DVD drive) in the command and input pathnames for the 
executable file, the atmosdat_E10.txt file, and the NCEP path name.
See later discussion on how to set file path names if using the
input trajectory option, the Range Reference Atmosphere input option,
the auxiliary profile input option, and/or the multiple-profile Monte
Carlo option.  Output file path names can be any name specified by the 
user.  All path names are limited to 99 characters in length.


Running in Command-Line Mode on a UNIX Machine

UNIX executable files must be compiled and created, as discussed above.  
For example, suppose the Earth-GRAM 2010 UNIX executable is created and 
named gram_E10.x in directory /usr/EarthGRAM10/UNIX_Executables.  If the 
setup procedure discussed in README0.txt is used, example input/output 
files are placed in folder /usr/EarthGRAM10/UNIX_IOfiles.  Earth-GRAM 2010 
can easily be run from any directory.  For example, to run from a directory 
named /usr/MyTest, do the following -

(1) Copy (and edit and/or rename, if desired) the file 
    /usr/EarthGRAM10/UNIX_IOfiles/NameRef.txt to the /usr/MyTest folder.
    Note that NAMELIST input file names are limited to 99  characters.
(2) Change directory by entering the command cd /usr/MyTest
(3) Execute Earth-GRAM 2010 by entering the command
     /usrEarthGRAM10/UNIX_Executables/gram_E10.x
(4) Enter the name of the input file (from step 1) when prompted
(5) Output files will appear in the /usr/MyTest folder

If file open errors are encountered, make sure that file pathnames are 
properly set in the NAMELIST input file.  For files residing in folders 
as described in README0.txt, set the following path names in the NAMELIST 
input file:
  atmpath = '/usr/EarthGRAM10/UNIX_IOfiles/atmosdat_E10.txt'
  NCEPpath = '/usr/EarthGRAM10/NCEPdata/'
    if the NCEP data are to be read from the /usr/EarthGRAM10 area, or
  NCEPpath = '/DVDROM/NCEPdata/'
    if NCEP data are to be read from the DVD on the /DVDROM/ drive

See later discussion on how to set file path names if using the input 
trajectory option, the Range Reference Atmosphere input option, the auxiliary 
profile input option, and/or the multiple-profile Monte Carlo option. Output 
file path names can be any name specified by the user.  All path names are 
limited to 99 characters in length.


Reference Input/Output Files for Testing

To facilitate runtime testing, a NAMELIST-formatted reference input file, 
called NameRef.txt, is provided on the distribution DVD (in folder 
EarthGRAM10\IOfiles on the distribution DVD).  Reference output files are also 
provided (OutputRef.txt and SpeciesRef.txt on the distribution DVD).  If a test 
run (as described in the above sections) is done using input file NameRef.txt, 
to produce output files named output.txt and species.txt, these output files 
can be tested against the reference output files provided.  On a PC (from the 
directory where the test output files reside), to compare against reference 
files on the DVD (E: drive), enter PC commands

             fc output.txt E:\PC_IOfiles\OutputRef.txt
             
             fc special.txt E:\PC_IOfiles\SpecialRef.txt

             fc species.txt E:\PC_IOfiles\SpeciesRef.txt
and             
             fc BLTest.txt E:\PC_IOfiles\BLTestRef.txt

Any file differences will be noted in the output from the fc (file compare)
command.  To do the same tests on a UNIX machine, the commands are

               diff output.txt /DVDROM/UNIX_IOfiles/OutputRef.txt

               diff special.txt /DVDROM/UNIX_IOfiles/SpecialRef.txt

               diff species.txt /DVDROM/UNIX_IOfiles/SpeciesRef.txt
and
               diff BLTest.txt /DVDROM/UNIX_IOfiles/BLTestRef.txt
    
Be sure to compare against files from the appropriate directory on the
distribution DVD, since there are differences in the PC and UNIX version
reference output files.  In UNIX format, numbers that are less than 1 (in 
magnitude) are written  (0.xxx -0.xxx etc.). For some compilers in the PC 
environment, such numbers  are written without the leading zeroes (i.e. as 
.xxx  -.xxx etc.). Because of variations in handling of round-off by 
different operating systems and/or compilers, a few numbers in the reference 
output files may differ from those in the user-generated test output files.  
Such differences should be no larger than about 1 in the last significant 
digit of the output.


Running with an Input Trajectory File

To run the stand-alone Earth-GRAM 2010 program with a pre-computed input 
trajectory file of positions and times, set the trajectory file pathname 
(trapath) in the NAMELIST input file to whatever file name is desired, set
the number of computed positions (nmax) to zero, and set the trajectory 
file unit number (iopt) to any non-conflicting value. The trajectory file 
contains one time and position per line, having values for elapsed time 
(seconds), height (km), geocentric latitude (degrees, North positive), and 
longitude (degrees, East positive).  Trajectory input heights greater than 
6000 km are treated as geocentric radius values, rather than as altitudes 
above reference ellipsoid.  Example programs illustrating how to drop 
Earth-GRAM into trajectory code cannot be run with trajectory file input, 
since trajectory positions are generated "on the fly" in these programs.  
For additional discussion of input options, see file README2.txt.


Running Multiple Profiles/Trajectories in Monte-Carlo Mode

Multiple profiles or trajectories can be computed in one run of the program,
in a Monte-Carlo simulation of various perturbation profiles, by providing
input for any number of random seed values.  Input random seeds are provided,
with one seed value per line, in a file whose path name is specified by
input parameter rndpath, the file unit number for which is specified by
input parameter iun. Integer random seed values must be in the range
0 <= seed <= 900,000,000.  Other aspects of the perturbation model can
be controlled by input parameters such as rpscale, initpert, patchy, rdinit,  
rtinit, ruinit, rvinit, rwinit, etc. For additional discussion of input 
options, see file README2.txt.


Running with Range Reference Atmosphere (RRA) Input

A major feature is the (optional) ability to use data (in the form of vertical 
profiles) from a set of Range Reference Atmospheres (RRA), as an alternate to 
the usual GRAM climatology, at a set of RRA site locations.  With this feature 
it is possible, for example, to simulate a flight profile that takes off from 
the location of one RRA site (e.g. Edwards AFB, using the Edwards RRA data), to 
smoothly transition into an atmosphere characterized by the GRAM climatology, 
then smoothly transition into an atmosphere characterized by a different RRA 
site (e.g. White Sands, NM), to be used as the landing site in the simulation.
Use of the RRA option is controlled by setting input parameters rrapath, 
rralist, iurra, sitelim, and sitenear.  For further discussion of input
options, see file README2.txt.  For further discussion of RRA data, see 
file README5.txt


Running with Auxiliary Profile Input

As an option, data read from an auxiliary profile may be used to replace
data from the conventional (NCEP/MAP/etc.) climatology. This option is 
controlled by setting parameters profile, sitenear, and sitelim in the 
NAMELIST input file. Each line of the auxiliary profile input file consists 
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
profile. Regular climatological values of wind components are used if BOTH wind 
components are zero in the auxiliary profile file. Note that RRA input option 
and auxiliary profile input option cannot be used simultaneously.  A sample 
auxiliary profile file (named KSCannRRA.txt) is provided.  See additional 
discussion in file README7.txt.


Compiling and Running on Other Platforms

MSFC Natural Environments Group does not have adequate resources to be able to
provide program versions and installation setups for all possible user
platforms.  Hopefully, the guidance provided above, and in README0.txt, will 
be adequate to allow users to set up, compile and run Earth-GRAM 2010 on 
whatever platform they desire to use.  For questions on compiling or running 
Earth-GRAM 2010, see contact information below.



If you have any questions or problems, please contact:

Patrick White
e-mail: Patrick.W.White@nasa.gov
phone:  (256)-961-1623

Mail Code EV44
NASA Marshall Space Flight Center
Huntsville, AL 35812


PDF versions of the background GRAM-95 and GRAM-99 reports [NASA Technical 
Memorandum 4715, "The NASA/MSFC Global Reference Atmospheric Model - 1995 
Version (GRAM-95)", August, 1995 and NASA/TM-1999-209630 "The NASA/MSFC Global 
Reference Atmospheric Model - 1999 Version (GRAM-99)", May 1999] are on the 
Earth-GRAM 2010 distribution DVD. 