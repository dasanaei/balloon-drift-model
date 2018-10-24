                             EARTH-GRAM 2010 OVERVIEW


General Introduction

    This is an overview of the program and data files provided with the 2010 
version of the Earth Global Reference Atmospheric Model (Earth-GRAM 2010).
An Earth-GRAM 2007 report and users guide (NASA/TM-2008-215581) is also provided on 
the distribution DVD as file GRAM07.pdf. 

    The user is assumed to be generally familiar with GRAM-95, as described in 
the report NASA Technical Memorandum 4715 "The NASA/MSFC Global Reference 
Atmospheric Model - 1995 Version (GRAM-95)", August 1995, also referred to here 
as the "GRAM-95 Report", or as NASA TM-4715. This report is provided as file 
GRAM95.pdf on the Earth-GRAM 2010 distribution DVD.

    The user is also assumed to be generally familiar with GRAM-99, as described 
in the report NASA/TM-1999-209630 "The NASA/MSFC Global Reference Atmospheric
Model - 1999 Version (GRAM-99)", May 1999, also referred to here as the
"GRAM-99 Report", or as NASA/TM-1999-209630. This report is provided as file
GRAM99.pdf on the Earth-GRAM 2010 distribution DVD.

    The Earth-GRAM 2010 distribution DVD contains:
      (1) A general introduction file, README.txt.
      (2) Several folders containing individual text and other files, suitable
          for browsing, or for "manual" set-up.  Contents of these folders
          are as follows -

Contents of AuxProfiles Folder

    Annual average Range Reference Atmosphere data from all 21 available sites,
in the form of Earth-GRAM "auxiliary profile" input files (see file README7.txt).
File names are RRAsssAnn.txt, where sss is the RRA site code (see file 
README5.txt).

Contents of DOCUMENTATION Folder

    The various README and other information files, and their general 
subject matter are:

ASR2004.pdf	         Portable Data Format (PDF) of Justus et al. (2004) paper on constituents
ASR2006.pdf	         Portable Data Format (PDF) of Justus et al. (2006) on MET model updates
GRAM90Part1.pdf          Part 1 of PDF file for GRAM-90 report
GRAM90Part2.pdf          Part 2 of PDF file for GRAM-90 report
GRAM95.pdf               PDF file for GRAM-95 report (referenced above)
GRAM99.pdf               PDF file for GRAM-99 report (referenced above)
GRAM07.pdf               GRAM 2007 Users Guide Tech Memo (NASA-TM-2008-215581)
GRAM2010.pdf             GRAM 2010 Users Guide Tech Memo (NASA/TM-2011-216467)
GRAM_VV&A.pdf            Earth-GRAM 2007 Verification, Validation and Accreditation document
gramfixV1.0.html         Description of fixes made between original Beta2 code and Version 1.0
gramfixV2.0.html         Description of fixes made between V1.0 and V2.0
gramhist.txt             Brief history of various GRAM versions
HWMfix.txt               Description of minor changes made to the Harmonic Wind Model (93)
MSISfix.txt              Description of minor changes made to MSIS (00) thermosphere model
NASA-CR-2062.pdf         PDF file of early upper-atmosphere study
NASA-CR-2203.pdf         PDF file of another early upper-atmosphere study
NASA-CR-181384.pdf       PDF file of report describing earlier GRAM perturbation updates
NASA-TM-4168.pdf         PDF file of NASA Atmospheric Turbulence model 
NASA-TM-100697.pdf       Fleming et al. global climatology report (MAP data)
NASA-TM-2008-215581.pdf  GRAM 2007 Users Guide Tech Memo
NCEPmods.pdf             Description of method used to adjust NCEP surface and near-surface
                         data, to achieve better correspondence with observations.
README0.txt              Description of files provided, and program setup procedure
README1.txt	         Instructions for compiling and running the program
README2.txt	         Description of Earth-GRAM 2010 input files
README3.txt	         Parameters available for the "Special Output" file
README4.txt	         National Centers for Environmental Prediction (NCEP) Climatology
README5.txt	         Range Reference Atmosphere (RRA) option and data files
README6.txt              MET-07, MSIS, and JB2008 thermosphere models
README7.txt              Other special GRAM features
README8.txt              Utility Programs for use with Earth-GRAM 2010
TurbTable.pdf            Turbulence statistics, based on NASA TM-4168, excerpted from NASA TM-4511
VerticalWindModel.pdf    Description and validation of revised vertical wind model,
                         added for Earth-GRAM 2007 Version 1.2, and revised for
                         Earth-GRAM 2010.


Contents of GRAM_code Folder

Earth-GRAM 2010 source code files in this folder are:

Cfiles_E10_C.f90  Module file for all common block variables
gram_E10.f90      Main Earth-GRAM 2010 program
gramsubs_E10.f90  General Earth-GRAM 2010 subroutines
HWMsubs_E10.f90   Harmonic Wind Model (used with MSIS thermosphere model)
Ifile_E10_I.f90   Interface module file for all subroutines
initial_E10.f90   Reads atmosdat file and initializes data
JB2008_E10.f90    Jacchia-Bowman 2008 model
MET07prg_E10.f90  The Marshall Engineering Thermosphere (MET-07) model
models_E10.f90    Other Earth-GRAM 2010 subroutines, not in gramsubs
MSISsubs_E10.f90  Naval Research Labs MSIS (00) thermosphere model
NCEPsubs_E10.f90  Reads and prepares NCEP climatology data
random_E10.f90    Random number generators
rramods_E10.f90   Reads the Range Reference Atmosphere (RRA) data
speconc_E10.f90   Species concentration subroutines
trajcalc_E10.f90  Subroutines used in trajectory demonstration programs.
gramtraj_E10.f90  Subroutine for use in user-provided trajectory program
trajdemo_E10.f90  Main driver (replaces gram_E10.f90) in simple trajectory
                  demonstration program
trajopts_E10.f90  Main driver (replaces gram_E10.f90) in example trajectory
                  program that demonstrates several special options
multtraj_E10.f90  Main driver (replaces gram_E10.f90) program illustrating 
                  how to call Earth-GRAM 2010 to evaluate several different
                  (multiple) trajectories in one program run. 
corrtraj_E10.f90  Main driver (replaces gram_E10.f90) program illustrating 
                  how to call Earth-GRAM 2010 to evaluate multiple profiles
                  in one program run, with small-scale and large-scale 
                  correlations preserved between the profiles.
multbody_E10.f90  Main driver illustrating how to call GRAM to evaluate
                  perturbations along trajectories of multiple bodies which
                  start out together and at various times separate from 
                  each other. The perturbations for each body start out the 
                  same, and then remain correlated over both position and 
                  time along the body trajectories, and across the spatial 
                  separation between the bodies.


Contents of NCEPdata Folder

  In subfolder FixedASCII:

  Climatology for NCEP monthly means and standard deviations, in ASCII
  format.  
  NCEPbinF.f90   Program to read the ASCII format NCEP files (once) and write 
                 them in binary, for use by Earth GRAM-2010. These are the NCEP
                 data files after having been adjusted by procedures described
                 in NCEPmods.pdf (in the DOCUMENTATION folder).
  NCEPbinF.exe   PC executable for the NCEPbinF.f90 utility program                      
                 See files README4.txt and README8.txt.

  In subfolder FixedBin:

  Climatology for adjusted NCEP monthly means and standard deviations, in PC 
  binary format.  See file README4.txt 


Contents of PC_Executables Folder

Executable files for stand-alone Earth-GRAM 2010, four versions of
programs to demonstrate GRAM used as a subroutine in a trajectory program
driver, and three utility programs provided.  UNIX executable files are
not provided, since not all UNIX machines have compatible executables.

gram_E10.exe        PC executable for normal (stand-alone) Earth-GRAM 2010
trajdemo_E10.exe    PC executable using trajdemo_E10.f90 as main driver for
                      example of using Earth-GRAM 2010 in trajectory 
                      calculations
trajopts_E10.exe    PC executable using trajopts_E10.f90 as main driver for
                      example of using Earth-GRAM 2010 in trajectory 
                      calculations
corrtraj_E10.exe    PC executable using corrtraj_E10.f90 as main driver for
                      example of using Earth-GRAM 2010 in trajectory 
                      calculations
multtraj_E10.exe    PC executable using multtraj_E10.f90 as main driver for
                      example of using Earth-GRAM 2010 in trajectory 
                      calculations
multbody_E10.exe    PC executable using multbody_E10.f90 as main driver for
                      example of using Earth-GRAM 2010 in trajectory 
                      calculations


Contents of IOfiles Folder

PC version of various input and output files from Earth-GRAM 2010. The only
difference in format between PC and UNIX I/O files are that OutputRef.txt 
and SpeciesRef.txt use the PC format, having no leading zero digit on
numbers that are less than one in magnitude.  Due to differences in
machine rounding between PC and UNIX platforms, some output numbers may 
differ by magnitude one in the last significant digit, between PC-version and
UNIX-version output files.

atmosdat_E10.txt   Atmospheric data file, always required by Earth-GRAM 2010
BLTestRef.txt      Boundary Layer model output file from reference case
NameRef.txt        Reference (test) NAMELIST input file
OutputRef.txt      Output file generated from NameRef.txt (PC format)
SpecialRef.txt     "Special-formatted" output file for reference case
SpeciesRef.txt     Species  file generated from NameRef.txt (PC format)
rand1000.txt       Sample random seed file for a 1000-profile Monte Carlo run
RRAtraj.txt        "Trajectory" input file, containing height, latitude,
                   longitude information for all available Range Reference
                   Atmosphere (RRA) sites (for comparing RRA option output
                   with NCEP climatology output)


Contents of RRAdata Folder

Folder containing 1983 and 2006 Range Reference Atmosphere (RRA) data files, 
for use with Earth-GRAM 2010, using RRA input option iyrrra.  File rrasites.txt
is a list of available RRA sites, and is read in by Earth-GRAM 2010.  To
use other sites or combinations of sites, build and read a site-list
file with similar information.  After discovery of an error in standard 
deviations of pressure and density in RRA data from the Roosevelt Roads 
Puerto Rico site (code rrd), data from this RRA site are no longer included.
For further details, see file README5.txt.


Contents of Utilities_Code Folder

bldtraj.f90   Program to build pseudo-trajectory file for using in Earth-
              GRAM 2010, to compute output for maps or cross-sections.
              See file README8.txt.
bldtraj.exe   PC executable for the bldtraj.f90 utility program
PathNames.txt Sample file with input path names, for use in bldtraj
makerand.f90  Program to generate random seeds for Monte-Carlo runs
makerand.exe  PC executable for the makerand.f90 utility program
timetraj.f90  Program to generate a trajectory file consisting of any
              number of profiles separated by a user-selected time step.
              See file README8.txt.
timetraj.exe  PC executable for the timetraj.f90 utility program
              

                  HOW TO SET UP EARTH-GRAM 2010 PROGRAM AND DATA FILES
                  
Manual Setup on a PC

The following example of manual setup on a PC assumes that the user wants
to set up Earth-GRAM 2010 in a directory called EarthGRAM2010 on the C: drive.
The distribution DVD is assumed to be on the E: drive

(1) Copy the E:\AuxProfiles folder to C:\EarthGRAM2010\AuxProfiles
(2) Copy the E:\DOCUMENTATION folder to C:\EarthGRAM2010\DOCUMENTATION
(3) Copy the E:\GRAM_code folder to C:\EarthGRAM2010\GRAM_code
(4) Copy the E:\NCEPdata folder to C:\EarthGRAM2010\NCEPdata
(5) Copy the E:\PC_Executables folder to C:\EarthGRAM2010\PC_Executables
(6) Copy the E:\IOfiles folder to C:\EarthGRAM2010\PC_IOfiles
(7) Copy the E:\RRAdata folder to C:\EarthGRAM2010\RRAdata
(8) Copy the E:\Utilities_Code folder to C:\EarthGRAM2010\Utilities_Code

NOTE: In Windows Explorer, all of these operations can be done at one time,
by creating a folder C:\EarthGRAM2010, then highlighting all folders on the 
E: drive (AuxProfiles, DOCUMENTATION, GRAM_code, NCEPdata, PC_Executables, 
IOfiles, RRAdata, and Utilities_Code), then dragging and dropping these 
folders to folder C:\EarthGRAM2010.

Note that PC executable files are provided.  Unless the user wishes to
make code changes (e.g. to modify the "special-formatted" output), then
it is not necessary to compile the GRAM code. 

Note that NCEP binary data files are provided.  These should work as-is
with the PC executable programs provided.  If you change any of the
source code provided and re-compile, then it is possible that you may 
need to re-build the NCEP binary files from the NCEP ASCII files provided.
In that case:
    Copy the E:\EarthGRAM2010\NCEPdata\ASCIIdata folder to folder 
    C:\EarthGRAM2010\NCEPdata\ASCIIdata; Compile and run utility program 
    NCEPbinF.f90 to generate binary format NCEP data file.  See further 
    details in README4.txt and README8.txt.


Manual Setup on a UNIX Platform

If the user has access to a UNIX machine with a DVD drive, called /DVDROM, for 
example, then a setup procedure similar to that for PC installation, given 
above, can be used. Note that source code and other text files may still have 
PC end-of-line markers, rather than UNIX end-of-line markers, when copied 
directly from the DVD.  The following assumes that the user wants to set up 
Earth-GRAM 2010 in a directory /usr/EarthGRAM2010 on the UNIX machine

(1) Copy the /DVDROM/AuxProfiles folder to /usr/EarthGRAM2010/AuxProfiles
(2) Copy the /DVDROM/DOCUMENTATION folder to /usr/EarthGRAM2010/DOCUMENTATION
(3) Copy the /DVDROM/GRAM_code folder to /usr/EarthGRAM2010/GRAM_code
(4) Copy the /DVDROM/NCEPdata/ASCIIdata folder to folder
    /usr/EarthGRAM2010/NCEPdata/ASCIIdata; Compile and run utility program 
    NCEPbinF.f90 to generate binary format NCEP data file.  See further details
    in README4.txt and README8.txt.
(5) Copy the /DVDROM/RRAdata folder to /usr/EarthGRAM2010/RRAdata
(6) Copy the /DVDROM/IOfiles folder to /usr/EarthGRAM2010/UNIX_IOfiles
(7) Copy the /DVDROM/Utilities_Code folder to /usr/EarthGRAM2010/Utilities_Code

Note that UNIX executable files are not provided.  Executable files must be
compiled on UNIX platforms, according to discussion given in file README1.txt.
For setup on UNIX machines for which a DVD drive is not easily accessible,
the distribution DVD may be loaded on a PC and File Transfer Protocol (FTP)
software may be used to transfer files, rather than using the copy commands
as outlined above.  If files are transferred to a UNIX platform from
a PC platform by using ASCII-mode FTP, conversion of end-of-line markers
is handled automatically.


Setup on Other Platforms

MSFC Natural Environments Group does not have adequate resources to provide 
program versions and installation setups for all possible user platforms.  
Hopefully, the guidance provided above will be adequate to allow users to set 
up Earth-GRAM 2010 on whatever platform they desire to use.  For questions on 
installation or operation of Earth-GRAM 2010, see contact information below.


How to Run Earth-GRAM 2010

For instructions on how to compile (if necessary) the Earth-GRAM 2010 program, 
and how to run it, see file README1.txt



CONTACT INFORMATION

If you have any questions or problems, please contact:

Patrick White
e-mail: Patrick.W.White@nasa.gov
phone:  (256)-961-1623

Mail Code EV44
NASA Marshall Space Flight Center
Huntsville, AL 35812





