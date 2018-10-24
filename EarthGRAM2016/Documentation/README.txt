             Earth Global Reference Atmospheric Model (Earth-GRAM 2010 Ver4.0)

                        NASA Marshall Space Flight Center
                        Natural Environments Branch, EV44
                      Marshall Space Flight Center, AL 35812


   This DVD contains:
   
File README.txt         - This general introduction file.

Folder DOCUMENTATION    - Folder containing various README and other information 
                          files, including PDF GRAM-95 and GRAM-99 and Earth-GRAM
                          2007 reports.

Folder GRAM_code        - Folder containing Earth-GRAM 2010 source code files.  

Folder Utilities_Code   - Folder containing source code for four utility
                          programs (bldtraj.f90, makerand.f90, NCEPbinF.f90,
                          and timetraj.f90)

Folder PC_Executables   - Folder containing PC executable files for stand-alone 
                          Earth-GRAM 2010, and various versions of programs to 
                          demonstrate GRAM used as a subroutine in a trajectory 
                          program driver.  UNIX executable files are not provided, 
                          since not all UNIX machines have compatible exectuables.

Folder IOfiles          - Folder containing various input and output files from 
                          Earth-GRAM 2010, including "reference" output files: 
			  OutputRef.txt, SpecialRef.txt, SpeciesRef.txt, and 
			  BLTestRef.txt. Due to differences in machine rounding 
			  between PC and UNIX or Linux platforms, some output  
			  numbers may differ by magnitude one in the last 
			  significant digit, between PC-version and UNIX-
			  version output files.

Folder NCEPdata         - Folder containing all ASCII and PC-binary data for
                          National Centers for Environmental Prediction (NCEP)
                          climatology.  Data are provided for all months, at
                          four UT times per day, plus daily statistics.  NCEP
                          period-of-record provided is 1990 through 2008. More
                          information on the NCEP climatology, including how
                          to make Unix binary files from the ASCII data provided,
                          is in file README4.txt.
                        
Folder RRAdata          - Folder containing 1983 and 2006 Range Reference Atmosphere
                          (RRA) data files, for use with the Earth-GRAM 2010 RRA 
                          input option.

For a more complete description of all files provided, and how to set up the 
data and program files, see file README0.txt, in the DOCUMENTATION folder.  For 
instructions on compiling and running Earth-GRAM 2010, see file README1.txt, in 
the DOCUMENTATION folder.

NOTE: Text files provided on the distribution DVD have end-of-line marks in PC
format.  These may need to be converted to end-of-line marks on your non-PC
platform.  FTP file transfer to your target platform from the DVD, read on a PC, 
can make the transform, if the FTP file transfer is done in ASCII mode.  On-line 
or GNU utilities are also available that can perform this end-of-line transform 
for you.






