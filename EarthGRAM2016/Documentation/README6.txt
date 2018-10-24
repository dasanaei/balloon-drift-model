MET-07, MSIS AND JB2008 THERMOSPHERIC MODELS AND PERTURBATION MODEL FEATURES
 

The Marshall Engineering Thermosphere (MET-07) Model

    The 2007 version of Marshall Engineering Thermosphere (MET-07) forms the 
default upper altitude portion (> 90 km) of Earth-GRAM 2007.  

    As with the version of MET in GRAM-99, MET-07 uses the Equinox of epoch 
J2000 to compute solar positions.  Other features of MET-07 include:

(1) Correction of number density and molecular weight, according to discussion
    in Justus et al. "Earth GRAM-99 and Trace Constituents", COSPAR, 2004.
(2) Change from spherical-Earth approximation to latitude-dependent surface
    gravity and effective Earth radius.
(3) Change from time resolution only to the nearest integer minute to (real)
    seconds time resolution.
(4) Correction of small discontinuities in the semi-annual variation term by
    converting day-of-year to real instead of integer, and treating each
    year as having either 365 or 366 days (as appropriate), rather than all
    years being treated as of length 365.2422 days.
(5) Additional output from MET07_TME subroutine of modified Julian Day, right 
    ascension of Sun, and right ascension at local lat-lon (used for input 
    to new JB2008 thermosphere model)


Mass Spectrometer, Incoherent Scatter (MSIS) Radar Extended Model  

   The option is now provided to use the 2000 version Naval Research Laboratory 
(NRL) MSIS Extended model, NRLMSISE-00, for Earth-GRAM thermospheric conditions.
If this option is selected, thermospheric winds are evaluated using the NRL
1993 Harmonic Wind Model, HWM-93.  If the MET option is selected, winds are 
computed from a geostrophic wind model, with modifications for thermospheric
effects of molecular viscosity.  This wind model has been used in GRAM since 
the 1990 version (Justus et al., NASA TM-4268, 1991).
   Information on the MSIS and HWM models is available at the following URLs:
   
http://www.nrl.navy.mil/content.php?P=03REVIEW105   
   
http://modelweb.gsfc.nasa.gov/atmos/nrlmsise00.html

http://www.answers.com/topic/nrlmsise-00

http://fact-archive.com/encyclopedia/NRLMSISE-00

http://www.eiscat.rl.ac.uk/svn/guisdap.svn/trunk/dist/g85/models/nrlmsise00/readme.txt

http://modelweb.gsfc.nasa.gov/atmos/hwm.html

http://adsabs.harvard.edu/abs/2006AGUFMSA11A..07D

http://nssdcftp.gsfc.nasa.gov/models/atmospheric/hwm07/readme.txt



   Selection of either MSIS thermospheric model option is controlled by input 
parameter itherm (see file README2.txt).  Minor corrections in MSIS and HWM (along 
the lines of corrections 3 and 4 for MET, above) have been made. Therefore, MSIS/HWM 
output from Earth-GRAM 2010 will not agree totally with output from the original 
NRLMSISE-00 version.  See file MSISfix.txt for more complete discussion of these 
corrections.


Jacchia-Bowman 2008 (JB2008) Thermosphere Model


!     This is the CIRA "Integration Form" of a Jacchia Model.                     !JB08 10
!     There are no tabular values of density.  Instead, the barometric            !JB08 11
!     equation and diffusion equation are integrated numerically using            !JB08 12
!     the Newton-Coates method to produce the density profile up to the           !JB08 13
!     input position.                                                             !JB08 14
!                                                                                 !JB08 15
!     INPUT:                                                                      !JB08 16
!                                                                                 !JB08 17
!           AMJD   : Date and Time, in modified Julian Days                       !JB08 18
!                    and Fraction (MJD = JD-2400000.5)                            !JB08 19
!           SUN(1) : Right Ascension of Sun (radians)                             !JB08 20
!           SUN(2) : Declination of Sun (radians)                                 !JB08 21
!           SAT(1) : Right Ascension of Position (radians)                        !JB08 22
!           SAT(2) : Geocentric Latitude of Position (radians)                    !JB08 23
!           SAT(3) : Height of Position (km)                                      !JB08 24
!           F10    : 10.7-cm Solar Flux (1.0E-22*Watt/(M**2*Hertz))               !JB08 25
!                    (Tabular time 1.0 day earlier)                               !JB08 26
!           F10B   : 10.7-cm Solar Flux, ave.                                     !JB08 27
!                    81-day centered on the input time                            !JB08 28
!                    (Tabular time 1.0 day earlier)                               !JB08 29
!           S10    : EUV index (26-34 nm) scaled to F10                           !JB08 30
!                    (Tabular time 1.0 day earlier)                               !JB08 31
!           S10B   : EUV 81-day ave. centered index                               !JB08 32
!                    (Tabular time 1.0 day earlier)                               !JB08 33
!           XM10   : MG2 index scaled to F10                                      !JB08 34
!                    (Tabular time 2.0 days earlier)                              !JB08 35
!           XM10B  : MG2 81-day ave. centered index                               !JB08 36
!                    (Tabular time 2.0 days earlier)                              !JB08 37
!           Y10    : Solar X-Ray & Lya index scaled to F10                        !JB08 38
!                    (Tabular time 5.0 days earlier)                              !JB08 39
!           Y10B   : Solar X-Ray & Lya 81-day ave. centered index                 !JB08 40
!                    (Tabular time 5.0 days earlier)                              !JB08 41
!           DSTDTC : Temperature change computed from Dst index                   !JB08 42
!                                                                                 !JB08 43
!     OUTPUT:                                                                     !JB08 44
!                                                                                 !JB08 45
!           TEMP(1): Exospheric Temperature above Input Position (deg K)          !JB08 46
!           TEMP(2): Temperature at Input Position (deg K)                        !JB08 47
!           RHO    : Total Mass-Density at Input Position (kg/m**3)               !JB08 48
!           pres   : Pressure (N/m**2)                                            !JB08 49
!           avgMW  : Mean molecular weight (kg/kmole)                             !JB08 50
!           aNd(1) : N2 number density (#/m**3)                                   !JB08 51
!           aNd(2) : O2 number density (#/m**3)                                   !JB08 52
!           aNd(3) : O  number density (#/m**3)                                   !JB08 53
!           aNd(4) : Ar number density (#/m**3)                                   !JB08 54
!           aNd(5) : He number density (#/m**3)                                   !JB08 55
!           aNd(6) : H  number density (#/m**3)                                   !JB08 56
!           SUMN   : total number density (#/m**3)                                !JB08 57
!                                                                                 !JB08 58
!                                                                                 !JB08 59
!     JB2008 Model Development: (Ref. 7)                                          !JB08 60
!                                                                                 !JB08 61
!                                                                                 !JB08 62
!     A. Development of the JB2006 model:                                         !JB08 63
!                                                                                 !JB08 64
!       1. Started with the CIRA72 model (Jacchia 71).                            !JB08 65
!                                                                                 !JB08 66
!       2. Converted to CIRA70 model replacing equations from Jacchia 70          !JB08 67
!          model (Ref. 5)                                                         !JB08 68
!                                                                                 !JB08 69
!       3. Replaced Tc equation using new solar indices (Ref. 1 and 2)            !JB08 70
!                                                                                 !JB08 71
!       4. Replaced semiannual equation with new global model based               !JB08 72
!          on F10B (Ref. 1 and 3)                                                 !JB08 73
!                                                                                 !JB08 74
!       5. Added correction for local solar time and latitude errors              !JB08 75
!          (Ref. 1)                                                               !JB08 76
!          Added smooth transition between altitude bands                         !JB08 77
!                                                                                 !JB08 78
!       6. Added high altitude ( z > 1500 km ) correction                         !JB08 79
!          (Ref. 1 and 4)                                                         !JB08 80
!                                                                                 !JB08 81
!       7. REV A of JB2006 - Oct 2006                                             !JB08 82
!                Smoothing of density corrections and scale height                !JB08 83
!                through different altitude bands in the latitude-                !JB08 84
!                local time correction subroutine DTSUB                           !JB08 85
!                dTx correction replaced with dTc correction                      !JB08 86
!                                                                                 !JB08 87
!     B. Modification to develop JB2008 model:                                    !JB08 88
!                                                                                 !JB08 89
!       1. Replaced Tc equation in JB2006 using new solar indices                 !JB08 90
!          (Ref. 7)                                                               !JB08 91
!                                                                                 !JB08 92
!       2. Replaced semiannual equation with new global model based               !JB08 93
!          on F10B and S10B (Ref. 6)                                              !JB08 94
!                                                                                 !JB08 95
!       3. Use dTc value based on Dst geomagnetic storm index                     !JB08 96
!          (This replaces ap use) (Ref. 7)                                        !JB08 97
!                                                                                 !JB08 98
!                                                                                 !JB08 99
!         All equation references below refer to the original                     !JB08100
!         Jacchia 1971 (CIRA 1972) model papers.                                  !JB08101
!                                                                                 !JB08102
!                                                                                 !JB08103
!     References:                                                                 !JB08104
!                                                                                 !JB08105
!      1. Bowman, Bruce R., etc. : "A New Empirical Thermospheric                 !JB08106
!         Density Model JB2006 Using New Solar Indices",                          !JB08107
!         AIAA/AAS Astrodynamics Specialists Conference, Keystone, CO,            !JB08108
!         21-24 Aug 2006, (Paper AIAA 2006-6166).                                 !JB08109
!                                                                                 !JB08110
!      2. Bowman, Bruce R., etc. : "Improvements in Modeling                      !JB08111
!         Thermospheric Densities Using New EUV and FUV Solar Indices",           !JB08112
!         AAS/AIAA Space Flight Mechanics Meeting, Tampa, FL,                     !JB08113
!         23-26 Jan 2006, (Paper AAS 06-237).                                     !JB08114
!                                                                                 !JB08115
!      3. Bowman, Bruce R.: "The Semiannual Thermospheric Density                 !JB08116
!         Variation From 1970 to 2002 Between 200-1100 km",                       !JB08117
!         AAS/AIAA Space Flight Mechanics Meeting, Maui, HI,                      !JB08118
!         8-12 Feb 2004, (Paper AAS 04-174).                                      !JB08119
!                                                                                 !JB08120
!      4. Bowman, Bruce R.; "Atmospheric Density Variations at                    !JB08121
!         1500 km to 4000 km Height Determined from Long Term                     !JB08122
!         Orbit Perturbation Analysis", AAS/AIAA Space Flight                     !JB08123
!         Mechanics Meeting, Santa Barbara, CA, 11-14 Feb 2001,                   !JB08124
!         (Paper AAS 01-132).                                                     !JB08125
!                                                                                 !JB08126
!      5. Jacchia, Luigi G.; "New Static Models of the                            !JB08127
!         Thermosphere and Exosphere with Empirical Temperature                   !JB08128
!         Profiles", (Smithsonian Astrophysical Observatory                       !JB08129
!         Special Report 313), 6 May 1970.                                        !JB08130
!                                                                                 !JB08131
!      6. Bowman, Bruce R., etc. : "The Thermospheric Semiannual Density          !JB08132
!         Response to Solar EUV Heating," JASTP, 2008                             !JB08133
!                                                                                 !JB08134
!      7. Bowman, Bruce R., etc. : "A New Empirical Thermospheric                 !JB08135
!         Density Model JB2008 Using New Solar and Geomagnetic Indices",          !JB08136
!         AIAA/AAS 2008, COSPAR CIRA 2008 Model                                   !JB08137


Other information and references to developmental papers for JB2008 are given at 
the JB2008 web site 

http://sol.spacenvironment.net/~JB2008/

http://sol.spacenvironment.net/~JB2008/code.html

http://adsabs.harvard.edu/abs/2008cosp...37..367B
      
These sites have links to solar indices required by JB2008 (s10 and xm10), JB2008 
source code, publications, contacts, figures, and the Space Environment Technologies 
Space Weather site.
   If JB2008 is selected for calculation of thermospheric density and temperature, 
winds are computed with the Harmonic Wind Model (HWM 93), used in conjunction with 
the MSIS model.


Thermospheric Perturbation Model

   A two-scale perturbation model is used in Earth-GRAM.  Small-scale 
perturbations are computed by a one-step Markov process (a 1st order auto-
regressive approach, equivalent to the 1st order autoregressive model of 
Hickey in NASA Contractor Report 4605).  A wave model for treating large-
scale perturbations was introduced in GRAM-99.  This model uses a cosine 
wave, with both horizontal and vertical wavelengths. For use in Monte-Carlo 
simulations, a degree of randomness is introduced into the large-scale wave 
model by randomly selecting the phase of the cosine wave (under control of 
the same random number seed values as used for the small-scale perturbations).
   Parameters for the small-scale perturbation model in GRAM-99 (and Earth-GRAM 
2010) were re-calculated, to produce good agreement with data given in Table 2 
of Hickey (NASA CR-4605 and NASA CR-201140).  The one-step correlation over 
15 seconds of movement for the Atmospheric Explorer (AE) satellites is 
equivalent to an average value of 0.846 (with a standard deviation of 0.040).  
The length scales used in the perturbation model of GRAM-99 and Earth-GRAM 2007
yield a correlation value of 0.870 over the distance which the AE satellite 
moves in 15 seconds, well within the range of variability of the correlation 
data from Hickey.  See further discussion of the thermospheric perturbation 
model in Section 3.3 of the GRAM-99 report.
   The large-scale wave model for perturbations in Earth-GRAM 2010 is used 
at all altitudes, not just in the MET or MSIS model altitude range.  
