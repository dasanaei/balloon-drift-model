                         SPECIAL OUTPUT FILE OPTIONS

Following is a listing of the section of GRAM where the "Special Output" can be prepared 
and generated.  This identifies logical places to do computations of special variables 
(e.g. the examples wind speed, wndspd, at line ATMD753, pressure scale height, Hgtp, and 
density scale height, Hgtd, at lines ATMD758 and ATMD759) or to do units conversions 
(e.g. multiplying by conversion factors, such as after lines ATMD886).  Comments 
in this code section also provide lists of the names of the variables that are available 
for writing to the Special Output. Examples of write statements and user-provided 
format for writing to the Special Output are given on lines ATMD905 through ATMD928.


!......................................................................           !ATMD747
!-----"Special" output option section                                             !ATMD748
      if (iopp /= 0) then                                                         !ATMD749
!       EXAMPLES OF SPECIALLY COMPUTED OUTPUT:                                    !ATMD750
!-----  Wind speed (m/s) and direction (meteorological convention) from           !ATMD751
!       mean winds                                                                !ATMD752
         wndspd = dsqrt(ugh**2 + vgh**2)                                          !ATMD753
         wnddir = 180.0D0                                                         !ATMD754
         if (dabs(ugh)>0.0D0 .and. dabs(vgh)>0.0D0) wnddir = 180.0D0*(1.0D0 +  &  !ATMD755
            datan2(ugh,vgh)/pi)                                                   !ATMD756
!-----  Pressure scale height (m), density scale height (m)                       !ATMD757
         hgtp = pgh/(dgh*g)                                                       !ATMD758
         if (ch <= hj1) then                                                      !ATMD759
            hgtd = hgtp/(1.0D0 + hgtp*dtz/tgh)                                    !ATMD760
         else                                                                     !ATMD761
            hgtd = hgtp/(1.0D0 + hgtp*dtz/tgh - hgtp*dmdz/wtmol)                  !ATMD762
         endif                                                                    !ATMD763
!-----  Mean free path (m)  (assume mean values of pressure and temp.)            !ATMD764
         mfpath = 7.071D17/(0.424D0*totnd)                                        !ATMD765
!-----  Other variables that are functions of pressure, density,                  !ATMD766
!       temperature, wind and/or moisture, can be calculated here.                !ATMD767
!       Some such variables are: coefficient of viscosity, kinematic              !ATMD768
!       viscosity, thermal conduction coefficient, potential                      !ATMD769
!       temperature, equivalent potential temperature, etc.                       !ATMD770
!                                                                                 !ATMD771
!----   **** To print out header information, refer to the section                !ATMD772
!       **** near format label 200 in init subroutine of initial_E10.f.           !ATMD773
!......................................................................           !ATMD774
!                                                                                 !ATMD775
!     As an aid to the user, the following tables give the names of               !ATMD776
!     the variables that are available for output:                                !ATMD777
!                                                                                 !ATMD778
!---- Position and time parameters                                                !ATMD779
!     ----------------------------                                                !ATMD780
!      ch -     Geocentric Height (km) above WGS84 Reference Ellipsoid            !ATMD781
!      Rlocal - Radius from Earth center (km) (Earth radius plus ch)              !ATMD782
!      phi -   Geocentric Latitude (deg)                                          !ATMD783
!      GdLat - Geodetic Latitude (deg)                                            !ATMD784
!      thet -  Longitude (deg), East(+) West(-)                                   !ATMD785
!      elt -   Elapsed Time (sec)                                                 !ATMD786
!                                                                                 !ATMD787
!---- Thermodynamic, wind and moisture parameters (on standard output)            !ATMD788
!     ----------------------------------------------------------------            !ATMD789
!                                                                                 !ATMD790
!                                                E-W   N-S   Vert.                !ATMD791
!           Pressure    Density/  Temperature    Wind/ Wind/ Wind                 !ATMD792
!           /Vap.Pr.    Vap.Dens.   /Dewpt.      SoSm  SoSp  (m/s)                !ATMD793
!           (Nt/m**2)   (kg/m**3)     (K)        (m/s) (m/s) RH(%)                !ATMD794
!           ---------   ---------  -----------   ----- ----- -----                !ATMD795
!*Mean         pgh         dgh        tgh         ugh   vgh   wgh                 !ATMD796
!                                                                                 !ATMD797
!*Mean-76      pghp (%)    dghp (%)   tghp (%)    n/a   n/a   n/a                 !ATMD798
!                                                                                 !ATMD799
!*Small-Scale                                                                     !ATMD800
! Perturbation prhs (%)    drhs (%)   trhs (%)    urhs  vrhs  n/a                 !ATMD801
!                                                                                 !ATMD802
!*Small-Scale                                                                     !ATMD803
! Stand. Dev.  sphs (%)    sdhs (%)   sths (%)    suhs  svhs  n/a                 !ATMD804
!                                                                                 !ATMD805
!*Large-Scale                                                                     !ATMD806
! Perturbation prhl (%)    drhl (%)   trhl (%)    urhl  vrhl  n/a                 !ATMD807
!                                                                                 !ATMD808
!*Large-Scale                                                                     !ATMD809
! Stand. Dev.  sphl (%)    sdhl (%)   sthl (%)    suhl  svhl  n/a                 !ATMD810
!                                                                                 !ATMD811
!*Total                                                                           !ATMD812
! Perturbation prh  (%)    drh  (%)   trh  (%)    urh   vrh   wrh                 !ATMD813
!                                                                                 !ATMD814
!*Total Stand.                                                                    !ATMD815
! Deviation    sph  (%)    sdh  (%)   sth  (%)    suh   svh   swh                 !ATMD816
!                                                                                 !ATMD817
!*Total=Mean+                                                                     !ATMD818
! Perturbation ph          dh         th          uh    vh    wh                  !ATMD819
!                                                                                 !ATMD820
!*Total-US76   php  (%)    dhp  (%)   thp  (%)   csp0   csp   n/a                 !ATMD821
!                                                                                 !ATMD822
!*Mean H2O     eofT        rhov       tdgh        n/a   n/a   rhp                 !ATMD823
!                                                                                 !ATMD824
!*Stand. Dev.                                                                     !ATMD825
! H2O          seofT       srhov      stdgh       n/a   n/a  srhp                 !ATMD826
!                                                                                 !ATMD827
!----   Species concentration parameters (on species output)                      !ATMD828
!       ----------------------------------------------------                      !ATMD829
!                                                                                 !ATMD830
!           H2O  O3   N2O  CO   CH4  CO2   N2  O2   O  Ar  He  H   N              !ATMD831
!           ---  --   ---  --   ---  ---   --  --   -  --  --  -   -              !ATMD832
! Concen-   ppmh2o    ppmn2o    ppmch4     ppmn2    ppmo   ppmhe  ppmn            !ATMD833
! tration        ppmo3     ppmco     ppmco2    ppmo2   ppmar   ppmh               !ATMD834
! -----------------------------------------------------------------               !ATMD835
! Number    h2ond     n2ond     ch4nd      n2nd     ond    hend   nnd             !ATMD836
! Density        o3nd      cond      co2nd     o2nd    arnd    hnd                !ATMD837
!                                                                                 !ATMD838
!       Mean molecular weight=mwnd        total number densty=totnd               !ATMD839
!                                                                                 !ATMD840
!                                                                                 !ATMD841
!----   Surface data (passed from Subroutine NCEPmod, via Common                  !ATMD842
!       Block srfdat)                                                             !ATMD843
!       -------------------------------------------------------                   !ATMD844
!       psrf     = average surface pressure (N/m**2)                              !ATMD845
!       dsrf     = average surface density (kg/m**3)                              !ATMD846
!       tsrf     = average surface temperature (K)                                !ATMD847
!       usrf     = average surface Eastward wind component (m/s)                  !ATMD848
!       vsrf     = average surface Northward wind component (m/s)                 !ATMD849
!       hsrf     = height of surface (km, above sea level)                        !ATMD850
!       tdsrf    = average surface dewpoint temperature (K)                       !ATMD851
!       spsrf    = standard deviation of surface pressure (N/m**2)                !ATMD852
!       sdsrf    = standard deviation of surface density (kg/m**3)                !ATMD853
!       stsrf    = standard deviation of surface temperature (K)                  !ATMD854
!       susrf    = standard deviation of surface Eastward wind (m/s)              !ATMD855
!       svsrf    = standard deviation of surface Northward wind (m/s)             !ATMD856
!       shsrf    = std. dev. (uncertainty) of surf. hgt (m, from spsrf)           !ATMD857
!       stdsrf   = standard deviation of surface dewpoint temp. (K)               !ATMD858
!       spdavsrf = average wind speed at surface (m/s)                            !ATMD859
!       spdsdsrf = standard deviation of surface wind speed (m/s)                 !ATMD860
!       uvtsrf   = U-V cross correlation at surface (unitless)                    !ATMD861
!                                                                                 !ATMD862
!...    Surface roughness information (passed from Subroutine pertrb,             !ATMD863
!       via Common Block vert) for sigma-w model                                  !ATMD864
!       -------------------------------------------------------------             !ATMD865
!       lc     = land surface code (0-13 if from 1-by-1 deg data base             !ATMD866
!                or 99 if z0 is from user input)                                  !ATMD867
!       z0     = surface roughness (m), determined from 1-by-1 deg                !ATMD868
!                land surface type, or as specified by user via z0in)             !ATMD869
!                                                                                 !ATMD870
!-----At this point, the user is invited to insert whatever output                !ATMD871
!     parameters (in whatever format) are desired for the "special"               !ATMD872
!     output instead of the Write and Format statements below. Any new            !ATMD873
!     variables introduced must be declared at the beginning of the               !ATMD874
!     atmod subroutine.                                                           !ATMD875
!                                                                                 !ATMD876
         isev = 0                                                                 !ATMD877
         if (densfact > 1.0D0) isev = 1                                           !ATMD878
!---  Release code version of write to special output file                        !ATMD879
         write (iopp, 9000) elt, ch, phi, thet, dgh, pgh, tgh, ugh, vgh, dh,  &   !ATMD880
            ph, th, uh, vh, drh, prh, sdh, sph, 0.01D0*sth*tgh, suh, svh,     &   !ATMD881
            swh, wh, csp0, csp, isev, lc, z0                                      !ATMD882
 9000    format(f10.2,f9.3,f10.5,f11.5,1p,2e11.4,0p,f8.2,2f8.2,1p,2e11.4,0p,  &   !ATMD883
            f8.2,2f8.2,9f7.2,2f8.2,i3,i4,f8.5)                                    !ATMD884
!                                                                                 !ATMD885
!-----  Any change of units for writing to "special" output file can              !ATMD886
!       be done here [e.g. International Standard (SI) to English                 !ATMD887
!       units, etc.]                                                              !ATMD888
!                                                                                 !ATMD889
!---  Samples of how to change units:                                             !ATMD890
!                                                                                 !ATMD891
!     Change pressure and standard deviation to millibars                         !ATMD892
!     pgh = pgh/100.0d0                                                           !ATMD893
!     sph = pgh*sph/100.0d0                                                       !ATMD894
!     ph =  ph/100.0d0                                                            !ATMD895
!---  Change density and standard deviation to gm/m**3                            !ATMD896
!     dgh = dgh*1000.0d0                                                          !ATMD897
!     sdh = dgh*sdh/100.                                                          !ATMD898
!     dh = dh*1000.0d0                                                            !ATMD899
!---  Change temperature and standard deviation to deg C                          !ATMD900
!     sth = tgh*sth/100.0d0                                                       !ATMD901
!     tgh = tgh - 273.15d0                                                        !ATMD902
!     th = th - 273.15d0                                                          !ATMD903
!                                                                                 !ATMD904
!---  Samples of some other write statements and formats for special              !ATMD905
!     output file, including some changes of units (e.g. dividing                 !ATMD906
!     height by 0.3048 to convert to k-feet), dividing mission                    !ATMD907
!     elapsed time by 3600. to convert to hours)                                  !ATMD908
!                                                                                 !ATMD909
!     Write(iopp,9000)ch,mfpath,dgh,tgh                                           !ATMD910
!9000 Format(F6.1,1p,2E11.3,0p,F7.1)                                              !ATMD911
!     Write(iopp,9000)ch,suh,svh,sdh,urh,vrh,drh                                  !ATMD912
!9000 Format(F6.2,6F8.2)                                                          !ATMD913
!---  Write(iopp,9000)elt/3600.,ch/0.3048,prhs,prhl,prh,drhs,drhl,drh,    &       !ATMD914
!---    trhs,trhl,trh,prhs/sphs,drhs/sdhs,trhs/sths                               !ATMD915
!9000 Format(2F8.1,12F8.3)                                                        !ATMD916
!---  Write(iopp,9000)elt/3600.,ch/0.3048,pgh,sph,ph,dgh,sdh,dh,tgh,      &       !ATMD917
!---    sth,th,prhs/sphs,drhs/sdhs,trhs/sths                                      !ATMD918
!9000 Format(2F8.1,1p,6E10.3,0p,3F10.2,3F8.3)                                     !ATMD919
!---  Write(iopp,9000)ch,sphs/sph,sdhs/sdh,sths/sth,suhs/suh,svhs/svh     &       !ATMD920
!---   sphl/sph,sdhl/sdh,sthl/sth,suhl/suh,svhl/svh,sph,sdh,sth,suh,      &       !ATMD921
!---   svh,suhs                                                                   !ATMD922
!9000 Format(F7.2,10F6.4,3F7.3,F6.2,2F7.2)                                        !ATMD923
!---  Write(iopp,9000)elt/3600.,ch,uh,vh                                          !ATMD924
!9000   Format(F8.1,F8.3,2F8.2)                                                   !ATMD925
!---  Write(iopp,9000)elt,ch,phi,thet,dgh,dghp,sdh,drh,tgh,sth,trh,ugh,   &       !ATMD926
!       suh,urh,vgh,svh,vrh                                                       !ATMD927
!9000 Format(F9.1,F9.3,F8.3,F9.3,1p,E10.3,0p,3F6.1,F7.1,2F6.1,6F7.1)              !ATMD928
!                                                                                 !ATMD929
!-----  The "special" output option Write and Format section ends here.           !ATMD930
      endif                                                                       !ATMD931
!......................................................................           !ATMD932



Note that the header for the special output is written in the subroutine
initial_E10.f.  The following code section shows the header that is written
for the special output list above (at line ATMD880).


!......................................................................           !INIT593
      if (iopp /= 0) then                                                         !INIT594
!-----  This is the header for the "special" output file.  The user               !INIT595
!       should insure that the 200 format is compatible with the output           !INIT596
!       format in subroutine atmod near label number 9000 in the                  !INIT597
!       models_E10.f90 file.                                                      !INIT598
         write (iopp, 200)                                                        !INIT599
!                                                                                 !INIT600
!---  Header format for write statement in normal release code                    !INIT601
!                                                                                 !INIT602
  200    format('      Time    Hgtkm GeocenLat  Lon(East)   DensMean',       &    !INIT603
            '   PresMean   Tmean  EWmean  NSmean   DensPert   PresPert ',    &    !INIT604
            '  Tpert  EWpert  NSpert Dpert% Ppert% SDden% SDprs%',           &    !INIT605
            ' SDtemK SDuwnd SDvwnd SDwwnd  Wpert  SpdAvg  SpdStd SoSmean',   &    !INIT606
            ' SoSpert Sev LC   z0(m)')                                            !INIT607
!---  Header formats corresponding to other example write statements              !INIT608
!     given in file models_E10.f90                                                !INIT609
!                                                                                 !INIT610
! 200 Format(' Hgtkm  SigmaU  SigmaV  SigmaD   Upert   Vpert   Dpert')            !INIT611
! 200 Format(' Timehrs  Hgtkft    PAvgmb    PStdmb   Ppertmb',                    !INIT612
!    &   '   DAvggm3   DStdgm3  Dpertgm3  TAvgdegC  TstddegC TpertdegC'           !INIT613
!    &   '  dPsmn%  dDsmn%  dTsmn%')                                              !INIT614
! 200 Format(' Timehrs  Hgtkft  dPsml%  dPlrg%  dPtot%  dDsml%',                  !INIT615
!    &  '  dDlrg%  dDtot%  dTsml%  dTlrg%  dTtot%  dPsmn%  dDsmn%',               !INIT616
!    &  '  dTsmn%')                                                               !INIT617
! 200   Format('  Hgtkm   spS   sdS   stS   suS   svS   spL   sdL',               !INIT618
!    &   '   stL   suL   svL    spT    sdT    stT   suT    svT ')                 !INIT619
! 200   Format(' Timehrs   Hgtkm Uwndmps Vwndmps')                                !INIT620
! 200   Format(' Time_sec   Hgt_km     Lat    LonE  Dens_kgm3 D-76%',             !INIT621
!    &    ' sigD% dDen% Temp_K sigT% dTem% Uwnm/s sigUms dU_m/s',                 !INIT622
!    &    ' Vwnm/s sigVms dV_m/s')                                                !INIT623
      endif                                                                       !INIT624
!-----  The "special" output header Write and Format section ends here.           !INIT625
!......................................................................           !INIT626



As configured above, (ATMD880-ATMD884 and INIT599-INIT607) the following
data are output to the "special-formatted" output file:

 Header     Description
---------  -------------------------------------------------------------
Time       Mission elapsed time (sec) from initial input date/time
Hgtkm      Current height above sea level (km)
GeocenLat  Current geocentric latitude (degrees)
Lon(East)  Current longitude (+ East or - West), degrees
DensMean   Mean atmospheric density (kg/m**3)
PresMean   Mean atmospheric pressure (N/m**2)
Tmean      Mean atmospheric temperature (K)
EWmean     Mean eastward wind component (m/s) (+ toward East)
NSmean     Mean northward wind component (m/s) (+ toward North)
DensPert   Perturbed (mean-plus-perturbation) density (kg/m**3)
PresPert   Perturbed pressure (mean-plus-perturbation) (N/m**2)
Tpert      Perturbed temperature (mean-plus-perturbation) (K)
EWpert     Perturbed eastward wind (mean-plus-perturbation) (m/s)
NSpert     Perturbed northward wind (mean-plus-perturbation) (m/s)
Dpert%     Density perturbation (percent of mean value)
Ppert%     Pressure perturbation (percent of mean value)
SDden%     Total density standard deviation (percent of mean)
SDprs%     Total pressure standard deviation (percent of mean)
SDtemK     Total temperature standard deviation (K)
SDuwnd     Total eastward wind standard deviation (m/s)
SDvwnd     Total northward wind standard deviation (m/s)
SDwwnd     Standard deviation of vertical wind component (m/s)
Wpert      Vertical wind perturbation (m/s)
SpdAvg     Mean wind speed (m/s)
SpdStd     Standard deviation of wind speed (m/s)
SoSmean    Mean sound speed (m/s)
SoSpert    Perturbed sound speed (m/s)
Sev        Turbulence severity level (0 = normal, 1 = severe)
LC         Surface type code used in vertical wind model (See README7.txt)
z0(m)      Surface roughness length (m) used in vertical wind model



By changing this header (at lines INIT603-INIT607, above) and the write statement
and format (at lines ATMD880-ATMD884, above), the user can output any of the 
variables described in the following table (and listed in abbreviated form in
code lines ATMD779-ATMD869.


               Variables Available For Special-Formatted Output

Position and Time:
-ch     Geocentric height above sea-level (km), based on WGS 84 reference ellipsoid
-Rlocal Radius from center of Earth (km) (WGS 84 Earth radius plus height h)
-phi    Geocentric Latitude (deg) 
-GdLat  Geodetic Latitude (deg) 
-thet   Longitude (deg), East(+) West(-) 
-elt    Elapsed Time (sec) from initial time in input file


Pressure or Vapor Pressure:
-pgh  	Mean pressure in N/m**2
-pghp  	Deviation of mean pressure from 1976 US Standard Atmosphere (%)
-prhs  	Small scale perturbation of pressure in % of mean
-sphs	Standard deviation of the small scale pressure perturbations in % of mean
-prhl	Large scale perturbation of pressure in % of mean
-sphl	Standard deviation of the large scale pressure perturbations in % of mean
-prh	Total of large and small perturbations of pressure in % of mean
-sph	Total standard deviation of the large and small pressure perturbations in % of mean
-ph	Total pressure (mean PLUS perturbations) in N/m**2
-php	Deviation of total pressure (mean PLUS perturbations) from 1976 US Standard Atmosphere (%)
-eofT	Mean vapor pressure of water in N/m**2
-seofT	Standard deviation of water vapor pressure in N/m**2

Density or Vapor Density:
-dgh	Mean density in kg/m**3
-dghp	Deviation of mean density from 1976 US Standard Atmosphere (%)
-drhs	Small scale perturbation of density in % of mean
-sdhs	Standard deviation of the small scale density perturbations in % of mean
-drhl	Large scale perturbation of density in % of mean
-sdhl	Standard deviation of the large scale density perturbations in % of mean
-drh	Total of large and small perturbations of density in % of mean
-sdh	Total standard deviation of the large and small density perturbations in % of mean	
-dh	Total density (mean PLUS perturbations) in kg/m**3
-dhp	Deviation of total density (mean PLUS perturbations) from 1976 US Standard Atmosphere (%)
-rhov	Mean water vapor density in kg/m**3
-srhov	Standard deviation of the water vapor density in kg/m**3

Temperature or Dewpoint Temperature:
-tgh	Mean temperature in K
-tghp	Deviation of mean temperature from 1976 US Standard Atmosphere (%)
-trhs	Small scale perturbation of temperature in % of mean
-sths	Standard deviation of the small scale temperature perturbations in % of mean
-trhl	Large scale perturbation of temperature in % of mean
-sthl	Standard deviation of the large scale temperature perturbations in % of mean
-trh	Total of large and small perturbations of temperature in % of mean
-sth	Total standard deviation of the large and small temperature perturbations in % of mean
-th	Total temperature (mean PLUS perturbations) in K
-thp	Deviation of total temperature (mean PLUS perturbations) from 1976 US Standard Atmosphere (%)
-tdgh	Mean dewpoint temperature in K
-stdgh	Standard deviation of the dewpoint temperature in K

Eastward (u) Wind Component (positive toward East):
-ugh	Mean u component (+ toward East) in m/s
-urhs	Small scale perturbation in u component in m/s
-suhs	Standard deviation of the small scale perturbation of the u component in m/s
-urhl	Large scale perturbation in u component in m/s
-suhl	Standard deviation of the large scale perturbation of the u component in m/s
-urh	Total perturbation (large and small) of the u component in m/s
-suh	Standard deviation of the total perturbations (large and small) of the u component in m/s
-uh	Total u component (mean plus large and small perturbations) in m/s

Northward (v) Wind Component (positive toward North):
-vgh	Mean v component (+ toward North) in m/s
-vrhs	Small scale perturbation in v component in m/s
-svhs	Standard deviation of the small scale perturbation of the v component in m/s
-vrhl	Large scale perturbation in v component in m/s
-svhl	Standard deviation of the large scale perturbation of the v component in m/s
-vrh	Total perturbation (large and small) of the v component in m/s
-svh	Standard deviation of the total perturbations (large and small) of the v component in m/s
-vh	Total v component (mean plus large and small perturbations) in m/s

Other Wind Statistics:
-uvt2   Cross correlation between Eastward and Northward wind components (unitless)
-spdgh  Monthly mean speed (m/s)
-sdsph  Monthly standard deviation of wind speed (m/s)

Vertical (w) Wind Component (positive upward):
-wgh	Mean w component in m/s	
-wrh	Total perturbation (small scale, there are no large scale) of w component in m/s
-swh	Standard deviation of the total perturbation (small scale, there are no large scale) of w 
          component in m/s
-wh	Total w component (mean plus perturbations) in m/s

Surface roughness information for vertical wind model:
-lc     Land surface code (0-13 if from 1-by-1 deg data base or 99 if z0 is from user input)  
-z0     Surface roughness length (m), determined from 1-by-1 deg land surface type, or as specified 
          by user via input z0in

Speed of Sound:
-csp0   Mean sound speed (m/s)
-csp    Perturbed sound speed (m/s)

Relative Humidity:
-rhp	Mean relative humidity in %
-srhp	Standard deviation of the relative humidity in %

Surface Data:
-psrf	Average Surface pressure in N/m**2
-dsrf	Average Surface density in kg/m**3
-tsrf	Average Surface temperature in K
-usrf	Average Surface eastward wind component in m/s
-vsrf	Average Surface northward wind component in m/s
-hsrf	Height of surface in m above sea level (elevation)
-tdsrf	Average Surface dewpoint in K
-spsrf	Standard deviation of surface pressure in N/m**2
-sdsrf	Standard deviation of surface density in kg/m**3
-stsrf	Standard deviation of surface temperature in K
-susrf	Standard deviation of surface eastward component in m/s
-svsrf	Standard deviation of surface northward component in m/s
-shsrf	Standard deviation (uncertainty) of the surface height (elevation) in m
-stdsrf	Standard deviation of the surface dewpoint in K
-spdavsrf Monthly average wind speed at surface (m/s)                   
-spdsdsrf Standard deviation of surface wind speed (m/s)         
-uvtsrf   U-V cross correlation at surface (unitless)              

Species Concentrations and Number Densities:
-ppmh2o Mole fraction (volume mixing ratio) for water vapor (H2O) in parts per million
-ppmo3  Mole fraction (volume mixing ratio) for ozone (O3) in parts per million
-ppmn2o Mole fraction (volume mixing ratio) for nitrous oxide (N2O) in parts per million
-ppmco  Mole fraction (volume mixing ratio) for carbon monoxide (CO) in parts per million
-ppmch4 Mole fraction (volume mixing ratio) for methane (CH4) in parts per million
-ppmco2 Mole fraction (volume mixing ratio) for carbon dioxide (CO2) in parts per million
-ppmn2  Mole fraction (volume mixing ratio) for molecular nitrogen (N2) in parts per million
-ppmo2  Mole fraction (volume mixing ratio) for molecular oxygen (O2) in parts per million
-ppmo   Mole fraction (volume mixing ratio) for atomic oxygen (O) in parts per million
-ppmar  Mole fraction (volume mixing ratio) for argon (Ar) in parts per million
-ppmhe  Mole fraction (volume mixing ratio) for helium (He) in parts per million
-ppmh   Mole fraction (volume mixing ratio) for atomic hydrogen (H) in parts per million
-ppmn   Mole fraction (volume mixing ratio) for atomic nitrogen (N) in parts per million
-h2ond  Number density for water vapor (H2O) in #/m**3
-o3nd   Number density for ozone (O3) in #/m**3
-n2ond  Number density for nitrous oxide (N2O) in #/m**3
-cond   Number density for carbon monoxide (CO) in #/m**3
-ch4nd  Number density for methane (CH4) in #/m**3
-co2nd  Number density for carbon dioxide (CO2) in #/m**3
-n2nd   Number density for molecular nitrogen (N2) in #/m**3
-o2nd   Number density for molecular oxygen (O2) in #/m**3
-ond    Number density for atomic oxygen (O) in #/m**3
-arnd   Number density for argon (Ar) in #/m**3
-hend   Number density for helium (He) in #/m**3
-hnd    Number density for atomic hydrogen (H) in #/m**3
-nnd    Number density for atomic nitrogen (N) in #/m**3
