!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_rad_data.f90
!!
!! \brief
!! mo_ham_rad_data holds the parameters and routines for the calculation of the 
!! optical parameters for the aerosol distribution simulated by the ECHAM/HAM aerosol module. 
!!
!! \author Philip Stier (MPI-Met)
!!
!! \responsible_coder
!! Philip Stier, philip.stier@physics.ox.ac.uk
!!
!! \revision_history
!!   -# P. Stier (MPI-Met) - original code (2003-05-08)
!!   -# P. Stier (Caltech)  (2006-02)
!!   -# S. Rast (MPI-Met) - adaption to echam5.3.2 (2007-03)
!!   -# D. O'Donnell (MPI-Met) - changed to flexible number of species and added SOA (XXXX)
!!   -# M.G. Schultz (FZ Juelich) - added irad... parameters and removed pseudo-flexibility
!!                                 (this module requires manual editing anyhow if number of species is changed) (XXXX)
!!   -# P. Stier (Uni Oxford) - adaptation to RRTM-SW (2010-02)
!!   -# H. Kokkola (FMI) - modified to include SALSA (2013-06) 
!!
!! \limitations
!! None
!!
!! \details
!! None
!!
!! \bibliographic_references
!!    - M. Hess, P. Koepke, and I. Schult (1998): Optical Properties of Aerosols and clouds:
!!      The software package OPAC, Bull. Am. Met. Soc., 79, 831-844.
!!    - Nilsson, B. (1979): Meteorological influence on aerosol extinction in the 0.2-40 um  
!!      wavelength range, Appl. Opt. 18, 3457-3473.
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!! Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!! licencing agreement to be found at:
!! https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!! The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_ham_rad_data

  USE mo_kind,          ONLY: dp

  IMPLICIT NONE

  PRIVATE

  !---public data
  PUBLIC :: iradso4, iradbc, iradoc, iradss, iraddu, iradwat, iradsoa

  PUBLIC :: Nwv_sw, Nwv_sw_opt, Nwv_lw,       &
            Nwv_tot, Nwv_swlw, Nwv_sw_tot,    &
            cnr, cni,                         &
            log_x0_min, log_x0_max,           &
            log_ni_min,                       &
            x0_min, x0_max,                   &
            nnrmax, nnimax, ndismax,          &
            nr_min, nr_max, ni_min, ni_max,   &
            inc_nr, inc_ni,                   &
            lambda, lambda_sw_opt,            &
            nradang


  !---public member functions

  PUBLIC :: ham_rad_data_initialize, nraddiagwv

  !--- species indices for radiation tables

  INTEGER, PARAMETER :: iradso4 = 1
  INTEGER, PARAMETER :: iradbc = 2
  INTEGER, PARAMETER :: iradoc = 3
  INTEGER, PARAMETER :: iradss = 4
  INTEGER, PARAMETER :: iraddu = 5
  INTEGER, PARAMETER :: iradwat = 6
  INTEGER, PARAMETER :: iradsoa = 7

  INTEGER, PARAMETER :: naeroradspec = 7    ! for dimensioning of tables

  !--- Parameter space for the 4 lookup tables:
  !
  !    Table 1 : SW fine    std. dev: sigma_fine (1.59)
  !          2 : SW coarse            sigma_coarse (2.00)
  !          3 : LW fine              sigma_fine (1.59)
  !          4 : LW coarse            sigma_coarse (2.00)

  INTEGER, PARAMETER  :: Ndismax(4)=(/ 100, 100, 100, 100 /) 
  INTEGER, PARAMETER  :: Nnrmax(4) =(/ 100, 100, 100, 100 /)
  INTEGER, PARAMETER  :: Nnimax(4) =(/ 200, 200, 200, 200 /)

  REAL(dp), PARAMETER :: nr_min(4)=(/  1.33_dp,   1.33_dp,  1.0_dp,   1.0_dp    /) !--min ref. ind. real
  REAL(dp), PARAMETER :: nr_max(4)=(/  2.00_dp,   2.00_dp,  3.0_dp,   3.0_dp    /) !--max ref. ind. real
  REAL(dp), PARAMETER :: ni_min(4)=(/  1.E-9_dp,  1.E-9_dp, 1.E-9_dp, 1.E-9_dp  /) !--min ref. ind. imag.
  REAL(dp), PARAMETER :: ni_max(4)=(/  1.00_dp,   1.00_dp,  2.0_dp,   2.0_dp    /) !--max ref. ind. imag.

  REAL(dp) :: x0_min(4)
  REAL(dp) :: x0_max(4)
  REAL(dp) :: log_x0_min(4) 
  REAL(dp) :: log_x0_max(4) 
  REAL(dp) :: log_nr_min(4) 
  REAL(dp) :: log_nr_max(4) 
  REAL(dp) :: log_ni_min(4) 
  REAL(dp) :: log_ni_max(4) 

  REAL(dp) :: inc_nr(4)
  REAL(dp) :: inc_ni(4)

  !--- Indices and dimensions:

  INTEGER, PARAMETER :: Nwv_sw     = 14, &   !-- RRTM-SW GCM wavebands
                        Nwv_sw_opt = 2,  &   !--Optional SW wavebands (550nm+865nm)
!                       Nwv_sw_opt = 2,  &   !--Optional SW wavebands (550nm+865nm+440nm)
                        Nwv_sw_tot = Nwv_sw+Nwv_sw_opt, &
                        Nwv_lw     = 16, &   !--RRTM-LW GCM wavebands
                        Nwv_swlw   = Nwv_sw+Nwv_lw, &
                        Nwv_tot    = Nwv_sw+Nwv_sw_opt+Nwv_lw

  REAL(dp), PARAMETER :: lambda_sw_opt(2)=(/ 0.550E-6_dp, 0.865E-6_dp /)
!  REAL(dp), PARAMETER :: lambda_sw_opt(2)=(/ 0.550E-6_dp, 0.865E-6_dp , 0.440E-6_dp/)

  !--- Define mask for output of wavelengths: 
  !
  !    nraddiagwv(jwv)=0   Off
  !                   =1   AOD
  !                   =2   AOD+AAOD
  !                   =3   AOD+AAOD+...

  INTEGER :: nraddiagwv(Nwv_tot)

  INTEGER :: nradang(2)

  !--- 1) Mid-band wavelengths:

  REAL(dp)     :: lambda(1:Nwv_tot) 

  !--- 2) LW refractive indices:
  
  REAL(dp)     :: cnr(Nwv_tot, naeroradspec),    &   ! real part of refractive index
                  cni(Nwv_tot, naeroradspec)         ! imaginary part of refractive index


CONTAINS

  SUBROUTINE ham_rad_data_initialize

    !---inherited types, data and functions
    !   -

    USE mo_ham,  ONLY: naerorad, nrad, nham_subm, HAM_M7, HAM_SALSA

    IMPLICIT NONE

    !---subroutine interface
    !   -

    !---local variables
    INTEGER :: itable

    !--- Initialize wavelength settings


    !--- Default output masks for wavelengths, currently only for the optional wavelengths:
    !ham_ps:radiation This could/should become namelist controlled 

    IF (naerorad>0) THEN

       nraddiagwv(1:NWv_tot)=0                          !--- Default: no diagnostics

       nraddiagwv(Nwv_sw+1:Nwv_sw+Nwv_sw_opt)=1         ! AOD only for all optional wavelengths
       nraddiagwv(Nwv_sw+1)=2                           ! Additional diagnostics for 550nm

       IF (ANY(nrad(:)==2) .OR. ANY(nrad(:)==3)) THEN
          nraddiagwv(Nwv_sw+Nwv_sw_opt+1)=1             ! AOD for first LW band
       END IF

       !--- Define wavelengths for Angstroem diagnostics (currently first two optional wavelengths):

       nradang(1)=Nwv_sw+1
       nradang(2)=nradang(1)+1

    END IF

    !--- 1.0) SW Refractive Indices:

    !    Note the RRTM-SW specific band ordering:
    !    3.46, 2.79, 2.33, 2.05, 1.78, 1.46, 1.27, 1.01, 0.70, 0.53, 0.39, 0.30, 0.23, 8.02 [um]

    !--- Sulfate (Palmer & Williams, Appl.Opt., 1975; provided by Stefan Kinne):

    cnr(1:Nwv_sw,iradso4) = &
         (/ 1.361_dp,      1.295_dp,      1.364_dp,      1.382_dp,      1.393_dp,      1.406_dp,      &
            1.413_dp,      1.422_dp,      1.427_dp,      1.432_dp,      1.445_dp,      1.450_dp,      &
            1.450_dp,      1.400_dp                                                                   /)

   cni(1:Nwv_sw,iradso4) = &
        (/ 1.400E-01_dp,   5.500E-02_dp,  2.100E-03_dp,  1.300E-03_dp,  5.100E-04_dp,  9.000E-05_dp,  &
           7.900E-06_dp,   1.300E-06_dp,  5.200E-08_dp,  1.000E-09_dp,  1.000E-09_dp,  1.000E-09_dp,  &
           1.000E-09_dp,   2.600E-01_dp                                                               /)

   !--- Optional Wavelengths at 550 and 865nm:

   cnr(Nwv_sw+1:Nwv_sw+Nwv_sw_opt,iradso4) =                                                          &
        (/ 1.432_dp,       1.424_dp                                                                   /)

   cni(Nwv_sw+1:Nwv_sw+Nwv_sw_opt,iradso4) =                                                          &
        (/ 1.000E-09_dp,   7.384E-07_dp                                                               /)

   !--- Black Carbon (Medium-absorbing values from Bond & Bergstrom, 2006):
   !
   !    - at 550 nm: n=1.85-0.71i
   !    - wavelength dependency scaled from OPAC (Hess et al., 1998)

   cnr(1:Nwv_sw,iradbc) = &
        (/1.984_dp,       1.936_dp,      1.917_dp,      1.905_dp,      1.894_dp,      1.869_dp,      &
          1.861_dp,       1.861_dp,      1.850_dp,      1.850_dp,      1.839_dp,      1.839_dp,      &
          1.713_dp,       2.245_dp                                                                   /)

   cni(1:Nwv_sw,iradbc) = &
        (/ 8.975E-01_dp,  8.510E-01_dp,  8.120E-01_dp,  7.939E-01_dp,  7.765E-01_dp,  7.397E-01_dp,  &
           7.274E-01_dp,  7.106E-01_dp,  6.939E-01_dp,  7.213E-01_dp,  7.294E-01_dp,  7.584E-01_dp,  &
           7.261E-01_dp,  1.088E+00_dp                                                               /)

   !--- Optional Wavelengths at 550nm and 865nm:

   cnr(Nwv_sw+1:Nwv_sw+Nwv_sw_opt,iradbc) =                                                          &
        (/ 1.85_dp,       1.85_dp                                                                    /)

   cni(Nwv_sw+1:Nwv_sw+Nwv_sw_opt,iradbc) =                                                          &
        (/ 7.10E-01_dp,   6.99E-01_dp                                                                /)

   !--- Organic Carbon (Medium-absorbing values from Bond & Bergstrom, 2006):
   !
   !    - wavelength dependency scaled from OPAC (Hess et al., 1998)

   cnr(1:Nwv_sw,iradoc) = &
        (/ 1.530_dp,      1.510_dp,      1.510_dp,      1.420_dp,      1.464_dp,      1.520_dp,      &
           1.420_dp,      1.420_dp,      1.530_dp,      1.530_dp,      1.530_dp,      1.443_dp,      &
           1.530_dp,      1.124_dp                                                                   /)

   cni(1:Nwv_sw,iradoc) = &
        (/ 2.75E-02_dp,   7.33E-03_dp,   7.33E-03_dp,   4.58E-03_dp,   6.42E-03_dp,   1.43E-02_dp,    &
           1.77E-02_dp,   2.01E-02_dp,   1.50E-02_dp,   7.70E-03_dp,   9.75E-03_dp,   1.63E-02_dp,    &
           5.27E-03_dp,   7.24E-02_dp                                                                /)

   !--- Optional Wavelengths at 550nm and 865nm:

   cnr(Nwv_sw+1:Nwv_sw+Nwv_sw_opt,iradoc) =                                                          &
        (/ 1.53_dp,       1.52_dp                                                                    /)

   cni(Nwv_sw+1:Nwv_sw+Nwv_sw_opt,iradoc) =                                                          &
        (/ 5.50E-03_dp,   1.10E-02_dp                                                                /)

   !--- Sea Salt (Nilsson, 1979):

   cnr(1:Nwv_sw,iradss) = &
        (/ 1.480_dp,      1.400_dp,      1.440_dp,      1.450_dp,      1.450_dp,      1.460_dp,      &
           1.470_dp,      1.470_dp,      1.480_dp,      1.490_dp,      1.500_dp,      1.510_dp,      &
           1.510_dp,      1.400_dp                                                                   /)

   cni(1:Nwv_sw,iradss) = &
        (/ 1.300E-02_dp,  8.000E-03_dp,  2.500E-03_dp,  1.500E-03_dp,  1.000E-03_dp,  5.500E-04_dp,  &
           3.300E-04_dp,  1.000E-04_dp,  1.000E-07_dp,  1.000E-08_dp,  2.000E-08_dp,  1.000E-06_dp,  &
           1.000E-05_dp,  1.400E-02_dp                                                               /)

   !--- Optional Wavelengths at 550nm and 865nm:

   cnr(Nwv_sw+1:Nwv_sw+Nwv_sw_opt,iradss) =                                                          &
        (/ 1.450_dp,      1.470_dp                                                                   /)

   cni(Nwv_sw+1:Nwv_sw+Nwv_sw_opt,iradss) =                                                          &
        (/ 1.000E-08_dp,  1.000E-8_dp                                                                /)

  !--- Dust (provided by Stefan Kinne, MPI-MET:
  !    mainly based on Sokolik & Toon, JGR, 1999)
  !    imaginary parts in the visible modified according to AERONET statistics
  !    much less absorption in the VIS (mainly based on Kinne et.al, JGR, 2003):

   cnr(1:Nwv_sw,iraddu) = &
        (/ 1.460_dp,      1.460_dp,      1.460_dp,      1.450_dp,      1.450_dp,      1.450_dp,      &
           1.450_dp,      1.450_dp,      1.450_dp,      1.450_dp,      1.450_dp,      1.450_dp,      &
           1.450_dp,      1.170_dp                                                                   /)

   cni(1:Nwv_sw,iraddu) = &
        (/ 1.180E-02_dp,  6.000E-03_dp,  2.500E-03_dp,  1.500E-03_dp,  1.000E-03_dp,  8.000E-04_dp,  &
           6.000E-04_dp,  7.500E-04_dp,  9.500E-04_dp,  1.000E-03_dp,  2.500E-03_dp,  2.000E-02_dp,  &
           2.500E-02_dp,  1.000E-01_dp                                                               /)

   !--- Optional Wavelengths at 550nm and 865nm:

   cnr(Nwv_sw+1:Nwv_sw+Nwv_sw_opt,iraddu) =                                                          &
        (/ 1.450_dp,      1.450_dp                                                                   /)

   cni(Nwv_sw+1:Nwv_sw+Nwv_sw_opt,iraddu) =                                                          &
        (/ 1.000E-03_dp,  8.400E-4_dp                                                                /)

   !--- Water (provided by Stefan Kinne interpolated with code by Andy Lacis (NASA-GISS):
   !          (Hale and Querry [1973] for 0.2 to 0.7 mm, 
   !          Palmer and Williams [1974] for 0.7 to 2.0 mm)

   cnr(1:Nwv_sw,iradwat) = &
        (/ 1.423_dp,      1.244_dp,      1.283_dp,      1.300_dp,      1.312_dp,      1.319_dp,      &
           1.324_dp,      1.328_dp,      1.331_dp,      1.335_dp,      1.341_dp,      1.350_dp,      &
           1.377_dp,      1.300_dp                                                                   /)

   cni(1:Nwv_sw,iradwat) = &
        (/ 5.000E-02_dp,  1.300E-01_dp,  6.500E-04_dp,  6.700E-04_dp,  1.200E-04_dp,  1.100E-04_dp,  &
           1.200E-05_dp,  2.100E-06_dp,  6.800E-08_dp,  2.800E-09_dp,  3.900E-09_dp,  1.700E-08_dp,  &
           6.400E-08_dp,  4.000E-02_dp                                                               /)

   !--- Optional Wavelengths at 550nm and 865nm:

   cnr(Nwv_sw+1:Nwv_sw+Nwv_sw_opt,iradwat) =                                                         &
        (/ 1.335_dp,      1.329_dp                                                                   /)

   cni(Nwv_sw+1:Nwv_sw+Nwv_sw_opt,iradwat) =                                                         &
        (/ 2.800E-09,     1.186E-06                                                                  /)


   !--- 2) LW refractive indices:
   !    Note that the RRTM spectral bands are ordered with increasing wavenumber,
   !    i.e. decreasing wavelengths.

   !--- Sulfate (Ammonium Sulfate from Toon and Pollack, 1976)
   !             Interpolated with subroutine of Andy Lacis (NASA-GISS)

   cnr(Nwv_sw+Nwv_sw_opt+1:Nwv_tot,iradso4) = &
        (/ 1.889_dp,      1.588_dp,      1.804_dp,      1.537_dp,      1.709_dp,      1.879_dp,      &
           2.469_dp,      0.685_dp,      1.427_dp,      0.956_dp,      1.336_dp,      1.450_dp,      &
           1.489_dp,      1.512_dp,      1.541_dp,      1.602_dp                               /)

   cni(Nwv_sw+Nwv_sw_opt+1:Nwv_tot,iradso4) = &
        (/ 0.967E-01_dp,  0.380E-01_dp,  0.287E-01_dp,  0.225E-01_dp,  0.200E-01_dp,  0.396E-01_dp,  &
           0.269E+00_dp,  0.111E+01_dp,  0.705E-01_dp,  0.678E+00_dp,  0.143E-01_dp,  0.664E-02_dp,  &
           0.657E-02_dp,  0.944E-02_dp,  0.148E-01_dp,  0.156E+00_dp                           /)

   !--- Black Carbon (Medium-absorbing values from Bond & Bergstrom, 2006):
   !
   !         - at 550 nm: n=1.85-0.71i
   !         - wavelength dependency scaled from OPAC (Hess et al., 1998)

   cnr(Nwv_sw+Nwv_sw_opt+1:Nwv_tot,iradbc) =  &
        (/ 2.84_dp,       2.63_dp,       2.53_dp,       2.46_dp,       2.42_dp,       2.36_dp,       &
           2.33_dp,       2.30_dp,       2.23_dp,       2.17_dp,       2.14_dp,       2.09_dp,       &
           2.06_dp,       2.04_dp,       2.03_dp,       1.98_dp                                /)

   cni(Nwv_sw+Nwv_sw_opt+1:Nwv_tot,iradbc) = &
        (/ 1.61E+00_dp,   1.42E+00_dp,   1.33E+00_dp,   1.28E+00_dp,   1.23E+00_dp,   1.18E+00_dp,   &
           1.16E+00_dp,   1.14E+00_dp,   1.08E+00_dp,   1.04E+00_dp,   1.00E+00_dp,   9.73E-01_dp,   &
           9.56E-01_dp,   9.46E-01_dp,   9.37E-01_dp,   8.91E-01_dp                            /)

   !--- Organic Carbon:
   !    OPAC WaSo Category (M. Hess, P. Koepke, and I. Schult, 1998):

   cnr(Nwv_sw+Nwv_sw_opt+1:Nwv_tot,iradoc) = &
        (/ 1.86_dp,       1.95_dp,       2.02_dp,       1.43_dp,       1.61_dp,       1.71_dp,       &
           1.81_dp,       2.64_dp,       1.23_dp,       1.42_dp,       1.42_dp,       1.45_dp,       &
           1.46_dp,       1.46_dp,       1.46_dp,       1.44_dp                                /)

   cni(Nwv_sw+Nwv_sw_opt+1:Nwv_tot,iradoc) = &
        (/ 4.58E-01_dp,   2.35E-01_dp,   1.86E-01_dp,   1.82E-01_dp,   5.31E-02_dp,   4.52E-02_dp,   &
           4.54E-02_dp,   3.76E-01_dp,   6.04E-02_dp,   5.30E-02_dp,   2.29E-02_dp,   1.27E-02_dp,   &
           1.17E-02_dp,   9.28E-03_dp,   4.88E-03_dp,   5.95E-03_dp                            /)

   !--- Sea Salt:
   !    (Shettle and Fenn, 1979; Nilsson, 1979; both based on Volz, 1972 measurements)
   !    Inter with subroutine of Andy Lacis (NASA-GISS)

   cnr(Nwv_sw+Nwv_sw_opt+1:Nwv_tot,iradss) = &
        (/ 1.668_dp,      1.749_dp,      1.763_dp,      1.447_dp,      1.408_dp,      1.485_dp,      &
           1.563_dp,      1.638_dp,      1.401_dp,      1.450_dp,      1.505_dp,      1.459_dp,      &
           1.483_dp,      1.488_dp,      1.478_dp,      1.484_dp                               /)

   cni(Nwv_sw+Nwv_sw_opt+1:Nwv_tot,iradss) = &
        (/ 0.981E+00_dp,  0.193E+00_dp,  0.111E+00_dp,  0.344E-01_dp,  0.192E-01_dp,  0.140E-01_dp,  &
           0.179E-01_dp,  0.293E-01_dp,  0.138E-01_dp,  0.543E-02_dp,  0.180E-01_dp,  0.288E-02_dp,  &
           0.251E-02_dp,  0.246E-02_dp,  0.175E-02_dp,  0.206E-02_dp                           /)

  !--- Dust (Irina Sokolik, personal communication, 2006):

   cnr(Nwv_sw+Nwv_sw_opt+1:Nwv_tot,iraddu) = &
        (/ 2.552_dp,      2.552_dp,      1.865_dp,      1.518_dp,      1.697_dp,      1.816_dp,      &
           2.739_dp,      1.613_dp,      1.248_dp,      1.439_dp,      1.423_dp,      1.526_dp,      &
           1.502_dp,      1.487_dp,      1.480_dp,      1.468_dp                               /)

   cni(Nwv_sw+Nwv_sw_opt+1:Nwv_tot,iraddu) = &
        (/ 0.7412_dp,     0.7412_dp,     0.5456_dp,     0.2309_dp,     0.1885_dp,     0.2993_dp,     &
           0.7829_dp,     0.4393_dp,     0.1050_dp,     0.0976_dp,     0.0540_dp,     0.0228_dp,     &
           0.0092_dp,     0.0053_dp,     0.0044_dp,     0.0101_dp                              /)

  !--- Water:
  !    (Downing and Williams, 1975)
  !     Interpolated with subroutine by Andy Lacis (NASA-GISS)

   cnr(Nwv_sw+Nwv_sw_opt+1:Nwv_tot,iradwat) = &
        (/ 1.689_dp,      1.524_dp,      1.401_dp,      1.283_dp,      1.171_dp,      1.149_dp,      &
           1.230_dp,      1.264_dp,      1.295_dp,      1.314_dp,      1.312_dp,      1.316_dp,      &
           1.327_dp,      1.333_dp,      1.348_dp,      1.416_dp                               /)

   cni(Nwv_sw+Nwv_sw_opt+1:Nwv_tot,iradwat) = &
        (/ 0.618E+00_dp,  0.392E+00_dp,  0.428E+00_dp,  0.395E+00_dp,  0.317E+00_dp,  0.107E+00_dp,  &
           0.481E-01_dp,  0.392E-01_dp,  0.347E-01_dp,  0.348E-01_dp,  0.132E+00_dp,  0.106E-01_dp,  &
           0.151E-01_dp,  0.881E-02_dp,  0.483E-02_dp,  0.169E-01_dp                           /)

   !--- Secondary Organic Carbons:
   !    Assume values of primary organics, i.e. OPAC WaSo Category (M. Hess, P. Koepke, and I. Schult, 1998):

   cni(:,iradsoa) = cni(:,iradoc)
   cnr(:,iradsoa) = cnr(:,iradoc)

    !---calculations for optimisation of reading of lookup tables
   
   SELECT CASE(nham_subm)
      
       CASE(HAM_M7)
    
          x0_min=(/  0.001_dp,  0.4_dp,   5.E-6_dp, 0.0015_dp /) !--min Mie parameter
          x0_max=(/ 25.0_dp,   40.0_dp,   3._dp,    4._dp     /) !--max Mie parameter
          
       CASE(HAM_SALSA)
    
          x0_min=(/  0.001_dp,  0.16_dp,   5.E-6_dp, 0.0015_dp /) !--min Mie parameter
    
          x0_max=(/ 25.0_dp,   210.0_dp,   3._dp,    17._dp     /) !--max Mie parameter
          
   END SELECT

    DO itable=1,4

       log_x0_min(itable) = LOG(x0_min(itable))
       log_x0_max(itable) = LOG(x0_max(itable))

       log_nr_min(itable) = LOG(nr_min(itable))
       log_nr_max(itable) = LOG(nr_max(itable))
       
       log_ni_min(itable) = LOG(ni_min(itable))
       log_ni_max(itable) = LOG(ni_max(itable))
       
       !--- Real part:
       inc_nr(itable)=(nr_max(itable)-nr_min(itable))/REAL(Nnrmax(itable),dp)

       !--- Imaginary part in log-space:
       inc_ni(itable)=(log_ni_max(itable)-log_ni_min(itable))/REAL(Nnimax(itable),dp)

    END DO

    
 END SUBROUTINE ham_rad_data_initialize

END MODULE mo_ham_rad_data
