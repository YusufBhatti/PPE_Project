!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_m7_emi_seasalt.f90
!!
!! \brief
!! This module contains several sea salt emission schemes.
!!
!! \author Kai Zhang (MPI-Met) (see also the individual subroutines)
!!
!! \responsible_coder
!! Kai Zhang, kai.zhang@pnnl.gov
!!
!! \revision_history
!!   -# Kai Zhang (MPI-Met) - original code (2009)
!!
!! \limitations
!! [ Start an optional warning here ]
!!
!! \details
!! Use nseasalt in namelist hamctl to select the scheme you want:
!!   - seasalt_emissions_monahan:   nseasalt=1
!!   - seasalt_emissions_lsce:      nseasalt=2
!!   - seasalt_emissions_mh:        nseasalt=4
!!   - seasalt_emissions_guelle:    nseasalt=5
!!   - seasalt_emissions_gong:      nseasalt=6
!!   - seasalt_emissions_long       nseasalt=7
!!   - seasalt_emissions_gong_SST   nseasalt=8
!!
!! \bibliographic_references
!!   - see individual seasalt scheme routines
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

!----

MODULE mo_ham_m7_emi_seasalt

  !---inherited types, data and functions
  USE mo_kind,           ONLY: dp

  IMPLICIT NONE

  !---public member functions
  PUBLIC :: start_emi_seasalt
  PUBLIC :: seasalt_emissions_monahan     ! nseasalt=1
  PUBLIC :: seasalt_emissions_lsce        ! nseasalt=2
  !PUBLIC :: seasalt_emissions_martensson  !nseasalt=3
  PUBLIC :: seasalt_emissions_mh          !nseasalt=4
  PUBLIC :: seasalt_emissions_guelle      !nseasalt=5
  PUBLIC :: seasalt_emissions_gong        !nseasalt=6 
  PUBLIC :: seasalt_emissions_long        !nseasalt=7
  PUBLIC :: seasalt_emissions_gong_SST    !nseasalt=8

  !---module data
  REAL(dp), PARAMETER, PRIVATE :: ppww = 3.41_dp      ! exponent of wind speed |u| (|u|**ppww) 

  INTEGER,  PARAMETER, PRIVATE :: nbin = 300          ! number of bins for the bin schemes
                                                      ! (Monahan (nseasalt=4), Guelle or Gong, Long)
  REAL(dp), PARAMETER, PRIVATE :: dmta = 0.100E-06_dp ! lower diameter [m], bin schemes 
  REAL(dp), PARAMETER, PRIVATE :: dmtd = 1.000E-05_dp ! upper diameter [m], bin schemes 

  REAL(dp), PARAMETER, PRIVATE :: dmtb_gong = 0.221E-06_dp ! diameter limit, [m], corresponding to 
                                                           ! 0.2um wet radius at 80% RH (G03), 0.2*2/1.814 
  REAL(dp), PARAMETER, PRIVATE :: dmtb_guelle = 8.000E-06_dp ! dry diameter limit, [m], corresonding 
                                                             ! to 4um radius in G01 

  REAL(dp), PARAMETER, PRIVATE :: dbeg(3) = (/0.050E-6_dp, 0.100E-6_dp, 1.000E-6_dp/)  ! ait, acc, coa
  REAL(dp), PARAMETER, PRIVATE :: dend(3) = (/0.100E-6_dp, 1.000E-6_dp, 1.000E-5_dp/)  ! ait, acc, coa

  REAL(dp), PRIVATE :: dmt(nbin)   

  REAL(dp), PRIVATE :: rm(nbin)
  REAL(dp), PRIVATE :: rd(nbin)
  REAL(dp), PRIVATE :: bmn(nbin) 

  REAL(dp), PRIVATE  :: ss1_mon                   ! sea salt flux factor 1, Monahan (nseasalt=1) scheme    
  REAL(dp), PRIVATE  :: ss2_mon                   ! sea salt flux factor 2, Monahan (nseasalt=1) scheme    


CONTAINS

  SUBROUTINE start_emi_seasalt

    USE mo_kind,         ONLY: dp
    USE mo_math_constants, ONLY: pi
    USE mo_ham,          ONLY: nseasalt
    USE mo_exception,    ONLY: message, message_text, em_param
    USE mo_util_string,  ONLY: separator

    IMPLICIT NONE

    !---local variables
    !   intermediate variables in calculating flux factors ss1 and ss2 for the monahan (nseasalt=1) scheme
    REAL(dp) :: zr1, zr2, zb1, zb2, zx1, zx2, zdr1, zdr2
    REAL(dp) :: zfact

    !   variables for constructing bins in the bin schemes
    REAL(dp) :: zdx, zdd
    INTEGER :: m

    !---executable procedure

    !---initialize the monahan (nseasalt=1) scheme (copy/paste ECHAM5/HAM2 code)
    zr1=0.416_dp
    zr2=3.49_dp
    zb1=0.58_dp-1.54_dp*LOG10(zr1)
    zb2=0.58_dp-1.54_dp*LOG10(zr2)
    zx1=10._dp**(1.19_dp*EXP(-zb1*zb1))
    zx2=10._dp**(1.19_dp*EXP(-zb2*zb2))
    zdr1=0.5_dp
    zdr2=4.5_dp
    zfact=1.373_dp*4._dp/3._dp*pi*1.15e3_dp
    ss1_mon=zfact*(1._dp+0.057_dp*zr1**1.05_dp)*zx1*zdr1*1.e-18_dp
    ss2_mon=zfact*(1._dp+0.057_dp*zr2**1.05_dp)*zx2*zdr2*1.e-18_dp

    IF (nseasalt==1) THEN
       CALL message('', separator)
       CALL message('', 'Monahan seasalt emissions:', level=em_param)
       WRITE(message_text, '(a,e25.15,a,e25.15)') 'Factors for sea salt fluxes: 1st = ',ss1_mon, '2nd = ', ss2_mon
       CALL message('', message_text, level=em_param)
       CALL message('', separator)
    END IF

    !---construction of bins
    ! dDp, take LOG scale 

    zdx = (LOG(dmtd) - LOG(dmta) ) / REAL(nbin,dp)

    zdd = 0._dp

    DO m = 1, nbin
       dmt(m) = EXP(LOG(dmta) + zdd)
       zdd = zdd + zdx
    END DO

    ! dry radius (m)  
  
    rd(:) = dmt(:) * 0.5_dp 

    ! wet radius (um) at RH=80%

    rm(:) = 1.814_dp*rd(:)*1.E+06_dp 

    ! B: monahan and guelle schemes, also for larger particle in Gong scheme
    bmn(:) = ( 0.380_dp - LOG10( rm(:) ) ) / 0.650_dp

    ! B: overwrite for smaller particles in the Gong scheme
    IF (nseasalt == 6) THEN
       DO m = 2,nbin
          IF (dmt(m).GT.dmta .and. dmt(m).le.dmtb_gong) THEN 
             bmn(m) = ( 0.433_dp - LOG10( rm(m) ) ) / 0.433_dp
          END IF
       END DO
    END IF

  END SUBROUTINE start_emi_seasalt
    
  SUBROUTINE seasalt_emissions_monahan(kproma, kbdim, krow, pmassf_as, pmassf_cs, pnumf_as, pnumf_cs)
    
    USE mo_kind,         ONLY: dp
    USE mo_memory_g3b,   ONLY: slm , seaice
    USE mo_vphysc,       ONLY: vphysc

    IMPLICIT NONE
 
    ! arguments
    INTEGER, INTENT(in)    :: kproma          ! geographic block number of locations
    INTEGER, INTENT(in)    :: kbdim           ! geographic block maximum number of locations 
    INTEGER, INTENT(in)    :: krow            ! geographic block number 
    REAL(dp), INTENT(out)  :: pmassf_as(kbdim)! seasalt emission mass flux in soluble accumulation mode
    REAL(dp), INTENT(out)  :: pmassf_cs(kbdim)! ...                        in soluble coarse mode
    REAL(dp), INTENT(out)  :: pnumf_as(kbdim) ! sesalt particle number flux in soluble accumulation mode
    REAL(dp), INTENT(out)  :: pnumf_cs(kbdim) ! ...                        in soluble coarse mode

    !--- Local Variables:
    REAL(dp)               :: zzzspeed(kbdim)
    LOGICAL                :: loocean(kbdim)  ! ocean mask
    REAL(dp)               :: zseaice(kbdim)  ! sea ice fraction
    REAL(dp)               :: ztmp1(kbdim)    !SF #458 temporary variable

    !--- Constants:
    REAL(dp), PARAMETER    :: zssmassa= 0.37E-15_dp,   & ! Mass of accu mode ss particle [kg]
                              zssmassc= 0.37E-12_dp      ! Mass of coar mode ss particle [kg]

    !--- 0) Initialisation
    loocean(1:kproma)  = ( slm(1:kproma,krow) < 1.e-2_dp )
    zseaice(1:kproma) = seaice(1:kproma, krow)
    pmassf_as(:) = 0._dp
    pmassf_cs(:) = 0._dp
    pnumf_as(:)  = 0._dp
    pnumf_cs(:)  = 0._dp

    !>>dod import of seasalt emissions schemes from HAM2:
    !      initialization moved to start_emi_seasalt
    !<<dod

    !--- 2) Calculate the emitted mass-flux: ------------------------------------------
    zzzspeed(1:kproma) = MIN(vphysc%velo10m(1:kproma,krow),20._dp)
    !>>SF #458 (replacing WHERE statements)
    ztmp1(1:kproma) = zzzspeed(1:kproma)**3.41_dp*(1._dp-zseaice(1:kproma))
    !SF Note that the fact to use the above variable ztmp1 in the following lines, instead of 
    !   writing 2 times the same expression yields non bit-id results as compared to previous
    !   version. Inevertheless prefer to keep this formulation, as it saves some computing.
    pmassf_as(1:kproma)  = MERGE( &
                                ss1_mon*ztmp1(1:kproma), &
                                pmassf_as(1:kproma), &
                                loocean(1:kproma))

     pmassf_cs(1:kproma) = MERGE( &
                                ss2_mon*ztmp1(1:kproma), &
                                pmassf_cs(1:kproma), &
                                loocean(1:kproma))
    !<<SF #458 (replacing WHERE statements)

    !--- 3) Calculate the emitted number-flux: ----------------------------------------
    pnumf_as(1:kproma) = pmassf_as(1:kproma) / zssmassa
    pnumf_cs(1:kproma) = pmassf_cs(1:kproma) / zssmassc

  END SUBROUTINE seasalt_emissions_monahan


  SUBROUTINE seasalt_emissions_lsce(kproma, kbdim, krow, pmassf_as, pmassf_cs, pnumf_as, pnumf_cs)
    !     -----------------------------------------------------------------------
    !     
    !   Author:
    !   -------
    !   Michael Schulz 
    !   Laboratoire des Sciences du Climat et de l'Environnement / Saclay
    !   10.1.2002
    !
    !   Modifications:
    !   --------------
    !   Philip Stier, MPI-MET  (Adaption to the ECHAM/HAM structure)       2002
    !   Michael Schulz, LSCE   (Modified source coefficients)        02/08/2002
    !   Martin Schultz, FZ Juelich (2010-02-11): adapted interface for ECHAM6-HAMMOZ
    !
    !   Purpose:
    !   --------
    !   Describe source flux of sea salt aerosol mass and number flux 
    !   as a function of wind speed
    !   for two aerosol modes: coarse soluble and accumulation soluble
    !
    !   Interface:
    !   ----------
    !      input 
    !       delt           time step duration       [sec]
    !       zspeed         wind speed at 10 m       [m/s]
    !       slf            sea land fraction        [ % 0.0 (sea) - 1. (land)]
    !       alake          lake fraction            [%]
    !
    !      output 
    !       pxtems         emission flux of sea salt mass and number [kg m-2 s-1] / [# m-2 s-1]
    !
    !
    !   Method:
    !   -------
    !     Tabulated mass and number fluxes following Monahan 86 and Smith&Harrison 98
    !     for two aerosol modes; interpolated according to actual wind speed;
    !     sea salt mass and number fluxes are fitted to two lognormal distributions
    !     following work published by Guelle et al. (2001) JGR 106, pp. 27509-27524
    !
    !

    USE mo_memory_g3b,   ONLY: slm , seaice, slf, alake
    USE mo_vphysc,       ONLY: vphysc

    IMPLICIT NONE
 
    ! arguments
    INTEGER, INTENT(in)    :: kproma          ! geographic block number of locations
    INTEGER, INTENT(in)    :: kbdim           ! geographic block maximum number of locations 
    INTEGER, INTENT(in)    :: krow            ! geographic block number 
    REAL(dp), INTENT(out)  :: pmassf_as(kbdim)! seasalt emission mass flux in soluble accumulation mode
    REAL(dp), INTENT(out)  :: pmassf_cs(kbdim)! ...                        in soluble coarse mode
    REAL(dp), INTENT(out)  :: pnumf_as(kbdim) ! sesalt particle number flux in soluble accumulation mode
    REAL(dp), INTENT(out)  :: pnumf_cs(kbdim) ! ...                        in soluble coarse mode

    !--- Local Variables:
    LOGICAL                :: loocean(kbdim)            ! ocean mask
    LOGICAL                :: ll1(kbdim)                !SF #458 temporary logical
    REAL(dp)               :: zseaice(kbdim)            ! sea ice fraction
    INTEGER                :: wcl(kbdim)                ! windclass in 1 m/s steps
    REAL(dp)               :: dzspeed(kbdim)
    REAL(dp)               :: zseafrac(kbdim)           ! fraction of the gridcell covered by

! ###mgs### I fear that the following is resolution-dependent (?)

    !--- Precalculated mass fluxes per wind class [kg m-2 s-1]:

    REAL(dp) , PARAMETER :: mass1flux(0:40) = (/                           &
         0.000E+00_dp, 2.483E-15_dp, 2.591E-14_dp, 1.022E-13_dp, 2.707E-13_dp, 5.761E-13_dp, &
         1.068E-12_dp, 1.800E-12_dp, 2.829E-12_dp, 4.215E-12_dp, 6.023E-12_dp, 8.317E-12_dp, &
         1.117E-11_dp, 1.464E-11_dp, 1.882E-11_dp, 2.378E-11_dp, 2.959E-11_dp, 3.633E-11_dp, &
         4.409E-11_dp, 5.296E-11_dp, 6.301E-11_dp, 7.433E-11_dp, 8.693E-11_dp, 1.012E-10_dp, &
         1.168E-10_dp, 1.342E-10_dp, 1.532E-10_dp, 1.741E-10_dp, 1.970E-10_dp, 2.219E-10_dp, &
         2.489E-10_dp, 2.781E-10_dp, 3.097E-10_dp, 3.437E-10_dp, 3.803E-10_dp, 4.195E-10_dp, &
         4.616E-10_dp, 5.065E-10_dp, 5.544E-10_dp, 6.054E-10_dp, 6.711E-10_dp             /)

    REAL(dp) , PARAMETER :: mass2flux(0:40) = (/                               &
         0.000E+00_dp, 2.319E-13_dp, 2.411E-12_dp, 9.481E-12_dp, 2.505E-11_dp, 5.321E-11_dp, &
         9.850E-11_dp, 1.658E-10_dp, 2.602E-10_dp, 3.874E-10_dp, 5.529E-10_dp, 7.628E-10_dp, &
         1.023E-09_dp, 1.341E-09_dp, 1.722E-09_dp, 2.175E-09_dp, 2.704E-09_dp, 3.319E-09_dp, &
         4.026E-09_dp, 4.832E-09_dp, 5.746E-09_dp, 6.776E-09_dp, 7.925E-09_dp, 9.214E-09_dp, &
         1.064E-08_dp, 1.221E-08_dp, 1.394E-08_dp, 1.584E-08_dp, 1.791E-08_dp, 2.016E-08_dp, &
         2.261E-08_dp, 2.526E-08_dp, 2.812E-08_dp, 3.120E-08_dp, 3.451E-08_dp, 3.806E-08_dp, &
         4.186E-08_dp, 4.592E-08_dp, 5.025E-08_dp, 5.486E-08_dp, 6.014E-08_dp             /)

    !--- Precalculated number fluxes per wind class [m-2 s-1]:

    REAL(dp) , PARAMETER :: numb1flux(0:40) = (/                               &
         0.000E+00_dp, 3.004E+01_dp, 3.245E+02_dp, 1.306E+03_dp, 3.505E+03_dp, 7.542E+03_dp, &
         1.410E+04_dp, 2.394E+04_dp, 3.787E+04_dp, 5.674E+04_dp, 8.147E+04_dp, 1.130E+05_dp, &
         1.523E+05_dp, 2.005E+05_dp, 2.586E+05_dp, 3.278E+05_dp, 4.091E+05_dp, 5.037E+05_dp, &
         6.129E+05_dp, 7.379E+05_dp, 8.800E+05_dp, 1.041E+06_dp, 1.220E+06_dp, 1.422E+06_dp, &
         1.646E+06_dp, 1.893E+06_dp, 2.166E+06_dp, 2.466E+06_dp, 2.794E+06_dp, 3.152E+06_dp, &
         3.541E+06_dp, 3.962E+06_dp, 4.419E+06_dp, 4.911E+06_dp, 5.441E+06_dp, 6.011E+06_dp, &
         6.621E+06_dp, 7.274E+06_dp, 7.972E+06_dp, 8.716E+06_dp, 8.801E+06_dp             /)

    REAL(dp) , PARAMETER :: numb2flux(0:40) = (/                               &
         0.000E+00_dp, 1.934E+01_dp, 2.068E+02_dp, 8.271E+02_dp, 2.211E+03_dp, 4.741E+03_dp, &
         8.841E+03_dp, 1.497E+04_dp, 2.363E+04_dp, 3.534E+04_dp, 5.066E+04_dp, 7.017E+04_dp, &
         9.447E+04_dp, 1.242E+05_dp, 1.600E+05_dp, 2.025E+05_dp, 2.525E+05_dp, 3.106E+05_dp, &
         3.776E+05_dp, 4.542E+05_dp, 5.413E+05_dp, 6.395E+05_dp, 7.501E+05_dp, 8.726E+05_dp, &
         1.009E+06_dp, 1.160E+06_dp, 1.327E+06_dp, 1.509E+06_dp, 1.709E+06_dp, 1.927E+06_dp, &
         2.163E+06_dp, 2.420E+06_dp, 2.697E+06_dp, 2.996E+06_dp, 3.318E+06_dp, 3.664E+06_dp, &
         4.034E+06_dp, 4.430E+06_dp, 4.852E+06_dp, 5.303E+06_dp, 5.740E+06_dp             /)

    !--- Precalculated mass flux gradient for each wind class for interpolation [kg m-3]:
    !    (dm/dv(i) where m = m(i) + dm/dv(i) * dv(i) )

    REAL(dp) , PARAMETER :: dmass1flux(0:40) = (/                              &
         2.483E-15_dp, 2.343E-14_dp, 7.630E-14_dp, 1.684E-13_dp, 3.054E-13_dp, 4.919E-13_dp, &
         7.319E-13_dp, 1.029E-12_dp, 1.386E-12_dp, 1.807E-12_dp, 2.294E-12_dp, 2.850E-12_dp, &
         3.477E-12_dp, 4.174E-12_dp, 4.960E-12_dp, 5.810E-12_dp, 6.745E-12_dp, 7.762E-12_dp, &
         8.863E-12_dp, 1.005E-11_dp, 1.132E-11_dp, 1.260E-11_dp, 1.423E-11_dp, 1.569E-11_dp, &
         1.733E-11_dp, 1.907E-11_dp, 2.090E-11_dp, 2.284E-11_dp, 2.487E-11_dp, 2.700E-11_dp, &
         2.924E-11_dp, 3.158E-11_dp, 3.403E-11_dp, 3.659E-11_dp, 3.925E-11_dp, 4.202E-11_dp, &
         4.490E-11_dp, 4.790E-11_dp, 5.100E-11_dp, 6.578E-11_dp, 6.578E-11_dp             /)

    REAL(dp) , PARAMETER :: dmass2flux(0:40) = (/                              &
         2.319E-13_dp, 2.179E-12_dp, 7.070E-12_dp, 1.557E-11_dp, 2.817E-11_dp, 4.528E-11_dp, &
         6.727E-11_dp, 9.446E-11_dp, 1.271E-10_dp, 1.655E-10_dp, 2.100E-10_dp, 2.606E-10_dp, &
         3.177E-10_dp, 3.812E-10_dp, 4.522E-10_dp, 5.297E-10_dp, 6.145E-10_dp, 7.068E-10_dp, &
         8.066E-10_dp, 9.142E-10_dp, 1.030E-09_dp, 1.149E-09_dp, 1.289E-09_dp, 1.425E-09_dp, &
         1.573E-09_dp, 1.730E-09_dp, 1.896E-09_dp, 2.070E-09_dp, 2.254E-09_dp, 2.447E-09_dp, &
         2.649E-09_dp, 2.860E-09_dp, 3.080E-09_dp, 3.311E-09_dp, 3.550E-09_dp, 3.800E-09_dp, &
         4.060E-09_dp, 4.329E-09_dp, 4.609E-09_dp, 5.279E-09_dp, 5.279E-09_dp             /)

    !--- Precalculated number flux gradient for each wind class for interpolation [m-3]:
    !    (dn/dv(i) where n = n(i) + dn/dv(i) * dv(i) )

    REAL(dp) , PARAMETER :: dnumb1flux(0:40) = (/                              &
         3.004E+01_dp, 2.945E+02_dp, 9.811E+02_dp, 2.200E+03_dp, 4.036E+03_dp, 6.562E+03_dp, &
         9.839E+03_dp, 1.393E+04_dp, 1.887E+04_dp, 2.473E+04_dp, 3.153E+04_dp, 3.935E+04_dp, &
         4.818E+04_dp, 5.808E+04_dp, 6.914E+04_dp, 8.130E+04_dp, 9.465E+04_dp, 1.092E+05_dp, &
         1.250E+05_dp, 1.421E+05_dp, 1.605E+05_dp, 1.798E+05_dp, 2.017E+05_dp, 2.237E+05_dp, &
         2.476E+05_dp, 2.729E+05_dp, 2.997E+05_dp, 3.280E+05_dp, 3.577E+05_dp, 3.890E+05_dp, &
         4.219E+05_dp, 4.564E+05_dp, 4.923E+05_dp, 5.301E+05_dp, 5.695E+05_dp, 6.105E+05_dp, &
         6.531E+05_dp, 6.976E+05_dp, 7.437E+05_dp, 8.550E+04_dp, 8.550E+04_dp             /)

    REAL(dp) , PARAMETER :: dnumb2flux(0:40) = (/                              &
         1.934E+01_dp, 1.875E+02_dp, 6.203E+02_dp, 1.384E+03_dp, 2.530E+03_dp, 4.100E+03_dp, &
         6.132E+03_dp, 8.659E+03_dp, 1.171E+04_dp, 1.532E+04_dp, 1.951E+04_dp, 2.430E+04_dp, &
         2.972E+04_dp, 3.582E+04_dp, 4.251E+04_dp, 4.997E+04_dp, 5.812E+04_dp, 6.700E+04_dp, &
         7.663E+04_dp, 8.702E+04_dp, 9.820E+04_dp, 1.106E+05_dp, 1.225E+05_dp, 1.366E+05_dp, &
         1.511E+05_dp, 1.664E+05_dp, 1.826E+05_dp, 1.997E+05_dp, 2.177E+05_dp, 2.366E+05_dp, &
         2.565E+05_dp, 2.773E+05_dp, 2.990E+05_dp, 3.218E+05_dp, 3.455E+05_dp, 3.702E+05_dp, &
         3.959E+05_dp, 4.226E+05_dp, 4.504E+05_dp, 4.376E+05_dp, 4.376E+05_dp             /)


    !--- 0) Initialisation
    loocean(1:kproma)  = ( slm(1:kproma,krow) < 1.e-2_dp )
    zseaice(1:kproma) = seaice(1:kproma, krow)
    pmassf_as(:) = 0._dp
    pmassf_cs(:) = 0._dp
    pnumf_as(:)  = 0._dp
    pnumf_cs(:)  = 0._dp

    !--- Windclass as a function of wind speed (limited to 40 m/s):

    wcl(1:kproma)=INT(vphysc%velo10m(1:kproma,krow))
    wcl(1:kproma)=MAX(0,MIN(wcl(1:kproma),40))

    !--- Calculate dv for interpolation in the windclas interval:

    dzspeed(1:kproma)=vphysc%velo10m(1:kproma,krow)-REAL(wcl(1:kproma),dp)

    !--- Calculate fraction of the gridcell of non ice-covered water:

    zseafrac(1:kproma)=(1._dp-slf(1:kproma,krow)-alake(1:kproma,krow))*(1._dp-zseaice(1:kproma))
    zseafrac(1:kproma)=MAX(0._dp,MIN(zseafrac(1:kproma),1._dp))
    !corinna: avoid emission from land !SF -> + #458 (replacing WHERE statements)
    ll1(1:kproma) = (slf(1:kproma,krow)+alake(1:kproma,krow) .GT. 0.95_dp)
    zseafrac(1:kproma) = MERGE(0._dp, zseafrac(1:kproma), ll1(1:kproma))
    !end corinna

    !--- Accumulation mode soluble:
    !    Mass flux [kg m-2 s-1]:
    pmassf_as(1:kproma) = (mass1flux(wcl(1:kproma))+dmass1flux(wcl(1:kproma))   &
                          * dzspeed(1:kproma))*zseafrac(1:kproma)
    !    Number flux [m-2 s-1]:
    pnumf_as(1:kproma)  = (numb1flux(wcl(1:kproma))+dnumb1flux(wcl(1:kproma))   &
                          * dzspeed(1:kproma))*zseafrac(1:kproma)

    !--- Coarse Mode Soluble:
    !    Mass flux [kg m-2 s-1]:
    pmassf_cs(1:kproma) = (mass2flux(wcl(1:kproma))+dmass2flux(wcl(1:kproma))   &
                          * dzspeed(1:kproma))*zseafrac(1:kproma)
    !    Number flux [m-2 s-1]:
    pnumf_cs(1:kproma)  = (numb2flux(wcl(1:kproma))+dnumb2flux(wcl(1:kproma))   &
                          * dzspeed(1:kproma))*zseafrac(1:kproma)

  END SUBROUTINE seasalt_emissions_lsce

  SUBROUTINE seasalt_emissions_mh(kproma, kbdim, krow, pmassf_as, pmassf_cs, pnumf_as, pnumf_cs)


    !  
    ! Description:
    ! ------------
    ! Calculates the emitted sea salt flux from the 10m wind speed following
    ! Monahan et al., 1986 
    !
    !   method: whitecap/lab
    !   size: 0.3um < r80 < 20um
    !   wind speed: N.A.
    !   SST: ~20 degree  
    ! 
    ! currently Aitken mode particles are negelected. 
    !  
    !
    ! Authors:
    ! ------------ 
    ! Kai Zhang, MPI-Met, 2009 
    !  
    !
    ! References:
    ! ------------
    ! 1. Monahan, E. C., D. E. Spiel, and K. L. Davidson, 
    !    A model of marine aerosol generation via whitecaps and wave disruption, 
    !    in Oceanic Whitecaps, edited by E. C. Monahan and G. MacNiochaill, 
    !    pp. 167–193, D. Reidel, Norwell, Mass., 1986. (M86) 
    ! 

    USE mo_kind,         ONLY: dp
    USE mo_math_constants, ONLY: pi
    USE mo_species,      ONLY: speclist
    USE mo_ham_species,  ONLY: id_ss
    USE mo_memory_g3b,   ONLY: slf, alake, seaice
    USE mo_vphysc,       ONLY: vphysc

    IMPLICIT NONE

    !--- Parameters:
    !    -

    !--- I/O:

    INTEGER, INTENT(in)    :: kproma
    INTEGER, INTENT(in)    :: kbdim
    INTEGER, INTENT(in)    :: krow
    REAL(dp),INTENT(out)   :: pmassf_as(kbdim)    ! mass flux of ss acc particles
    REAL(dp),INTENT(out)   :: pmassf_cs(kbdim)    ! mass flux of ss coa particles
    REAL(dp),INTENT(out)   :: pnumf_as(kbdim)    ! number flux of ss acc particles
    REAL(dp),INTENT(out)   :: pnumf_cs(kbdim)    ! number flux of ss coa particles

    !--- Local:

    REAL(dp):: zseafrac(kbdim)         ! fraction of the gridcell covered by
                                       ! non-iced sea water [0.-1.]
    REAL(dp):: zmassf_ks(kbdim)      ! mass   flux of ss ait particles  (not currently supported)
    REAL(dp):: znumf_ks(kbdim)      ! number flux of ss ait particles  (not currently supported)

    REAL(dp):: fi(kbdim,nbin) 
    REAL(dp):: p1,p2,p3,dr
    REAL(dp):: zav !particle volumn  

    INTEGER :: m 

    !----------------------------------------------------------------------------------------------


    ! initialize number flux for each bin (#/m2/s) 

    fi = 0._dp

    !--- Calculate fraction of the gridcell of non ice-covered water:

    zseafrac(1:kproma) = (1._dp-slf(1:kproma,krow)-alake(1:kproma,krow))*(1._dp-seaice(1:kproma,krow))
    zseafrac(1:kproma) = MAX(0._dp,MIN(zseafrac(1:kproma),1._dp))

    !>>SF #458 (replacing where statements)
    zseafrac(1:kproma) = MERGE( &
                               0._dp, &
                               zseafrac(1:kproma), &
                               (slf(1:kproma,krow) > 0.5_dp))
    !<<SF #458 (replacing where statements)

    zmassf_ks(1:kproma) = 0._dp 
    pmassf_as(1:kproma) = 0._dp 
    pmassf_cs(1:kproma) = 0._dp 

    znumf_ks(1:kproma) = 0._dp
    pnumf_as(1:kproma) = 0._dp
    pnumf_cs(1:kproma) = 0._dp

    ! formula (5) in M03 

    DO m = 2,nbin

       dr = rm(m) - rm(m-1) 

       IF (dmt(m).GT.dmta .and. dmt(m).le.dmtd) THEN 

         p1 = 1.373_dp*rm(m)**(-3) 
         p2 = 1._dp + 0.057_dp*rm(m)**1.05_dp
         p3 = 10**(1.19_dp*EXP(-bmn(m)**2)) 
 
         fi(1:kproma,m) = p1*p2*p3*dr*vphysc%velo10m(1:kproma,krow)**ppww

       END IF

       zav = speclist(id_ss)%density*4._dp/3._dp*pi*rd(m)**3

       ! mass flux [kg m-2 s-1] and number flux [m-2 s-1]

       IF (dmt(m).GT.dbeg(1) .AND. dmt(m).LE.dend(1) ) THEN
          znumf_ks(1:kproma) = znumf_ks(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)
          zmassf_ks(1:kproma) = zmassf_ks(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)*zav 
       END IF

       IF (dmt(m).GT.dbeg(2) .AND. dmt(m).LE.dend(2) ) THEN
          pnumf_as(1:kproma) = pnumf_as(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)
          pmassf_as(1:kproma) = pmassf_as(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)*zav
       END IF

       IF (dmt(m).GT.dbeg(3) .AND. dmt(m).LE.dend(3) ) THEN
          pnumf_cs(1:kproma) = pnumf_cs(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)
          pmassf_cs(1:kproma) = pmassf_cs(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)*zav
       END IF
 
    END DO
 

  END SUBROUTINE seasalt_emissions_mh


  SUBROUTINE seasalt_emissions_guelle(kproma, kbdim, krow, pmassf_as, pmassf_cs, pnumf_as, pnumf_cs)


    !  
    ! Description:
    ! ------------
    ! Calculates the emitted sea salt flux from the 10m wind speed following
    ! Guelle et al., 2001 
    !
    !   method: M86+S98  
    !   size: 0.3um < r80 < 20um
    !   wind speed: N.A.
    !   SST: ~20 degree  
    ! 
    ! currently Aitken mode particles are negelected. 
    !  
    ! Authors:
    ! ------------ 
    ! Kai Zhang, MPI-Met, 2009  
    !
    ! References:
    ! ------------
    ! 1. Monahan, E. C., D. E. Spiel, and K. L. Davidson, 
    !    A model of marine aerosol generation via whitecaps and wave disruption, 
    !    in Oceanic Whitecaps, edited by E. C. Monahan and G. MacNiochaill, 
    !    pp. 167–193, D. Reidel, Norwell, Mass., 1986. (M86) 
    ! 
    ! 2. M.H. Smith and N.M. Harrison, 
    !    The sea spray generation function, 
    !    Journal of Aerosol Science 29 (1998) (Suppl. 1), pp. S189–S190. (SM98) 
    ! 
    ! 3. Guelle, W., M. Schulz, Y. Balkanski, and F. Dentener (2001), 
    !    Influence of the source formulation on modeling the atmospheric global distribution of sea salt aerosol, 
    !    J. Geophys. Res., 106(D21), 27,509–27,524. (G01) 
    ! 

    USE mo_kind,         ONLY: dp
    USE mo_math_constants, ONLY: pi
    USE mo_species,      ONLY: speclist
    USE mo_ham_species,  ONLY: id_ss
    USE mo_memory_g3b,   ONLY: slf, alake, seaice
    USE mo_vphysc,       ONLY: vphysc

    IMPLICIT NONE


    !--- Parameters:
    REAL(dp),PARAMETER :: pr01 = 3.00_dp      !r01 in SM98 
    REAL(dp),PARAMETER :: pr02 = 30.0_dp      !r02 in SM98 
    REAL(dp),PARAMETER :: pf01 = 1.5_dp       !f01 in SM98 
    REAL(dp),PARAMETER :: pf02 = 1.0_dp       !f02 in SM98 
    REAL(dp),PARAMETER :: pa01 = 0.2_dp       !a01 in SM98 
    REAL(dp),PARAMETER :: pa02 = 6.8E-03_dp   !a02 in SM98 
    REAL(dp),PARAMETER :: pp01 = 3.5_dp       !power of wind 01 in SM98 
    REAL(dp),PARAMETER :: pp02 = 3.0_dp       !power of wind 02 in SM98   

    !--- I/O:

    INTEGER, INTENT(in)    :: kproma
    INTEGER, INTENT(in)    :: kbdim
    INTEGER, INTENT(in)    :: krow
    REAL(dp),INTENT(out)   :: pmassf_as(kbdim)    ! mass flux of ss acc particles
    REAL(dp),INTENT(out)   :: pmassf_cs(kbdim)    ! mass flux of ss coa particles
    REAL(dp),INTENT(out)   :: pnumf_as(kbdim)    ! number flux of ss acc particles
    REAL(dp),INTENT(out)   :: pnumf_cs(kbdim)    ! number flux of ss coa particles

    !--- Local:

    REAL(dp):: zseafrac(kbdim)         ! fraction of the gridcell covered by
                                       ! non-iced sea water [0.-1.]
    REAL(dp):: zmassf_ks(kbdim)      ! mass   flux of ss acc particles  (currently not supported)
    REAL(dp):: znumf_ks(kbdim)      ! number flux of ss acc particles  (currently not supported)

    REAL(dp):: fi(kbdim,nbin)
    REAL(dp):: p1,p2,p3,dr
    REAL(dp):: zav !particle volumn  

    INTEGER :: m

 
    ! initialize number flux for each bin (#/m2/s) 

    fi = 0._dp

    !--- Calculate fraction of the gridcell of non ice-covered water:

    zseafrac(1:kproma) = (1._dp-slf(1:kproma,krow)-alake(1:kproma,krow))*(1._dp-seaice(1:kproma,krow))
    zseafrac(1:kproma) = MAX(0._dp,MIN(zseafrac(1:kproma),1._dp))

    !>>SF #458 (replacing where statements)
    zseafrac(1:kproma) = MERGE( &
                               0._dp, &
                               zseafrac(1:kproma), &
                               (slf(1:kproma,krow) > 0.5_dp))
    !<<SF #458 (replacing where statements)

    zmassf_ks(1:kproma) = 0._dp 
    pmassf_as(1:kproma) = 0._dp 
    pmassf_cs(1:kproma) = 0._dp 

    znumf_ks(1:kproma) = 0._dp
    pnumf_as(1:kproma) = 0._dp
    pnumf_cs(1:kproma) = 0._dp

    ! formula (5) in M03 

    DO m = 2,nbin

       !dLOGdrm  

       dr = rm(m) - rm(m-1) 

       IF (dmt(m).GT.dmta .and. dmt(m).le.dmtb_guelle) THEN 

         p1 = 1.373_dp*rm(m)**(-3) 
         p2 = 1._dp + 0.057_dp*rm(m)**1.05_dp
         p3 = 10**(1.19_dp*EXP(-bmn(m)**2)) 

         fi(1:kproma,m) = p1*p2*p3*dr*vphysc%velo10m(1:kproma,krow)**ppww

       ELSEIF (dmt(m).GT.dmtb_guelle .and. dmt(m).le.dmtd) THEN 

         !pf01 = 1.5_dp       !f01 in SM98 
         !pf02 = 1.0_dp       !f02 in SM98 
         !pr01 = 3.0_dp       !r01 in SM98 
         !pr02 = 30.0_dp      !r02 in SM98 
         !pa01 = 0.2_dp       !a01 in SM98 
         !pa02 = 6.8E-03_dp   !a02 in SM98 
         !pp01 = 3.5_dp       !w01 in SM98 
         !pp02 = 3.0_dp       !w02 in SM98   

         p1 = EXP( -pf01* LOG( rm(m)/pr01 ) ) 
         p2 = EXP( -pf02* LOG( rm(m)/pr02 ) ) 

         fi(1:kproma,m) = pa01*p1*dr*vphysc%velo10m(1:kproma,krow)**pp01 + & 
                          pa02*p2*dr*vphysc%velo10m(1:kproma,krow)**pp02

       END IF

       zav = speclist(id_ss)%density*4._dp/3._dp*pi*rd(m)**3

       ! mass flux [kg m-2 s-1] and number flux [m-2 s-1]

       IF (dmt(m).GT.dbeg(1) .AND. dmt(m).LE.dend(1) ) THEN
          znumf_ks(1:kproma) = znumf_ks(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)
          zmassf_ks(1:kproma) = zmassf_ks(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)*zav 
       END IF

       IF (dmt(m).GT.dbeg(2) .AND. dmt(m).LE.dend(2) ) THEN
          pnumf_as(1:kproma) = pnumf_as(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)
          pmassf_as(1:kproma) = pmassf_as(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)*zav
       END IF

       IF (dmt(m).GT.dbeg(3) .AND. dmt(m).LE.dend(3) ) THEN
          pnumf_cs(1:kproma) = pnumf_cs(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)
          pmassf_cs(1:kproma) = pmassf_cs(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)*zav
       END IF
 
    END DO
 
  END SUBROUTINE seasalt_emissions_guelle


  SUBROUTINE seasalt_emissions_gong(kproma, kbdim, krow, pmassf_as, pmassf_cs, pnumf_as, pnumf_cs)

    !  
    ! Description:
    ! ------------
    ! Calculates the emitted sea salt flux from the 10m wind speed following
    ! Gong, 2003 
    !
    !   method: M86/lab
    !   size: 0.07um < r80 < 20um
    !   wind speed: N.A.
    !   SST: ~20 degree
    ! 
    ! currently Aitken mode particles are negelected. 
    !   
    ! Authors:
    ! ------------ 
    ! Kai Zhang, MPI-Met, 2009 
    !   
    ! References:
    ! ------------
    ! 1. Monahan, E. C., D. E. Spiel, and K. L. Davidson, 
    !    A model of marine aerosol generation via whitecaps and wave disruption, 
    !    in Oceanic Whitecaps, edited by E. C. Monahan and G. MacNiochaill, 
    !    pp. 167–193, D. Reidel, Norwell, Mass., 1986. (M86)  
    ! 
    ! 2. S.L. Gong, 
    !    A parameterization of sea-salt aerosol source function for sub- and super-micron particles, 
    !    Global Biogeochemical Cycles 17 (2003) (4), p. 1097.   (G03) 
    ! 

    USE mo_kind,         ONLY: dp
    USE mo_math_constants, ONLY: pi
    USE mo_species,      ONLY: speclist
    USE mo_ham_species,  ONLY: id_ss
    USE mo_memory_g3b,   ONLY: slf, alake, seaice
    USE mo_vphysc,       ONLY: vphysc

    IMPLICIT NONE

    !--- Parameters:
    !    -


    !--- I/O:

    INTEGER, INTENT(in)    :: kproma               !kproma
    INTEGER, INTENT(in)    :: kbdim                !column  
    INTEGER, INTENT(in)    :: krow                 !chunk 
    REAL(dp),INTENT(out)   :: pmassf_as(kbdim)    ! mass flux of ss acc particles
    REAL(dp),INTENT(out)   :: pmassf_cs(kbdim)    ! mass flux of ss coa particles
    REAL(dp),INTENT(out)   :: pnumf_as(kbdim)    ! number flux of ss acc particles
    REAL(dp),INTENT(out)   :: pnumf_cs(kbdim)    ! number flux of ss coa particles

    !--- Local:

    REAL(dp):: zseafrac(kbdim)         ! fraction of the gridcell covered by
                                       ! non-iced sea water [0.-1.]
    REAL(dp):: zmassf_ks(kbdim)        ! mass   flux of ss ait particles  (currently not supported)
    REAL(dp):: znumf_ks(kbdim)         ! number flux of ss ait particles  (currently not supported)

    REAL(dp):: fi(kbdim,nbin) 
    REAL(dp):: p0,p1,p2,p3,dr
    REAL(dp):: zav !particle volumn  

    INTEGER :: m

    ! initialize number flux for each bin (#/m2/s) 

    fi = 0._dp

    ! calculate fraction of the gridcell of non ice-covered water

    zseafrac(1:kproma) = (1._dp-slf(1:kproma,krow)-alake(1:kproma,krow))*(1._dp-seaice(1:kproma,krow))
    zseafrac(1:kproma) = MAX(0._dp,MIN(zseafrac(1:kproma),1._dp))

    !>>SF #458 (replacing where statements)
    zseafrac(1:kproma) = MERGE( &
                               0._dp, &
                               zseafrac(1:kproma), &
                               (slf(1:kproma,krow) > 0.5_dp))
    !<<SF #458 (replacing where statements)

    ! initialization 

    zmassf_ks(1:kproma) = 0._dp 
    pmassf_as(1:kproma) = 0._dp 
    pmassf_cs(1:kproma) = 0._dp 

    znumf_ks(1:kproma) = 0._dp
    pnumf_as(1:kproma) = 0._dp
    pnumf_cs(1:kproma) = 0._dp

    ! loop over bins 

    DO m = 2,nbin

       !dLOGdrm 

       dr = rm(m) - rm(m-1) 

       IF (dmt(m).GT.dmta .and. dmt(m).le.dmtb_gong) THEN 

         p0 = 4.7_dp*(1._dp+30._dp*rm(m))**(-0.017_dp*rm(m)**(-1.44_dp)) 
         p1 = 1.373_dp*rm(m)**(-p0) 
         p2 = 1._dp + 0.057_dp*rm(m)**3.45_dp
         p3 = 10**(1.607_dp*EXP(-bmn(m)**2)) 
 
         fi(1:kproma,m) = p1*p2*p3*dr*vphysc%velo10m(1:kproma,krow)**ppww

       ELSEIF (dmt(m).GT.dmtb_gong .and. dmt(m).le.dmtd) THEN 

         p1 = 1.373_dp*rm(m)**(-3) 
         p2 = 1._dp + 0.057_dp*rm(m)**1.05_dp
         p3 = 10**(1.19_dp*EXP(-bmn(m)**2))  

         fi(1:kproma,m) = p1*p2*p3*dr*vphysc%velo10m(1:kproma,krow)**ppww

       END IF

       zav = speclist(id_ss)%density*4._dp/3._dp*pi*rd(m)**3

       ! mass flux (kg/m2/s) and number flux (#/m2/s) 

       IF (dmt(m).GT.dbeg(1) .AND. dmt(m).LE.dend(1) ) THEN
          znumf_ks(1:kproma) = znumf_ks(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)
          zmassf_ks(1:kproma) = zmassf_ks(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)*zav 
       END IF

       IF (dmt(m).GT.dbeg(2) .AND. dmt(m).LE.dend(2) ) THEN
          pnumf_as(1:kproma) = pnumf_as(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)
          pmassf_as(1:kproma) = pmassf_as(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)*zav
       END IF

       IF (dmt(m).GT.dbeg(3) .AND. dmt(m).LE.dend(3) ) THEN
          pnumf_cs(1:kproma) = pnumf_cs(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)
          pmassf_cs(1:kproma) = pmassf_cs(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)*zav
       END IF
 
    END DO
 
  END SUBROUTINE seasalt_emissions_gong


  SUBROUTINE seasalt_emissions_long(kproma, kbdim, krow, pmassf_as, pmassf_cs, pnumf_as, pnumf_cs)

    !  
    ! Description:
    ! ------------
    ! Calculates the emitted sea salt flux from the 10m wind speed following Long et al., 2011
    !
    !
    !   method: Based on Laboratory data (Keene et al., 2007, JGR) 
    !   size: 0.03um < d80 < 20um
    !   wind speed: N.A.
    !   SST: Not given, SST correction according to Sofiev et al 2011 
    !   Expansion to primary organics possible
    !   
    ! Authors:
    ! ------------ 
    ! I. Tegen, S. Barthel, 2016 (TROPOS)
    !   
    ! References:
    ! ------------
    ! Long M.S., Keene W.D., Kieber D.J., Erickson D.J., Maring H., A sea-state
    ! based source function for size- and composition-resolved marine aerosol
    ! production, Atmospheric Chemistry and Physics, 11, 1203-1216, 2011
    !
    ! Sofiev, M., Soares,J.,Prank, M., deLeeuw, G., Kukkonen, J., A
    ! regional-to-global model of emission and transport of sea salt particles
    ! in the atmosphere,JGR (doi:0148‐0227/11/2010JD014713), 2011


!--------------------------------------------------------------------------------------------


    USE mo_kind,         ONLY: dp
    USE mo_math_constants, ONLY: pi
    USE mo_species,      ONLY: speclist
    USE mo_ham_species,  ONLY: id_ss
    USE mo_memory_g3b,   ONLY: slf, alake, seaice
    USE mo_vphysc,       ONLY: vphysc
    USE mo_exception,     ONLY: message, message_text, finish

    IMPLICIT NONE

    !--- Parameters:
    !    -
    REAL(dp), PARAMETER :: dmtb_long = 0.551E-06_dp ! as dry diameter(m), limit wet diameter 1um at 80% RH, 
                                                    ! according to Long et al 2011 (1./1.814)
    REAL(dp), PARAMETER :: ppww_long=3.74_dp        !

    !--- I/O:

    INTEGER, INTENT(in)    :: kproma               ! kproma
    INTEGER, INTENT(in)    :: kbdim                ! column  
    INTEGER, INTENT(in)    :: krow                 ! chunk 
    REAL(dp),INTENT(out)   :: pmassf_as(kbdim)     ! mass flux of ss acc particles
    REAL(dp),INTENT(out)   :: pmassf_cs(kbdim)     ! mass flux of ss coa particles
    REAL(dp),INTENT(out)   :: pnumf_as(kbdim)      ! number flux of ss acc particles
    REAL(dp),INTENT(out)   :: pnumf_cs(kbdim)      ! number flux of ss coa particles

    !--- Local:

    REAL(dp):: zseafrac(kbdim)         ! fraction of the gridcell covered by
                                       ! non-iced sea water [0.-1.]
    REAL(dp):: zmassf_ks(kbdim)        ! mass   flux of ss ait particles  (currently not supported)
    REAL(dp):: znumf_ks(kbdim)         ! number flux of ss ait particles  (currently not supported)

    REAL(dp):: fi(kbdim,nbin) 
    REAL(dp):: p0,p1,p2,p11,p12,p13,p14,p21,p22,p23,p24,dr,logdp
    REAL(dp):: zav !particle volumn  
    REAL(dp):: SST_corr(1:kproma,krow),SST_corr_1,SST_corr_2,dmtum(nbin)
    REAL(dp):: SST_mask1(1:kproma,krow),SST_mask2(1:kproma,krow), &
               SST_mask3(1:kproma,krow),SST_corr_all(1:kproma,krow)

    INTEGER :: m

    ! initialize number flux for each bin (#/m2/s) 

    fi = 0._dp

    ! calculate fraction of the gridcell of non ice-covered water

    zseafrac(1:kproma) = (1._dp-slf(1:kproma,krow)-alake(1:kproma,krow))*(1._dp-seaice(1:kproma,krow))
    zseafrac(1:kproma) = MAX(0._dp,MIN(zseafrac(1:kproma),1._dp))

    !>>SF #458 (replacing where statements)
    zseafrac(1:kproma) = MERGE( &
                               0._dp, &
                               zseafrac(1:kproma), &
                               (slf(1:kproma,krow) > 0.5_dp))
    !<<SF #458 (replacing where statements)

    ! initialization 

    zmassf_ks(1:kproma) = 0._dp 
    pmassf_as(1:kproma) = 0._dp 
    pmassf_cs(1:kproma) = 0._dp 

    znumf_ks(1:kproma) = 0._dp
    pnumf_as(1:kproma) = 0._dp
    pnumf_cs(1:kproma) = 0._dp
    
    p0=2.E-8_dp
    p11=1.46E0_dp
    p12=1.33E0_dp
    p13=-1.82E0_dp
    p14=8.83E0_dp
    p21=-1.53E0_dp
    p22=-8.1E-2_dp
    p23=-4.26E-1_dp
    p24=8.84E0_dp

! loop over bins

    DO m = 2,nbin

       !dLOGdrm 

       dr = 2*(rm(m) - rm(m-1))
       logdp=log10(2*rm(m))-log10(2*rm(m-1))

       IF (dmt(m).GT.dmta .and. dmt(m).le.dmtb_long) THEN
         
         p1=p11*(Log10(2*rm(m)))**3+p12*(Log10(2*rm(m)))**2+p13*(Log10(2*rm(m)))+p14 
         fi(1:kproma,m) = 10.**p1*(p0*vphysc%velo10m(1:kproma,krow)**ppww_long)*logdp

       ELSEIF (dmt(m).GT.dmtb_long .and. dmt(m).le.dmtd) THEN
           
         p2=p21*(Log10(2*rm(m)))**3+p22*(Log10(2*rm(m)))**2+p23*(Log10(2*rm(m)))+p24 
         fi(1:kproma,m) = 10.**p2*(p0*vphysc%velo10m(1:kproma,krow)**ppww_long)*logdp

       END IF

       ! SST correction according to Sofiev et al 2011,
   
       dmtum(m)=dmt(m)*1.e06  ! dmt -> micrometers
   
       !SST_corr_1 = 0.092e0*dmtum(m)**(-0.96e0) ! valid for Tw=271.15K ; -2°C      ! T base (Long)  25 deg
       !SST_corr_2 = 0.15e0*dmtum(m)**(-0.88e0) ! valid for Tw=278.15K ; 5°C
   
       ! SST_corr_1 = 0.19e0*dmtum(m)**(-0.60e0) ! valid for Tw=271.15K ; -2°CC      T base  (Long) 15 deg
       ! SST_corr_2 = 0.31e0*dmtum(m)**(-0.56e0) ! valid for Tw=278.15K ; 5°C
   
       SST_corr_1 = 0.13e0*dmtum(m)**(-0.78e0) ! valid for Tw=271.15K ; -2°CC        T base (Long) 20 deg
       SST_corr_2 = 0.22e0*dmtum(m)**(-0.70e0) ! valid for Tw=278.15K ; 5°C
   
       SST_corr(1:kproma,krow) = (SST_corr_1*(278.15e0-vphysc%tsw(1:kproma,krow)) &
                                 +SST_corr_2*(vphysc%tsw(1:kproma,krow)-271.15e0))/7.e0
   
       SST_mask1(1:kproma,krow) = MERGE (SST_corr(1:kproma,krow), 0._dp, &
                                        (vphysc%tsw(1:kproma,krow) .LE. 278.15 ))
   
       ! SST_corr_1 = 0.15e0*dmtum(m)**(-0.88e0) ! valid for Tw=278.15K ; 5°C
       ! SST_corr_2 = 0.48e0*dmtum(m)**(-0.36e0) ! valid for Tw=288.15K ; 15°C
   
       ! SST_corr_1 = 0.31e0*dmtum(m)**(-0.56e0) ! valid for Tw=278.15K ; 5°C
       ! SST_corr_2 = 1.e0 ! valid for Tw=288.15K ; 15°
   
       SST_corr_1 = 0.22e0*dmtum(m)**(-0.70e0) ! valid for Tw=278.15K ; 5°C
       SST_corr_2 = 0.70e0*dmtum(m)**(-0.18e0)  ! valid for Tw=288.15K ; 15°C
   
   
       SST_corr(1:kproma,krow) = (SST_corr_1*(288.15e0-vphysc%tsw(1:kproma,krow)) &
                                 +SST_corr_2*(vphysc%tsw(1:kproma,krow)-278.15e0))/1.e1
   
       SST_mask2(1:kproma,krow) = MERGE(SST_corr(1:kproma,krow), 0._dp,         &
                                        (vphysc%tsw(1:kproma,krow) .GT. 278.15  &
                                   .AND. vphysc%tsw(1:kproma,krow) .LE. 288.15))
   
   
       !SST_corr_1 = 0.48e0*dmtum(m)**(-0.36e0) ! valid for Tw=288.15K ; 15°C
       !SST_corr_2 = 1.e0 ! valid for Tw=298.15K ; 25°C
   
       ! SST_corr_1 = 1.e0! valid for Tw=288.15K ; 15°C
       ! SST_corr_2 = 2.08e0*dmtum(m)**(0.36e0) ! valid for Tw=298.15K ; 25°C
   
       SST_corr_1 = 0.70e0*dmtum(m)**(-0.18e0)  ! valid for Tw=288.15K ; 15°C
       SST_corr_2 = 1.45e0*dmtum(m)**(0.18e0) ! valid for Tw=298.15K ; 25°C
   
       SST_corr(1:kproma,krow) = (SST_corr_1*(298.15e0-vphysc%tsw(1:kproma,krow)) &
                                 +SST_corr_2*(vphysc%tsw(1:kproma,krow)-288.15e0))/1.e1
   
       SST_mask3(1:kproma,krow) = MERGE(SST_corr(1:kproma,krow), 0._dp, &
                                       (vphysc%tsw(1:kproma,krow) .GT. 288.15))
   
       SST_corr_all(1:kproma,krow) = SST_mask1(1:kproma,krow) &
                                    +SST_mask2(1:kproma,krow) &
                                    +SST_mask3(1:kproma,krow)

       !SST_corr_all(1:kproma,krow) = MERGE(1._dp, SST_corr_all(1:kproma,krow), &
       !                                    (vphysc%tsw(1:kproma,krow) .GT. 298.15)) ! optional, limit T dependence to <25 Deg C

      ! >>>> OPTIONAL: To exclude SST dependence comment out next line 
      fi(1:kproma,m)=fi(1:kproma,m)*SST_corr_all(1:kproma,krow)

       zav = speclist(id_ss)%density*4._dp/3._dp*pi*rd(m)**3

       ! mass flux (kg/m2/s) and number flux (#/m2/s) 

       IF (dmt(m).GT.dbeg(1) .AND. dmt(m).LE.dend(1) ) THEN
          znumf_ks(1:kproma) = znumf_ks(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)
          zmassf_ks(1:kproma) = zmassf_ks(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)*zav 
       END IF

       IF (dmt(m).GT.dbeg(2) .AND. dmt(m).LE.dend(2) ) THEN
          pnumf_as(1:kproma) = pnumf_as(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)
          pmassf_as(1:kproma) = pmassf_as(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)*zav
       END IF

       IF (dmt(m).GT.dbeg(3) .AND. dmt(m).LE.dend(3) ) THEN
          pnumf_cs(1:kproma) = pnumf_cs(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)
          pmassf_cs(1:kproma) = pmassf_cs(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)*zav
       END IF
 
    END DO

  END SUBROUTINE seasalt_emissions_long


  SUBROUTINE seasalt_emissions_gong_SST(kproma, kbdim, krow, pmassf_as, pmassf_cs, pnumf_as, pnumf_cs)

    !  
    ! Description:
    ! ------------
    ! Calculates the emitted sea salt flux from the 10m wind speed following
    ! Gong, 2003 
    !
    !   method: M86/lab
    !   size: 0.07um < r80 < 20um
    !   wind speed: N.A.
    !   SST: Explicit SST dependence, according to Sofiev et al 2011
    ! 
    ! currently Aitken mode particles are negelected. 
    !   
    ! Authors:
    ! ------------ 
    ! Kai Zhang, MPI-Met, 2009, modified by I. Tegen, 2016
    !   
    ! References:
    ! ------------
    ! 1. Monahan, E. C., D. E. Spiel, and K. L. Davidson, 
    !    A model of marine aerosol generation via whitecaps and wave disruption, 
    !    in Oceanic Whitecaps, edited by E. C. Monahan and G. MacNiochaill, 
    !    pp. 167–193, D. Reidel, Norwell, Mass., 1986. (M86)  
    ! 
    ! 2. S.L. Gong, 
    !    A parameterization of sea-salt aerosol source function for sub- and super-micron particles, 
    !    Global Biogeochemical Cycles 17 (2003) (4), p. 1097.   (G03) 
    !
    ! 3. Sofiev, M., Soares,J.,Prank, M., deLeeuw, G., Kukkonen, J., A
    ! regional-to-global model of emission and transport of sea salt particles
    ! in the atmosphere, JGR (doi:0148‐0227/11/2010JD014713)

    !

    USE mo_kind,         ONLY: dp
    USE mo_math_constants, ONLY: pi
    USE mo_species,      ONLY: speclist
    USE mo_ham_species,  ONLY: id_ss
    USE mo_memory_g3b,   ONLY: slf, alake, seaice
    USE mo_vphysc,       ONLY: vphysc

    IMPLICIT NONE

    !--- Parameters:
    !    -


    !--- I/O:

    INTEGER, INTENT(in)    :: kproma               !kproma
    INTEGER, INTENT(in)    :: kbdim                !column  
    INTEGER, INTENT(in)    :: krow                 !chunk 
    REAL(dp),INTENT(out)   :: pmassf_as(kbdim)    ! mass flux of ss acc particles
    REAL(dp),INTENT(out)   :: pmassf_cs(kbdim)    ! mass flux of ss coa particles
    REAL(dp),INTENT(out)   :: pnumf_as(kbdim)    ! number flux of ss acc particles
    REAL(dp),INTENT(out)   :: pnumf_cs(kbdim)    ! number flux of ss coa particles

    !--- Local:

    REAL(dp):: zseafrac(kbdim)         ! fraction of the gridcell covered by
                                       ! non-iced sea water [0.-1.]
    REAL(dp):: zmassf_ks(kbdim)        ! mass   flux of ss ait particles  (currently not supported)
    REAL(dp):: znumf_ks(kbdim)         ! number flux of ss ait particles  (currently not supported)

    REAL(dp):: fi(kbdim,nbin) 
    REAL(dp):: p0,p1,p2,p3,dr
    REAL(dp):: zav !particle volumn  
    REAL(dp):: SST_corr(1:kproma,krow),SST_corr_1,SST_corr_2,dmtum(nbin)
    REAL(dp):: SST_mask1(1:kproma,krow),SST_mask2(1:kproma,krow),SST_mask3(1:kproma,krow),SST_corr_all(1:kproma,krow)

    INTEGER :: m

    ! initialize number flux for each bin (#/m2/s) 

    fi = 0._dp

   ! calculate fraction of the gridcell of non ice-covered water

    zseafrac(1:kproma) = (1._dp-slf(1:kproma,krow)-alake(1:kproma,krow))*(1._dp-seaice(1:kproma,krow))
    zseafrac(1:kproma) = MAX(0._dp,MIN(zseafrac(1:kproma),1._dp))

    !>>SF #458 (replacing where statements)
    zseafrac(1:kproma) = MERGE( &
                               0._dp, &
                               zseafrac(1:kproma), &
                               (slf(1:kproma,krow) > 0.5_dp))
    !<<SF #458 (replacing where statements)

    ! initialization 

    zmassf_ks(1:kproma) = 0._dp 
    pmassf_as(1:kproma) = 0._dp 
    pmassf_cs(1:kproma) = 0._dp 

    znumf_ks(1:kproma) = 0._dp
    pnumf_as(1:kproma) = 0._dp
    pnumf_cs(1:kproma) = 0._dp

    ! loop over bins 

    DO m = 2,nbin

       !dLOGdrm 

       dr = rm(m) - rm(m-1) 

       IF (dmt(m).GT.dmta .and. dmt(m).le.dmtb_gong) THEN 

         p0 = 4.7_dp*(1._dp+30._dp*rm(m))**(-0.017_dp*rm(m)**(-1.44_dp)) 
         p1 = 1.373_dp*rm(m)**(-p0) 
         p2 = 1._dp + 0.057_dp*rm(m)**3.45_dp
         p3 = 10**(1.607_dp*EXP(-bmn(m)**2)) 
 
         fi(1:kproma,m) = p1*p2*p3*dr*vphysc%velo10m(1:kproma,krow)**ppww

       ELSEIF (dmt(m).GT.dmtb_gong .and. dmt(m).le.dmtd) THEN 

         p1 = 1.373_dp*rm(m)**(-3) 
         p2 = 1._dp + 0.057_dp*rm(m)**1.05_dp
         p3 = 10**(1.19_dp*EXP(-bmn(m)**2))  

         fi(1:kproma,m) = p1*p2*p3*dr*vphysc%velo10m(1:kproma,krow)**ppww

       END IF

       ! SST correction according to Sofiev et al 2011
   
       dmtum(m)=dmt(m)*1.e06  ! dmt -> micrometers
   
        SST_corr_1 = 0.092e0*dmtum(m)**(-0.96e0) ! valid for Tw=271.15K ; -2°C      ! T base (Long)  25 deg
        SST_corr_2 = 0.15e0*dmtum(m)**(-0.88e0) ! valid for Tw=278.15K ; 5°C
   
       ! SST_corr_1 = 0.19e0*dmtum(m)**(-0.60e0) ! valid for Tw=271.15K ; -2°CC      T base  (Long) 15 deg
       ! SST_corr_2 = 0.31e0*dmtum(m)**(-0.56e0) ! valid for Tw=278.15K ; 5°C
   
       !SST_corr_1 = 0.13e0*dmtum(m)**(-0.78e0) ! valid for Tw=271.15K ; -2°CC        T base (Long) 20 deg
       !SST_corr_2 = 0.22e0*dmtum(m)**(-0.70e0) ! valid for Tw=278.15K ; 5°C
   
   
   
       SST_corr(1:kproma,krow) = (SST_corr_1*(278.15e0-vphysc%tsw(1:kproma,krow)) &
                                 +SST_corr_2*(vphysc%tsw(1:kproma,krow)-271.15e0))/7.e0
   
       SST_mask1(1:kproma,krow) = MERGE (SST_corr(1:kproma,krow),0._dp,(vphysc%tsw(1:kproma,krow) .LE. 278.15 ))
   
        SST_corr_1 = 0.15e0*dmtum(m)**(-0.88e0) ! valid for Tw=278.15K ; 5°C
        SST_corr_2 = 0.48e0*dmtum(m)**(-0.36e0) ! valid for Tw=288.15K ; 15°C
   
       ! SST_corr_1 = 0.31e0*dmtum(m)**(-0.56e0) ! valid for Tw=278.15K ; 5°C
       ! SST_corr_2 = 1.e0 ! valid for Tw=288.15K ; 15°
   
       !SST_corr_1 = 0.22e0*dmtum(m)**(-0.70e0) ! valid for Tw=278.15K ; 5°C
       !SST_corr_2 = 0.70e0*dmtum(m)**(-0.18e0)  ! valid for Tw=288.15K ; 15°C
   
   
       SST_corr(1:kproma,krow) = (SST_corr_1*(288.15e0-vphysc%tsw(1:kproma,krow)) &
                               +  SST_corr_2*(vphysc%tsw(1:kproma,krow)-278.15e0))/1.e1
   
       SST_mask2(1:kproma,krow) = MERGE (SST_corr(1:kproma,krow), 0._dp,        &
                                         (vphysc%tsw(1:kproma,krow) .GT. 278.15 &
                                    .AND. vphysc%tsw(1:kproma,krow) .LE. 288.15))
   
   
        SST_corr_1 = 0.48e0*dmtum(m)**(-0.36e0) ! valid for Tw=288.15K ; 15°C
        SST_corr_2 = 1.e0 ! valid for Tw=298.15K ; 25°C
   
       ! SST_corr_1 = 1.e0! valid for Tw=288.15K ; 15°C
       ! SST_corr_2 = 2.08e0*dmtum(m)**(0.36e0) ! valid for Tw=298.15K ; 25°C
   
   
       !SST_corr_1 = 0.70e0*dmtum(m)**(-0.18e0)  ! valid for Tw=288.15K ; 15°C
       !SST_corr_2 = 1.45e0*dmtum(m)**(0.18e0) ! valid for Tw=298.15K ; 25°C
   
       SST_corr(1:kproma,krow) = (SST_corr_1*(298.15e0-vphysc%tsw(1:kproma,krow)) &
                                 +SST_corr_2*(vphysc%tsw(1:kproma,krow)-288.15e0))/1.e1
   
       SST_mask3(1:kproma,krow) = MERGE (SST_corr(1:kproma,krow), 0._dp, &
                                         (vphysc%tsw(1:kproma,krow) .GT. 288.15))
   
       SST_corr_all(1:kproma,krow)=SST_mask1(1:kproma,krow)+SST_mask2(1:kproma,krow)+SST_mask3(1:kproma,krow)
       !SST_corr_all(1:kproma,krow) = MERGE (1._dp, SST_corr_all(1:kproma,krow), &
            !(vphysc%tsw(1:kproma,krow) .GT. 298.15))  ! optional, limit T dependence to <25 deg C
   
       fi(1:kproma,m)=fi(1:kproma,m)*SST_corr_all(1:kproma,krow)

       zav = speclist(id_ss)%density*4._dp/3._dp*pi*rd(m)**3

       ! mass flux (kg/m2/s) and number flux (#/m2/s) 

       IF (dmt(m).GT.dbeg(1) .AND. dmt(m).LE.dend(1) ) THEN
          znumf_ks(1:kproma) = znumf_ks(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)
          zmassf_ks(1:kproma) = zmassf_ks(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)*zav 
       END IF

       IF (dmt(m).GT.dbeg(2) .AND. dmt(m).LE.dend(2) ) THEN
          pnumf_as(1:kproma) = pnumf_as(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)
          pmassf_as(1:kproma) = pmassf_as(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)*zav
       END IF

       IF (dmt(m).GT.dbeg(3) .AND. dmt(m).LE.dend(3) ) THEN
          pnumf_cs(1:kproma) = pnumf_cs(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)
          pmassf_cs(1:kproma) = pmassf_cs(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)*zav
       END IF
 
    END DO

    END SUBROUTINE seasalt_emissions_gong_SST


END MODULE mo_ham_m7_emi_seasalt
