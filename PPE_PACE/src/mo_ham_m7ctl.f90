!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_m7ctl.f90
!!
!! \brief
!! mo_ham_m7ctl contains parameters, switches and initialization routines for the m7 aerosol module.
!!
!! \author Elisabetta Vignatti (JRC/EI)
!! \author Philip Stier (MPI-Met)
!!
!! \responsible_coder
!! Martin G. Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# E. Vignati and J. Wilson (JRC/EI) - original code (2000)
!!   -# P. Stier (MPI-Met) (2001/2002)
!!   -# J. Kazil (MPI-Met) (2008)
!!   -# D. O'Donnell (MPI-Met) (2007-2007)
!!   -# M.G. Schultz (FZ Juelich) - new module struture (2009)
!! 
!! \limitations
!! Currently, there are two index lists for aerosol species: aero_idx in mo_species
!! and subm_aerospec in this module. I hope these are identical for the current model set-up 
!! in preparation for CMIP5. Later, one may wish to distinguish between the two: aero_idx
!! could contain additional aerosol species (e.g. from MOZART or climatologies), and this could
!! mess up the M7 code. If this can be generalized: fine. if not we should keep the two 
!! lists separate. mo_ham_rad (for example) works on aero_idx to be independent of M7 specifics.
!!
!! \details
!! None
!!
!! \bibliographic_references
!! None
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

MODULE mo_ham_m7ctl

  USE mo_kind,             ONLY: dp
  USE mo_math_constants,   ONLY: pi
  USE mo_physical_constants, ONLY: avo
  USE mo_species,          ONLY: nmaxspec
  USE mo_ham,              ONLY: HAM_M7, naeroclass, sigma_fine, sigma_coarse

  IMPLICIT NONE

  PRIVATE

  ! -- subroutines
  PUBLIC :: sethamM7, m7_initialize

  ! -- variables
  PUBLIC :: nwater, nsnucl, nonucl
  PUBLIC :: lnucl_stat
  PUBLIC :: inucs, iaits, iaccs, icoas, iaiti, iacci, icoai

  PUBLIC :: critn, cmin_aernl, cmin_aerml, cminvol, cminrad, cminrho, cdconv
  PUBLIC :: wna2so4, wh2so4, wnacl, wnahso4
  PUBLIC :: dna2so4, dh2so4, dnacl, dnahso4
  PUBLIC :: dbc, doc, ddust
  PUBLIC :: crh

  PUBLIC :: sigma, sigmaln
  PUBLIC :: crdiv, caccso4, gmb, wvb, dh2o
  PUBLIC :: bk, rerg, r_kcal

  PUBLIC :: so4_coating_threshold
  PUBLIC :: cmr2ras, cmr2mmr, cmedr2mmedr, cmr2ram, ram2cmr

  !--- 1) Define and pre-set switches for the processes of M7: -----------------------

  !--- Physical:

  INTEGER :: nwater     = 1         ! Aerosol water uptake scheme:
                                    !
                                    ! nwater = 0 Jacobson et al., JGR 1996
                                    !        = 1 Kappa-Koehler theory based approach (Petters and Kreidenweis, ACP 2007)
  
  INTEGER :: nsnucl     = 2         ! Choice of the H2SO4/H2O nucleation scheme:
                                    ! 
                                    ! nsnucl = 0 off
                                    !        = 1 Vehkamaeki et al., JGR 2002
                                    !        = 2 Kazil and Lovejoy, ACP 2007
  
  INTEGER :: nonucl     = 1         ! Choice of the organic nucleation scheme:
                                    ! 
                                    ! nonucl = 0 off
                                    !        = 1 Activation nucleation, Kulmala et al., ACP 2006
                                    !        = 2 Activation nucleation, Laakso et al., ACP 2004
  
  LOGICAL :: lnucl_stat = .FALSE.   ! Sample the cloud-free volume as function of T, RH, [H2SO4(g)],
                                    ! H2SO4 condensation sink, and ionization rate (memory intensive)

  !--- Mass index (in array aerml and ttn): 
  !
  !         Attention:
  !         The mass of sulfate compounds is always given in [molec. cm-3] 
  !         whilst the mass of other compounds is given in [ug cm-3].
  !
  !         Compounds:
  !
  !           so4 = sulphate
  !           bc  = black carbon
  !           oc  = organic carbon, 
  !           ss  = sea salt
  !           du  = dust 
  !
  !         Modes:
  !
  !           n   = nucleation mode
  !           k   = aitken mode 
  !           a   = accumulation mode
  !           c   = coarse mode
  !
  !         Type:
  !
  !           s   = soluble mode
  !           i   = insoluble mode
  
  ! M7 aerosol mode indices
  ! i*[n|k|a|c][i|s]   :   n = nucleation mode
  !                        k = aitken mode
  !                        a = accumulation mode
  !                        c = coarse mode
  !                        i = insoluble
  !                        s = soluble
  ! empty matrix entries are not populated.
  INTEGER, PUBLIC :: iso4ns, iso4ks, iso4as, iso4cs, &  
                             ibcks,  ibcas,  ibccs,  &
                             ibcki,                  &
                             iocks,  iocas,  ioccs,  &
                             iocki,                  &
                                     issas,  isscs,  &
                                     iduas,  iducs,  &
                                     iduai,  iduci

  !--- Number index (in array aernl):
  !

  INTEGER, PARAMETER ::                                      &
           inucs=1,  iaits=2,  iaccs=3,  icoas=4,  iaiti=5,  iacci=6,  icoai=7    
  ! MODE:           |         |         |         |         |
  !         nucl.   | aitk.   | acc.    | coar.   | aitk.   | acc.    | coar.   |
  !         soluble | soluble | soluble | soluble | insol.  | insol.  | insol.  |


  !--- 4) Definition of the modes of M7: ------------------------------------------------------

  !--- 4.1) Threshold radii between the different modes [cm]:
  !         Used for the repartititioning in m7_dconc.
  !         crdiv(jclass) is the lower bound and crdiv(jclass+1) is 
  !         the upper bound of the respective geometric mode.

  REAL(dp) :: crdiv(4)=(/ 0.0005E-4_dp, 0.005E-4_dp, 0.05E-4_dp, 0.5E-4_dp /)
  !                                   |            |           |          |     
  !                                   |            |           |          |
  !                              nucleation     aitken       accum     coarse mode

  !--- 4.2) Standard deviation for the modes:

  REAL(dp), PARAMETER :: sigma(naeroclass(HAM_M7))=(/ sigma_fine, sigma_fine, sigma_fine, sigma_coarse, &
                                                                  sigma_fine, sigma_fine, sigma_coarse /)

  !--- Natural logarithm of the standard deviation of each mode:
  !    Calulated in m7_initialize. 

  REAL(dp) :: sigmaln(naeroclass(HAM_M7))

  !--- 5) Conversion factors for lognormal particle size distributions: -------------
  !       Calulated in m7_initialize. 

  REAL(dp) :: cmr2ras(naeroclass(HAM_M7)) ! Conversion factor: count median radius to radius of average surface

  REAL(dp) :: cmr2mmr(naeroclass(HAM_M7)) ! Conversion factor: count median radius to mass mean radius

  REAL(dp) :: cmedr2mmedr(naeroclass(HAM_M7)) ! Conversion factor: count median radius to mass median radius

  REAL(dp) :: cmr2ram(naeroclass(HAM_M7)) ! Conversion factor: count median radius to radius of average mass

  REAL(dp) :: ram2cmr(naeroclass(HAM_M7)) ! Conversion factor: radius of average mass to count median radius


  !--- 6) Assumed thresholds for occurence of specific quantities: -------------

  REAL(dp), PARAMETER :: cmin_aerml = 1.E-15_dp , &   ! threshold for aerosol mass
                         cmin_aernl = 1.E-10_dp , &   ! threshold for aerosol number
!>>gf #240
                         cminvol    = 1.E-23_dp, & ! threshold for aerosol volume 
                         cminrad    = 1.E-8_dp,  & ! threshold for aerosol radius 
                         cminrho    = 1.E-10_dp    ! threshold for aerosol density

  !--- 6.1) Density conversion

  REAL(dp), PARAMETER :: cdconv = 1.E-3_dp     ! density conversion from kg m-3 to g cm-3

!<<gf
  
  !--- 7) Chemical constants: ----------------------------------------------------
  !
  !--- Accomodation coefficient of H2SO4 on aerosols:
  !    (reduced for insoluble modes)

  REAL(dp), PARAMETER :: caccso4(naeroclass(HAM_M7)) = (/ 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.3_dp, 0.3_dp, 0.3_dp /)

  !--- Critical relative humidity:

  REAL(dp), PARAMETER :: crh    = 0.45_dp           ! Assumed relative humidity for the 
                                             ! Na2SO4 / NaCl system below which 
                                             ! crystalization occurs.
                                             ! (estimated from Tang, I.N.; JGR 102, D2 1883-1893)

  !--- 8) Physical constants: ----------------------------------------------------
  !
  !--- 8.1) General physical constants: 

  REAL(dp), PARAMETER :: bk      = 1.38e-16_dp,   & ! Bolzman constant []   ! ### use ak from mo_constants and scale!
                     rerg    = 8.314E+7_dp,   & ! Ideal gas constant [erg.K-1.mole-1]  ! ### use ar from mo_constants and scale
                     r_kcal  = 1.986E-3_dp      ! Ideal gas constant [kcal K-1.mole-1] ! ### scale from rerg or ar
  
  !--- 8.2) Type specific physical constants:
  !
  REAL(dp), PARAMETER :: dh2so4  = 1.841_dp,      & ! Density          H2SO4  [g cm-3]
                     ddust   = 2.650_dp,      & ! Density          du     [g cm-3]
                     dbc     = 2._dp,         & ! Density          bc     [g cm-3]
                     doc     = 2._dp,         & ! Density          oc     [g cm-3]
                     dnacl   = 2.165_dp,      & ! Density          NaCl   [g cm-3]
                     dna2so4 = 2.68_dp,       & ! Density          Na2SO4 [g cm-3]
                     dnahso4 = 2.435_dp,      & ! Density          NaHSO4 [g cm-3]
                     dh2o    = 1.0_dp,        & ! Density          H2O    [g cm-3]

!### mgs: preferable to make use of mw constants in mo_aero ! ###
                     wh2so4  = 98.0734_dp,    & ! Molecular weight H2SO4  [g mol-1]
                     wh2o    = 18.0_dp,       & ! Molecular weight H2O    [g mol-1]
!!mgs!!                     wso4    = 96.0576_dp,    & ! Molecular weight SO4    [g mol-1]
!!mgs!!                     wso2    = 64.0_dp,       & ! Molecular weight SO2    [g mol-1]
                     wna     = 22.99_dp,      & ! Atomic    weight Na     [g mol-1]
                     wcl     = 35.453_dp,     & ! Atomic    weight Cl     [g mol-1]
                     wnacl   = 58.443_dp,     & ! Molecular weight NaCl   [g mol-1]
                     wna2so4 = 142.0376_dp,   & ! Molecular weight Na2SO4 [g mol-1]
                     wnahso4 = 120.0555_dp      ! Molecular weight NaHSO4 [g mol-1]

!!mgs!!  !>>dod soa
!!mgs!!  REAL(dp), PARAMETER :: ws = 32._dp            ! atomic weight of sulphur
!!mgs!!  !<<dod

  !--- 9) Assumed parameters: ------------------------------------------------------

!++mgs: renamed! was cLayerThickness
  REAL(dp), PARAMETER :: so4_coating_threshold = 1.0_dp   ! Assumed required layer thickness of
                                                    ! sulfate to transfer an insoluble 
                                                    ! particle to a soluble mode. It is 
                                                    ! given in units of layers of 
                                                    ! monomolecular sulfate. Determines the
                                                    ! transfer rate from insoluble to 
                                                    ! soluble modes. 
  
  !--- 10) Nucleation constants: ---------------------------------------------------
  
  REAL(dp) :: critn ! Smallest possible number of H2SO4 molecules in a nucleation mode particle
  
  !--- 11) Data used for the calculation of the aerosol properties -----------------
  !        under ambient conditions:
  !        (Included the conversion from Pa to hPa in the first parameter.)

  REAL(dp), PARAMETER :: wvb(17)=                                                   &
                     (/   95.80188_dp,     -28.5257_dp,     -1.082153_dp,     0.1466501_dp, &
                         -20627.51_dp,    0.0461242_dp,     -0.003935_dp,      -3.36115_dp, &
                       -0.00024137_dp,  0.067938345_dp, 0.00000649899_dp,   8616124.373_dp, &
                       1.168155578_dp, -0.021317481_dp,   0.000270358_dp, -1353332314.0_dp, &
                      -0.002403805_dp                                              /)

  REAL(dp), PARAMETER :: gmb(9)=                                                 &
                     (/ 1.036391467_dp, 0.00728531_dp, -0.011013887_dp, -0.068887407_dp, &
                        0.001047842_dp, 0.001049607_dp, 0.000740534_dp, -1.081202685_dp, &
                       -0.0000029113_dp                                         /)

  !--- 4) Logical mask for coagulation kernel: -------------------------------------
  !       (The coagulation kernel mask is symmetric and not all 
  !       values are used for physical considerations. As its 
  !       calculation is very expensive, a mask is used to 
  !       calculate only the necessarey elements.)

  !>>dod changed handling of coagulation
  TYPE, PUBLIC :: t_coag
     INTEGER :: mode1
     INTEGER :: mode2
  END TYPE t_coag

  INTEGER, PARAMETER, PUBLIC :: ncoag = 16

  TYPE(t_coag), PUBLIC :: coag_modes(ncoag)
  !<<dod

  !--- 12) Service routines for initialization and auxiliary computations ----------

CONTAINS

  SUBROUTINE m7_initialize

    ! Purpose:
    ! ---------
    ! Initializes constants and parameters 
    ! used in the m7 aerosol model.
    !
    ! Author:
    ! ---------
    ! Philip Stier, MPI                          may 2001
    ! Declan O'Donnell, MPI-M, 2008
    !
    ! Interface:
    ! ---------
    ! *m7_initialize*  is called from *start_ham* in mo_ham_init
    !
    USE mo_ham,              ONLY: sizeclass, nclass
    USE mo_ham_wetdep_data,  ONLY: csr_strat_wat, csr_strat_mix, csr_strat_ice, csr_conv, &
                                   cbcr, cbcs

    IMPLICIT NONE

    INTEGER :: jclass
 
    LOGICAL :: lsedi(naeroclass(HAM_M7)) = (/.FALSE., .FALSE., .TRUE., .TRUE.,  & ! soluble modes
                                        .FALSE., .TRUE., .TRUE. /)  ! insoluble modes
 

    !>>dod soa 
    !---executable procedure

    sizeclass(1)%classname   = "Nucleation soluble"
    sizeclass(1)%shortname  = "NS"
    sizeclass(1)%self       = 1
    sizeclass(1)%lsoluble   = .TRUE.
    sizeclass(1)%lsed       = .FALSE.
    sizeclass(1)%lsoainclass = .FALSE.
    sizeclass(1)%lactivation = .FALSE.

    sizeclass(2)%classname   = "Aitken soluble"
    sizeclass(2)%shortname  = "KS"
    sizeclass(2)%self       = 2
    sizeclass(2)%lsoluble   = .TRUE.
    sizeclass(2)%lsed       = .FALSE.                   
    sizeclass(2)%lsoainclass = .TRUE.
    sizeclass(2)%lactivation = .TRUE.

    sizeclass(3)%classname   = "Accumulation soluble"
    sizeclass(3)%shortname  = "AS"
    sizeclass(3)%self       = 3
    sizeclass(3)%lsoluble   = .TRUE.
    sizeclass(3)%lsed       = .TRUE.                   
    sizeclass(3)%lsoainclass = .TRUE.
    sizeclass(3)%lactivation = .TRUE.

    sizeclass(4)%classname   = "Coarse soluble"
    sizeclass(4)%shortname  = "CS"
    sizeclass(4)%self       = 4
    sizeclass(4)%lsoluble   = .TRUE.
    sizeclass(4)%lsed       = .TRUE.                   
    sizeclass(4)%lsoainclass = .TRUE.
    sizeclass(4)%lactivation = .TRUE.

    sizeclass(5)%classname   = "Aitken insoluble"
    sizeclass(5)%shortname  = "KI"
    sizeclass(5)%self       = 5
    sizeclass(5)%lsoluble   = .FALSE.
    sizeclass(5)%lsed       = .FALSE.                   
    sizeclass(5)%lsoainclass = .TRUE.
    sizeclass(5)%lactivation = .FALSE.

    sizeclass(6)%classname   = "Accumulation insoluble"
    sizeclass(6)%shortname  = "AI"
    sizeclass(6)%self       = 6
    sizeclass(6)%lsoluble   = .FALSE.
    sizeclass(6)%lsed       = .TRUE.                   
    sizeclass(6)%lsoainclass = .FALSE.
    sizeclass(6)%lactivation = .FALSE.

    sizeclass(7)%classname   = "Coarse insoluble"
    sizeclass(7)%shortname  = "CI"
    sizeclass(7)%self       = 7
    sizeclass(7)%lsoluble   = .FALSE.
    sizeclass(7)%lsed       = .TRUE.                   
    sizeclass(7)%lsoainclass = .FALSE.
    sizeclass(7)%lactivation = .FALSE.
    !<<dod

    !---backward compatibility stuff
    DO jclass=1,nclass
       sizeclass(jclass)%lsed = lsedi(jclass) 
    END DO
    
    
    ! The following properties could be incorporated into sizeclass but it means updating
    ! many impacted modules and subroutines...

    DO jclass=1, nclass

       !--- 1) Calculate conversion factors for lognormal distributions:----
       !       Radius of average mass (ram) to count median radius (cmr) and 
       !       vice versa. Count median radius to radius of average 
       !       mass (ram).
       !       These factors depend on the standard deviation (sigma)
       !       of the lognormal distribution.
       !       (Based on the Hatch-Choate Conversins Equations;
       !        see Hinds, Chapter 4.5, 4.6 for more details.
       !        In particular equation 4.53.)

       !--- Count Median Radius to Mass Median Radius:

       cmedr2mmedr(jclass) = EXP(3.0_dp*(LOG(sigma(jclass)))**2)

       !--- Count Median Radius to Mass Mean Radius:

       cmr2mmr(jclass) = EXP(3.5_dp*(LOG(sigma(jclass)))**2)

       !--- Count Median Radius to Radius of Average Mass:

       cmr2ram(jclass) = EXP(1.5_dp*(LOG(sigma(jclass)))**2)

       !--- Radius of Average Mass to Count Median Radius:

       ram2cmr(jclass) = 1._dp / cmr2ram(jclass)

       !--- Count Median Radius to Radius of Average Surface:

       cmr2ras(jclass) = EXP(1.0_dp*(LOG(sigma(jclass)))**2)


       !--- 2) Calculate the natural logarithm of the standard deviation:

       sigmaln(jclass) = LOG(sigma(jclass))

    END DO

    !--- 3) Nucleation mode constants:
    !
    !    3.1) Set the lower mode boundary particle dry radius for the nucleation
    !         mode (does not depend on the choice of the nucleation scheme, as
    !         we use different ones which produce particles of different sizes):
    
    crdiv(1) = 0.5E-7_dp ! cm
    
    !    3.2) Smallest possible number of H2SO4 molecules in a nucleation mode
    !         particle:
    
    critn = crdiv(1)**(3.0_dp)*pi*avo*dh2so4/wh2so4/0.75_dp
    
    !--------------------------------------------------------------------

    !>>dod optimisation of coagulation
    ! nucleation mode 
    coag_modes(1)%mode1 = inucs
    coag_modes(1)%mode2 = inucs
    coag_modes(2)%mode1 = inucs
    coag_modes(2)%mode2 = iaits
    coag_modes(3)%mode1 = inucs
    coag_modes(3)%mode2 = iaccs
    coag_modes(4)%mode1 = inucs
    coag_modes(4)%mode2 = icoas
    coag_modes(5)%mode1 = inucs
    coag_modes(5)%mode2 = iaiti
    coag_modes(6)%mode1 = inucs
    coag_modes(6)%mode2 = iacci
    coag_modes(7)%mode1 = inucs
    coag_modes(7)%mode2 = icoai

    ! aitken soluble mode
    coag_modes(8)%mode1 = iaits
    coag_modes(8)%mode2 = iaits
    coag_modes(9)%mode1 = iaits
    coag_modes(9)%mode2 = iaccs
    coag_modes(10)%mode1 = iaits
    coag_modes(10)%mode2 = icoas
    coag_modes(11)%mode1 = iaits
    coag_modes(11)%mode2 = iaiti
    coag_modes(12)%mode1 = iaits
    coag_modes(12)%mode2 = iacci
    coag_modes(13)%mode1 = iaits
!gf(#135)    coag_modes(13)%mode2 = icoas
    coag_modes(13)%mode2 = icoai

    ! accumulation soluble mode
    coag_modes(14)%mode1 = iaccs
    coag_modes(14)%mode2 = iaccs
    coag_modes(15)%mode1 = iaccs
!gf #158    coag_modes(15)%mode2 = iacci
    coag_modes(15)%mode2 = iaiti

    ! aitken insoluble mode
    coag_modes(16)%mode1 = iaiti
    coag_modes(16)%mode2 = iaiti

    !--- 4) Set prescribed scavenging ratios
    !
    !SF Note: this was formerly hardcoded in mo_ham_m7_wetdep_data.
    !         In order to handle more transparently M7 and SALSA, it is
    !         necessary to make this initialization dynamic.

    IF (.NOT. ALLOCATED(csr_strat_wat)) ALLOCATE(csr_strat_wat(naeroclass(HAM_M7)))
    IF (.NOT. ALLOCATED(csr_strat_mix)) ALLOCATE(csr_strat_mix(naeroclass(HAM_M7)))
    IF (.NOT. ALLOCATED(csr_strat_ice)) ALLOCATE(csr_strat_ice(naeroclass(HAM_M7)))
    IF (.NOT. ALLOCATED(csr_conv))      ALLOCATE(csr_conv(naeroclass(HAM_M7)))
    IF (.NOT. ALLOCATED(cbcr))          ALLOCATE(cbcr(naeroclass(HAM_M7)))
    IF (.NOT. ALLOCATED(cbcs))          ALLOCATE(cbcs(naeroclass(HAM_M7)))

    csr_strat_wat(1:naeroclass(HAM_M7)) = (/0.10_dp, 0.25_dp, 0.85_dp, 0.99_dp, 0.20_dp, 0.40_dp, 0.40_dp/)
    csr_strat_mix(1:naeroclass(HAM_M7)) = (/0.10_dp, 0.40_dp, 0.75_dp, 0.75_dp, 0.10_dp, 0.40_dp, 0.40_dp/)
    csr_strat_ice(1:naeroclass(HAM_M7)) = (/0.10_dp, 0.10_dp, 0.10_dp, 0.10_dp, 0.10_dp, 0.10_dp, 0.10_dp/)
    csr_conv(1:naeroclass(HAM_M7))      = (/0.20_dp, 0.60_dp, 0.99_dp, 0.99_dp, 0.20_dp, 0.40_dp, 0.40_dp/)
    !--- Mean mass scavenging coefficients normalized by rain-rate [kg m-2]:
    ! Rain: Seinfeld & Pandis, Fig 20.15:
    cbcr(1:naeroclass(HAM_M7))          = &
                        (/ 5.0E-4_dp, 1.0E-4_dp, 1.0E-3_dp, 1.0E-1_dp, 1.0E-4_dp, 1.0E-3_dp, 1.0E-1_dp /)
    !    Snow: Little available, graphs in Prupbacher show similar order of magnitude in
    !          collection efficiency as rain,therefore assume typical mean value as above:
    cbcs(1:naeroclass(HAM_M7))          = &
                        (/ 5.0E-3_dp, 5.0E-3_dp, 5.0E-3_dp, 5.0E-3_dp, 5.0E-3_dp, 5.0E-3_dp, 5.0E-3_dp /)

  END SUBROUTINE m7_initialize


  SUBROUTINE sethamM7
    
    ! *sethamM7* modifies pre-set switches of the aeroM7ctl
    !             namelist for the configuration of the 
    !             M7 component of the ECHAM/HAM aerosol model
    ! 
    ! Authors:
    ! --------
    ! Philip Stier, MPI-M                                                12/2002
    ! Jan Kazil, MPI-M                                       2008-03-03 20:34:52
    !
    ! *sethamM7* is called from *init_subm* in mo_submodel_interface
    !
    
    USE mo_mpi,         ONLY: p_parallel_io, p_bcast, p_io
    USE mo_namelist,    ONLY: open_nml, position_nml, POSITIONED
    USE mo_exception,   ONLY: message, em_warn, em_info
    USE mo_util_string, ONLY: separator
    USE mo_ham,         ONLY: lgcr
    
    IMPLICIT NONE
    
    INCLUDE 'ham_m7ctl.inc'
    
    ! Local variables:
    
    INTEGER :: ierr, inml, iunit

    ! Read the namelist with the switches:
    
    CALL message('',separator)
    CALL message('sethamM7', 'Reading namelist ham_m7ctl...', level=em_info)
    IF (p_parallel_io) THEN

      inml = open_nml('namelist.echam') 
      iunit = position_nml ('HAM_M7CTL', inml, status=ierr)
      SELECT CASE (ierr)
      CASE (POSITIONED)
      READ (iunit, ham_m7ctl)
      END SELECT
      
    ENDIF
    
    ! Broadcast the switches over the processors:
    
    CALL p_bcast (nwater,     p_io)
    CALL p_bcast (nsnucl,     p_io)
    CALL p_bcast (nonucl,     p_io)
    CALL p_bcast (lnucl_stat, p_io)

    !--- error checking
    IF (nsnucl > 1 .AND. .NOT. lgcr) THEN
      CALL message('sethamM7', 'nsnucl > 1 requires lgcr=.TRUE.! Setting lgcr=.TRUE. now.', &
                   level=em_warn)
      lgcr = .TRUE.
    END IF
    
    !--- write the values of the switches:
    CALL sethamM7_log(nwater,nsnucl,nonucl,lnucl_stat)
        
  END SUBROUTINE sethamM7
  
  SUBROUTINE sethamM7_log(nwater,nsnucl,nonucl,lnucl_stat)
    
    ! *sethamM7_log* writes the values of the given switches in a given output
    ! unit and conducts selceted consistency checks.
    !
    ! Authors:
    ! --------
    !
    ! 2008 Jan Kazil, MPI-M
    
    USE mo_exception,    ONLY: message, message_text, em_param, em_error, em_info
    USE mo_submodel,     ONLY: print_value
    USE mo_util_string,  ONLY: separator 
    USE mo_ham,          ONLY: lgcr
    
    IMPLICIT NONE
    
    !
    ! Input variables:
    !
    
    LOGICAL :: lnucl_stat
    INTEGER :: nwater, nsnucl, nonucl
  
    CALL message('','',level=em_param) 
    CALL message('sethamM7','Initialization of the M7 aerosol module', level=em_info) 
   
    CALL print_value('nwater', nwater) 
    SELECT CASE(nwater)
    CASE (0)
      CALL message('', ' --> Jacobson, Tabazadeh, and Turco (Jacobson et al., JGR 1996)',  &
                   level=em_param)
    CASE (1)
      CALL message('', ' --> Kappa-Koehler theory based approach (Petters and Kreidenweis, ACP 2007)', &
                   level=em_param)
    CASE DEFAULT
      WRITE(message_text,'(a,i0,a)') 'nwater must be 0 or 1 (present value = ',nwater,')'
      CALL message('sethamM7', message_text, level=em_error)
    END SELECT
    
   
    CALL print_value('nsnucl', nsnucl)
    SELECT CASE(nsnucl)
    CASE (0) 
      CALL message('', 'H2SO4/H2O nucleation off',  &
                   level=em_param)
    CASE (1) 
      CALL message('', ' --> Neutral H2SO4/H2O nucleation (Vehkamaeki et al., JGR 2002)',  &
                   level=em_param)
    CASE (2) 
      CALL message('', ' --> Neutral and charged H2SO4/H2O nucleation (Kazil and Lovejoy, ACP 2007)',  &
                   level=em_param)
    CASE DEFAULT
      WRITE(message_text,'(a,i0,a)') 'nsnucl must be 0, 1 or 2 (present value = ',nsnucl,')'
      CALL message('sethamM7', message_text, level=em_error)
    END SELECT
    
    
    CALL print_value('nonucl', nonucl)
    SELECT CASE(nonucl)
    CASE (0)
      CALL message('', ' --> Organic aerosol nucleation off',  &
                   level=em_param)
    CASE (1)
      CALL message('', ' --> Organic aerosol (activation) nucleation (Kulmala et al., ACP 2006)',  &
                   level=em_param)
    CASE (2)
      CALL message('', ' --> Organic aerosol (kinetic) nucleation after (Laakso et al., ACP 2004)',  &
                   level=em_param)
    CASE DEFAULT
      WRITE(message_text,'(a,i0,a)') 'nonucl must be 0, 1 or 2 (present value = ',nonucl,')'
      CALL message('sethamM7', message_text, level=em_error)
    END SELECT
    
    
    IF (lnucl_stat)  THEN
      CALL message('', 'lnucl_stat=.TRUE. : Sampling the cloud-free volume as function of',  &
                   level=em_param)
      CALL message('', 'T, RH, [H2SO4(g)], H2SO4 condensation sink, and ionization rate',    &
                   level=em_param)
    END IF
    
    CALL message('', separator) 
    
  END SUBROUTINE sethamM7_log

END MODULE mo_ham_m7ctl
