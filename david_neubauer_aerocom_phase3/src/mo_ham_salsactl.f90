!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_salsactl.f90
!!
!! \brief
!! mo_ham_salsactl contains parameters, switches and initialization routines for the salsa aerosol module.
!!
!! \responsible_coder
!! Harri Kokkola, harri.kokkola@fmi.fi
!!
!! \revision_history
!!   -# T. Bergman (FMI) - nmod->nclass to facilitate new aerosol models (2013-02-05)
!!   -# A. Laakso (FMI) - salsa_initialize (2013-02-25)
!!   -# H. Kokkola (FMI) - Adopted common subroutines for M7 and SALSA from mo_ham_m7.f90 (2014)
!! 
!! \limitations
!! Currently, there are two index lists for aerosol species: aero_idx in mo_species
!! and subm_aerospec in this module. I hope these are identical for the current model set-up 
!! in preparation for CMIP5. Later, one may wish to distinguish between the two: aero_idx
!! could contain additional aerosol species (e.g. from MOZART or climatologies), and this could
!! mess up the SALSA code. If this can be generalized: fine. if not we should keep the two 
!! lists separate. mo_ham_rad (for example) works on aero_idx to be independent of SALSA specifics.
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

MODULE mo_ham_salsactl

  USE mo_kind,             ONLY: dp
  USE mo_species,          ONLY: nmaxspec


  IMPLICIT NONE

  PRIVATE

  ! -- subroutines
  PUBLIC :: setham_salsa


  ! -- variables

  PUBLIC :: nsnucl
  PUBLIC :: locgas, lsol2b, act_coeff, nj3


  !++alaak 
  PUBLIC :: in1a,in2a,in2b,fn1a,fn2a,fn2b,nbins
  PUBLIC :: nbin, nbin2, nbin3,reglim,nlim,nreg
  PUBLIC :: epsv,vhilim,vlolim,vratiohi,vratiolo,dpmid

  PUBLIC :: n3,massacc,d_sa,epsoc,slim,ions,surfw0
  PUBLIC :: recalc
  !--alaak 


  !--- 1) Define and pre-set switches for the processes of M7: -----------------------

  !--- Physical:
  INTEGER :: nsnucl = 2 ! Choice of the sulfate aerosol nucleation scheme:
                        !  nsnucl = 1  Binary
                        !         = 2  activation type nucleation
                        !         = 3  Kinetic
                        !         = 4  Ternary
                        !         = 5  nucleation with ORGANICs
                        !         = 6  activation type of nucleation with H2SO4+ORG
                        !         = 7  heteromolecular nucleation with H2SO4*ORG
                        !         = 8  homomolecular nucleation of  H2SO4 + 
                        !              heteromolecular nucleation with H2SO4*ORG
                        !         = 9  homomolecular nucleation of  H2SO4 and ORG + 
                        !              heteromolecular nucleation with H2SO4*ORG

  LOGICAL :: locgas = .FALSE.,&   ! emission of organic carbon in gas phase
             lsol2b = .FALSE.     ! repartitioning of insoluble material in 
                                  ! case of increase in solubility 

  LOGICAL :: recalc   = .FALSE.   ! recalculation of wet diameter between
                                  ! calculation of microphysical processes

  INTEGER ::                    & ! J3 parametrization
             nj3 = 1              ! 1 = condensational sink (Kerminen&Kulmala, 2002)
                                  ! 2 = coagulational sink (Lehtinen et al. 2007)
                                  ! 3 = coagS+self-coagulation (Anttila et al. 2010)
  REAL(dp) :: act_coeff=1.e-7_dp  ! activation coefficient

!---------------------------------------------------------------------------------------

  !--Indices corresponding to number concentration, radius
  !    and chemical components in each subregime

  INTEGER, PARAMETER ::            &
   nreg = 2                          ! number of main size regimes


  REAL(dp), PARAMETER ::                       &
   reglim(nreg+2) =                            & ! low/high diameter limits
    (/ 3.e-9_dp, 5.e-8_dp, 7.e-7_dp, 1.e-5_dp /) ! of main size regimes [m]

   INTEGER, PARAMETER :: &
   nbin(nreg)   = (/ 3, 7 /)   ! number of bins in each main regime

  INTEGER, PARAMETER ::      &
   nbin2 = 4,                & ! number of bins in former 2-region
   nbin3 = nbin(2) - nbin2     ! number of bins in former 3-region

  INTEGER, PARAMETER ::      & ! number/radius: start index
   in1a = 1,                 & ! regime 1a
   in2a = in1a + nbin(1),    & ! regime 2a
   in2b = in2a + nbin(2),    & ! regime 2b

                               ! number/radius: last index
   fn1a = in2a - 1,          & ! regime 1a
   fn2a = fn1a + nbin(2),    & ! regime 2a
   fn2b = fn2a + nbin(2),    & ! regime 2b

   nbins = fn2b                ! total number of size bins


  REAL(dp) ::               &
       epsv(nbins),         &
       vhilim(fn2b),        &
       vlolim(fn2b),        &
       vratiohi(fn2b),      &
       vratiolo(fn2b),      &
       dpmid(fn2b),         &
       csr_strat_wat(fn2b), &
       csr_strat_mix(fn2b), &
       csr_strat_ice(fn2b), &
       csr_conv(fn2b),      &
       zbcr(fn2b)

  REAL(dp), PARAMETER ::     &
   deltav = 1.096e-7_dp,     & ! vapor jump length (m)
   deltaT = 2.16e-7_dp,      & ! thermal jump length (m)
   alphaT = 0.96_dp,         & ! thermal accomodation coefficient
   alphac = 1.0_dp             ! condensation coefficient

  REAL(dp), PARAMETER ::     & !
   n3 = 158.79_dp               ! number of H2SO4 molecules in 3 nm cluster 
                               !  assuming d_sa = 5.54 ???     
  !-- 4.3) Properties of condensing vapours

  REAL(dp), PARAMETER ::                               & ! diameter of condensing molecule [m]
      d_sa   = 5.539376964394570e-10_dp,               &

      d_oc   = 6.195906936656752e-10_dp,               &

      d_h2o  = 3.851565216195334e-10_dp

  REAL(dp), PARAMETER :: &
       slim = 1.005_dp,  & ! water saturation used as limit
       ions = 3.0_dp,    & ! van't Hoff factor (ions produced upon dissociation)
       surfw0 = 0.073_dp,& ! surface tension of pure water @ ~ 293 K [J/m2]
       epsoc = 0.15_dp     ! water uptake of organic material

  !-- 7) Parameters for cloud activation

  REAL(dp), PARAMETER :: crcut=0.035*1E-6_dp ! Assumed lower cut-off of the
                                             ! aerosol size distribution [m]

  !--- Ulrike: included for activation in convective clouds
  REAL(dp), PARAMETER :: crcut_cv=0.025*1E-6_dp ! Assumed lower cut-off of the
  
  
  REAL(dp) :: cfracn(fn2b)
  
  REAL(dp) :: zfracn(fn2b),   &
              zfracn_cv(fn2b)
  REAL(dp), PARAMETER :: &
   massacc(nbins) = 1._dp

  REAL(dp), PARAMETER :: &
   nlim = 1._dp,         & ! number conc. limit below which bin empty  [#/m3] 
   m3_2_um3 = 1.e+18_dp    ! conversion factor for volume from m3 to um3

  INTEGER, PUBLIC :: iso4b(fn2b), iocb(fn2b), ibcb(fn2b), idub(fn2b), issb(fn2b)

  !--- 12) Service routines for initialization and auxiliary computations ----------

CONTAINS





  SUBROUTINE setham_salsa
    
    ! *setham_salsa* modifies pre-set switches of the hamsalsactl
    !             namelist for the configuration of the 
    !             SALSA component of the ECHAM/HAM aerosol model
    ! 
    ! Authors:
    ! --------
    ! Philip Stier, MPI-M                                                12/2002
    ! Jan Kazil, MPI-M                                       2008-03-03 20:34:52
    ! Anton Laakso, FMI  (SALSA                              6/2013
    ! *setham_salsa* is called from *start_ham* in mo_ham_init
    !
    
    USE mo_mpi,         ONLY: p_parallel, p_parallel_io, p_bcast, p_io
    USE mo_namelist,    ONLY: open_nml, position_nml, POSITIONED
    USE mo_exception,   ONLY: message, message_text, em_error, em_warn, em_info
    USE mo_util_string, ONLY: separator
    
    IMPLICIT NONE
    
    INCLUDE 'ham_salsactl.inc'
    
    ! Local variables:
    
    INTEGER :: ierr, inml, iunit

    ! Read the namelist with the switches:
    
    CALL message('',separator)
    CALL message('setham_salsa', 'Reading namelist ham_salsactl...', level=em_info)
    IF (p_parallel_io) THEN

      inml = open_nml('namelist.echam') 
      iunit = position_nml ('HAM_SALSACTL', inml, status=ierr)
      SELECT CASE (ierr)
      CASE (POSITIONED)
      READ (iunit, ham_salsactl)
      END SELECT
      
    ENDIF
    
    ! Broadcast the switches over the processors:
    CALL p_bcast (nsnucl,     p_io)
    CALL p_bcast (lsol2b, p_io)
    CALL p_bcast (locgas, p_io)
    CALL p_bcast (nj3,  p_io)

    !--- write the values of the switches:
    CALL sethamsalsa_log(locgas,lsol2b,nsnucl)

  END SUBROUTINE setham_salsa

  SUBROUTINE sethamsalsa_log(locgas,lsol2b,nsnucl)
    USE mo_exception,    ONLY: finish, message, message_text, em_param, em_error
    USE mo_submodel,     ONLY: print_value
    USE mo_util_string,  ONLY: separator 
    
    IMPLICIT NONE

    LOGICAL :: locgas,lsol2b
    INTEGER :: nsnucl

    CALL message('', separator) 
    CALL message('sethamsalsa_log','Initialization of settings for aerosol module SALSA', level=em_param) 

    CALL print_value('locgas', locgas)
    CALL print_value('lsol2b', lsol2b)

    SELECT CASE(nsnucl)
    CASE (0)
    CALL message('', 'Nucleation = OFF',  &
                   level=em_param)
    CASE (1)
      CALL message('', 'Binary nucleation',  &
                   level=em_param)
    CASE (2)
      CALL message('', 'Activation type nucleation', &
                   level=em_param)
    CASE (3)
      CALL message('', 'Kinetic nucleation', &
                   level=em_param)
    CASE (4)
      CALL message('', 'Ternary nucleation', &
                   level=em_param)
    CASE (5)
      CALL message('', 'nucleation with organics', &
                   level=em_param)
    CASE (6)
      CALL message('', 'activation type of nucleation with H2SO4+ORG', &
                   level=em_param)
    CASE (7)
      CALL message('', 'heteromolecular nucleation with H2SO4*ORG', &
                   level=em_param)
    CASE (8)
      CALL message('', 'homomolecular nucleation of  H2SO4 + heteromolecular nucleation with H2SO4*ORG', &
                   level=em_param)
    CASE (9)
      CALL message('', 'homomolecular nucleation of  H2SO4 and ORG +heteromolecular nucleation with H2SO4*ORG', &
                   level=em_param)
    CASE DEFAULT
      WRITE(message_text,*) 'nucleation requested but no scheme selected (present value = ',nsnucl,')'
      CALL message('sethamsalsa_log', message_text, level=em_error)
    END SELECT

    ! >> thk: nuclation that includes organics is currently not working:
    IF (nsnucl > 4) THEN
       CALL message(&
            '',&
            'ERROR: Nucleation that includes organics is currently not working.'&
            )

       CALL finish(&
            'sethamsalsa_log',&
            'run terminated'&
            )
    END IF


    CALL message('', separator)

  END SUBROUTINE sethamsalsa_log

END MODULE mo_ham_salsactl
