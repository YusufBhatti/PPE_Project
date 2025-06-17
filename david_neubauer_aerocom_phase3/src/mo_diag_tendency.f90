!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_diag_tendency
!BOP
  ! !MODULE: mo_diag_tendency

  ! !DESCRIPTION:
  !\begin{verbatim}
  !-------------------------------------------------------------------
  ! overview about the separation of terms in diagnostic equations
  !
  !  LTDIAG      additional switch in &RUNCTL
  !              .true. ==> addional diagnostics
  !
  ! ***************** Interface to ECHAM *************************
  !
  ! *ECHAM*       *DIAGNOSTICS*
  !
  ! CONTROL -+
  !          +--> DIAG_Init(1)                general initialization
  !                  |
  !                  +---> DIAG_Rerun('R')    read rerun data
  !
  ! STEPON --+
  !          +--> DIAG_SumAll                 count total tendency
  !          |
  !          +--> DIAG_Write                  store data in model output
  !          |       |
  !          |       +---> DIAG_Check
  !          |       |
  !          |       +---> DiagWriteSP
  !          |
  !          +--> DIAG_Init(-1)               clear diagnostic memory
  !                  |
  !                  +---> DIAG_Rerun('W')    store rerun data
  !
  ! SCAN1SL -+
  !          +--> DIAG_fftd
  !          |
  !          SI2
  !          |
  !          +--> DIAG_SpecTrans
  !                  |
  !                  +--> DIAG_sym1
  !                  |
  !                  +--> DIAG_ltd
  !
  !---------------------------------------------------------------------
  ! count terms in additional arrays
  !   *spectral*    *Gaussian*
  !
  !                 DYN            adiabatic terms
  !                  |
  !                 TF2            time filter
  !                  |
  !                 TF1            time filter
  !                  |
  !                 GPC
  !                  +--> PHYSC    diabatic terms
  !                  |      |
  !                  |      +---> M_RADHEAT long/short wave
  !                  |
  !                  +--> SI1      semi-implicit terms
  !                  |
  !                 DIAG_fftd
  !                 |
  !               SI2                semi-implicit terms
  !               |
  !             DIAG_SpecTrans
  !                     |
  !                   DIAG_sym1
  !                     |
  !                   DIAG_ltd
  ! SCCD                           semi-implicit terms
  !   |
  ! SCCTP                          semi-implicit terms
  !   |
  ! HDIFF                          horizontal diffusion
  !   |
  ! DIAG_SumAll
  !   |
  ! DIAG_Write
  !
  !**************************************************************
  !\end{verbatim}

  ! !REVISION HISTORY: 
  ! I. Kirchner, MPI, October 1998, basic version
  ! I. Kirchner, MPI, May-2000, patch E4v3.22p1
  ! I. Kirchner, MPI, May-2002, revision E5R1.04

  ! !USES:
  USE mo_exception,     ONLY: finish, message
  USE mo_linked_list,   ONLY: t_stream
  USE mo_kind,          ONLY: dp
!BOX
  IMPLICIT NONE

  PRIVATE
!EOX
  PUBLIC  :: DIAG_Init      ! allocate/deallocate memory, initialize diagnostic arrays
  PUBLIC  :: DIAG_SpecTrans ! second part of spectral transform
  PUBLIC  :: DIAG_fftd      ! Fourier transformation of diagnostics arrays
  PRIVATE :: DIAG_sym1      ! composition of symetric/antisymetric part
  PRIVATE :: DIAG_ltd       ! legendre transform of diagnostic arrays
  PUBLIC  :: DIAG_Write     ! store diagnostic fields
  PRIVATE :: DIAG_Check     ! global check of diagnostic arrays
!BOX
  ! special diagnostic fields

  REAL(dp), ALLOCATABLE, TARGET, PUBLIC :: &    ! ARRAYS in spectral space
       pdvor (:,:,:,:), &  !   terms for vorticity equation
       pddiv (:,:,:,:), &  !   terms for divergence equation
       pdtem (:,:,:,:)     !   terms for temperature equation
  REAL(dp), POINTER, PUBLIC :: &
       pdprs (:,:,:), &    !   terms for surface pressure equation
       pdprl (:,:,:), &    !   surface pressure term vertical integration
       p4prs (:,:,:,:)

  REAL(dp), ALLOCATABLE, TARGET :: &    ! old value for accumulation of tendency
       pdovor (:,:,:,:), pdodiv (:,:,:,:), pdotep (:,:,:,:)


  INTEGER, PARAMETER :: NO_PTFH = 3
  REAL(dp), POINTER, PUBLIC ::   &  ! ARRAYS in grid point space
       ptfh1 (:,:,:,:), &  !   time integration parts
       ptfh2 (:,:)
  LOGICAL, SAVE, PUBLIC :: lset_fh1    = .FALSE.

  INTEGER, PARAMETER :: NO_PDIGA=25  ! no. of 3-dim grid point space fields
  INTEGER, PARAMETER :: NO_PDSGA=2   ! no. of 2-dim grid point space fields
  INTEGER, PARAMETER :: NO_PDIGB=10  ! no. of 3-dim fourier space fields
  INTEGER, PARAMETER :: NO_PDIGS=4   ! no. of 2-dim mixed grid point and fourier space fields
  INTEGER, PARAMETER :: NO_PDIGS_FRAC=3   ! no. of 2-dim grid point space fields fraction

  REAL(dp), POINTER, PUBLIC :: pdiga(:,:,:,:), pdsga (:,:,:)
  REAL(dp), ALLOCATABLE, PUBLIC, TARGET :: &    ! ARRAYS for accumulation
                        pdigaa(:,:,:,:), pdigas(:,:,:,:), &
                        pdsgaa(:,:,:),   pdsgas(:,:,:),   &
       pdigb (:,:,:,:), pdigba(:,:,:,:), pdigbs(:,:,:,:), &
       pdigs (:,:,:,:), pdigsa(:,:,:,:), pdigss(:,:,:,:), &
       pdsgs (:,:),     pdsgsa(:,:),     pdsgss(:,:)

!EOX
  TYPE(t_stream), POINTER :: tdiag
!BOX
  INTEGER, PUBLIC         :: dio_index ! event index of the TDIAG stream

  INTEGER, PARAMETER, PUBLIC :: &  ! number of terms in different equations
       NDVOR = 9, &! vorticity
       NDDIV = 9, &! divergence
       NDTEM =14, &! temperature
       NDPRS = 4   ! surface pressure

  INTEGER, PARAMETER, PUBLIC :: &  ! function separation
       IDIAG_INIT     = 1, &
       IDIAG_INI_GBUF = 2, &
       IDIAG_INI_PREV = 3, &
       IDIAG_FREE     = 4

  LOGICAL, SAVE :: &
       laccu_ok    = .FALSE., &! mean of total tendency is fine
       ldiag_start = .TRUE.

  CHARACTER(len=256) :: diag_mess = ''

!EOX
!EOP
CONTAINS
!======================================================================
!BOP
  ! !IROUTINE:  DIAG_Init
  ! !INTERFACE:

  SUBROUTINE DIAG_Init(itype)

    ! !DESCRIPTION: 
    ! initialisation of tendency diagnostics

    ! !USES:
    USE mo_memory_base, ONLY: new_stream, default_stream_setting, &
                              add_stream_element,                 &
                              SPECTRAL, GAUSSIAN, GRIB,     &
                              SURFACE
    USE mo_mpi,         ONLY: p_parallel
    USE mo_output,      ONLY: nudging_table
    USE mo_memory_sp,   ONLY: sd, svo, stp
    USE mo_control,     ONLY: nlev, nlevp1, nsp, ngl, nmp1, nlon

    INTEGER, INTENT(in) :: itype       ! function separator
!EOP
!BOC
!BOX
    REAL(dp), POINTER  :: pdum(:,:,:), p4(:,:,:,:)
    LOGICAL            :: lpdo1 = .FALSE., lpdo2 = .FALSE.
!EOX
    IF (p_parallel) &
         CALL finish('mo_diag_tendency:DIAG_Init','not prepared for parallel mode')
    SELECT CASE(itype)

    CASE(IDIAG_INIT)             ! startup initialization
!BOX
      CALL message('mo_diag_tendency:DIAG_Init','  ------- start initialization -----')
!EOX
      CALL message('','Tendency Diagnostics (E5R1.04) 22-Aug-2002 (kirchner@dkrz.de)')
!BOX
      ! define tendency diagnostic output stream

      CALL new_stream ( tdiag, 'tdiag', lrerun=.TRUE., filetype=GRIB)
      dio_index = tdiag%post_idx

      CALL default_stream_setting ( tdiag, &
           table=nudging_table, bits=24, units='(X)/s', contnorest=.TRUE., &
           lpost  = .TRUE., lrerun = .TRUE., laccu  = .TRUE.,              &
           repr=SPECTRAL, ldims=(/nlev, 2, nsp/), gdims=(/nlev, 2, nsp/)   )
!EOX
      ! -----------------------------------------------------------
      ! allocate spectral space arrays, these fields are stored

      ALLOCATE (pdvor (nlev,2,nsp,NDVOR)); pdvor(:,:,:,:) = 0.0_dp !  vorticity equation

      p4 => pdvor(:,:,:,1:1)
      CALL add_stream_element(tdiag,'VOR1',  pdum,p4=p4, code=41, &
           longname = 'VorEQ horizontal advec + press.grad. + cori.term (DYN)')
      p4 => pdvor(:,:,:,2:2)
      CALL add_stream_element(tdiag,'VOR2',  pdum,p4=p4, code=42, &
           longname = 'VorEQ vertical advection (DYN)')
      p4 => pdvor(:,:,:,3:3)
      CALL add_stream_element(tdiag,'VOR3',  pdum,p4=p4, code=43, &
           longname = 'VorEQ vertical diffusion due to impuls (VDIFF)')
      p4 => pdvor(:,:,:,4:4)
      CALL add_stream_element(tdiag,'VOR4',  pdum,p4=p4, code=44, &
           longname = 'VorEQ gravity wave drag (GWDRAG)')
      p4 => pdvor(:,:,:,5:5)
      CALL add_stream_element(tdiag,'VOR5',  pdum,p4=p4, code=45, &
           longname = 'VorEQ moisture mass flux (CUCALL)')
      p4 => pdvor(:,:,:,6:6)
      CALL add_stream_element(tdiag,'VOR6',  pdum,p4=p4, code=46, &
           longname = 'VorEQ timefilter')
      p4 => pdvor(:,:,:,7:7)
      CALL add_stream_element(tdiag,'VOR7',  pdum,p4=p4, code=47, &
           longname = 'VorEQ semi-implicit part of time integration')
      p4 => pdvor(:,:,:,8:8)
      CALL add_stream_element(tdiag,'VOR8',  pdum,p4=p4, code=48, &
           longname = 'VorEQ horizontal diffusion')
      p4 => pdvor(:,:,:,9:9)
      CALL add_stream_element(tdiag,'VORSUM',pdum,p4=p4, code=49, &
           longname = 'VorEQ total tendency')
      
      ALLOCATE (pddiv (nlev,2,nsp,NDDIV)); pddiv(:,:,:,:) = 0.0_dp !  divergence equation

      p4 => pddiv(:,:,:,1:1)
      CALL add_stream_element(tdiag,'DIV1',  pdum,p4=p4, code=61, &
           longname = 'DivEQ horizontal advec + cori. + press.grad. + G-term (DYN)')
      p4 => pddiv(:,:,:,2:2)
      CALL add_stream_element(tdiag,'DIV2',  pdum,p4=p4, code=62, &
           longname = 'DivEQ vertical advection')
      p4 => pddiv(:,:,:,3:3)
      CALL add_stream_element(tdiag,'DIV3',  pdum,p4=p4, code=63, &
           longname = 'DivEQ vertical diffusion due to impuls (VDIFF)')
      p4 => pddiv(:,:,:,4:4)
      CALL add_stream_element(tdiag,'DIV4',  pdum,p4=p4, code=64, &
           longname = 'DivEQ gravity wave drag (GWDRAG)')
      p4 => pddiv(:,:,:,5:5)
      CALL add_stream_element(tdiag,'DIV5',  pdum,p4=p4, code=65, &
           longname = 'DivEQ moisture mass flux (CUCALL)')
      p4 => pddiv(:,:,:,6:6)
      CALL add_stream_element(tdiag,'DIV6',  pdum,p4=p4, code=66, &
           longname = 'DivEQ timefilter')
      p4 => pddiv(:,:,:,7:7)
      CALL add_stream_element(tdiag,'DIV7',  pdum,p4=p4, code=67, &
           longname = 'DivEQ semi-implicit part of time integration')
      p4 => pddiv(:,:,:,8:8)
      CALL add_stream_element(tdiag,'DIV8',  pdum,p4=p4, code=68, &
           longname = 'DivEQ horizontal diffusion')
      p4 => pddiv(:,:,:,9:9)
      CALL add_stream_element(tdiag,'DIVSUM',pdum,p4=p4, code=69, &
           longname = 'DivEQ total tendency')

      ALLOCATE (pdtem (nlev,2,nsp,NDTEM)); pdtem(:,:,:,:) = 0.0_dp ! temperature

      p4 => pdtem(:,:,:,1:1)
      CALL add_stream_element(tdiag,'TEM01', pdum,p4=p4,  code=81, &
           longname = 'TempEQ horizontal advection (DYN)')
      p4 => pdtem(:,:,:,2:2)
      CALL add_stream_element(tdiag,'TEM02', pdum,p4=p4,  code=82, &
           longname = 'TempEQ vertical advection (DYN)')
      p4 => pdtem(:,:,:,3:3)
      CALL add_stream_element(tdiag,'TEM03', pdum,p4=p4,  code=83, &
           longname = 'TempEQ energy conversion (DYN)')
      p4 => pdtem(:,:,:,4:4)
      CALL add_stream_element(tdiag,'TEM04', pdum,p4=p4,  code=84, &
           longname = 'TempEQ radiation (RADHEAT)')
      p4 => pdtem(:,:,:,5:5)
      CALL add_stream_element(tdiag,'TEM05', pdum,p4=p4,  code=85, &
           longname = 'TempEQ vertical diffusion due to turbulence (VDIFF)')
      p4 => pdtem(:,:,:,6:6)
      CALL add_stream_element(tdiag,'TEM06', pdum,p4=p4,  code=86, &
           longname = 'TempEQ gravity wave drag (GWDRAG)')
      p4 => pdtem(:,:,:,7:7)
      CALL add_stream_element(tdiag,'TEM07', pdum,p4=p4,  code=87, &
           longname = 'TempEQ convection (CUCALL)')
      p4 => pdtem(:,:,:,8:8)
      CALL add_stream_element(tdiag,'TEM08', pdum,p4=p4,  code=88, &
           longname = 'TempEQ large scale cloud processes (COND)')
      p4 => pdtem(:,:,:,9:9)
      CALL add_stream_element(tdiag,'TEM09', pdum,p4=p4,  code=89, &
           longname = 'TempEQ timefilter')
      p4 => pdtem(:,:,:,10:10)
      CALL add_stream_element(tdiag,'TEM10', pdum,p4=p4, code=90, &
           longname = 'TempEQ semi-implicit part of time integration')
      p4 => pdtem(:,:,:,11:11)
      CALL add_stream_element(tdiag,'TEM11', pdum,p4=p4, code=91, &
           longname = 'TempEQ horizontal diffusion')
      p4 => pdtem(:,:,:,12:12)
      CALL add_stream_element(tdiag,'TEM12', pdum,p4=p4, code=92, &
           longname = 'TempEQ longwave radiation')
      p4 => pdtem(:,:,:,13:13)
      CALL add_stream_element(tdiag,'TEM13', pdum,p4=p4, code=93, &
           longname = 'TempEQ shortwave radiation')
      p4 => pdtem(:,:,:,14:14)
      CALL add_stream_element(tdiag,'TEMSUM',pdum,p4=p4, code=94, &
           longname = 'TempEQ total tendency')
      
      ! divergence for each layer
      CALL add_stream_element(tdiag,'PRL',   pdprl, code=100, &
           longname = 'PresEQ convergence in each layer')

      ALLOCATE (p4prs (1,2,nsp,NDPRS)); p4prs(:,:,:,:) = 0.0_dp !  pressure equation
!BOX
      pdprs => p4prs(1,:,:,:)
      CALL default_stream_setting ( tdiag, &
            table=nudging_table, bits=24, &
            lpost  = .TRUE., lrerun = .TRUE., laccu  = .TRUE., &
            repr=SPECTRAL, leveltype=SURFACE, &
            ldims=(/1, 2, nsp/), gdims=(/1, 2, nsp/),&
            units='(X)/s', no_default=.TRUE.)
!EOX
      p4 => p4prs(1:1,:,:,1:1)
      CALL add_stream_element(tdiag,'PRS1',  pdum, code=101, klev=1, p4=p4, &
           ldims=(/1, 2, nsp/) ,gdims=(/1, 2, nsp/), &
           longname = 'PresEQ vertical integrated convergence')
      p4 => p4prs(1:1,:,:,2:2)
      CALL add_stream_element(tdiag,'PRS2',  pdum, code=102, klev=1, p4=p4, &
           ldims=(/1, 2, nsp/) ,gdims=(/1, 2, nsp/), &
           longname = 'PresEQ timefilter')
      p4 => p4prs(1:1,:,:,3:3)
      CALL add_stream_element(tdiag,'PRS3',  pdum, code=103, klev=1, p4=p4, &
           ldims=(/1, 2, nsp/) ,gdims=(/1, 2, nsp/), &
           longname = 'PresEQ semi-implicit part of time integration')
      p4 => p4prs(1:1,:,:,4:4)
      CALL add_stream_element(tdiag,'PRSSUM',pdum, code=104, klev=1, p4=p4, &
           ldims=(/1, 2, nsp/) ,gdims=(/1, 2, nsp/), &
           longname = 'PresEQ total step to step')
!BOX
      !-------------------------------------------------------------------------
      ! fields used for restart

      CALL default_stream_setting ( tdiag, &
            lrerun=.TRUE., lpost = .FALSE., laccu = .FALSE., &
            repr=SPECTRAL, &
            ldims=(/nlev, 2, nsp/), gdims=(/nlev, 2, nsp/), &
            units='(X)', no_default=.TRUE.)

      ALLOCATE (pdovor(nlev,2,nsp,2)); pdovor(:,:,:,:) = 0.0_dp
      p4 => pdovor(:,:,:,1:1)
      CALL add_stream_element(tdiag,'PDOV1',  pdum,p4=p4, &
           longname = 'PDO Vor time level 1')
      p4 => pdovor(:,:,:,2:2)
      CALL add_stream_element(tdiag,'PDOV2',  pdum,p4=p4, &
           longname = 'PDO Vor time level 2')

      ALLOCATE (pdodiv(nlev,2,nsp,2)); pdodiv(:,:,:,:) = 0.0_dp
      p4 => pdodiv(:,:,:,1:1)
      CALL add_stream_element(tdiag,'PDOD1',  pdum,p4=p4, &
           longname = 'PDO Div time level 1')
      p4 => pdodiv(:,:,:,2:2)
      CALL add_stream_element(tdiag,'PDOD2',  pdum,p4=p4, &
           longname = 'PDO Div time level 2')

      ALLOCATE (pdotep(nlevp1,2,nsp,2)); pdotep(:,:,:,:) = 0.0_dp
      p4 => pdotep(:,:,:,1:1)
      CALL add_stream_element(tdiag,'PDOTP1', pdum,p4=p4, &
           longname = 'PDO STP time level 1',klev=nlevp1)
      p4 => pdotep(:,:,:,2:2)
      CALL add_stream_element(tdiag,'PDOTP2', pdum,p4=p4, &
           longname = 'PDO STP time level 2',klev=nlevp1)

!EOX       
      ! -----------------------------------------------------------
      ! allocate gridpoint space arrays
!BOX
       CALL default_stream_setting ( tdiag, &
            lrerun=.TRUE., repr=GAUSSIAN, &
            ldims=(/nlon, nlev, ngl/), gdims=(/nlon, nlev, ngl/), &
            no_default=.TRUE.)
!EOX
      !
      !   PTFH1   memory of timefilter
      ! updated in TF1
!BOX
      CALL add_stream_element(tdiag,'PTFH1', ptfh1, &
           longname = 'time filter buffer 1 ', &
           ldims=(/nlon, nlev,NO_PTFH,ngl/) ,gdims=(/nlon, nlev,NO_PTFH,ngl/) )
!EOX
      !           1 ... vorticity
      !           2 ... divergence
      !           3 ... temperature
      !   PTFH2   pressure
      ! updated in TF1
!BOX
      CALL add_stream_element(tdiag,'PTFH2', ptfh2, &
           longname = 'time filter buffer 2', &
           ldims=(/nlon,ngl/) ,gdims=(/nlon, ngl/))
!EOX
      !
      ! -----------------------------------------------------------
      !
      !     g .... accumulated in grid point space
      !     gs ... parts from gridpoint space, accumulated in spectral space
      !     s .... accumulated in spectral space
      !
      ! workspace for accumulation of terms in grid point space
      !
      ! **** PDIGA
      !
      !     Index
      !     ............. vorticity and divergenc equation
      !     VOM = Fu, VOL = Fv
      !
      !     Fu (VOM parts in GPC)
      ! DYN    1   coriolisterm and pressure tendency of Fu
      ! DYN    2   vertical advection of Fu
      ! PHYSC  3   diffusion due to impuls Pu VDIFF
      ! PHYSC  4   diffusion due to gravity wave drag Pu GWDRAG
      ! PHYSC  5   diffusion due to mass flux Pu CUCALL
      !
      !     Fv (VOL parts in GPC)
      ! DYN    6   coriolisterm and pressure tendency of Fv
      ! DYN    7   vertical advection of Fv
      ! PHYSC  8   diffusion due to impuls Pv VDIFF
      ! PHYSC  9   diffusion due to gravity wave drag Pv GWDRAG
      ! PHYSC  10  diffusion due to mass flux Pv CUCALL
      !
      ! DYN    11  potential and kinetic energy term G
      !
      !     ............. temperature tendency equation
      ! DYN    12  horizontal advection term
      ! DYN    13  vertical advection term
      ! DYN    14  energy conversion
      ! PHYSC  15  RADHEAT radiation tendency
      ! PHYSC  16  VDIFF turbulence
      ! PHYSC  17  GWDRAG gravity wave drag
      ! PHYSC  18  CUCALL convective condensation
      ! PHYSC  19  COND large+scale condensation
      !
      !     ............. pressure tendency equation
      ! DYN    20 level dependend divergence part, 
      !
      ! TF2    21  timefilter of vorticity
      ! TF2    22  timefilter of divergence
      ! TF2    23  timefilter of temperature
      !
      ! RADHEAT 24 longwave radiation
      ! RADHEAT 25 shortwave radiation
      !
!BOX
      CALL add_stream_element(tdiag,'PDIGA', pdiga, &
           longname = 'diagnostic term DIGA', &
           ldims=(/nlon, nlev,NO_PDIGA,ngl/) ,gdims=(/nlon, nlev,NO_PDIGA,ngl/))
      ALLOCATE (pdigaa (nlev,2,nmp1,NO_PDIGA))
      ALLOCATE (pdigas (nlev,2,nmp1,NO_PDIGA))
!EOX
      !
      ! **** PDSGA
      !
      ! one level arrays of surface pressure
      !      1     last level total integral
      !      2     timefilter for pressure
      !
!BOX
      CALL add_stream_element(tdiag,'PDSGA', pdsga, &
           longname = 'diagnostic term DSGA', &
           ldims=(/nlon,NO_PDSGA,ngl/) ,gdims=(/nlon,NO_PDSGA,ngl/))
      ALLOCATE (pdsgaa (2,nmp1,NO_PDSGA))
      ALLOCATE (pdsgas (2,nmp1,NO_PDSGA))
!EOX
      !
      ! ***** PDIGB
      !
      !     array for d/dlambda derivatives calculated in fourierspace
      !     corresponding to pdiga(*,1...10,*)
      !
!BOX
      ALLOCATE (pdigb    (nlon,nlev,NO_PDIGB,ngl))
      ALLOCATE (pdigba (nlev,2,nmp1,NO_PDIGB))
      ALLOCATE (pdigbs (nlev,2,nmp1,NO_PDIGB))
!EOX
      !
      ! ***** PDIGS
      !
      !     local memory buffer for accumulation of semi-implicit parts
      !     fraction solved in grid point space
      !
      !    1 vorticity and implicit and explicit part (L)
      !    2 divergence
      !    3 temperatur and pressure
      !    4 implicit part of vorticity (M) used in Fourierspace
!BOX
      ALLOCATE (pdigs    (nlon,nlev,NO_PDIGS,ngl))
      ALLOCATE (pdigsa (nlev,2,nmp1,NO_PDIGS))
      ALLOCATE (pdigss (nlev,2,nmp1,NO_PDIGS))
!EOX
      !
      ! ****** PDSGS
      !
      !     semi-implicit part of log surface pressure
!BOX
      ALLOCATE (pdsgs    (nlon,ngl))
      ALLOCATE (pdsgsa (2,nmp1))
      ALLOCATE (pdsgss (2,nmp1))
      
      CALL message('mo_diag_tendency:DIAG_Init','  ------- end of initialization -----')
!EOX
    CASE (IDIAG_FREE)        ! get memory free
!BOX
      DEALLOCATE (pdvor, pddiv, pdtem, p4prs)
      DEALLOCATE (pdovor, pdodiv, pdotep)
      DEALLOCATE        (pdigaa, pdigas)
      DEALLOCATE        (pdsgaa, pdsgas)
      DEALLOCATE (pdigb, pdigba, pdigbs)
      DEALLOCATE (pdigs, pdigsa, pdigss)
      DEALLOCATE (pdsgs, pdsgsa, pdsgss)
!EOX
    CASE (IDIAG_INI_GBUF)    ! reset local accumulation buffer
!BOX
      pdigb (:,:,:,:) = 0.0_dp
      pdigs (:,:,:,:) = 0.0_dp
      pdsgs (:,:)     = 0.0_dp
!EOX
    CASE (IDIAG_INI_PREV)    ! prepare reference values
!BOX
      ! detect restart fields using the global mean of temperature
      lpdo1 = pdotep(1,1,1,1) > 2._dp*EPSILON(1.0_dp)
      lpdo2 = pdotep(1,1,1,2) > 2._dp*EPSILON(1.0_dp)

      ! calculate tendencies
      IF (lpdo1) THEN
        ! mean is correct, also for the first call
        IF (ldiag_start) laccu_ok = .TRUE.

        ! memory is filled, the tendency can be calculated as [X(t+1) - X(t-1)]
        pdvor(:,:,:,NDVOR) = pdvor(:,:,:,NDVOR)+(svo(:,:,:)      - pdovor(:,:,:,1))
        pddiv(:,:,:,NDDIV) = pddiv(:,:,:,NDDIV)+(sd (:,:,:)      - pdodiv(:,:,:,1))
        pdtem(:,:,:,NDTEM) = pdtem(:,:,:,NDTEM)+(stp(1:nlev,:,:) - pdotep(1:nlev,:,:,1))
        pdprs(  :,:,NDPRS) = pdprs(  :,:,NDPRS)+(stp(nlevp1,:,:) - pdotep(nlevp1,:,:,1))
      END IF

      IF (lpdo2) THEN
        ! move (t) into (t-1) memory
        pdovor(:,:,:,1) = pdovor(:,:,:,2)
        pdodiv(:,:,:,1) = pdodiv(:,:,:,2)
        pdotep(:,:,:,1) = pdotep(:,:,:,2)
      END IF

      ! store (t+1) for next integration loop
      ! it is the unfiltered valu ein the spectral space
      ! in TF1 the same field is available at grid points
      pdovor(:,:,:,2) = svo(:,:,:)
      pdodiv(:,:,:,2) = sd (:,:,:)
      pdotep(:,:,:,2) = stp(:,:,:)

      ldiag_start = .FALSE.  ! switch after the first pass

    CASE default
      WRITE (diag_mess,*) 'type not implemented, ITYPE= ',itype
      CALL finish('mo_diag_tendency:DIAG_Init',diag_mess)

    END SELECT
    
  END SUBROUTINE DIAG_Init
!EOX
!EOC
!======================================================================
!BOP
  ! !IROUTINE:  DIAG_SpecTrans
  ! !INTERFACE:

  SUBROUTINE DIAG_SpecTrans

    ! !DESCRIPTION: 
    ! second part of spectral transform
    !
    ! insert in SCAN1SL after CALL SI2

    ! !USES:
    USE mo_legendre,      ONLY: legmod
    USE mo_control,       ONLY: nhgl
!EOP
!BOC
!BOX
    INTEGER :: i

    north_south_loop: DO i  = 1, nhgl

      CALL DIAG_sym1(i)

      CALL legmod   (i)
      CALL DIAG_ltd

    END DO north_south_loop

  END SUBROUTINE DIAG_SpecTrans
!EOX
!EOC
!======================================================================
!BOP
  ! !IROUTINE:  DIAG_fftd 
  ! !INTERFACE:
  SUBROUTINE DIAG_fftd

    ! !DESCRIPTION: 
    ! transform the diagnostic arrays into fourierspace
    !
    ! insert in SCAN1SL after CALL FFTD

    ! !USES:
    USE mo_time_control,  ONLY: l_putdata
#ifdef MKL_DFT
    USE mo_mkl_dft,       ONLY: mkl_dft
#elif ESSL_DFT
    USE mo_essl_dft,      ONLY: essl_dft   
#else
    USE mo_fft992,        ONLY: fft992
#endif
    USE mo_control,       ONLY: nlp2, nlon, nlev, ngl
    USE mo_spectral,      ONLY: compress, expand
!EOP
!BOC
!BOX
    INTEGER, PARAMETER :: inc = 1, isign = -1
    INTEGER            :: iamount
    REAL(dp), TARGET :: zinp(nlp2*(nlev*NO_PDIGA+NO_PDSGA)*ngl) ! work space for transformation
    REAL(dp) , POINTER :: f1d(:), f4d(:,:,:,:)

    ! transform 3-dim fields semi-implicit parts counted in grid space
    iamount = nlev*NO_PDIGS_FRAC*ngl
    f1d => zinp( (iamount*nlp2+1): )
    f4d => pdigs(:,:,1:NO_PDIGS_FRAC,:)
    CALL expand  (zinp, f4d  ,2)
    CALL expand  (f1d,  pdsgs,2)
#ifdef MKL_DFT
    CALL mkl_dft(zinp,inc,nlp2,nlon,iamount+ngl,isign)
#elif ESSL_DFT
    CALL essl_dft(zinp,nlon,nlp2,iamount+ngl,isign)
#else
    CALL fft992(zinp,inc,nlp2,nlon,iamount+ngl,isign)
#endif
    CALL compress(zinp, f4d  ,2)
    CALL compress(f1d,  pdsgs,2)

    ! transform other terms during output time step
    IF (l_putdata(dio_index)) THEN
      iamount = nlev*NO_PDIGA*ngl
      f1d => zinp( (iamount*nlp2+1): )
      CALL expand  (zinp, pdiga ,2)
      CALL expand  (f1d,  pdsga ,2)
#ifdef MKL_DFT
      CALL mkl_dft(zinp,inc,nlp2,nlon,iamount+NO_PDSGA*ngl,isign)
#elif ESSL_DFT
      CALL essl_dft(zinp,nlon,nlp2,iamount+NO_PDSGA*ngl,isign)
#else
      CALL fft992(zinp,inc,nlp2,nlon,iamount+NO_PDSGA*ngl,isign)
#endif
      CALL compress(zinp, pdiga ,2)
      CALL compress(f1d,  pdsga ,2)
    ENDIF

  END SUBROUTINE DIAG_fftd
!EOX
!EOC
!======================================================================
!BOP
  ! !IROUTINE:  DIAG_sym1
  ! !INTERFACE:

  SUBROUTINE DIAG_sym1(ihrow)

    ! !DESCRIPTION: 
    ! separate the diagnostics arrays into symetric and asymetric part
    !
    ! insert in SCAN1SL after CALL SYM1

    ! !USES:
    USE mo_time_control,  ONLY: l_putdata
    USE mo_control,       ONLY: nmp1, nlev, ngl

    INTEGER, INTENT(in) :: ihrow   ! latitude index
!EOP
!BOC
!BOX
    INTEGER :: jl, jm, irow_n, irow_s

    !     even and odd components
    irow_n = ihrow
    irow_s = ngl + 1 - ihrow

    spec_loop1: DO jm = 1,nmp1

      level_loop1: DO jl = 1,nlev

        pdigss(jl,:,jm,:) = 0.5_dp*(pdigs(2*jm-1:2*jm,jl,:,irow_n) + pdigs(2*jm-1:2*jm,jl,:,irow_s))
        pdigsa(jl,:,jm,:) = 0.5_dp*(pdigs(2*jm-1:2*jm,jl,:,irow_n) - pdigs(2*jm-1:2*jm,jl,:,irow_s))

      ENDDO level_loop1

      pdsgss(:,jm) = 0.5_dp*(pdsgs(2*jm-1:2*jm,irow_n) + pdsgs(2*jm-1:2*jm,irow_s))
      pdsgsa(:,jm) = 0.5_dp*(pdsgs(2*jm-1:2*jm,irow_n) - pdsgs(2*jm-1:2*jm,irow_s))

    ENDDO spec_loop1

    ! transform only during output step
    IF (l_putdata(dio_index)) THEN

      spec_loop2: DO jm = 1,nmp1

        level_loop2: DO jl = 1,nlev

          pdigas(jl,:,jm,:) = 0.5_dp*(pdiga(2*jm-1:2*jm,jl,:,irow_n) + pdiga(2*jm-1:2*jm,jl,:,irow_s))
          pdigaa(jl,:,jm,:) = 0.5_dp*(pdiga(2*jm-1:2*jm,jl,:,irow_n) - pdiga(2*jm-1:2*jm,jl,:,irow_s))
          pdigbs(jl,:,jm,:) = 0.5_dp*(pdigb(2*jm-1:2*jm,jl,:,irow_n) + pdigb(2*jm-1:2*jm,jl,:,irow_s))
          pdigba(jl,:,jm,:) = 0.5_dp*(pdigb(2*jm-1:2*jm,jl,:,irow_n) - pdigb(2*jm-1:2*jm,jl,:,irow_s))

        ENDDO level_loop2

        pdsgas(:,jm,:) = 0.5_dp*(pdsga(2*jm-1:2*jm,:,irow_n) + pdsga(2*jm-1:2*jm,:,irow_s))
        pdsgaa(:,jm,:) = 0.5_dp*(pdsga(2*jm-1:2*jm,:,irow_n) - pdsga(2*jm-1:2*jm,:,irow_s))

      ENDDO spec_loop2

    ENDIF

  END SUBROUTINE DIAG_sym1
!EOX
!EOC
!======================================================================
!BOP
  ! !IROUTINE:  DIAG_ltd 
  ! !INTERFACE:
  SUBROUTINE DIAG_ltd

    ! !DESCRIPTION: 
    ! perform legendre transform for diagnostic arrays
    !
    ! insert in SCAN1SL after CALL LTD

    ! !USES:
    USE mo_time_control,  ONLY: l_putdata
    USE mo_truncation,    ONLY: nmp, nnp
    USE mo_legendre,      ONLY: pnmd, anmd, rnmd
    USE mo_control,       ONLY: nmp1
!EOP
!BOC
!BOX
    INTEGER :: jm, ims, ins, is, jn, jh, iu

    north_south_loop: DO jh = 1,2 ! 1: north, 2:south
      iu = 2-jh

      spec_loop: DO jm = 1,nmp1
        ims = nmp(jm)-iu
        ins = nnp(jm)+iu

        sym_loop: DO jn=2,ins,2
          is = ims+jn
          IF (jh == 1) THEN   !     calculations for northern hemisphere
            !
            ! semi-implicit parts of vorticity
            pdvor(:,:,is, 7) = pdvor(:,:,is, 7) + pdigss(:,:,jm,1)*pnmd(is) &
                                                - pdigsa(:,:,jm,4)*anmd(is)

            ! explicit part of divergence, temperature and pressure
            pddiv(:,:,is, 7) = pddiv(:,:,is, 7) + pdigss(:,:,jm,2)*rnmd(is)
            pdtem(:,:,is,10) = pdtem(:,:,is,10) + pdigss(:,:,jm,3)*pnmd(is)
            pdprs  (:,is, 3) = pdprs  (:,is, 3) + pdsgss  (:,jm  )*pnmd(is)
            !
            IF (l_putdata(dio_index)) THEN
              ! dynamic and physical tendencies
              ! remark: ANMD contains the minus sign
              ! vorticity
              pdvor(:,:,is,1:5) = pdvor (:,:,is,1:5) + pdigbs(:,:,jm,6:10)*pnmd(is) &
                                                     + pdigaa(:,:,jm,1:5 )*anmd(is)
              pdvor(:,:,is,  6) = pdvor (:,:,is,  6) + pdigas(:,:,jm,21  )*pnmd(is)
              !
              ! divergence
              pddiv(:,:,is,1:5) = pddiv (:,:,is,1:5) + pdigbs(:,:,jm,1:5 )*pnmd(is) &
                                                     - pdigaa(:,:,jm,6:10)*anmd(is)
              pddiv(:,:,is,  6) = pddiv (:,:,is,  6) + pdigas(:,:,jm,22  )*pnmd(is)
              ! laplacian operation with G-term
              pddiv(:,:,is,  1) = pddiv (:,:,is,  1) + pdigas(:,:,jm,11  )*rnmd(is)
              !
              ! temperature
              pdtem(:,:,is, 1: 8) = pdtem(:,:,is, 1: 8) + pdigas(:,:,jm,12:19)*pnmd(is)
              pdtem(:,:,is,    9) = pdtem(:,:,is,    9) + pdigas(:,:,jm,23   )*pnmd(is)
              pdtem(:,:,is,12:13) = pdtem(:,:,is,12:13) + pdigas(:,:,jm,24:25)*pnmd(is)
              !
              ! pressure
              pdprl(:,:,is)     = pdprl(:,:,is)      + pdigas(:,:,jm,20)*pnmd(is)
              pdprs  (:,is,1:2) = pdprs  (:,is,1:2)  + pdsgas  (:,jm, :)*pnmd(is)
            ENDIF
            !
          ELSE                ! calculations for southern hemisphere
            !
            ! semi-implicit parts of vorticity
            pdvor(:,:,is, 7) = pdvor (:,:,is, 7) + pdigsa(:,:,jm,1)*pnmd(is) &
                                                 - pdigss(:,:,jm,4)*anmd(is)
            ! explicit part divergence, temperature and pressure
            pddiv(:,:,is, 7) = pddiv(:,:,is, 7) + pdigsa(:,:,jm,2)*rnmd(is)
            pdtem(:,:,is,10) = pdtem(:,:,is,10) + pdigsa(:,:,jm,3)*pnmd(is)
            pdprs  (:,is, 3) = pdprs  (:,is, 3) + pdsgsa  (:,jm  )*pnmd(is)
            !
            IF (l_putdata(dio_index)) THEN
              ! dynamic and physical tendencies
              ! vorticity
              pdvor(:,:,is,1:5) = pdvor(:,:,is,1:5) + pdigba(:,:,jm,6:10)*pnmd(is) &
                                                    + pdigas(:,:,jm,1:5 )*anmd(is)
              pdvor(:,:,is,  6) = pdvor(:,:,is,  6) + pdigaa(:,:,jm,21  )*pnmd(is)
              !
              ! divergence
              pddiv(:,:,is,1:5) = pddiv(:,:,is,1:5) + pdigba(:,:,jm,1:5 )*pnmd(is) &
                                                    - pdigas(:,:,jm,6:10)*anmd(is)
              pddiv(:,:,is,  6) = pddiv(:,:,is,  6) + pdigaa(:,:,jm,22  )*pnmd(is)
              pddiv(:,:,is,  1) = pddiv(:,:,is,  1) + pdigaa(:,:,jm,11  )*rnmd(is)
              !
              ! temperature
              pdtem(:,:,is, 1: 8) = pdtem(:,:,is, 1: 8) + pdigaa(:,:,jm,12:19)*pnmd(is)
              pdtem(:,:,is,    9) = pdtem(:,:,is,    9) + pdigaa(:,:,jm,23   )*pnmd(is)
              pdtem(:,:,is,12:13) = pdtem(:,:,is,12:13) + pdigaa(:,:,jm,24:25)*pnmd(is)
              !
              ! pressure
              pdprl(:,:,is)     = pdprl(:,:,is)     + pdigaa(:,:,jm,20)*pnmd(is)
              pdprs  (:,is,1:2) = pdprs  (:,is,1:2) + pdsgaa  (:,jm, :)*pnmd(is)
            ENDIF
            !
          ENDIF

        ENDDO sym_loop

      ENDDO spec_loop

     ENDDO north_south_loop

  END SUBROUTINE DIAG_ltd
!EOX
!EOC
!======================================================================
!BOP
  ! !IROUTINE:  DIAG_Write
  ! !INTERFACE:

  SUBROUTINE DIAG_Write

    ! !DESCRIPTION: 
    ! correction of output buffer, reset local buffer
!EOP
!BOC
!BOX

    IF (laccu_ok) THEN
      CALL DIAG_Check
    ELSE
      pdvor(:,:,:,NDVOR) = 0.0_dp
      pddiv(:,:,:,NDDIV) = 0.0_dp
      pdtem(:,:,:,NDTEM) = 0.0_dp
      pdprs(  :,:,NDPRS) = 0.0_dp
    END IF

    ! correction of all terms due to the leap frog scheme
    pdvor(:,:,:,:) = 0.5_dp*pdvor(:,:,:,:)
    pddiv(:,:,:,:) = 0.5_dp*pddiv(:,:,:,:)
    pdtem(:,:,:,:) = 0.5_dp*pdtem(:,:,:,:)
    pdprl(:,:,:)   = 0.5_dp*pdprl(:,:,:)
    p4prs(:,:,:,:) = 0.5_dp*p4prs(:,:,:,:)

    ! reset accumulated grid point arrays
    pdiga(:,:,:,:) = 0.0_dp
    pdsga(:,:,:)   = 0.0_dp

    ! for the next postprocessing cycle the accumulation should be fine
    IF ( pdotep(1,1,1,1) > 2._dp*EPSILON(1.0_dp) ) laccu_ok = .TRUE.

  END SUBROUTINE DIAG_Write
!EOX
!EOC
!======================================================================
!BOP
  ! !IROUTINE:  DIAG_Check
  ! !INTERFACE:

  SUBROUTINE DIAG_Check
    ! !DESCRIPTION: 
    ! The procedure diagnoses the tendency calculation. The sum of all terms
    ! is compared with the total tendency. For the first accumulation interval
    ! the correlation can not be used. But for all following intervals all 
    ! correlations must be 1.00, except for nudging mode. In nudging mode
    ! the nudging term will not be put into the account, therefore the
    ! diagnostics are not complete.
    !
    
    ! !USES:
    USE mo_control,  ONLY: nsp, nlev
    USE mo_spectral, ONLY: corrsp
!EOP
!BOC
!BOX
    REAL(dp)          :: ccv, ccd, cct, ccp
    REAL(dp), POINTER :: p3(:,:,:), p2a(:,:), p2b(:,:)
    REAL(dp), TARGET  :: wrk(nlev,2,nsp)
    INTEGER           :: i

    wrk = 0.0_dp
    DO i=1,NDVOR-1
       wrk(:,:,:) = wrk(:,:,:) + pdvor(:,:,:,i)
    END DO
    p3 => pdvor(:,:,:,NDVOR)
    ccv = corrsp (wrk,p3)

    wrk = 0.0_dp
    DO i=1,NDDIV-1
       wrk(:,:,:) = wrk(:,:,:) + pddiv(:,:,:,i)
    END DO
    p3 => pddiv(:,:,:,NDDIV)
    ccd = corrsp (wrk,p3)

    wrk = 0.0_dp
    DO i=1,NDTEM-3
       wrk(:,:,:) = wrk(:,:,:) + pdtem(:,:,:,i)
    END DO
    p3 => pdtem(:,:,:,NDTEM)
    cct = corrsp (wrk,p3)

    wrk = 0.0_dp
    DO i=1,NDPRS-1
       wrk(1,:,:) = wrk(1,:,:) + pdprs(:,:,i)
    END DO
    p2a => wrk(1,:,:)
    p2b => pdprs(:,:,NDPRS)
    ccp = corrsp (p2a,p2b)

    WRITE(diag_mess,'(4(a,f8.5))') 'VOR ',ccv,' DIV ',ccd,' TEM ',cct,' PRS ',ccp
    CALL message('mo_diag_tendency:DIAG_Check',diag_mess)

  END SUBROUTINE DIAG_Check
!EOX
!EOC
END MODULE mo_diag_tendency
