!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_nmi
!BOP
  ! !MODULE: mo_nmi

  ! !DESCRIPTION: 
  ! nonlinear normal mode functions, in nudging mode the reference data will be
  ! optional filtered using normal modes, without nudging the initialization
  ! procedure is acting
  !
  !\begin{verbatim}
  ! ******************** Interface to ECHAM **********************************
  !
  ! *ECHAM*       *NMI*
  !
  ! CONTROL ----> NMI_Init
  !                 +---> NMI_Normal_Modes
  !                         +---> NMI_Horizontal_Modes
  !
  ! STEPON -----> NMI_Make(NMI_MAKE_NMI) first call at beginning of time step
  !
  !               modify PHYSC
  !
  !        -----> NMI_Make(NMI_MAKE_NMI) second call after hdiff
  !
  ! Nudging ----> NMI_Make(NMI_MAKE_FOD)
  !                 +---> NMI_LoadBuffer
  !                 +---> NMI_Filter
  !\end{verbatim}
!BOX
  ! development steps
  ! 4.0     implementation of the initial part of the nmi
  ! 3.0     revision march 2000
  ! 1.1     revision
  ! 1.0     final version with explcit filtering only
  ! 0.10    last checks global average
  ! 0.9     new methode for separation, Lagrangian multiplier
  ! 0.8     check of separation of H in T and lnPS
  ! 0.4     filtering in MNI_Filter
  ! 0.3     calculation of modes checked
  !         no tides, no initialisation, no diabtic tendencies
!EOX
  ! !REVISION HISTORY: 

  ! W. Wergen, ECMWF, July 1982 - December 1982
  ! W. Heckley, ECMWF, July 1986 - November 1988
  ! C. Temperton, ECMWF, January 1989 - May 1991
  ! L. Kornblueh, MPI, June 1998
  ! I. Kirchner, MPI, August 1998 - January 1999
  ! I. Kirchner, MPI, October 1999, March 2000, March 2002
  ! R. Johanni, IPP Garching, May-2002, parallel version
  ! I. Kirchner, MPI, August 2002, revision

  ! !USES:
  USE mo_kind,            ONLY: dp
  USE mo_exception,       ONLY: message_text, message, finish
  USE mo_time_conversion, ONLY: time_days
  USE mo_linked_list,     ONLY: t_stream
  USE mo_decomposition, ONLY: ldc => local_decomposition

!BOX
  IMPLICIT NONE

  PRIVATE
!EOX
  ! !PUBLIC MEMBER FUNCTIONS:
  PUBLIC  :: NMI_Init
  PUBLIC  :: NMI_Close
  PUBLIC  :: NMI_Make
!BOX
  PRIVATE :: NMI_Normal_Modes
  PRIVATE :: NMI_Horizontal_Modes
  PRIVATE :: NMI_Filter
  PRIVATE :: unitmx
  PRIVATE :: df, bf, cf, ef, hf, of, pf, qf 
  PRIVATE :: NMI_LoadBuffer
!EOX
!EOP

  ! define the filter characteristics
  ! normal filter, used for nudging
  INTEGER, ALLOCATABLE :: iplim(:,:,:) ! cut off limit for fast modes
  INTEGER              :: nvm          ! number of vertical modes included
  ! diabatic tendencies, used for initialization
  INTEGER, ALLOCATABLE :: iplimd(:,:,:)! cut off limit for fast modes
  INTEGER              :: nvmd         ! number of vertical modes included

  TYPE(time_days),SAVE :: nmi_start_date

  ! accumulation of physical tendencies and iteration
  ! NMI_Make is called twice, the first call is at the beginning
  ! of the time step and the second call is after HDIFF, both implemented
  ! in STEPON, the phase is controlled by the following variables
  !
  INTEGER, PARAMETER         :: NMI_NULL       =  0 ! inactive phase
  INTEGER, PARAMETER         :: NMI_PRE_OPEN   =  1 ! pre integration phase starts
  INTEGER, PARAMETER         :: NMI_PRE        =  2 ! pre integration phase, no accumlation of diabatic heat source
  INTEGER, PARAMETER         :: NMI_ACCU_OPEN  =  3 ! start accumulation of heat source
  INTEGER, PARAMETER, PUBLIC :: NMI_ACCU       =  4 ! during accumulation of heat source
  INTEGER, PARAMETER         :: NMI_ACCU_CLOSE =  5 ! close heat source accumulation and average heat source
  INTEGER, PARAMETER, PUBLIC :: NMI_USE_AVG    =  6 ! use average heat source
  INTEGER, PARAMETER         :: NMI_RESTART    =  7 ! restart model with balanced fields
  INTEGER, PARAMETER         :: NMI_PRINT_NEXT =  8 ! prepare for print out
  INTEGER, PARAMETER         :: NMI_PRINT      =  9 ! nmi diagnostics print out and close nmi

  INTEGER, PARAMETER :: NMI_NDGFILT    = 20 ! use nmi only as filter

  INTEGER, PUBLIC  :: nmi_phase      = NMI_NULL  ! hold NMI state during run
  INTEGER  :: nmi_counter    = 0

  TYPE(t_stream), POINTER :: stream_nmi

  ! memory for diabatic accumulation
  REAL(dp), POINTER, DIMENSION(:,:,:), PUBLIC :: dh_t, dh_m, dh_l ! grid fields: long, lev, lat
  REAL(dp), POINTER, DIMENSION(:,:,:)         :: buf_vo, buf_d    ! spectral fields: nlev, 2, nsp
  REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET :: buf_tp
  REAL(dp), POINTER, DIMENSION(:,:), PUBLIC   :: buf_t, buf_m, buf_l ! workspace grid buffer

  ! internal used parameters and arrays
  ! vertical modes
  REAL(dp), ALLOCATABLE :: zdb(:), &! delpr/b
                      phib  (:),   &! phibar, equivalent depths
                      vm_ev (:,:), &! m (vertical modes) eigenvectors
                      vm_evi(:,:)   ! inverse(m)

  ! horizontal modes
  INTEGER  :: ing, &                      ! number of gravity modes
              inr                         ! number of rossby modes
  REAL(dp), ALLOCATABLE :: zhmod(:)       ! eigenvectors and eigenvalues
  INTEGER, ALLOCATABLE  :: ihmod(:,:,:,:) ! pointer array for horizontal modes
                                          !  (1,.,.,.) block starting pointer
                                          !  (2,.,.,.) number of Gravity modes (east+west)
                                          !  (3,.,.,.) number of Rossby modes
                                          !  (.,X,.,.) 1 .. symmetric, 2 .. antisymmetric part
                                          !  (.,.,X,.) zonal wave number
                                          !  (.,.,.,X) vertical mode

  ! work space for NMI calculation
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: &
       svo_nmi, sd_nmi, stp_nmi, &  ! buffer A absolut
       svo_dt,  sd_dt,  stp_dt      ! buffer B tendencies

  CHARACTER (len=256) :: nmi_mess

  ! options for NMI_make
  INTEGER, PARAMETER         :: NMI_MAKE_TEST=0
  INTEGER, PARAMETER, PUBLIC :: NMI_MAKE_FOD =1   ! Filter Observed Data
  INTEGER, PARAMETER, PUBLIC :: NMI_MAKE_FMMO=3   ! Filter Model Minus Observed
  INTEGER, PARAMETER, PUBLIC :: NMI_MAKE_NMI =5   ! Normal Mode Initialization


  ! options for data buffer handling
  INTEGER, PARAMETER :: NMI_LOAD_CLEAR   =0  ! clear buffer
  INTEGER, PARAMETER :: NMI_LOAD_ECH_BUFA=1  ! echam data into buffer A
  INTEGER, PARAMETER :: NMI_LOAD_ECH_BUFB=2  ! echam data into buffer B
  INTEGER, PARAMETER :: NMI_LOAD_ORD_BUFA=3  ! observed read buffer into buffer A
  INTEGER, PARAMETER :: NMI_LOAD_ORD_BUFB=4  ! observed read buffer into buffer B
  INTEGER, PARAMETER :: NMI_LOAD_OIN_BUFA=5  ! interpolated data into buffer A
  INTEGER, PARAMETER :: NMI_LOAD_OIN_BUFB=6  ! interpolated data into buffer B
  INTEGER, PARAMETER :: NMI_CORREL_BUFFER=7  ! correlate buffer

  INTEGER, PARAMETER :: NMI_COPY_AB   =10  ! exchange buffer
  INTEGER, PARAMETER :: NMI_COPY_BA   =15
  INTEGER, PARAMETER :: NMI_CALC_TEND =20  ! calculate tendencies
  INTEGER, PARAMETER :: NMI_CALC_DIFF =25  ! difference A minus B
  INTEGER, PARAMETER :: NMI_CALC_ADD  =30  ! add both buffer

  INTEGER, PARAMETER :: NMI_STORE_FAST=50  ! store the fast mode part


  ! filter options, type of mode separation
  INTEGER, PARAMETER :: NMI_SEL_ALL     = 0   ! use all modes without global average
  INTEGER, PARAMETER :: NMI_SEL_ALL_W0  = 1   ! use all modes, result includes global average
  INTEGER, PARAMETER :: NMI_SEL_GM      = 2   ! use GW modes without global average
  INTEGER, PARAMETER :: NMI_SEL_GM_W0   = 3   ! use GW modes, result includes global average
  INTEGER, PARAMETER :: NMI_SEL_RM      = 4   ! use RW modes without global average
  INTEGER, PARAMETER :: NMI_SEL_RM_IMPL = 5   ! subtract GW modes, result includes global average
  INTEGER, PARAMETER :: NMI_SEL_RM_W0   = 6   ! use RW modes, result includes global average
  INTEGER, PARAMETER :: NMI_SEL_FILT    = 7   ! GW modes without global average inside window
  INTEGER, PARAMETER :: NMI_SEL_FILT_W0 = 8   ! GW modes partly, result includes global average, inside window
  INTEGER, PARAMETER :: NMI_SEL_FILTD   = 9   ! GW modes without global average, inside second window
  INTEGER, PARAMETER :: NMI_SEL_FILTD_W0=10   ! GW modes partly, result includes global average

  LOGICAL, SAVE :: nmi_mem = .FALSE.

  ! min and max values are experimental (I.K. 4-mar-02)
  INTEGER, PARAMETER :: NTPRE_MIN  = 2, NTPRE_MAX  = 20  ! 2 ... 20 valid
  INTEGER, PARAMETER :: NTDIA_MIN  = 8, NTDIA_MAX  = 20  ! 8 ... 20 valid
  INTEGER, PARAMETER :: NTITER_MIN = 2, NTITER_MAX = 10  ! 2 ... 10 valid

  INTEGER  :: ntpre  = NTPRE_MIN   ! length of pre-integration interval
  INTEGER  :: ntdia  = NTDIA_MIN   ! length of accumulation interval (diabatic tendencies)
  INTEGER  :: ntiter = NTITER_MIN  ! length of iteration interval
  
  ! separate between fast and slow modes using a period
  REAL(dp) :: pcut  = 12.0_dp ! cut off period in hours, used for nudging
  REAL(dp) :: pcutd = 6.0_dp  ! cut off period in hours, used for initialization

  LOGICAL, PUBLIC  :: lnmi_cloud = .TRUE.  ! run initialization with clouds

  LOGICAL, PUBLIC  :: lnmi_run   = .FALSE. ! detect initialisation mode

  ! To be more memory friendly set reduce_memory_usage to .true.
  LOGICAL, PARAMETER :: reduce_memory_usage = .TRUE.

  ! Decide if we have to keep all log-messagess from root (p_pe==0) in the same way
  ! as with reduce_memory_usage=.FALSE., even if there is no work for root:
  LOGICAL, PARAMETER :: keep_root_messages = .TRUE.

  ! management of memory requirements:
  LOGICAL :: zero_nsm = .FALSE.
  LOGICAL, ALLOCATABLE :: have_wave(:)
  INTEGER, ALLOCATABLE :: zhmod_wave_mem(:)
 
CONTAINS

!======================================================================
!BOP
  ! !IROUTINE: NMI_Init
  ! !INTERFACE:

  SUBROUTINE NMI_Init(filter_mode)

    ! !DESCRIPTION:
    ! initialization of NMI for ECHAM

    ! !USES:
    USE mo_namelist,        ONLY: open_nml, position_nml, POSITIONED, MISSING
    USE mo_output,          ONLY: nudging_table
    USE mo_memory_base,     ONLY: new_stream, default_stream_setting, &
                                  add_stream_element, SPECTRAL, GAUSSIAN, GRIB, SURFACE
    USE mo_time_control,    ONLY: current_date, next_date, delta_time, write_date
    USE mo_time_conversion, ONLY: TC_set, TC_convert, time_native, OPERATOR(<)
    USE mo_truncation,      ONLY: nnp
    USE mo_mpi,             ONLY: p_parallel, p_parallel_io, p_bcast, p_io, p_pe
    USE mo_control,         ONLY: nmp1, nk, nlev, nlevp1, nsp, nm
!BOX
    IMPLICIT NONE
!EOX
    ! !INPUT PARAMETERS:
    LOGICAL, INTENT(in) :: filter_mode   ! true for nudging mode
!BOX
    TYPE(time_native) :: my_date
    INTEGER :: imem, ii, imo, idim, ierr, snsp, idia, isec, inml, iunit, im, jm, ipos, ipos_extra
    REAL(dp), POINTER :: p3(:,:,:), p4(:,:,:,:)

    INTEGER  :: dt_nmi_start(6) = 0   ! start date of NMI procedure

!EOX
!EOP
!BOP
!EOP

    INCLUDE 'nmictl.inc'

    snsp = ldc%snsp

    CALL message('','NMI Version based on E5R1.04 (kirchner@dkrz.de)')

    IF (p_parallel_io) THEN
      inml = open_nml('namelist.echam')
      iunit = position_nml ('NMICTL', inml, status=ierr)
      SELECT CASE (ierr)
      CASE (POSITIONED) ; READ (iunit,nmictl)
      CASE (MISSING)    ; CALL message('','no namelist NMICTL found, use defaults')
      END SELECT
    END IF
    IF (p_parallel) THEN
      CALL p_bcast (ntpre,        p_io)
      CALL p_bcast (ntdia,        p_io)
      CALL p_bcast (ntiter,       p_io)
      CALL p_bcast (pcut,         p_io)
      CALL p_bcast (pcutd,        p_io)
      CALL p_bcast (dt_nmi_start, p_io)
      CALL p_bcast (lnmi_cloud,   p_io)
    END IF

    ! general filter characteristic
    nvm     = nlev      ! maximal number of vertical modes normal filter
    nvmd    = nlev      ! maximal number of vertical modes diabatic filter

    ! check switches for initialization process
    ntpre  = MIN(MAX(ntpre, NTPRE_MIN), NTPRE_MAX)
    ntdia  = MIN(MAX(ntdia, NTDIA_MIN), NTDIA_MAX)
    ntiter = MIN(MAX(ntiter,NTITER_MIN),NTITER_MAX)

    IF (nk >= 106) THEN
       ! For resolutions greater or equal to t106 make the default
       ! number of timesteps used in the accumulation of physical
       ! tendencies for the initialisation equal to 2 hours as a
       ! function of the timestep
       idia = NINT(REAL(7200,dp)/delta_time+0.499_dp)+2
       IF (ntdia < idia) &
            CALL message('','Accumulation of diabatic heating less than 2 hours')
    ENDIF

    ! correct negative values
    pcut   = MAX(pcut, 0.0_dp)
    pcutd  = MAX(pcutd,0.0_dp)

    IF (reduce_memory_usage) THEN
      zero_nsm = (ldc%nsm == 0)
      IF (p_pe == 0 .AND. keep_root_messages) zero_nsm = .FALSE.
    ENDIF
   
    IF (.NOT. nmi_mem) THEN

      ! allocate space for vertical modes
      ALLOCATE (zdb (nlev), phib(nlev), vm_ev (nlev,nlev), vm_evi(nlev,nlev))
      zdb (:)   = 0.0_dp
      phib(:)   = 0.0_dp
      vm_ev (:,:) = 0.0_dp
      vm_evi(:,:) = 0.0_dp

      ! allocate space for horizontal mode index table
      ALLOCATE (ihmod(3,2,nmp1,nlev))
      ihmod(:,:,:,:) = 0

      ! estimate the memory for horizontal modes
      imem = 0
      DO ii = 1,nmp1
        imo  = 3*nnp(ii)
        imem = imem + (1+imo)*imo
      ENDDO

      IF (reduce_memory_usage) THEN
        IF (.NOT. zero_nsm) THEN
          ALLOCATE(have_wave(0:nm))
          have_wave = .FALSE.
          DO im = 1,ldc%nsm
            jm = ldc%sm(im)
            have_wave(jm) = .TRUE.
          ENDDO
          ALLOCATE(zhmod_wave_mem(0:nm))
          ipos = 0
          ipos_extra = 0
          DO jm = 0,nm
            zhmod_wave_mem(jm) = jm_to_zhmod_memory_estimation(jm)
            IF (have_wave(jm)) THEN
              ipos = ipos + zhmod_wave_mem(jm)
            ELSE
              ipos_extra = MAX(ipos_extra, zhmod_wave_mem(jm))
            ENDIF
          ENDDO
          ipos = ipos + ipos_extra ! extra buffer space
          ALLOCATE (zhmod(ipos))
          zhmod(:) = 0.0_dp
        ENDIF
      ELSE
        ALLOCATE (zhmod(nlev*imem))
        zhmod(:) = 0.0_dp
      ENDIF
      
      ! allocate memory for cut off limit of NMI filter
      ALLOCATE (iplim(2,nmp1,nlev))
      iplim(:,:,:) = 0

      ! filter characteristic of diabatic tendencies
      ALLOCATE (iplimd(2,nmp1,nlev))
      iplimd(:,:,:) = 0

      ! allocate spectral buffer A and B
      idim = nlev  *2*snsp
      ALLOCATE (svo_nmi(idim), sd_nmi(idim), svo_dt(idim), sd_dt(idim))
      idim = nlevp1*2*snsp
      ALLOCATE (stp_nmi(idim), stp_dt(idim))
      svo_nmi(:) = 0.0_dp; svo_dt(:) = 0.0_dp
      sd_nmi (:) = 0.0_dp; sd_dt (:) = 0.0_dp
      stp_nmi(:) = 0.0_dp; stp_dt(:) = 0.0_dp

      nmi_mem = .TRUE.
      CALL message('','NMI memory allocated')

    END IF

    ! find the nmi state
    nmi_phase = NMI_NULL

    IF (filter_mode) THEN
      ! for nudging mode the NMI works as filter
      CALL message ('','NMI works as filter for nudging data only')
      nmi_phase = NMI_NDGFILT
      lnmi_run  = .FALSE.
    ELSE
      CALL message ('','NMI procedure initiated')
      lnmi_run  = .TRUE.
    END IF

    CALL message('','Optional parameters:')

    WRITE(nmi_mess,*) '   cut off period (hours)           PCUT  = ',pcut
    CALL message('',nmi_mess)

    IF (nmi_phase == NMI_NULL) THEN
      ! the following parameters are only used in initialization mode
      CALL message('','Control initialization mode with ...')
      WRITE(nmi_mess,*) '   cut off period (hours) diabatic  PCUTD = ',pcutd
      CALL message('',nmi_mess)

      IF (SUM(dt_nmi_start(:)) /= 0) THEN
        CALL TC_set(&
             dt_nmi_start(1), dt_nmi_start(2), &
             dt_nmi_start(3), dt_nmi_start(4), &
             dt_nmi_start(5), dt_nmi_start(6), my_date)
        CALL TC_convert(my_date, nmi_start_date)

      ELSE
        CALL finish('NMI_Init',&
             'start date (DT_NMI_START) for initialization procedure missing')
      END IF

      CALL message('','Start initialization at ...')
      CALL write_date(nmi_start_date,'NMI start date is ...')

      IF (nmi_start_date < current_date) THEN
        CALL write_date(next_date,'First prognostic date is ...')
        CALL finish('NMI_Init','NMI start date not in future')
      END IF

      isec = INT(ntpre*delta_time)
      WRITE(nmi_mess,*) '   start up period (steps)         NTPRE  = ',ntpre,' ',isec,' s'
      CALL message('',nmi_mess)

      isec = INT(ntdia*delta_time)
      WRITE(nmi_mess,*) '   accumulation period (steps)     NTDIA  = ',ntdia,' ',isec,' s'
      CALL message('',nmi_mess)

      isec = INT(ntiter*delta_time)
      WRITE(nmi_mess,*) '   iteration period (steps)        NTITER = ',ntiter,' ',isec,' s'
      CALL message('',nmi_mess)

      IF (lnmi_cloud) THEN
        CALL message('','   run NMI with clouds')
      ELSE
        CALL message('','   run NMI without clouds')
      END IF

    ENDIF

    CALL NMI_Normal_Modes(nlev, filter_mode)    ! Initialisation is required - establish modes

    ! allocate additional memory for initialization procedure
    IF (nmi_phase == NMI_NULL) THEN

      ! allocata an additional stream for NMI
      CALL new_stream( stream_nmi, 'nmi', lrerun=.FALSE., lpout=.TRUE., filetype=GRIB)
      CALL default_stream_setting(stream_nmi, &
           table=nudging_table, bits=24, &
           lpost  = .TRUE., lrerun = .FALSE., laccu  = .FALSE.,              &
           repr=GAUSSIAN )

      ! the diabatic foring terms
      CALL add_stream_element(stream_nmi,'DAC_T',dh_t,&
           code=201, longname='diabatic forcing temp (mean)',         units='K/s')

      CALL add_stream_element(stream_nmi,'DAC_M',dh_m,&
           code=202, longname='diabatic forcing meridional (mean)',   units='m/s^2')

      CALL add_stream_element(stream_nmi,'DAC_L',dh_l,&
           code=203, longname='diabatic forcing longitudinal (mean)', units='m/s^2')

      ALLOCATE(buf_t(ldc%nglon,nlev),buf_m(ldc%nglon,nlev), buf_l(ldc%nglon,nlev))

      ! the correction due to the initialization procedure
      CALL default_stream_setting(stream_nmi, &
           table=nudging_table, bits=24, &
           lpost  = .TRUE., lrerun = .FALSE., laccu  = .FALSE.,            &
           repr=SPECTRAL, ldims=(/nlev, 2, snsp/), gdims=(/nlev, 2, nsp/), &
           no_default=.TRUE. )

      ALLOCATE(buf_tp(nlevp1,2,snsp,1))
      p4 => buf_tp(:,:,:,1:1)
      CALL add_stream_element(stream_nmi,'COR_TP',p3,p4=p4,&
           code=205, longname='correction tp',          units='K,LN(Pa)/s',&
           ldims=(/nlevp1, 2, snsp/), gdims=(/nlevp1, 2, nsp/), klev=nlevp1 )

      p4 => buf_tp(1:nlev,:,:,1:1)
      CALL add_stream_element(stream_nmi,'COR_T',p3,p4=p4,&
           code=206, longname='correction temperature', units='K/s',&
           ldims=(/nlev, 2, snsp/), gdims=(/nlev, 2, nsp/) )

      p4 => buf_tp(nlevp1:nlevp1,:,:,1:1)
      CALL add_stream_element(stream_nmi,'COR_P',p3,p4=p4,&
           code=207, longname='correction pressure',    units='LN(Pa)/s',&
           leveltype=SURFACE, ldims=(/1, 2, snsp/), gdims=(/1, 2, nsp/) , klev=1)

      CALL add_stream_element(stream_nmi,'COR_VO',buf_vo,&
           code=208, longname='correction vorticity',   units='1/s',&
           ldims=(/nlev, 2, snsp/), gdims=(/nlev, 2, nsp/) )

      CALL add_stream_element(stream_nmi,'COR_D', buf_d, &
           code=209, longname='correction divergence',  units='1/s',&
           ldims=(/nlev, 2, snsp/), gdims=(/nlev, 2, nsp/) )

      CALL message('','NMI stream initialized')
    END IF

    CALL message('','NMI initialized successful')

  END SUBROUTINE NMI_Init

  INTEGER FUNCTION jm_to_zhmod_memory_estimation(jm) RESULT(ipos)
    USE mo_control,       ONLY: nn, nk, nlev
    INTEGER, INTENT(in) :: jm ! Wavenumber, 0 <= jm <= nm
    INTEGER :: imo, ing, inge, inr, ipk, ingt, km, jk, jsym

    ipos = 0

    ! from NMI_Horizontal_Modes:
    km = jm
    ipk = MIN(nn+km,nk)

    level_loop: DO jk = 1, nlev       ! all vertical modes are included
      sym_loop: DO jsym = 0, 1        ! separate loop for symmetric and antisymmetric case

        imo = 2*INT((ipk-km+jsym)/2)+INT((ipk-km+1-jsym)/2)
        IF (km == 0) THEN
        ELSE
          imo = imo + 2 - jsym
        ENDIF

        IF (jsym == 0) THEN
          ingt = (nk-km+2)/2
          IF (km == 0) ingt = ingt-1
        ELSE
          ingt = (nk-km+1)/2
        ENDIF
        inge = ingt

        !inge = imo ! == worst case assumption

        ing = 2*inge
        inr = imo-ing

        ! ---- store eigenvalues ZR
        ipos = ipos + ing
        ipos = ipos + inr

        ! ---- store eigenvectors ZA
        ipos = ipos + imo*inge
        ipos = ipos + imo*inge
        ipos = ipos + imo*inr
      ENDDO sym_loop
    ENDDO level_loop

  END FUNCTION jm_to_zhmod_memory_estimation
  
!======================================================================
!BOP
  ! !IROUTINE: NMI_Close
  ! !INTERFACE:

  SUBROUTINE NMI_Close

    ! !DESCRIPTION:
    ! get all NMI memory free

!EOP
    IMPLICIT NONE

    IF (nmi_mem) THEN
      
      IF (ALLOCATED(have_wave)) DEALLOCATE(have_wave)
      IF (ALLOCATED(zhmod_wave_mem)) DEALLOCATE(zhmod_wave_mem)
      IF (ALLOCATED(zhmod)) DEALLOCATE(zhmod)
      DEALLOCATE(zdb, phib)
      DEALLOCATE(vm_ev, vm_evi, ihmod, iplim, iplimd)
      DEALLOCATE(svo_nmi, sd_nmi, stp_nmi)
      DEALLOCATE(svo_dt,  sd_dt,  stp_dt)
      DEALLOCATE(buf_tp, buf_t, buf_m, buf_l)

      nmi_mem = .FALSE.

    ELSE
      CALL message('NMI_Close','no NMI-memory used before')

    END IF

  END SUBROUTINE NMI_Close

!======================================================================
!BOP
  ! !IROUTINE: NMI_Make
  ! !INTERFACE:

  SUBROUTINE NMI_Make(funcmode)
    ! !DESCRIPTION:
    ! Control the NMI procedure, the subroutine is called twice in 
    ! initialization mode. The first call is in scan1 before the physical
    ! forcings will be calculated. The second call is at the end
    ! of a time step.

    ! !USES:
    USE mo_memory_sp,       ONLY: svo, sd, stp
    USE mo_time_control,    ONLY: l_putrerun, l_putdata, next_date, &
                                  change_current_date, current_date
    USE mo_control,         ONLY: lnmi, nlevp1
    USE mo_time_conversion, ONLY: OPERATOR(<), OPERATOR(==)
    USE mo_memory_base,     ONLY: set_stream
    USE mo_output,          ONLY: open_output_streams, close_output_streams, &
                                  backup_output_streams
!BOX
    IMPLICIT NONE
!EOX
    ! !INPUT PARAMETERS:
    INTEGER, INTENT(in)   :: funcmode  ! specific function
!EOP
!BOC
    SELECT CASE (funcmode)
!BOX
    CASE(NMI_MAKE_TEST) ! that's only for testing
       CALL finish('TEST BENCH ENDE','Run aborted.')
!EOX
    CASE (NMI_MAKE_FOD)            ! *** filter ERA data
!BOX
       CALL NMI_LoadBuffer(NMI_LOAD_ORD_BUFA)    ! observed --> buf A
       CALL NMI_Filter(.TRUE.,NMI_SEL_FILT_W0)   ! filter ERA
       CALL NMI_LoadBuffer(-NMI_LOAD_ORD_BUFA)   ! buf A --> observed
!EOX
    CASE (NMI_MAKE_FMMO)           ! *** filter anomaly (model-obs)
!BOX
       CALL NMI_LoadBuffer(NMI_LOAD_ECH_BUFA)    ! model --> buf A

       CALL NMI_LoadBuffer(NMI_LOAD_OIN_BUFB)    ! interpolated --> buf B
       CALL NMI_LoadBuffer(NMI_CALC_DIFF)        ! buf A - buf B --> buf A
       CALL NMI_Filter(.TRUE.,NMI_SEL_FILT)      ! calculate fast part of model-era
       CALL NMI_LoadBuffer(NMI_STORE_FAST)       ! store fast part
       CALL NMI_LoadBuffer(NMI_CALC_ADD)         ! buf A + buf B --> buf B
       CALL NMI_LoadBuffer(-NMI_LOAD_OIN_BUFB)   ! buf B --> interpolated
!EOX
    CASE (NMI_MAKE_NMI)            ! *** initialsation loop
!BOX
      WRITE(nmi_mess,*) 'Run NMI initialization, start with phase = ',nmi_phase
      CALL message('',nmi_mess)

      IF ( (nmi_phase /= NMI_NULL) .AND. l_putrerun) &
           CALL finish('NMI_Make','rerun not allowed in NMI mode')

      SELECT CASE(nmi_phase)

      CASE(NMI_NULL)       ! nmi is active but not running
        ! check date
        IF ( ( (current_date < nmi_start_date) .OR. (current_date == nmi_start_date) ) &
             .AND. (nmi_start_date < next_date) )THEN
          nmi_start_date = current_date   ! store restart date
          nmi_phase = NMI_PRE_OPEN
        END IF
        WRITE(nmi_mess,*) 'NMI_NULL phase ',nmi_phase,' counter ',nmi_counter

      CASE(NMI_PRE_OPEN)        ! normal integration, before nmi start date 

        ! store sd, svo, stp in temporary buffer
        buf_tp(:,:,:,1) = stp(:,:,:)
        buf_vo(:,:,:)   = svo(:,:,:)
        buf_d (:,:,:)   = sd (:,:,:)
        ! unklar ob nicht diese felder fuer die iteration als ausgangspunkt
        ! dienen sollten
        !
        !CALL NMI_LoadBuffer(NMI_LOAD_ECH_BUFB)  ! atm -> DT
        nmi_counter = 2*ntpre ! NMI_Make is called twice at each time step
        nmi_phase = NMI_PRE
        WRITE(nmi_mess,*) 'NMI_PRE_OPEN phase ',nmi_phase,' counter ',nmi_counter

      CASE(NMI_PRE)
        IF (nmi_counter == 0) THEN
          dh_t(:,:,: ) = 0.0_dp  ! prepare fields for accumulation
          dh_m(:,:,:)  = 0.0_dp
          dh_l(:,:,:)  = 0.0_dp
          nmi_phase = NMI_ACCU_OPEN
        ELSE
          nmi_counter = nmi_counter - 1
        END IF
        WRITE(nmi_mess,*) 'NMI_PRE phase ',nmi_phase,' counter ',nmi_counter

      CASE(NMI_ACCU_OPEN)
        ! set accumulation counter in scan1
        nmi_phase = NMI_ACCU
        nmi_counter = 2*ntdia
        WRITE(nmi_mess,*) 'NMI_ACCU_OPEN phase ',nmi_phase,' counter ',nmi_counter

      CASE(NMI_ACCU)       ! accumulation phase
        IF(nmi_counter == 1) THEN
          ! average diabatic forcing
          dh_t = dh_t/ntdia
          dh_m = dh_m/ntdia
          dh_l = dh_l/ntdia
          nmi_phase = NMI_ACCU_CLOSE
          ! ?? unklar ob hier nicht das ausgangsfeld vom start der initialisierung
          ! benutzt werden sollte, wegen leap-frog jedoch inkonsistenz
          ! d.h. buf_tp, buf_vo, buf_d laden, bzw. oben laden
          ! at the beginning of the iteration the data from the last prognostic
          ! time step will be load into the DT buffer
          CALL NMI_LoadBuffer(NMI_LOAD_ECH_BUFB)  ! atm -> DT
          ! filter the diabatic forcing, not inserted here now
          ! first it must be transformed into spectral space
          ! notwendig ? oder wird das durch das Modell erledigt

        ELSE IF (nmi_counter == 0) THEN
          CALL finish('NMI_Make','no accumulation phase of diabatic heating defined') 
        END IF
        nmi_counter = nmi_counter - 1
        WRITE(nmi_mess,*) 'NMI_ACCU phase ',nmi_phase,' counter ',nmi_counter

      CASE(NMI_ACCU_CLOSE)
        nmi_phase = NMI_USE_AVG
        nmi_counter = 2*ntiter
        WRITE(nmi_mess,*) 'NMI_ACCU_CLOSE phase ',nmi_phase,' counter ',nmi_counter

      CASE(NMI_USE_AVG)

        IF(MOD(nmi_counter,2) == 0) THEN
          ! perform an iteration step
          ! buffer DT (B) contains the previous value
          ! buffer NMI (A) is load with the new forcing
          !
          ! ** prepare the forcing field
          CALL NMI_LoadBuffer(NMI_LOAD_ECH_BUFA)  ! atm -> NMI
          CALL NMI_LoadBuffer(NMI_CORREL_BUFFER)  ! check iteration
          ! the tendency is calculated and stored in NMI buffer
          CALL NMI_LoadBuffer(NMI_CALC_TEND)      ! NMI = (NMI - DT) / delta_time
          ! store the field for next iteration
          CALL NMI_LoadBuffer(NMI_LOAD_ECH_BUFB)  ! atm ->DT
          !
          ! ** fast mode part iteration
          ! calculates the fast mode part balance using second filter
          CALL NMI_Filter(.FALSE.,NMI_SEL_FILTD)  ! NMI modified
          !
          ! ** prepare the new field for iteration
          CALL NMI_LoadBuffer(NMI_CALC_DIFF)      ! subtract fast mode part
          ! restore atm fields
          CALL NMI_LoadBuffer(-NMI_LOAD_ECH_BUFA) ! NMI -> atm
          !
          CALL message('','NMI_USE_AVG run nmi iteration')
        ELSE
          CALL message('','NMI_USE_AVG skip NMI_Make')
        END IF

        IF (nmi_counter == 0) THEN
          ! archive the final correction into output
          buf_tp(:,:,:,1) = buf_tp(:,:,:,1) - stp(:,:,:)
          buf_vo(:,:,:)   = buf_vo(:,:,:)   - svo(:,:,:)
          buf_d (:,:,:)   = buf_d (:,:,:)   - sd (:,:,:)

          nmi_phase = NMI_RESTART
        END IF

        nmi_counter = nmi_counter - 1
        WRITE(nmi_mess,*) 'NMI_USE_AVG phase ',nmi_phase,' counter ',nmi_counter

      CASE(NMI_RESTART)

        CALL backup_output_streams
        CALL close_output_streams

        ! reset time stepping
        CALL change_current_date(nmi_start_date)

        CALL open_output_streams

        nmi_phase = NMI_PRINT_NEXT
        WRITE(nmi_mess,*) 'NMI_RESTART phase ',nmi_phase,' counter ',nmi_counter

      CASE(NMI_PRINT_NEXT)
        WRITE(nmi_mess,*) 'global averaged pressure = ',stp(nlevp1,1,1)
        CALL message('NMI_Make',nmi_mess)

        IF (l_putdata(1)) nmi_phase = NMI_PRINT
        WRITE(nmi_mess,*) 'NMI_PRINT_NEXT phase ',nmi_phase,' counter ',nmi_counter

      CASE(NMI_PRINT)
        nmi_phase = NMI_NULL
        lnmi = .FALSE.
        lnmi_run = .FALSE.
        CALL set_stream(stream_nmi, lpout=.FALSE., lrerun=.FALSE.)
        WRITE(nmi_mess,*) 'NMI_PRINT phase ',nmi_phase,' counter ',nmi_counter

      END SELECT

      CALL message('NMI_Make',nmi_mess)

    END SELECT

  END SUBROUTINE NMI_Make
!EOX
!EOC
!======================================================================
!BOP
  ! !IROUTINE: NMI_Normal_Modes
  ! !INTERFACE:

  SUBROUTINE NMI_Normal_Modes(nlev, filter_only)

    ! !DESCRIPTION:
    ! calculates the full set of normal modes

    ! !USES:
    USE mo_control,       ONLY: nm, vct, nvclev, earth_angular_velocity
    USE mo_semi_impl,     ONLY: apr, betadt, tr
    USE mo_math_constants,ONLY: pi
    USE mo_physical_constants,     ONLY: grav
    USE mo_time_control,  ONLY: delta_time
    USE mo_hyb,           ONLY: bb, delpr
!BOX
    IMPLICIT NONE
!EOX
    ! !INPUT PARAMETERS:
    INTEGER, INTENT(in) :: nlev          ! number of model levels
    LOGICAL, INTENT(in) :: filter_only   ! specify for nudging mode (filter only)

!EOP
    REAL(dp) :: &
         zevima, zevimi, zdtim2, renorm1,       &! time unit factor
         zwrk1(4*nlev), zwrk2(nlev,nlev),       &! workspace
         zwrk3(nlev,nlev),                      &
         zvr(nlev,nlev), zvi(nlev,nlev)          ! real/imaginary part of eigenvectors

    INTEGER :: &
         ipivot(nlev),             &! permutation index
         ifail,                    &! failure code for external procedures
         imar(1), imai(1),         &
         ipos,                     &! pointer in horizontal mode array
         ings,                     &! no. of modes symmetric part
         jk, js, jm, jj, jk1, jk2, &
         iposgw                     ! pointer for different wave types

    INTEGER :: ipos_prev
    
    CALL message('NMI_Normal_Modes','Vertical mode setup:')

    ! Find, sort and print eigenvalues and vectors of bb

    IF (betadt == 0.0_dp) THEN
      zdtim2 = (delta_time*0.5_dp)**2         ! neccessary for full explicit scheme
    ELSE
      zdtim2 = (delta_time*0.5_dp*betadt)**2  ! use this factor to get meaningfull data range
    ENDIF

    ! the vertival structure matrix is a function of the Basic State
    ! calculate eigenvalues and eigenvectors
    ! rescale structur matrix (gravitation-matrix)

    zwrk2(1:nlev,1:nlev) = TRANSPOSE(bb(1:nlev,1:nlev)) / zdtim2

    ! all other possibilities are not correct
    ! bb^T or bb --> eigenvalues larger then expected
    ! bb or bb/zdtim2 --> eigenvektor others then given in paper

    ifail = 0
    CALL dgeev (&    ! LAPACK eigenvalue solver A*X = lamdba*X
         'n',&          ! calculate no left eigenvectors
         'v',&          ! calculate right eigenvectors
         nlev,&         ! order of matrix A
         zwrk2,&        ! inp: structure matrix A, out: overwritten
         nlev,&         ! leading dimension of A
         zdb,&          ! real part of eigenvalues
         phib,&         ! imaginary part of eigenvalues
         zvr,&          ! contains nothing
         nlev,&         ! leading dimension of left eigenvectors
         zvi,&          ! right eigenvectors in columns
         nlev,&         ! leading dimension of right eigenvectors
         zwrk1,&        ! work space
         4*nlev,&       ! dimension of workspace
         ifail)         ! error return code
    IF (ifail /= 0) THEN
      CALL message('','Calculation of eigenvectors/values failed.')
      WRITE (nmi_mess,'(a,i4)') ' Failure in s/dgeev, ifail = ', ifail
      CALL finish ('NMI_Normal_Modes',nmi_mess)
    ENDIF

    imai   = MAXLOC(ABS(phib(1:nlev)))    ! find eigenvalue with largest imaginary part
    zevima = phib(imai(1))
    imai   = MINLOC(ABS(phib(1:nlev)))    ! find eigenvalue with smallest imaginary part
    zevimi = phib(imai(1))
    IF ((ABS(zevima) > 0.0_dp).OR.(ABS(zevimi) > 0.0_dp)) &
         CALL finish ('NMI_Normal_Modes',&
         'Found complex eigenvalues. Failure in vertical structure matrix')

    ! sort eigenvalues and eigenvectors
    jk = 1
    sort_loop: DO

      IF (ABS(phib(jk)) > 0.0_dp) THEN     ! found complex and complex conjugate eigenvalue
        vm_ev(:,jk)    =  zvi(:,jk)
        vm_evi(:,jk)   =  zvi(:,jk+1)
        vm_ev(:,jk+1)  =  zvi(:,jk)
        vm_evi(:,jk+1) = -zvi(:,jk+1)
        jk = jk+2
      ELSE                          ! real eigenvalue
        vm_ev(:,jk)  = zvi(:,jk)
        vm_evi(:,jk) = 0.0_dp
        jk = jk+1
      ENDIF
      IF (jk > nlev) EXIT

    ENDDO sort_loop

    DO jk=1,nlev
      imar         = MAXLOC(zdb(:))      ! find next largest real part of all eigenvalues
      phib(jk)     = zdb(imar(1))
      zdb(imar(1)) = 0.0_dp
      vm_evi(:,jk) = vm_ev(:,imar(1))       ! store real part of eigenvector
    ENDDO

    ! equivalent depths in PHIB
    ! corresponding eigenvectors in VM_EVI

    WRITE (nmi_mess,'(a,f6.1,a,f8.0,a)') &
         '  Vertical modes for tr = ', tr, ' k;   pr = ', apr,' pa'
    CALL message('',nmi_mess)

    WRITE (nmi_mess,'(a)') &
         '   level       a(k)        b(k)     eigenvalues    phase velocity   EV/g[gpm]'
    CALL message('',nmi_mess)

    DO jk = 1,nlev
      WRITE (nmi_mess,'(3x,i4,2x,f13.6,f12.8,f13.2,f14.2,f13.2)') &
           jk,vct(jk),vct(nvclev+jk), phib(jk), SQRT(phib(jk)), phib(jk)/grav
      CALL message('',nmi_mess)
    END DO

    ! eigenvectors normalized
    ! convert vertical modes to "universal" scaling

    zevima = SUM(delpr(:))                ! check the normalization
    IF (ABS(zevima-apr) > 1.e-5_dp) THEN
      CALL message('','Universal scaling failure.')
      WRITE (nmi_mess,'(2(a,e20.10))') &
           '  surface pressure ... expected = ',apr,' calculated = ',zevima
      CALL finish ('NMI_Normal_Modes',nmi_mess)
    ENDIF

    ! for rescaling and orthogonalization of vertical modes see paper (198?) by C.Temperton
    !
    !    WRITE (nout,'(a)') ' **NMI_Normal_Modes** Rescaling factor for: '
    !    DO  jk=1,nlev
    !       zsum      = SUM(vm_evi(:,jk)*vm_evi(:,jk)*delpr(:))
    !       zscale    = SQRT(apr/zsum)
    !       !       vm_ev(:,jk)  = zscale*vm_evi(:,jk)
    !       vm_ev(:,jk)  = vm_evi(:,jk)
    !       zrscale   = 1.0/zscale
    !       WRITE (nout,'(a,i3,a,f14.10,a,f14.10,a)') &
    !&            '   vertical mode', jk, ' = ', zscale,&
    !&            '   ( ~ for coefficients = ', zrscale, ')'
    !    ENDDO
    !    WRITE (nout,'(a)') '  do not use the universal rescaling now !!!'
    !
    ! scaled eigenvectors in ZM

    ! Compute h**(-1) and delpr*b**(-1)
    ! calculate inverse of ZM

    zwrk2(:,:) = vm_ev(:,:)
    CALL unitmx(zwrk3,nlev)
    ifail = 0
    CALL dgesv (&       ! BLAS linear equation solver
         nlev,&            ! number of linear equations NEQ [int]
         nlev,&            ! number of right hand sides NRHS [int]
         zwrk2,&           ! inp: coefficients matrix A(LDA,NEQ), out: l-u-factorization [real]
         nlev,&            ! leading dimension LDA [int]
         ipivot,&          ! out: permutation matrix IPIV(NEQ) [int]
         zwrk3,&           ! inp: right hand side X(LDX,RHS), out: solution matrix [real]
         nlev,&            ! leading dimension LDX [int]
         ifail)            ! error code [int]
    IF (ifail /= 0) THEN
      WRITE (nmi_mess,'(a,i4)') 'Inverting VM_EV, failure in s/dgesv, ifail = ', ifail
      CALL finish ('NMI_Normal_Modes',nmi_mess)
    ENDIF
    vm_evi(:,:) = zwrk3(:,:)    ! inverse eigenvector matrix

    ! calculate inverse of BB
    zwrk2(1:nlev,1:nlev) = TRANSPOSE(bb(1:nlev,1:nlev)) / zdtim2
    CALL unitmx(zwrk3,nlev)
    ifail = 0
    CALL dgesv (nlev,nlev,zwrk2,nlev,ipivot,zwrk3,nlev,ifail)
    IF (ifail /= 0) THEN
      WRITE (nmi_mess,'(a,i4)') 'Inverting BB, failure in s/dgesv, ifail = ', ifail
      CALL finish ('NMI_Normal_Modes',nmi_mess)
    ENDIF

    DO  jk=1,nlev                ! final calculation of ZDB (S^T*BB^-1)
      zdb(jk) = DOT_PRODUCT(delpr(:),zwrk3(:,jk))/apr
    ENDDO

    ! delpr/bb/apr stored in ZDB

    CALL message('NMI_Normal_Modes','Calculate horizontal modes:')
    ipos = 0   ! index in full array with modes
    IF (reduce_memory_usage) THEN
      IF (zero_nsm) RETURN
    ENDIF
    wave_loop: DO  jm = 0, nm                 ! calculation for all zonal wave numbers

      IF (reduce_memory_usage) THEN
        ipos_prev = ipos
        CALL assert(ipos+zhmod_wave_mem(jm) <= SIZE(zhmod), __LINE__)
      ENDIF

      CALL NMI_Horizontal_Modes (jm, ipos)
                                              ! Print frequencies of horizontal modes
      renorm1 = pi / (earth_angular_velocity * 3600.0_dp)        ! rescale factor to get the period in hours

      ! Set the filter limits for the NMI
      level_loop: DO jk = 1,nlev

        sym_loop: DO js = 1,2

          ings = ihmod(2,js,jm+1,jk)/2      ! get number of one set of gravity modes
          gwest_mode_loop: DO jj = 1,ings   ! loop over one set of gravity modes (westerly part)

            iposgw = ihmod(1,js,jm+1,jk) -1 +ings +jj  ! position of eigenvalue

            ! search periods smaller then the cut-off period

            IF (renorm1/zhmod(iposgw) < pcut) THEN  ! normal filter
              nvm = jk ; iplim(js,jm+1,jk) = iplim(js,jm+1,jk) + 1
            END IF

            IF (renorm1/zhmod(iposgw) < pcutd) THEN ! second filter for diabatic tendencies
              nvmd = jk ; iplimd(js,jm+1,jk) = iplimd(js,jm+1,jk) + 1
            END IF

          END DO gwest_mode_loop

        END DO sym_loop

      END DO level_loop

      IF (reduce_memory_usage) THEN
        IF (.NOT. have_wave(jm)) ipos = ipos_prev ! we don't need zhmod for this jm anymore:
      ENDIF

      WRITE (nmi_mess,'(a,i5,a)') &
           'No of horizontal modes for wavenumber ', jm, ' in cut-off window:'
      CALL message('',nmi_mess)

      ! control print out for normal filter range
      jk1 = 1
      DO
        jk2 = jk1 + 8
        jk2 = MIN(jk2,nvm)
        WRITE(nmi_mess,'(20(a,i3.3,a,i5,1x))') &
             ('VM_(',js,'): ',iplim(1,jm+1,js)+iplim(2,jm+1,js),js=jk1,jk2)
        CALL message('',nmi_mess)
        IF (jk2 == nvm) EXIT
        jk1 = jk2+1
      END DO

      IF (filter_only) CYCLE

      ! control print out for diabatic filter range
      jk1 = 1
      DO
        jk2 = jk1 + 8
        jk2 = MIN(jk2,nvmd)
        WRITE(nmi_mess,'(20(a,i3.3,a,i5,1x))') &
             ('VMD(',js,'): ',iplim(1,jm+1,js)+iplim(2,jm+1,js),js=jk1,jk2)
        CALL message('',nmi_mess)
        IF (jk2 == nvmd) EXIT
        jk1 = jk2+1
      END DO
      
    ENDDO wave_loop

  END SUBROUTINE NMI_Normal_Modes



!======================================================================
!BOP
  ! !IROUTINE: NMI_Horizontal_Modes
  ! !INTERFACE:

  SUBROUTINE NMI_Horizontal_Modes(km, ipos)

    ! !DESCRIPTION:
    ! calculate horizontal modes for one zonal wave number

    ! !USES:
    USE mo_control,       ONLY: nn, nnp1, nk, nlev

    ! !INPUT/OUTPUT PARAMETERS:
    INTEGER, INTENT(IN)    :: km     ! actual zonal wave numbers for which
                                     ! modes are to be computed
    INTEGER, INTENT(INOUT) :: ipos   ! record position in ihmod
!EOP

    REAL(dp), ALLOCATABLE :: &
         za(:,:), &         ! horizontal struktur matrix
         zr(:), &           ! eigenvalues of horizontal structur matrix
         work(:)

    INTEGER :: &
         ipk, &             ! maximal number of meridional index
         kmo, &             ! maximal field dimension
         jk, &              ! vertical loop index
         jsym, &            ! symmetry loop
         jr, il, &
         imo, &             ! number of modes
         jms, &
         jn, &
         ifail, inge, j, jj, ingt, imo2

    REAL(dp) :: zed        ! local equivalent depths

    INTEGER :: last_ipos
    last_ipos = ipos
    
    ! allocate local memory
    ipk = MIN(nn+km,nk)
    kmo = 3*nnp1 + 6
    ALLOCATE ( za(kmo,kmo), zr(kmo), work(3*kmo-1) )

    ! start loop on equivalent depths

    level_loop: DO jk = 1, nlev       ! all vertical modes are included

      zed = phib(jk)                  ! get the equivalent depths

      sym_loop: DO jsym = 0, 1        ! separate loop for symmetric and antisymmetric case
          za(:,:) = 0.0_dp
          zr(:)   = 0.0_dp
          work(:) = 0.0_dp

          ! Define matrix A (za)
          ! fill in upper triangel + diagonal of A
          ! IMO ... number of modes can be solved

          imo = 2*INT((ipk-km+jsym)/2)+INT((ipk-km+1-jsym)/2)
          IF (km == 0) THEN
             jms = 2
             il  = 2 + jsym
          ELSE
             imo = imo + 2 - jsym
             jms = km
             il  = 1
          ENDIF

          imo2 = 0
          sym_asym: IF (jsym == 0) THEN                ! symmetric flow

            DO jn = jms, ipk, 2             ! loop for velocity potential and mass
              za(il,il)   =  cf(km,jn)
              za(il,il+1) =  ef(jn,zed)
              za(il,il+2) = -bf(km,jn+1)
              il   = il   + 3
              imo2 = imo2 + 2
            ENDDO
            il = 3                          ! loop for streamfct.
            IF (km == 0) il = 1
            DO jn = km+1, ipk, 2
              za(il,il)   =  cf(km,jn)
              za(il,il+1) = -bf(km,jn+1)
              il   = il   + 3
              imo2 = imo2 + 1
            ENDDO

          ELSE                              ! antisymmetric flow

            DO jn = jms, ipk, 2             ! loop for streamfct.
              za(il,il)   =  cf(km,jn)
              za(il,il+1) = -bf(km,jn+1)
              il   = il   + 3
              imo2 = imo2 + 1
            ENDDO
            il = 2                          ! loop for velocity pot. and mass
            IF (km == 0) il = 1
            DO jn = km+1, ipk, 2
              za(il,il)   =  cf(km,jn)
              za(il,il+1) =  ef(jn,zed)
              za(il,il+2) = -bf(km,jn+1)
              il   = il   + 3
              imo2 = imo2 + 2
            ENDDO

          ENDIF sym_asym

          ! setup of matrix A following W.Wergen (OK 02-feb-99)
          IF (imo/=imo2) THEN
            WRITE(nmi_mess,*) ' ERROR: count modes ',imo2,' expected ',imo
            CALL message('',nmi_mess)
            CALL finish ('NMI_Horizontal_Modes', 'Run aborted.')
          ENDIF

          DO jr = 1, imo                   ! fill in lower triangal
            za( jr+1:imo, jr) = za( jr, jr+1:imo)
          ENDDO

          ! Solve eigenproblem
          ifail = 0
          CALL dsyev (&  ! LAPACK get eigenvalues and eigenvectors of real symmetric matrix
               'V',&        ! get eigenvalues and eigenvectors
               'L',&        ! lower triangle is in  A stored
               imo,&        ! order of matrix, N
               za,&         ! inp: symmetric matrix A(LDA,N), out: orthonormal eigenvectors
               kmo,&        ! leading dimension of matrix A LDA
               zr,&         ! out: eigenvalues in ascending order
               work,&       ! workspace
               3*kmo-1,&    ! dimension of workspace
               ifail)       ! return code
          IF (ifail /= 0) THEN
            WRITE (nmi_mess,'(a,i2,a,g13.5,a,i3,a,i2)') &
                 ' Error in s/dsyev - ifail = ', ifail, &
                 ' for phi = ', zed, ' m = ', km, ' and isym = ', jsym
            CALL finish ('NMI_Horizontal_Modes',nmi_mess)
          ENDIF

          ! Store sorted frequencies and modes in buffer ZHMOD
          !
          ! zr == free frequencies of normal modes
          ! order is (gravity east, rossby, gravity west)
          ! in ZR (GE[max..0] R[max..0] GW[0..max])
          ! following the frequency [meridionalindex]
          !
          ! reorder in ZHMOD
          ! (GE[0..max] GW[max..0] R[0..max]
          ! slow......fast ........ slow

          inge = 0
          DO  j = 1,imo      ! count eastward moving gravity waves
            IF (zr(j) > -1e-10_dp) EXIT
            inge = inge+1
          ENDDO
          ing = 2*inge       ! all gravity modes
          inr = imo-ing      ! residue is number of rossby modes

          ! check the number of modes
          IF (jsym == 0) THEN
            ingt = (nk-km+2)/2
            IF (km == 0) ingt = ingt-1
          ELSE
            ingt = (nk-km+1)/2
          ENDIF
          IF (inge /= ingt) THEN
            WRITE (nmi_mess,'(a,/,a,3i3,a,i4,a,i4)') &
                 '   Possible problem in counting eastward gravity waves', &
                 '   - for jm, jk, jsym (', km, jk, jsym, ') theory: ', &
                 inge, ' reality: ', ingt
            CALL finish ('NMI_Horizontal_Modes',nmi_mess)
          ENDIF

          ! set pointer array for horizontal modes
          ! ipos is the last index of the previous array
          ! ihmod(1,...) is the position of the first place in the actual array

          ihmod(1,jsym+1,km+1,jk) = ipos + 1    ! start position
          ihmod(2,jsym+1,km+1,jk) = ing         ! no. of gravity modes
          ihmod(3,jsym+1,km+1,jk) = inr         ! no. of rossby modes


          ! reset frequencies of rossby modes for wave number zero
          ! these Rossby modes stationary
          ! the eigenvectors are non zero 
          ! the mixed rossby-gravity mode and the kelvin mode are
          ! gravity modes here
          IF (km==0) THEN
            DO jj=inge+1,inge+inr
              IF (ABS(zr(jj)) > EPSILON(1.0_dp)) THEN
                WRITE(nmi_mess,*) 'Rossby mode frequency non zero ',zr(jj),&
                     ': set to zero; WM= ',km,' VM= ',jk
                CALL message('',nmi_mess)
              ELSE
                WRITE(nmi_mess,*) 'Rossby mode with zero frequency found;',&
                     ' WN= ',km,' VM= ',jk
                CALL message('',nmi_mess)
              END IF
              zr(jj) = 0.0_dp
            END DO
          END IF

          ! ---- store eigenvalues ZR
          DO jj=1,inge
            zhmod(ipos     +jj) = zr(inge+1-jj)  ! GW east
            zhmod(ipos+inge+jj) = zr(imo +1-jj)  ! GW west
          END DO
          ipos = ipos + ing
          DO jj= inge+inr, inge+1, -1
            zhmod(ipos+jj) = zr(jj)   ! Rossby modes
          END DO
          ipos = ipos + inr

          ! ---- store eigenvectors ZA
          DO jj= inge , 1, -1
            zhmod(ipos+1:ipos+imo) = za(1:imo, jj)
            ipos = ipos + imo
          END DO
          DO jj= imo , imo-inge+1, -1
            zhmod(ipos+1:ipos+imo) = za(1:imo, jj)
            ipos = ipos + imo
          END DO
          DO jj= inge+inr, inge+1, -1
            zhmod(ipos+1:ipos+imo) = za(1:imo, jj)
            ipos = ipos + imo
          END DO

          ! IPOS is now the last memory position in ZHMOD used
        ENDDO sym_loop

     ENDDO level_loop

     IF (reduce_memory_usage) THEN
       CALL assert(ipos - last_ipos == zhmod_wave_mem(km), __LINE__)
     ENDIF
     
     DEALLOCATE (work, zr, za)

  END SUBROUTINE NMI_Horizontal_Modes

!======================================================================
!BOP
  ! !IROUTINE: NMI_Filter
  ! !INTERFACE:

  SUBROUTINE NMI_Filter(lfilter,sel_modes)
    ! !DESCRIPTION:
    ! perform the NMI filter or the non-linear iteration

    ! !USES:
    USE mo_hyb,           ONLY: delpr, rdelpr, nlevm1
    ! rdtr /= rd*tr, see INHYSI for details
    USE mo_truncation,    ONLY: nnp
    !USE mo_semi_impl,     ONLY: apr, betadt, tr
    USE mo_semi_impl,     ONLY: apr, tr
    USE mo_control,       ONLY: nn, nnp1, nk, nlev, nlevp1, earth_angular_velocity
    !USE mo_time_control,  ONLY: delta_time
    USE mo_physical_constants,     ONLY: rd, earth_radius
!BOX
    IMPLICIT NONE
!EOX
    ! !INPUT PARAMETERS:
    LOGICAL, INTENT(in) :: lfilter
    ! switch for projection or initialization
    ! .true.    only projection of absolut values and filtering
    ! .false.   solve the tendecy equation for next iteration
    INTEGER, INTENT(in) :: sel_modes     ! type of projection, mode selection
!EOP
    INTEGER :: invm          ! number of vertical modes actual used

    REAL(dp) :: &
         zzdb   (nlev),              &
         zph    (nlevp1),            &! half level pressure
         zrlnpr (nlev),              &! log (delta p)
         zralphr(nlev),              &! alpha
         zralprr(nlev)                ! 1/alpha

    REAL(dp), ALLOCATABLE, DIMENSION(:) :: &
         zbs,     &! spectral coefficients T, lnPS
         zworks,  &! spectral coefficients VO, D, H
         zbb,     &! separation matrix
         zbv,     &! parts of different waves in vertical space
         zworkv,  &! subset of coefficients in vertical space
         zworkm    ! projection in horizontal mode space

    REAL(dp)    :: &
         zedra, &   ! sqrt(equivalent depth) / (radius of earth)
         zs, zo     ! helpvariables explicit solution

    INTEGER :: &
         jpfs3d, jpfsm, jpbnm, &   ! array dimensions
         isc1, isc2, isc3, &       ! index for scaling, rescaling
!         zdtim, &         ! scaling of equivalent depths
         ijump, &       ! record length of actual vertical mode data set
         isp, &         ! local pointer in spectral array (nlev)
         ismp, &        ! local pointer in spectral array (nlevp1)
         im, &          ! index of local wavenumber
         jm, &          ! index of global wavenumber
         jn, &          ! zonal index number
         ipk, &         ! maximal zonal index for local wave number
         innp2, &       ! actual number of one set of spectral coefficients real+imaginary
         ivmo, &        ! one set of spectral coefficients (VO+D+H)
         innp2nm, &     ! one set of components for all zonal index
         ilt, &         ! actual record length in VO and D
         iltp, &        ! actual record length in TP
         ist, &         ! pointer in T and lnPS
         ish, &         ! pointer in H
         jk, &          ! loop index of level and vertical modes
         izr, izi, idr, idi, ihr, ihi, &     ! pointer in vertical space arrays
         ipos, &        ! pointer in horizontal mode array
         imo, &         ! total number of horizontal modes
         ingrx, &       ! number of modes used
         iposx, &       ! local position in mode array
         iposv, &       ! local position in eigenvector array
         jsym, &        ! index for symmetry loop
         idx1, idx2, idx3, idx4, idx5, idx6 ! symmetry loop steering 

    REAL(dp)    :: oro

    IF (reduce_memory_usage) THEN
      IF (zero_nsm) RETURN
    ENDIF
    
    ! initialization and checks

    SELECT CASE (sel_modes)
    CASE (NMI_SEL_FILT_W0, NMI_SEL_FILT) ; invm = nvm  ! normal filter
    CASE (NMI_SEL_FILTD_W0,NMI_SEL_FILTD); invm = nvmd ! diabatic filter
    CASE default                         ; invm = nlev ! all vertical modes included
    END SELECT

!    IF (betadt == 0.0_dp) THEN
!      zdtim = INT(delta_time*0.5_dp)         ! neccessary for full explicit scheme
!    ELSE
!      zdtim = INT(delta_time*0.5_dp*betadt)  ! use this factor to get meaningfull data range
!    ENDIF

    ! rdtr is modified in model setup  (rdtr .ne. rd*tr) here

    ! ===========================================================================

    jpfs3d = 2 * nlev   * nnp1  ! max length of spectral array for D and VO
    jpfsm  = 2 * nlevp1 * nnp1  ! max length of spectral array for T and lnPS
    jpbnm  = 3 * nnp1           ! max number of horizontal modes

    ! workspace in spectral space
    ALLOCATE ( zbs    (2 * jpfsm) )        ! ZBS(TP)
    ALLOCATE ( zworks (3 * jpfs3d) )       ! ZWORKS(V+D+H)
    ALLOCATE ( zbb    (2 * nnp1) )         ! ZBB(SPECTCOEF)

    ! workspace in vertical mode space
    ALLOCATE ( zbv (invm * 2 * jpbnm) )    ! ZBV(invm,3*(REAL+IMAG)
    ALLOCATE ( zworkv     (2 * jpbnm) )    ! ZWORKV(3*(REAL+IMAG))

    ! workspace in horizontal mode space
    ALLOCATE ( zworkm (2*3*nnp1) )         ! ZWORKM(imo,2)

    ! ===========================================================================

    ! setup basic state, only one column
    CALL pres(zph,1,apr,1) ! pressure levels for basic state
    zrlnpr(:) = LOG(apr)
    CALL auxhyb(&
         delpr,&      ! OUT: pressure difference
         rdelpr,&     ! OUT: reziproce of DELPR
         zrlnpr,&     ! OUT: logarithm of pressure differences
         zralphr,&    ! OUT: ALPHAs for integration
         zph,&        ! INP: half-level pressure
         1,&          ! first dimension of arrays
         1)           ! number of points
    zralprr(:) = 1._dp/zralphr(:)

    ! initialize pointer in spectral arrays

    ijump = 2*invm
    isp   = 1
    ismp  = 1

    nmi_waveno_loop : DO im = 1,ldc%nsm    !>>>>>>>>>>>>>>>>>>> wavenumber loop START

      jm = ldc%sm(im) ! Wavenumber

      IF (reduce_memory_usage) THEN
        CALL assert(have_wave(jm), __LINE__)
      ENDIF

      ! Check if complete wavenumber resides on this PE

      IF(ldc%snn0(im) /= 0 .OR. ldc%snnp(im) /= nnp(jm+1) ) &
           CALL finish('NMI_Filter','NMI selected and lfull_m not set!')

      ! define local dimensions and pointer in work space

      ipk     = MIN(nn+jm,nk)
      innp2   = nnp(jm+1) * 2   ! spectral components of given wavenumber (R+I)
      ilt     = innp2 * nlev    ! offset in VO and D spectral field
      iltp    = innp2 * nlevp1  ! offset in T/lnPS spectral field
      ivmo    = innp2 * 3       ! normal mode components (VO+D+H)
      innp2nm = innp2 * invm    ! vertical space components

      zworks(    1:  ilt) = svo_nmi(isp :isp +ilt -1)   ! Fill in VO and D
      zworks(ilt+1:2*ilt) = sd_nmi (isp :isp +ilt -1)
      zbs(       1:iltp)  = stp_nmi(ismp:ismp+iltp-1)   ! Fill in T and lnPS

      ! Calculate H with T' and lnPS'
      !
      ! compute h for lowest level
      DO jn = 1, innp2
        ist =         jn*nlevp1                ! pointer  lnPS
        ish = 2*ilt + jn*nlev                  ! pointer  H
        
        ! orography correction not necessary ... test only
        oro = 0.0_dp
        zworks(ish) = &                                   ! ==> H(nlev)
               zralphr(nlev) * zbs(ist-1) &               ! g * T(nlev)
             + rd*tr       * zbs(ist)     &               ! R*Tr*lnPS
             + oro    ! orography correction
       END DO

       ! compute H for remaining levels
       DO jn = 1,innp2
         ist =         (jn-1)*nlevp1            ! pointer  T
         ish = 2*ilt + (jn-1)*nlev              ! pointer  H
         DO jk = nlevm1,1,-1
           zworks(ish+jk) = &                                 ! ==> H(level)
                  zworks(ish+jk+1) &                          ! H(level+1)
                + zbs(ist+jk)  *              zralphr(jk) &   ! T(level)
                + zbs(ist+jk+1)*(zrlnpr(jk+1)-zralphr(jk+1))  ! T(level+1)
          END DO
       END DO

       IF (jm == 0) THEN

         ! conserve wave 0 first meridional index (= global average unchanged)
         ! this should be the basic state
         zworks(      1:      2*nlev)   = 0.0_dp
         zworks(  ilt+1:  ilt+2*nlev)   = 0.0_dp
         zworks(2*ilt+1:2*ilt+2*nlev)   = 0.0_dp

       ENDIF

       ! ZWORKS complete defined (VO+D+H) structure is (NLEV,2,3*NNP)
       !
       ! ========================================================================
       ! Expansion into vertical mode space (page 12 W.Wergen)
       !
       ! M == VM_EV and invers(M) == vm_evi
       ! use vm_evi for projection, organizes as VM_EVI(nlev,invm+)
       !
       ! transpose[VM_EVI(nlev,invm)]*ZWORKS(nlev,ivmo) = ZBV(invm,ivmo)

       CALL dgemm &
            ('N', 'N', invm, ivmo, nlev, 1.0_dp, vm_evi, nlev, zworks, nlev, 0.0_dp, zbv,invm)

       ! ZBV(invm,ivmo) now defined
       !>>>>>>>>>>>>>>>>> vertical mode loop START
       ! Treat horizontal dependence for first invm vertical modes

       vmode_loop : DO jk = 1, invm

         zedra  = earth_radius/SQRT(phib(jk))

         ! phib is calculated using rescaled BB
         !zedra  = a/(SQRT(phib(jk))*zdtim)
         ! rescaling not necessary
         ! set pointer for accumulation of wavenumber fraction
         
         izr    = jk                      ! ZETA real part
         izi    = jk             + invm   !      imaginary part 
         idr    = jk +   innp2nm          ! DIV  real part 
         idi    = jk +   innp2nm + invm   !      imaginary part
         ihr    = jk + 2*innp2nm          ! H    real part
         ihi    = jk + 2*innp2nm + invm   !      imaginary part

         ! Scale fields in ZBV (will be multiplied by i later)
         ! the global average is skipped

         IF (jm == 0) THEN   ! set up indicies for scaling loop
           isc1 = 1
           isc2 = ijump
         ELSE
           isc1 = jm
           isc2 = 0
         ENDIF

         isc3 = isc2
         DO jn = isc1,ipk
           zbv(izr+isc3) = zbv(izr+isc3) * qf(jn)
           zbv(izi+isc3) = zbv(izi+isc3) * qf(jn)
           zbv(idr+isc3) = zbv(idr+isc3) * qf(jn)
           zbv(idi+isc3) = zbv(idi+isc3) * qf(jn)
           zbv(ihr+isc3) = zbv(ihr+isc3) * zedra
           zbv(ihi+isc3) = zbv(ihi+isc3) * zedra
           isc3 = isc3 + ijump
         ENDDO

         ! Separate the type of the solution for nonlinear terms
         ! no separation, implicit excluded
         ! in projection mode all times explicit
         
         IF ( lfilter .OR. &      ! separation during initialization
              ((iplim(1,jm+1,jk)+iplim(2,jm+1,jk)) < nnp(jm+1)) ) THEN

           ! case A: projection only, then always explicit method
           ! case B: filtering
           !         only if gravity waves of all meridional index
           !         inside the filter window the implicit method is applicable ?
           !
           ! Explicit solution
           !
           ! explicit nmi for this wavenumber & vertical mode            
           ! the following code is passed through twice,
           ! symmetric flow (jsym = 0), antisymmetric flow (jsym = 1)

           sym_loop : DO jsym = 0,1

             ! Get eigenvalues and eigenvectors of horizontal modes

             ipos = ihmod(1,jsym+1,jm+1,jk)
             ing  = ihmod(2,jsym+1,jm+1,jk)
             inr  = ihmod(3,jsym+1,jm+1,jk)
             imo  = ing + inr

             !======== load ZWORKV (REAL+IMAG) == (IMO,2)
             ! real part      == (VOr,Dr,Hr)*(No of spectral coefficients)
             ! imaginary part == (VOi,Di,Hi)*(No of spectral coefficients)

             zworkv(1:2*imo) = 0._dp
                   
             ! get number of horizontal normal modes used
             ! ZHMOD(ipos+) contains gravity and rossby eigenvectors ==> matrix A
             ! gravity part = ZHMOD(1:imo,    1:ing)
             ! rossby part  = ZHMOD(1:imo,ing+1:ing+inr)
             ! for projection ZHMOD should be transposed using xGEMM moduls

             SELECT CASE (sel_modes)
             CASE (NMI_SEL_ALL, NMI_SEL_ALL_W0)
               ingrx = imo                      ! number of modes used
               iposx = ipos+imo                 ! start pointer for eigenvectors
               iposv = ipos                     ! start pointer for eigenvalues

             CASE (NMI_SEL_GM, NMI_SEL_RM_IMPL, NMI_SEL_GM_W0)                   
               ingrx = ing
               iposx = ipos+imo
               iposv = ipos

             CASE (NMI_SEL_RM, NMI_SEL_RM_W0)
               ingrx = inr
               iposx = ipos+imo+imo*ing
               iposv = ipos+ing

             CASE (NMI_SEL_FILT, NMI_SEL_FILT_W0)
               ingrx = 2*iplim(jsym+1,jm+1,jk)
               iposx = INT(ipos+imo+imo*(0.5_dp*ing-iplim(jsym+1,jm+1,jk)))
               iposv = INT(ipos+(0.5_dp*ing-iplim(jsym+1,jm+1,jk)))

             CASE (NMI_SEL_FILTD, NMI_SEL_FILTD_W0)
               ingrx = 2*iplimd(jsym+1,jm+1,jk)
               iposx = INT(ipos+imo+imo*(0.5_dp*ing-iplimd(jsym+1,jm+1,jk)))
               iposv = INT(ipos+(0.5_dp*ing-iplimd(jsym+1,jm+1,jk)))

             CASE default
               WRITE (nmi_mess,'(a,i5)') 'not implemented SEL_MODES = ',sel_modes
               CALL finish ('NMI_Filter',nmi_mess)

             END SELECT

             IF (ingrx /= 0) THEN   ! perform filtering explicit

               ! Compose vector x
               CALL setidx ( jsym, ijump, jm,  idx1, idx2, idx3, idx4, idx5, idx6 )
               DO jn = idx5,ipk,2                       ! load ZETA
                 zworkv(idx2)     = -zbv(izr+idx1)
                 zworkv(idx2+imo) = -zbv(izi+idx1)
                 idx1 = idx1 + 2*ijump
                 idx2 = idx2 + 3
               ENDDO

               DO jn = idx6,ipk,2                       ! load DIV and H
                 zworkv(idx4)       = -zbv(idi+idx3)
                 zworkv(idx4+imo)   =  zbv(idr+idx3)
                 zworkv(idx4+1)     =  zbv(ihr+idx3)
                 zworkv(idx4+1+imo) =  zbv(ihi+idx3)
                 idx3 = idx3 + 2*ijump
                 idx4 = idx4 + 3
               ENDDO

               ! Project into horizontal normal mode space
               CALL dgemm &
                    ('T','N',ingrx,2,imo,1.0_dp,zhmod(iposx),imo,zworkv,imo,0.0_dp,zworkm,ingrx)

               ! ZWORKV(imo,2) ==>> ZWORKM(ingrx,2)
               
               IF (.NOT.lfilter) THEN  ! explicit solution for initialisation
                 DO jn = 1,ingrx
                   zs               =  zworkm(ingrx+jn)
                   zo               =  2._dp*earth_angular_velocity*zhmod(iposv+jn-1)
                   IF (ABS(zo) > EPSILON(1.0_dp)) THEN
                     zworkm(ingrx+jn) =  zworkm(jn) / zo
                     zworkm(jn)       =  - zs / zo
                   ELSE
                     WRITE(nmi_mess,*) 'eigenvalue zero, no iteration possible, JN= ',jn,' IPOSV= ',iposv
                     CALL message('',nmi_mess)
                   END IF

                 ENDDO

                 ! ZWORKM contains now the nonlinear forcing increment
               ENDIF

               ! Modify normal mode coefficients in normal mode space
               ! Exclude tidal signal from initialisation ... removed
               ! Project back into vertical space
               !
               !  ZWORKM(ingrx,2) ==>> ZWORKV(imo,2)

               CALL dgemm &
                    ('N','N',imo,2,ingrx,1.0_dp,zhmod(iposx),imo,zworkm,ingrx,0.0_dp,zworkv,imo)

             ENDIF

             ! Decompose vector x
             CALL setidx ( jsym, ijump, jm, idx1, idx2, idx3, idx4, idx5, idx6 )

             ! store explicit part in ZBV
             DO jn = idx5,ipk,2
               zbv(izr+idx1) = -zworkv(idx2)
               zbv(izi+idx1) = -zworkv(idx2+imo)
               idx1 = idx1 + 2*ijump
               idx2 = idx2 + 3
             ENDDO

             DO jn = idx6,ipk,2
               zbv(idi+idx3) = -zworkv(idx4)
               zbv(idr+idx3) =  zworkv(idx4+imo)
               zbv(ihr+idx3) =  zworkv(idx4+1)
               zbv(ihi+idx3) =  zworkv(idx4+1+imo)
               idx3 = idx3 + 2*ijump
               idx4 = idx4 + 3
             ENDDO

           END DO sym_loop

           ! end of explicit calculation for given zonal wavenumber
           !==========================
           ! no implicit calculation included
           ! Implicit nmi for this wavenumber & vertical mode
         ENDIF

         ! now ZBV new calculated
         isc3 = isc2                  ! Rescaling
         DO jn = isc1,ipk
           zbv(izr+isc3) = zbv(izr+isc3) / qf(jn)
           zbv(izi+isc3) = zbv(izi+isc3) / qf(jn)
           zbv(idr+isc3) = zbv(idr+isc3) / qf(jn)
           zbv(idi+isc3) = zbv(idi+isc3) / qf(jn)
           zbv(ihr+isc3) = zbv(ihr+isc3) / zedra
           zbv(ihi+isc3) = zbv(ihi+isc3) / zedra
           isc3 = isc3 + ijump
         ENDDO

       END DO vmode_loop

       ! Inverse vertical transform
       CALL dgemm &
            ('N', 'N', nlev, ivmo, invm, 1.0_dp, vm_ev, nlev, zbv, invm, 0.0_dp, zworks,nlev)

       ! correct global average calculation... removed
       ! alternative separation methode ... removed
       ! calculate changes of surface pressure --> eq. 61 in W.Wergen
       ! ZDB(1,nlev)*ZWORKS(nlev,innp2) = ZBB(1,innp2)
       ! ZBB(1,innp2) ... changes of log(surface pressure)
       ! correct H with orography ... removed
       !
       ! Split H into lnPS and T
       zzdb(:)=zdb(:)
       CALL dgemm &
            ('T', 'N', 1,innp2,nlev, 1.0_dp, zzdb, nlev, zworks(2*ilt+1), nlev, 0.0_dp, zbb,1)

       ! possible modifications
       ! 1. no pressure changes
       ! 2. calculate changes of lnPs from changes of H

       DO jn = 1,innp2
         !zbs(jn*nlevp1) = 0.0_dp
         zbs(jn*nlevp1) = zbb(jn)
       ENDDO

       IF (jm==0) THEN
         ! no changes in global averaged pressure
         zbs(  nlevp1) = 0.0_dp
         zbs(2*nlevp1) = 0.0_dp
       END IF

       ! Compute changes in T

       DO jn = 1,innp2
         !     lowest level
         ish = 2*ilt + jn*nlev           ! pointer  H
         ist =         jn*nlevp1         ! pointer  lnPS
         zbs(ist-1) = zralprr(nlev) * &  ! ==> T(nlev)
              ( zworks(ish) - &          ! H(nlev)
              rd*tr*zbs(ist) )           ! R*Tr*lnPS

         !     remaining levels
         ish = (jn-1)*nlev   + 2*ilt    ! pointer  H
         ist = (jn-1)*nlevp1            ! pointer  T
         DO  jk = nlevm1,1,-1
           zbs(ist+jk) =  zralprr(jk) * ( &                           ! ==> T(level)
                  zworks(ish+jk) &                                    ! H(level)
                - zworks(ish+jk+1) &                                  ! H(level+1)
                - zbs(ist+jk+1) * ( zrlnpr(jk+1) - zralphr(jk+1) ) &  ! T(level+1)
                )
         ENDDO
       ENDDO

       ! Compute new values

       IF (.NOT. lfilter) THEN
         !        zero surface pressure changes of all components
         !DO jn=1,innp2
         !  zbs(jn*nlevp1) = 0._dp
         !ENDDO
         ! zero change of global average surface pressure
         zbs(  nlevp1) = 0.0_dp
         zbs(2*nlevp1) = 0.0_dp
       ENDIF

       SELECT CASE(sel_modes)
       CASE(NMI_SEL_FILT, NMI_SEL_FILTD,&
            NMI_SEL_GM,   NMI_SEL_RM,   NMI_SEL_ALL,&
            NMI_SEL_GM_W0,NMI_SEL_RM_W0,NMI_SEL_ALL_W0)        ! store filtered fields
         svo_nmi(isp :isp +ilt -1) = zworks(    1:  ilt)
         sd_nmi (isp :isp +ilt -1) = zworks(ilt+1:2*ilt)
         SELECT CASE(sel_modes)
         CASE(NMI_SEL_GM_W0, NMI_SEL_RM_W0, NMI_SEL_ALL_W0)
           IF (jm==0) THEN
             ! global average T and lnPs unchanged
             stp_nmi(ismp+2*nlevp1:ismp+iltp-1) = zbs(1+2*nlevp1:iltp)
           ELSE
             stp_nmi(ismp:ismp+iltp-1) = zbs(1:iltp)
           END IF
         CASE default
           stp_nmi(ismp:ismp+iltp-1) = zbs(1:iltp)
         END SELECT

       CASE(NMI_SEL_FILT_W0, NMI_SEL_FILTD_W0, NMI_SEL_RM_IMPL) ! subtract the calculated fraction
         ! VO and D
         svo_nmi(isp :isp +ilt -1) = svo_nmi(isp :isp +ilt -1) - zworks(    1:ilt)
         sd_nmi (isp :isp +ilt -1) = sd_nmi (isp :isp +ilt -1) - zworks(ilt+1:2*ilt)
         ! ln(ps) and T
         stp_nmi(ismp:ismp+iltp-1) = stp_nmi(ismp:ismp+iltp-1) - zbs   (    1:iltp)
         
       END SELECT

       ! increment pointer in work space
       isp     = isp  + nnp(jm+1)*2*nlev
       ismp    = ismp + nnp(jm+1)*2*nlevp1

     END DO nmi_waveno_loop

    IF (isp /= ldc%snsp*2*nlev+1) &
      CALL finish('NMI_Filter','isp at end of nmi_waveno_loop inconsistent')

    ! clear work space
    DEALLOCATE (zbb, zbv, zworkv, zworkm, zbs, zworks)

  CONTAINS

    SUBROUTINE setidx ( jsym, jump, m, i1, i2, i3, i4, i5, i6 )
      ! set flow dependent index table
      INTEGER, INTENT(in)  :: jsym, jump, m
      INTEGER, INTENT(out) :: i1, i2, i3, i4, i5, i6

      i1 = (1-jsym)*jump      ! i1 ok
      i2 = 3 - 2*jsym         ! i2 ok
      i3 = jsym*jump          ! i3 ok
      i4 = 1 + jsym           ! i4 ok
      i5 = m + (1-jsym)       ! i5 ok
      i6 = m + jsym           ! i6 ok
      IF (m == 0) THEN
        i4 = 2 - jsym ! ok
        IF (jsym == 0) THEN  ! symmetric flow
          i2 = 1             ! i2 ok
          i3 = 2 * jump      ! i3 ok
          i6 = 2             ! i6 ok
        ELSE                 ! antisymmetric flow
          i1 = 2 * jump      ! i1 ok
          i2 = 3             ! i2 ok
          i5 = 2             ! i5 ok
        ENDIF
      ENDIF

    END SUBROUTINE setidx

  END SUBROUTINE NMI_Filter

!======================================================================
!BOP
  ! !IROUTINE: unitmx
  ! !INTERFACE:

  SUBROUTINE unitmx(array,dim)

    ! !DESCRIPTION:
    ! set a unit matrix

!BOX
    IMPLICIT NONE
!EOX
    ! !INPUT/OUTPUT PARAMETERS:
    REAL(dp), INTENT(inout) :: array(:,:)   ! out: unit matrix
    INTEGER, INTENT(in)     :: dim          ! order of matrix
!EOP
    INTEGER             :: j

    array(1:dim,1:dim) = 0.0_dp
    DO j=1,dim
      array(j,j) = 1.0_dp
    ENDDO

  END SUBROUTINE unitmx
!======================================================================
  FUNCTION df (m, n) RESULT (d)
    INTEGER, INTENT(IN) :: m, n
    REAL(dp) :: d

    d = SQRT(REAL((n**2-m**2),dp)/REAL((4*n**2-1),dp))

  END FUNCTION df
!======================================================================
  FUNCTION bf (m, n) RESULT (b)
    INTEGER, INTENT(IN) :: m, n
    REAL(dp) :: b

    b = df(m,n)*SQRT(REAL((n-1)*(n+1),dp)/REAL(n**2,dp))

  END FUNCTION bf
!======================================================================
  FUNCTION cf (m, n) RESULT (c)
    INTEGER, INTENT(IN) :: m, n
    REAL(dp) :: c

    c = REAL(m,dp)/REAL(n*(n+1),dp)

  END FUNCTION cf
!======================================================================
  FUNCTION ef (n, ped) RESULT (e)
    USE mo_physical_constants,     ONLY: earth_radius
    USE mo_control,                ONLY: earth_angular_velocity
    INTEGER, INTENT(IN)   :: n
    REAL(dp), INTENT(IN)  :: ped
    REAL(dp) :: e

    e = SQRT(ped*REAL(n*(n+1),dp))/(2.0_dp*earth_radius*earth_angular_velocity)

  END FUNCTION ef
!======================================================================
  FUNCTION hf (m, n) RESULT (h)
    USE mo_control,                ONLY: earth_angular_velocity
    INTEGER, INTENT(IN) :: n,m
    REAL(dp)  :: h

    h = 2.0_dp*earth_angular_velocity*REAL(m,dp)/REAL(n*(n+1),dp)

  END FUNCTION hf
!======================================================================
  FUNCTION of (n) RESULT (o)
    USE mo_control,                ONLY: earth_angular_velocity
    INTEGER, INTENT(IN) :: n
    REAL(dp) :: o

    o = 2.0_dp*earth_angular_velocity*SQRT(REAL(n**2,dp)-1._dp)/REAL(n,dp)

  END FUNCTION of
!======================================================================
  FUNCTION pf(n) RESULT (p)
    INTEGER, INTENT(IN) :: n
    REAL(dp) :: p

    p = SQRT(REAL(n*(n+1),dp))

  END FUNCTION pf
!======================================================================
  FUNCTION qf (n) RESULT (q)
    USE mo_physical_constants,     ONLY: earth_radius
    INTEGER, INTENT(IN)  :: n
    REAL(dp) :: q

    q = earth_radius**2/pf(n)

  END FUNCTION qf
!======================================================================
!BOP
  ! !IROUTINE: NMI_LoadBuffer
  ! !INTERFACE:

  SUBROUTINE NMI_LoadBuffer(buffer_action)
    ! !DESCRIPTION:
    ! load and reload different data sets into NMI work space
    ! !USES:
    USE mo_nudging_buffer,ONLY: &
         sdobs3,      svobs3,      stobs3, &
         sdobs,       svobs,       stobs,&
         sdfast_a,    svfast_a,    stfast_a, &
         sdfast_b,    svfast_b,    stfast_b, &
         sdfast_accu, svfast_accu, stfast_accu, &
         lfill_a, lfill_b, ifast_accu
    USE mo_memory_sp,     ONLY: svo, sd, stp
    USE mo_control,       ONLY: nlev, nlevp1
    USE mo_time_control,  ONLY: delta_time
    USE mo_spectral,      ONLY: corrsp
!BOX
    IMPLICIT NONE
!EOX
    ! !INPUT PARAMETERS:
    INTEGER :: buffer_action  ! define data storage direction
                              ! >0 load something into NMI arrays
                              ! <0 use NMI arrays as source
!EOP
    REAL(dp):: zrdt, corr_t, corr_p, corr_d, corr_v
    INTEGER :: ia1, ia2, snsp
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:), TARGET :: &
         wrk1_tp, wrk1_dvo, wrk2_tp, wrk2_dvo
    REAL(dp), POINTER :: p2a(:,:), p2b(:,:), p3a(:,:,:), p3b(:,:,:)

#if defined (__SX__) || defined (ES)
    EXTERNAL util_reshape
#else
    INTRINSIC RESHAPE
#endif

    snsp = ldc%snsp
    ia1  = nlev  *2*snsp
    ia2  = nlevp1*2*snsp

    SELECT CASE (buffer_action)

      !***** operations with buffers

    CASE(NMI_CORREL_BUFFER)          ! correlate buffer
      ALLOCATE(wrk1_tp(nlevp1,2,snsp), wrk1_dvo(nlev,2,snsp))
      ALLOCATE(wrk2_tp(nlevp1,2,snsp), wrk2_dvo(nlev,2,snsp))
#if defined (__SX__) || defined (ES)
      CALL util_reshape(wrk1_dvo, svo_nmi, ia1)
      CALL util_reshape(wrk2_dvo, svo_dt,  ia1)
#else
      wrk1_dvo(:,:,:) = RESHAPE(svo_nmi(:),(/nlev  ,2,snsp/))
      wrk2_dvo(:,:,:) = RESHAPE(svo_dt (:),(/nlev  ,2,snsp/))
#endif
      corr_v = corrsp(wrk1_dvo, wrk2_dvo)

#if defined (__SX__) || defined (ES)
      CALL util_reshape(wrk1_dvo, sd_nmi, ia1)
      CALL util_reshape(wrk2_dvo, sd_dt,  ia1)
#else
      wrk1_dvo(:,:,:) = RESHAPE(sd_nmi(:),(/nlev  ,2,snsp/))
      wrk2_dvo(:,:,:) = RESHAPE(sd_dt (:),(/nlev  ,2,snsp/))
#endif
      corr_d = corrsp(wrk1_dvo, wrk2_dvo)
#if defined (__SX__) || defined (ES)
      CALL util_reshape(wrk1_tp, stp_nmi, ia2)
      CALL util_reshape(wrk2_tp, stp_dt,  ia2)
#else
      wrk1_tp(:,:,:) = RESHAPE(stp_nmi(:),(/nlevp1,2,snsp/))
      wrk2_tp(:,:,:) = RESHAPE(stp_dt (:),(/nlevp1,2,snsp/))
#endif
      p3a => wrk1_tp(1:nlev,:,:)
      p3b => wrk2_tp(1:nlev,:,:)
      corr_t = corrsp(p3a,p3b)

      p2a => wrk1_tp(nlevp1,:,:)
      p2b => wrk2_tp(nlevp1,:,:)
      corr_p = corrsp(p2a,p2b)

      WRITE(nmi_mess,*) 'correlations NMI buffer: ',&
           'DIV= ',corr_d,' VOR= ',corr_v,' TEMP= ',corr_t, &
           ' LNPS= ',corr_p
      CALL message('NMI_LoadBuffer',nmi_mess)
      DEALLOCATE(wrk1_tp, wrk1_dvo)
      DEALLOCATE(wrk2_tp, wrk2_dvo)


    CASE (NMI_LOAD_CLEAR)            ! clear NMI data
      svo_nmi(:) = 0.0_dp
      sd_nmi (:) = 0.0_dp
      stp_nmi(:) = 0.0_dp


      !***** transfer   ATM <--> NMI (A)

    CASE (NMI_LOAD_ECH_BUFA)         ! load ECHAM spectral arrays
#if defined (__SX__) || defined (ES)
      CALL util_reshape(svo_nmi, svo,ia1)
      CALL util_reshape(sd_nmi,  sd, ia1)
      CALL util_reshape(stp_nmi, stp,ia2)
#else
      svo_nmi(:) = RESHAPE(svo(:,:,:),(/ia1/))
      sd_nmi (:) = RESHAPE(sd (:,:,:),(/ia1/))
      stp_nmi(:) = RESHAPE(stp(:,:,:),(/ia2/))
#endif

    CASE (-NMI_LOAD_ECH_BUFA)        ! reload ECHAM spectral arrays
#if defined (__SX__) || defined (ES)
      CALL util_reshape(svo, svo_nmi,ia1)
      CALL util_reshape(sd,  sd_nmi, ia1)
      CALL util_reshape(stp, stp_nmi,ia2)
#else
      svo(:,:,:) = RESHAPE(svo_nmi(:),(/nlev  ,2,snsp/))
      sd (:,:,:) = RESHAPE(sd_nmi (:),(/nlev  ,2,snsp/))
      stp(:,:,:) = RESHAPE(stp_nmi(:),(/nlevp1,2,snsp/))
#endif


      !***** transfer   ATM <--> DT (B)

    CASE (NMI_LOAD_ECH_BUFB)         ! load data in second buffer
#if defined (__SX__) || defined (ES)
      CALL util_reshape(svo_dt, svo,ia1)
      CALL util_reshape(sd_dt,  sd, ia1)
      CALL util_reshape(stp_dt, stp,ia2)
#else
      svo_dt(:) = RESHAPE(svo(:,:,:),(/ia1/))
      sd_dt (:) = RESHAPE(sd (:,:,:),(/ia1/))
      stp_dt(:) = RESHAPE(stp(:,:,:),(/ia2/))
#endif   

    CASE (-NMI_LOAD_ECH_BUFB)        ! reload ECHAM spectral arrays from second buffer
#if defined (__SX__) || defined (ES)
      CALL util_reshape(svo, svo_dt,ia1)
      CALL util_reshape(sd,  sd_dt, ia1)
      CALL util_reshape(stp, stp_dt,ia2)
#else
      svo(:,:,:) = RESHAPE(svo_dt(:),(/nlev  ,2,snsp/))
      sd (:,:,:) = RESHAPE(sd_dt (:),(/nlev  ,2,snsp/))
      stp(:,:,:) = RESHAPE(stp_dt(:),(/nlevp1,2,snsp/))
#endif


      !****** transfer  OBS <--> NMI (A)

    CASE (NMI_LOAD_ORD_BUFA)         ! load nudging data in first buffer
#if defined (__SX__) || defined (ES)
      CALL util_reshape(svo_nmi, svobs3,ia1)
      CALL util_reshape(sd_nmi,  sdobs3,ia1)
      CALL util_reshape(stp_nmi, stobs3,ia2)
#else
      svo_nmi(:) = RESHAPE(svobs3(:,:,:),(/ia1/))
      sd_nmi (:) = RESHAPE(sdobs3(:,:,:),(/ia1/))
      stp_nmi(:) = RESHAPE(stobs3(:,:,:),(/ia2/))
#endif

    CASE (-NMI_LOAD_ORD_BUFA)        ! reload nudging spectral arrays from first buffer
#if defined (__SX__) || defined (ES)
      CALL util_reshape(svobs3, svo_nmi,ia1)
      CALL util_reshape(sdobs3, sd_nmi, ia1)
      CALL util_reshape(stobs3, stp_nmi,ia2)
#else
      svobs3(:,:,:) = RESHAPE(svo_nmi(:),(/nlev  ,2,snsp/))
      sdobs3(:,:,:) = RESHAPE(sd_nmi (:),(/nlev  ,2,snsp/))
      stobs3(:,:,:) = RESHAPE(stp_nmi(:),(/nlevp1,2,snsp/))
#endif


      !****** transfer   OBS-INPUT <--> DT (B)

    CASE (NMI_LOAD_ORD_BUFB)         ! load nudging data in second buffer
#if defined (__SX__) || defined (ES)
      CALL util_reshape(svo_dt, svobs3,ia1)
      CALL util_reshape(sd_dt,  sdobs3,ia1)
      CALL util_reshape(stp_dt, stobs3,ia2)
#else
      svo_dt(:) = RESHAPE(svobs3(:,:,:),(/ia1/))
      sd_dt (:) = RESHAPE(sdobs3(:,:,:),(/ia1/))
      stp_dt(:) = RESHAPE(stobs3(:,:,:),(/ia2/))
#endif

    CASE (-NMI_LOAD_ORD_BUFB)        ! reload ECHAM spectral arrays from second buffer
#if defined (__SX__) || defined (ES)
      CALL util_reshape(svobs3, svo_dt,ia1)
      CALL util_reshape(sdobs3, sd_dt, ia1)
      CALL util_reshape(stobs3, stp_dt,ia2)
#else
      svobs3(:,:,:) = RESHAPE(svo_dt(:),(/nlev  ,2,snsp/))
      sdobs3(:,:,:) = RESHAPE(sd_dt (:),(/nlev  ,2,snsp/))
      stobs3(:,:,:) = RESHAPE(stp_dt(:),(/nlevp1,2,snsp/))
#endif


      !***** transfer   OBS-interpolated  <--> DT (B)

    CASE (NMI_LOAD_OIN_BUFB)         ! load nudging data in second buffer
#if defined (__SX__) || defined (ES)
      CALL util_reshape(svo_dt, svobs,ia1)
      CALL util_reshape(sd_dt,  sdobs,ia1)
      CALL util_reshape(stp_dt, stobs,ia2)
#else
      svo_dt(:) = RESHAPE(svobs(:,:,:),(/ia1/))
      sd_dt (:) = RESHAPE(sdobs(:,:,:),(/ia1/))
      stp_dt(:) = RESHAPE(stobs(:,:,:),(/ia2/))
#endif

    CASE (-NMI_LOAD_OIN_BUFB)        ! reload ECHAM spectral arrays from second buffer
#if defined (__SX__) || defined (ES)
      CALL util_reshape(svobs, svo_dt,ia1)
      CALL util_reshape(sdobs, sd_dt, ia1)
      CALL util_reshape(stobs, stp_dt,ia2)
#else
      svobs(:,:,:) = RESHAPE(svo_dt(:),(/nlev  ,2,snsp/))
      sdobs(:,:,:) = RESHAPE(sd_dt (:),(/nlev  ,2,snsp/))
      stobs(:,:,:) = RESHAPE(stp_dt(:),(/nlevp1,2,snsp/))
#endif


      !****** copy buffer   NMI (A) <--> DT (B)

    CASE (NMI_COPY_AB)        ! copy NMI buffer  TEN -->> NMI
      svo_nmi = svo_dt
      sd_nmi  = sd_dt 
      stp_nmi = stp_dt

    CASE (NMI_COPY_BA)        ! copy NMI buffer NMI -->> TEN
      svo_dt = svo_nmi
      sd_dt  = sd_nmi 
      stp_dt = stp_nmi


      !***** special operations

    CASE (NMI_CALC_TEND)      ! calculate tendencies
      zrdt = 1.0_dp / delta_time
      svo_nmi(:) = zrdt * (svo_nmi(:) - svo_dt(:))
      sd_nmi (:) = zrdt * (sd_nmi (:) - sd_dt (:))
      stp_nmi(:) = zrdt * (stp_nmi(:) - stp_dt(:))

    CASE (NMI_CALC_DIFF)       ! calculate difference
      svo_nmi(:) = (svo_dt(:) - svo_nmi(:))
      sd_nmi (:) = (sd_dt (:) - sd_nmi (:))
      stp_nmi(:) = (stp_dt(:) - stp_nmi(:))

    CASE (NMI_CALC_ADD)        ! add both NMI buffer, result in second buffer
      svo_dt(:) = svo_nmi(:) + svo_dt(:)
      sd_dt (:) = sd_nmi (:) + sd_dt (:)
      stp_dt(:) = stp_nmi(:) + stp_dt(:)


      !***** accumulation of filtered part
    CASE (NMI_STORE_FAST)      ! store filtered fields
      IF (lfill_a) THEN
        svfast_accu(:,:,:) = svfast_accu(:,:,:) - svfast_a(:,:,:)
        sdfast_accu(:,:,:) = sdfast_accu(:,:,:) - sdfast_a(:,:,:)
        stfast_accu(:,:,:) = stfast_accu(:,:,:) - stfast_a(:,:,:)
        ifast_accu = ifast_accu + 1
      END IF
      IF (lfill_b) THEN
        svfast_a(:,:,:) = svfast_b(:,:,:)
        sdfast_a(:,:,:) = sdfast_b(:,:,:)
        stfast_a(:,:,:) = stfast_b(:,:,:)
      END IF

#if defined (__SX__) || defined (ES)
      CALL util_reshape(svfast_b, svo_nmi,ia1)
      CALL util_reshape(sdfast_b, sd_nmi, ia1)
      CALL util_reshape(stfast_b, stp_nmi,ia2)
#else
      svfast_b(:,:,:) = RESHAPE(svo_nmi(:),(/nlev  ,2,snsp/))
      sdfast_b(:,:,:) = RESHAPE(sd_nmi (:),(/nlev  ,2,snsp/))
      stfast_b(:,:,:) = RESHAPE(stp_nmi(:),(/nlevp1,2,snsp/))
#endif
      IF (lfill_a) THEN
        svfast_accu(:,:,:) = svfast_accu(:,:,:) + svfast_b(:,:,:)
        sdfast_accu(:,:,:) = sdfast_accu(:,:,:) + sdfast_b(:,:,:)
        stfast_accu(:,:,:) = stfast_accu(:,:,:) + stfast_b(:,:,:)
      END IF
      IF (lfill_b) lfill_a = .TRUE.
      lfill_b = .TRUE.


    CASE default
      CALL message('NMI_LoadBuffer','illegal BUFFER_ACTION selection')
      CALL finish('NMI_LoadBuffer','Run aborted.')

    END SELECT

  END SUBROUTINE NMI_LoadBuffer

  SUBROUTINE assert(lvar, line)
    LOGICAL, INTENT(in) :: lvar
    INTEGER, INTENT(in) :: line
    IF (.NOT. lvar) THEN
      WRITE(message_text,'(a,i0,a)') 'assert at line ',line,' failed'
      CALL message('mo_nmi', message_text)
      CALL finish('mo_nmi','Internal error.')
    ENDIF
  END SUBROUTINE assert
  
END MODULE mo_nmi
