!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! This module contains all subroutines necessary to handle a CCN climatology
!!
!! @author 
!! <ol>
!! <li>S. Ferrachat (ETHZ)
!! </ol>
!!
!! @par Revision History
!! <ol>
!! <li>S. Ferrachat   (ETHZ) -  original code structure - (2010-03-xx) 
!!                                  
!! </ol>
!!
!! @par This module is used by
!! and to_be_added
!! 
!! @par Responsible coder
!! sylvaine.ferrachat@env.ethz.ch
!!
!! @par Copyright
!! 2010 by MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ECHAM is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!! violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!! copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!! an according license agreement with MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_ccnclim

  USE mo_kind,          ONLY: dp
  USE mo_linked_list,   ONLY: t_stream,  HYBRID, SURFACE

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_ccnclim_submodel
  PUBLIC :: ccnclim_define_tracer
  PUBLIC :: read_ccnclim
  PUBLIC :: ccn_3d
  PUBLIC :: ccnclim_avail_activ_lin_leaitch
  PUBLIC :: ccnclim_IN_setup

  INTEGER, PUBLIC :: idt_cdnc_ccnclim, idt_icnc_ccnclim

  REAL(dp), PUBLIC :: fracdusol     = 0.010_dp  !fraction of dust particules in all soluble modes
  REAL(dp), PUBLIC :: fracduai      = 0.094_dp  !fraction of dust particules in the insoluble accumulation mode
  REAL(dp), PUBLIC :: fracduci      = 0.229_dp  !fraction of dust particules in the insoluble coarse mode
  REAL(dp), PUBLIC :: fracbcsol     = 0.042_dp  !fraction of BC particules in all soluble modes
  REAL(dp), PUBLIC :: fracbcinsol   = 0._dp     !fraction of BC particules in all insoluble modes
  REAL(dp), PUBLIC :: rwetki        = 1._dp       !wet radius for the insoluble aitken mode (m) dummy value
  REAL(dp), PUBLIC :: rwetai        = 2.59e-08_dp !wet radius for the insoluble accumulation mode (m)
  REAL(dp), PUBLIC :: rwetci        = 7.97e-08_dp !wet radius for the insoluble coarse mode (m)
  REAL(dp), PUBLIC :: rwetas        = 1.09e-07_dp !wet radius for the soluble accumulation mode (m)

  TYPE (t_stream), PUBLIC, POINTER :: ccnclim

  REAL(dp), PUBLIC, POINTER :: ccn_an_strat_surf(:,:)
  REAL(dp), PUBLIC, POINTER :: ccn_an_strat_up(:,:)
  REAL(dp), PUBLIC, POINTER :: ccn_na_strat_surf(:,:)
  REAL(dp), PUBLIC, POINTER :: ccn_na_strat_up(:,:)
  REAL(dp), PUBLIC, POINTER :: ccn_an_conv_surf(:,:)
  REAL(dp), PUBLIC, POINTER :: ccn_an_conv_up(:,:)
  REAL(dp), PUBLIC, POINTER :: ccn_na_conv_surf(:,:)
  REAL(dp), PUBLIC, POINTER :: ccn_na_conv_up(:,:)
  REAL(dp), PUBLIC, POINTER :: ccn_tot_strat(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ccn_tot_conv(:,:,:)

CONTAINS

  !---------------------------------------------------------------------------
  !>
  !! @brief Initializes the CCN climatology submodel
  !!
  !! @remarks This subroutine reads into the ccnclim namelist
  !! and constructs the CCN climatology stream
  !! 

  SUBROUTINE init_ccnclim_submodel

    USE mo_memory_base,   ONLY: new_stream, add_stream_element,     &
                                default_stream_setting, add_stream_reference, AUTO
    USE mo_submodel,      ONLY: print_value
    USE mo_mpi,           ONLY: p_io, p_parallel_io, p_parallel, p_bcast
    USE mo_namelist,      ONLY: open_nml, position_nml, POSITIONED, &
                                 LENGTH_ERROR, READ_ERROR
    USE mo_filename,      ONLY: out_filetype
    USE mo_exception,     ONLY: finish, message, em_info
    USE mo_activ,         ONLY: nfrzmod
    USE mo_param_switches, ONLY: ncd_activ

    INTEGER :: ierr, inml, iunit

    INCLUDE 'ccnclimctl.inc'

    IF (p_parallel_io) THEN
       inml = open_nml('namelist.echam')
       iunit = position_nml ('CCNCLIMCTL', inml, status=ierr)
       SELECT CASE (ierr)
       CASE (POSITIONED)
          READ (iunit, ccnclimctl)
       CASE (LENGTH_ERROR)
          CALL finish ('init_ccnclim_submodel', &
                       'length error in namelist ccnclimctl')
       CASE (READ_ERROR)
          CALL finish ('init_ccnclim_submodel', &
                       'read error in namelist.echam')
       END SELECT
    END IF

    IF (p_parallel) THEN
       CALL p_bcast (fracdusol, p_io)
       CALL p_bcast (fracduai, p_io)
       CALL p_bcast (fracduci, p_io)
       CALL p_bcast (fracbcsol, p_io)
       CALL p_bcast (fracbcinsol, p_io)
       CALL p_bcast (rwetki, p_io)
       CALL p_bcast (rwetai, p_io)
       CALL p_bcast (rwetci, p_io)
       CALL p_bcast (rwetas, p_io)
    END IF

    IF (p_parallel_io) THEN
       CALL print_value('init_ccnclim_submodel, fract. of dust part. in sol modes', &
                       fracdusol)
       CALL print_value('init_ccnclim_submodel, fract. of dust part. in insol acc. mode', &
                       fracduai)
       CALL print_value('init_ccnclim_submodel, fract. of dust part. in insol coarse mode', &
                       fracduci)
       CALL print_value('init_ccnclim_submodel, fract. of BC part. in sol modes', &
                       fracbcsol)
       CALL print_value('init_ccnclim_submodel, fract. of BC part. in insol modes', &
                       fracbcinsol)
       CALL print_value('init_ccnclim_submodel, wet radius in aikten insoluble mode', &
                       rwetki)
       CALL print_value('init_ccnclim_submodel, wet radius in acc. insoluble mode', &
                       rwetai)
       CALL print_value('init_ccnclim_submodel, wet radius in coarse insoluble mode', &
                       rwetci)
       CALL print_value('init_ccnclim_submodel, wet radius in acc. soluble mode', &
                       rwetas)
    ENDIF

    !--- set the number of freezing modes:
    nfrzmod = 1

    !--- security checks
    IF (ncd_activ == 2) THEN
        ncd_activ = 1
        CALL message('init_ccn_submodel','ncd_activ reset to 1 because 2 is not supported with ccn climatology', &
                   level=em_info)
    ENDIF 

    !--- creates the ccnclim stream
    CALL new_stream (ccnclim, 'ccnclim', filetype = out_filetype, &
         post_suf = '_ccnclim', rest_suf = '_ccnclim')

    !--- add standard fields for post-processing:

    CALL add_stream_reference (ccnclim, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (ccnclim, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (ccnclim, 'aps'     ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (ccnclim, 'gboxarea','geoloc',lpost=.TRUE.)

    !--- add stream elements:

    !----- 2D fields: 
    CALL default_stream_setting (ccnclim, units     = 'number m-3',&
                                          lrerun    = .TRUE. ,     &
                                          laccu     = .FALSE. ,    &
                                          lpost     = .TRUE. ,     &
                                          leveltype = SURFACE,     &
                                          code      = AUTO         )

    CALL add_stream_element (ccnclim, 'CCN_AN_STRAT_SURF', ccn_an_strat_surf,  &
                            longname='anthrop. CCN for stratiform clouds (surface)')

    CALL add_stream_element (ccnclim, 'CCN_AN_STRAT_UP', ccn_an_strat_up,  &
                            longname='anthrop. CCN for stratiform clouds (300 hPa)')

    CALL add_stream_element (ccnclim, 'CCN_NA_STRAT_SURF', ccn_na_strat_surf,  &
                            longname='natural CCN for stratiform clouds (surface)')

    CALL add_stream_element (ccnclim, 'CCN_NA_STRAT_UP', ccn_na_strat_up,  &
                            longname='natural CCN for stratiform clouds (300 hPa)')

    CALL add_stream_element (ccnclim, 'CCN_AN_CONV_SURF', ccn_an_conv_surf,  &
                            longname='anthrop. CCN for convective clouds (surface)')

    CALL add_stream_element (ccnclim, 'CCN_AN_CONV_UP', ccn_an_conv_up,  &
                            longname='anthrop. CCN for convective clouds (300 hPa)')

    CALL add_stream_element (ccnclim, 'CCN_NA_CONV_SURF', ccn_na_conv_surf,  &
                            longname='natural CCN for convective clouds (surface)')

    CALL add_stream_element (ccnclim, 'CCN_NA_CONV_UP', ccn_na_conv_up,  &
                            longname='natural CCN for convective clouds (300 hPa)')

    !----- 3D fields:
    CALL default_stream_setting (ccnclim, leveltype = HYBRID)

    CALL add_stream_element (ccnclim, 'CCN_TOT_STRAT', ccn_tot_strat,    &
                            longname='total CCN for stratiform clouds')

    CALL add_stream_element (ccnclim, 'CCN_TOT_CONV', ccn_tot_conv,      &
                            longname='total CCN for convective clouds')   

  END SUBROUTINE init_ccnclim_submodel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! ccnclim_define_tracer: create ECHAM tracers for CDNC and ICNC
!!
!! @author see module info
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see module info
!!
!! @par This subroutine is called by
!! init_submodels
!!
!! @par Externals:
!! <ol>
!! <li>none
!! </ol>
!!
!! @par Responsible coder
!! sylvaine.ferrachat@env.ethz.ch
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ccnclim_define_tracer

    USE mo_exception,         ONLY: message, message_text, em_debug, em_error
    USE mo_tracer,            ONLY: new_tracer, new_diag_burden
    USE mo_tracdef,           ONLY: OFF, ON, OK

    IMPLICIT NONE

    INTEGER :: iburdenid, ierr

    !---executable procedure

    CALL message('ccnclim_define_tracer', 'Starting ...', level=em_debug)

    iburdenid = new_diag_burden('CDNC', itype=1, lclobber=.false.)

    CALL new_tracer('CDNC', 'CCNCLIM', ierr=ierr,      &
                     idx=idt_cdnc_ccnclim,         &
                     units='1 kg-1',               &
                     table=131,                    &
                     code=131,                     &
                     nwrite=1,                     &
                     burdenid=iburdenid,           &
                     nconv=OFF,                    &
              !!     nconvmassfix=OFF,             &
                     nvdiff=ON,                    &
                     nint=ON,                      &
                     longname='cloud droplet number concentration')

    IF (ierr /= OK) THEN
      WRITE(message_text,'(a,i0)') 'new_tracer CDNC returned error code ',ierr
      CALL message('ham_define_tracer', message_text, level=em_error)
    END IF

    iburdenid = new_diag_burden('ICNC', itype=1, lclobber=.false.)

    CALL new_tracer('ICNC', 'CCNCLIM', ierr=ierr,      &
                     idx=idt_icnc_ccnclim,         &
                     units='1 kg-1',               &
                     table=131,                    &
                     code=132,                     &
                     nwrite=1,                     &
                     burdenid=iburdenid,           &
                     nconv=OFF,                    &
              !!     nconvmassfix=OFF,             &
                     nvdiff=ON,                    &
                     nint=ON,                      &
                     longname='ice crystal number concentration')

    IF (ierr /= OK) THEN
      WRITE(message_text,'(a,i0)') 'new_tracer ICNC returned error code ',ierr
      CALL message('ham_define_tracer', message_text, level=em_error)
    END IF

  END SUBROUTINE ccnclim_define_tracer

  !---------------------------------------------------------------------------
  !>
  !! @brief Read monthly CCN values
  !! 
  !! @remarks This routine reads surface and upper troposphere
  !! CCN number concentrations
  !!
 
  SUBROUTINE read_ccnclim(kyear,kmonth)

    USE mo_io_units,          ONLY: nout
    USE mo_mpi,               ONLY: p_parallel_io
    USE mo_decomposition,     ONLY: gc => global_decomposition,   &
                                    dc => local_decomposition

    USE mo_transpose,         ONLY: scatter_gp
    USE mo_read_netcdf77,     ONLY: read_var_hs_nf77_2d

    INTEGER, INTENT(IN)             :: kyear
    INTEGER, INTENT(IN)             :: kmonth
    INTEGER, PARAMETER              :: nfield=8
    CHARACTER(len=10)     :: filename
    CHARACTER(len=256)    :: yvar_name(nfield)

    INTEGER                           :: i, ierr, jr
    INTEGER                           :: nproma

    REAL(dp) ,POINTER     :: field(:,:,:)
    REAL(dp) ,ALLOCATABLE :: lfield(:,:,:)


    IF (p_parallel_io) THEN
       ALLOCATE (field (dc% nlon,nfield,dc% nlat))
    ENDIF

    WRITE(filename,'(''ccn'',i4.4,''.nc'')') kyear

    !--- 1) Read monthly CCN field:

    yvar_name(1)   = 'CCN_an_strat_mo_surf'
    yvar_name(2)   = 'CCN_na_strat_mo_surf'
    yvar_name(3)   = 'CCN_an_strat_mo_up'
    yvar_name(4)   = 'CCN_na_strat_mo_up'
    yvar_name(5)   = 'CCN_an_conv_mo_surf'
    yvar_name(6)   = 'CCN_na_conv_mo_surf'
    yvar_name(7)   = 'CCN_an_conv_mo_up'
    yvar_name(8)   = 'CCN_na_conv_mo_up'

    IF (p_parallel_io) THEN
       DO i=1, nfield
         CALL read_var_hs_nf77_2d (filename, 'lon', 'lat', 'time', kmonth,     &
                   TRIM(yvar_name(i)), field(:,i,:), ierr)
       END DO
    WRITE (nout,*) '------------------------------------------------------------------------------'
    WRITE (nout,*) '--- Initialization of CCN -------------------------------------'
    WRITE (nout,*)
    WRITE (nout,*) '   read ',filename,': '
    WRITE (nout,*) '      reading CCN_an_strat_mo_surf'
    WRITE (nout,*) '      reading CCN_na_strat_mo_surf'
    WRITE (nout,*) '      reading CCN_an_strat_mo_up'
    WRITE (nout,*) '      reading CCN_na_strat_mo_up'
    WRITE (nout,*) '      reading CCN_an_conv_mo_surf'
    WRITE (nout,*) '      reading CCN_na_conv_mo_surf'
    WRITE (nout,*) '      reading CCN_an_conv_mo_up'
    WRITE (nout,*) '      reading CCN_na_conv_mo_up'
    WRITE (nout,*) '------------------------------------------------------------------------------'
    END IF
  ALLOCATE (lfield(dc%nproma, nfield, dc% ngpblks))
  CALL scatter_gp (field, lfield, gc)

  DO jr=1,dc%ngpblks

     IF ( jr == dc% ngpblks ) THEN
        nproma = dc% npromz
     ELSE
        nproma = dc% nproma
     END IF

     ccn_an_strat_surf(1:nproma,jr) = lfield(1:nproma,1,jr)
     ccn_na_strat_surf(1:nproma,jr) = lfield(1:nproma,2,jr)
     ccn_an_strat_up(1:nproma,jr)   = lfield(1:nproma,3,jr)
     ccn_na_strat_up(1:nproma,jr)   = lfield(1:nproma,4,jr)
     ccn_an_conv_surf(1:nproma,jr)  = lfield(1:nproma,5,jr)
     ccn_na_conv_surf(1:nproma,jr)  = lfield(1:nproma,6,jr)
     ccn_an_conv_up(1:nproma,jr)    = lfield(1:nproma,7,jr)
     ccn_na_conv_up(1:nproma,jr)    = lfield(1:nproma,8,jr)

  ENDDO

  DEALLOCATE (lfield)
  IF (p_parallel_io) DEALLOCATE (field)
  RETURN

  END SUBROUTINE read_ccnclim

  !---------------------------------------------------------------------------
  !>
  !! @brief Computes 3D CCN values
  !! 
  !! @remarks This routine takes surface and upper troposphere
  !! CCN number concentrations and creates a 3D field
  !! by computing a corresponding vertical profile for each location
  !!

  SUBROUTINE ccn_3d(kproma, kbdim, klev, krow, papm1)

    INTEGER,  INTENT(in) :: kproma, kbdim, klev, krow
    REAL(dp), INTENT(in) :: papm1(kbdim,klev)

    !--- Local vars:

    REAL(dp) :: zccnsurf(kbdim), zccnup(kbdim), zccntot(kbdim,klev)

    !--- strat case:

    zccnsurf(1:kproma) = ccn_an_strat_surf(1:kproma,krow) &
                       + ccn_na_strat_surf(1:kproma,krow)

    zccnup(1:kproma)   = ccn_an_strat_up(1:kproma,krow) &
                       + ccn_na_strat_up(1:kproma,krow)

    CALL ccn_profile(kproma, kbdim, klev, papm1, zccnsurf, zccnup, zccntot)

    ccn_tot_strat(1:kproma,:,krow) = zccntot(1:kproma,:)

    !--- conv case:

    zccnsurf(1:kproma) = ccn_an_conv_surf(1:kproma,krow) &
                       + ccn_na_conv_surf(1:kproma,krow)

    zccnup(1:kproma)   = ccn_an_conv_up(1:kproma,krow) &
                       + ccn_na_conv_up(1:kproma,krow)

    CALL ccn_profile(kproma, kbdim, klev, papm1, zccnsurf, zccnup, zccntot)

    ccn_tot_conv(1:kproma,:,krow) = zccntot(1:kproma,:)

  END SUBROUTINE ccn_3d
 
  !---------------------------------------------------------------------------
  !>
  !! @brief Utility routine to compute the vertical profiles
  !! 
  !! @remarks The calculation reproduces what is implemented for
  !! acdnc setup in case of fully prescribed CDNC

  SUBROUTINE ccn_profile(kproma, kbdim, klev, papm1, pccnsurf, pccnup, pccn)

    INTEGER,  INTENT(in)  :: kproma, kbdim, klev
    REAL(dp), INTENT(in)  :: papm1(kbdim,klev)
    REAL(dp), INTENT(in)  :: pccnsurf(kbdim), pccnup(kbdim)
    REAL(dp), INTENT(out) :: pccn(kbdim,klev)

    !--- Local variables:

    INTEGER  :: jexp, & !exponent
                jk

    REAL(dp) :: zpthresh,          & !pressure threshold (boundary layer)
                zprat(kbdim,klev), &
                ztmp1(kbdim,klev), &
                ztmp2(kbdim,klev)

    LOGICAL  :: lobl(kbdim,klev)

    jexp = 2
    zpthresh = 80000._dp

    zprat(1:kproma,:) = ( MIN(8._dp, zpthresh/papm1(1:kproma,:)) )**jexp

    lobl(1:kproma,:)  = (papm1(1:kproma,:) < zpthresh)

    DO jk=1,klev
       ztmp1(1:kproma,jk) = pccnup(1:kproma) &
                          + (pccnsurf(1:kproma)-pccnup(1:kproma))*(EXP(1._dp-zprat(1:kproma,jk)))

       ztmp2(1:kproma,jk) = pccnsurf(1:kproma)
    ENDDO

    pccn(1:kproma,:) = MERGE(ztmp1(1:kproma,:), ztmp2(1:kproma,:), lobl(1:kproma,:))

  END SUBROUTINE ccn_profile

  !---------------------------------------------------------------------------
  !>
  !! @brief Computes available particules for activation
  !! 
  !! @remarks Preparatory routine for Lin&Leaitch activation scheme

  SUBROUTINE ccnclim_avail_activ_lin_leaitch(kproma, kbdim, klev, krow)

    USE mo_activ, ONLY: na
    USE mo_conv,  ONLY: na_cv

    INTEGER, INTENT(IN) :: kproma, kbdim, klev, krow

    na(1:kproma,:,krow)    = ccn_tot_strat(1:kproma,:,krow) ![m-3]
    na_cv(1:kproma,:,krow) = ccn_tot_conv(1:kproma,:,krow)  ![m-3] 

  END SUBROUTINE ccnclim_avail_activ_lin_leaitch

  !---------------------------------------------------------------------------
  !>
  !! @brief Routine to set up necessary vars for mixed-phase and cirrus freezing calculations
  !! 
  !! @remarks Here, some quantities are surrogates for their m7 modes counterparts
  !! and for this reason they keep names refering to m7 modes.
  !! There is no dependency to HAM and m7.
  !! These quantities were diagnosed in a preliminary full HAM run
  !!

  SUBROUTINE ccnclim_IN_setup(kproma, kbdim, klev, krow,                                  &
                              prho, prho_rcp,                                             &
                              prwetki, prwetai, prwetci,                                  &
                              pfracdusol, pfracduai, pfracduci, pfracbcsol, pfracbcinsol, &
                              pascs, papnx, paprx, papsigx, ld_het                        )

    USE mo_activ, ONLY: nfrzmod

    INTEGER,  INTENT(in)  :: kproma, kbdim, klev, krow
    REAL(dp), INTENT(in)  :: prho(kbdim,klev)     ! air density
    REAL(dp), INTENT(in)  :: prho_rcp(kbdim,klev) ! reciprocal air density
    !-- for mixed-phase freezing:
    REAL(dp), INTENT(out) :: prwetki(kbdim,klev)  ! wet radius, aitken insoluble mode
    REAL(dp), INTENT(out) :: prwetai(kbdim,klev)  ! wet radius, accumulation insoluble mode
    REAL(dp), INTENT(out) :: prwetci(kbdim,klev)  ! wet radius, coarse insoluble mode
    REAL(dp), INTENT(out) :: pfracdusol(kbdim,klev)   ! total fraction of dust particules in all soluble modes
    REAL(dp), INTENT(out) :: pfracduai(kbdim,klev)    ! fraction of dust particules in the accum. soluble mode 
    REAL(dp), INTENT(out) :: pfracduci(kbdim,klev)    ! fraction of dust particules in the coarse soluble mode 
    REAL(dp), INTENT(out) :: pfracbcsol(kbdim,klev)   ! total fraction of BC particules in all soluble modes 
    REAL(dp), INTENT(out) :: pfracbcinsol(kbdim,klev) ! total fraction of BC particules in all insoluble modes 
    !-- for cirrus freezing:
    REAL(dp), INTENT(out) :: pascs(kbdim,klev)           ! soluble aerosol number conc. 
    REAL(dp), INTENT(out) :: papnx(kbdim,klev,nfrzmod)   ! aerosol number available for freezing [1/cm3] 
    REAL(dp), INTENT(out) :: paprx(kbdim,klev,nfrzmod)   ! radius of aerosols avail. for freezing  [cm] 
    REAL(dp), INTENT(out) :: papsigx(kbdim,klev,nfrzmod) ! std. dev. of aerosols available for freezing
    LOGICAL, INTENT(out)  :: ld_het  !switch to set heterogeneous freezing below 235K (cirrus scheme)

    !--- Mixed-phase freezing setup:

     !-- Wet radii:
     prwetki(1:kproma,:) = rwetki
     prwetai(1:kproma,:) = rwetai
     prwetci(1:kproma,:) = rwetci

     !-- Various dust and BC fractions:
     pfracdusol(1:kproma,:)   = fracdusol
     pfracduai(1:kproma,:)    = fracduai
     pfracduci(1:kproma,:)    = fracduci
     pfracbcsol(1:kproma,:)   = fracbcsol
     pfracbcinsol(1:kproma,:) = fracbcinsol

    !--- Cirrus freezing setup:

     !-- Soluble aerosol number concentration:
     pascs(1:kproma,:) = ccn_tot_strat(1:kproma,:,krow)*prho_rcp(1:kproma,:) ![kg-1]
     pascs(1:kproma,:) = MAX(pascs(1:kproma,:), 10.E6_dp)

     !-- Number, radius and std. dev. of aerosols available for cirrus freezing
     !   note: the number will be later reduced by the actual icnc in the cloud routine

     ld_het = .false.

     papnx(1:kproma,:,1) = prho(1:kproma,:) * pascs(1:kproma,:) ![1/m3]

     paprx(1:kproma,:,1) = 100._dp * rwetas ![cm]
     paprx(1:kproma,:,1) = MAX(paprx(1:kproma,:,1), 0.05E-4_dp)

     papsigx(1:kproma,:,1) = 1._dp

  END SUBROUTINE ccnclim_IN_setup

END MODULE mo_ccnclim
