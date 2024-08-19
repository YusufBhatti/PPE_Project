!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! mo_diag_tendency_new: diagnostic of various tendencies
!!
!! @author Sebastian Rast, MPI Met, Hamburg
!!
!! $ID: n/a$
!!
!! @par Revision History
!! original source by J.S.Rast (2010-12-16), partly based on code 
!! by I. Kirchner 
!! added hdiff diagnostic by J.S.Rast (2012-10-12)
!!
MODULE mo_diag_tendency_new
  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: finish, message
  USE mo_memory_base,             ONLY: new_stream, add_stream_element, &
                                        default_stream_setting, &
                                        set_stream_element_info, &
                                        add_stream_reference, AUTO, t_stream
  USE mo_decomposition,           ONLY: lc=>local_decomposition
  USE mo_control,                 ONLY: nsp
  USE mo_linked_list,             ONLY: SPECTRAL
  USE mo_time_event,              ONLY: io_time_event
  USE mo_time_control,            ONLY: p_bcast_event
  USE mo_mpi,                     ONLY: p_parallel_io, p_bcast, p_io
  USE mo_namelist,                ONLY: position_nml, open_nml, &
                                        MISSING, POSITIONED
  USE mo_util_string,             ONLY: tolower
  USE mo_string_utls,             ONLY: st1_in_st2_proof,st1_in_st2_idx

  IMPLICIT NONE

  TYPE t_tdiag_vars
    REAL(wp), POINTER :: dudt_hines(:,:,:),dvdt_hines(:,:,:),dtdt_hines(:,:,:)
    REAL(wp), POINTER :: dudt_sso(:,:,:),dvdt_sso(:,:,:),dtdt_sso(:,:,:)
    REAL(wp), POINTER :: dtdt_rheat_sw(:,:,:),dtdt_rheat_lw(:,:,:)
    REAL(wp), POINTER :: dudt_vdiff(:,:,:),dvdt_vdiff(:,:,:), &
                         dtdt_vdiff(:,:,:),dqdt_vdiff(:,:,:), &
                         dxldt_vdiff(:,:,:),dxidt_vdiff(:,:,:)
    REAL(wp), POINTER :: dudt_cucall(:,:,:),dvdt_cucall(:,:,:), &
                         dtdt_cucall(:,:,:),dqdt_cucall(:,:,:)
    REAL(wp), POINTER :: dtdt_cloud(:,:,:),dqdt_cloud(:,:,:), &
                         dxldt_cloud(:,:,:),dxidt_cloud(:,:,:)
    REAL(wp), POINTER :: dsddt_hdiff(:,:,:),dsvodt_hdiff(:,:,:), &
                         dstdt_hdiff(:,:,:)
  END TYPE t_tdiag_vars

  PRIVATE
  PUBLIC :: init_tdiag, set_tendency
  PUBLIC :: tdiag_vars, t_tdiag_vars, t3d_physc
  ! accumulate tendencies of specific process
  INTERFACE set_tendency
    MODULE PROCEDURE set_tendency_gp2d
    MODULE PROCEDURE set_tendency_sp3d
  END INTERFACE set_tendency
  ! tdiag stream
  INTEGER, PARAMETER  :: n_tdiag_vars=25
  CHARACTER (LEN=32)  :: c_tdiag_vars(1:n_tdiag_vars)=       &
                       (/'dudt_hines                      ', &
                         'dvdt_hines                      ', &
                         'dtdt_hines                      ', &
                         'dudt_sso                        ', &
                         'dvdt_sso                        ', &
                         'dtdt_sso                        ', &
                         'dtdt_rheat_sw                   ', &
                         'dtdt_rheat_lw                   ', &
                         'dudt_vdiff                      ', &
                         'dvdt_vdiff                      ', &
                         'dtdt_vdiff                      ', &
                         'dqdt_vdiff                      ', &
                         'dxldt_vdiff                     ', &
                         'dxidt_vdiff                     ', &
                         'dudt_cucall                     ', &
                         'dvdt_cucall                     ', &
                         'dtdt_cucall                     ', &
                         'dqdt_cucall                     ', &
                         'dtdt_cloud                      ', &
                         'dqdt_cloud                      ', &
                         'dxldt_cloud                     ', &
                         'dxidt_cloud                     ', &
                         'dsddt_hdiff                     ', &
                         'dsvodt_hdiff                    ', &
                         'dstdt_hdiff                     '  &
                       /)
  CHARACTER (LEN=32)  :: tdiagnam(n_tdiag_vars)
  TYPE t_varlist
    INTEGER           :: n
    CHARACTER(LEN=32) :: var(n_tdiag_vars)
  END TYPE t_varlist

  TYPE(io_time_event), SAVE     :: puttdiag

  TYPE (t_tdiag_vars)           :: tdiag_vars  !> contains all tendency 
                                               !> diagnostic variables
  REAL (wp), POINTER  :: t3d_physc !> 3d temperature from physc
  TYPE (t_varlist)    :: vdiff_list      ,hdiff_list      ,radheat_list       &
                        ,gwspectrum_list ,ssodrag_list    ,cucall_list        &
                        ,cloud_list      ,uwind_list      ,vwind_list         &
                        ,temp_list       ,qhum_list       ,xl_list            &
                        ,xi_list
  REAL(wp), PARAMETER :: daylength=86400._wp
  
CONTAINS

  SUBROUTINE init_tdiag

    CHARACTER(LEN=32)             :: ichar
    CHARACTER(LEN=32)             :: c_vars(1:n_tdiag_vars)
    INTEGER                       :: i, j, n, ierr, n_end
    LOGICAL                       :: lpost
    TYPE (t_stream), POINTER      :: tdiag

    INTEGER :: inml, iunit

    INCLUDE 'tdiagctl.inc'

    !> set default output interval
    puttdiag%counter      = 6
    puttdiag%unit         = 'hours'
    puttdiag%adjustment   = 'first'
    puttdiag%offset       = 0
    !> set default list of diagnostic variables
    tdiagnam(1)='all                             '
    DO i=2,n_tdiag_vars
      tdiagnam(i)='end                             '
    END DO
    !> set list associated with key words
    n=6
    vdiff_list%n=n
    vdiff_list%var(1:n)=(/'dudt_vdiff                      ', &
                          'dvdt_vdiff                      ', &
                          'dtdt_vdiff                      ', &
                          'dqdt_vdiff                      ', &
                          'dxldt_vdiff                     ', &
                          'dxidt_vdiff                     '  &
                        /)
    n=3
    hdiff_list%n=n
    hdiff_list%var(1:n)=(/'dsddt_hdiff                     ', &
                          'dsvodt_hdiff                    ', &
                          'dstdt_hdiff                     '  &
                        /)
    n=2
    radheat_list%n=n
    radheat_list%var(1:n)=(/'dtdt_rheat_sw                   ', &
                            'dtdt_rheat_lw                   '  &
                          /)
    n=3
    gwspectrum_list%n=n
    gwspectrum_list%var(1:n)=(/'dudt_hines                      ', &
                               'dvdt_hines                      ', &
                               'dtdt_hines                      '  &
                             /)
    n=3
    ssodrag_list%n=n
    ssodrag_list%var(1:n)=(/'dudt_sso                        ', &
                            'dvdt_sso                        ', &
                            'dtdt_sso                        '  &
                          /)
    n=4
    cucall_list%n=n
    cucall_list%var(1:n)=(/'dudt_cucall                     ', &
                           'dvdt_cucall                     ', &
                           'dtdt_cucall                     ', &
                           'dqdt_cucall                     '  &
                         /)
    n=4
    cloud_list%n=n
    cloud_list%var(1:n)=(/'dtdt_cloud                      ', &
                          'dqdt_cloud                      ', &
                          'dxldt_cloud                     ', &
                          'dxidt_cloud                     '  &
                        /)
    n=6
    uwind_list%n=n
    uwind_list%var(1:n)=(/'dudt_vdiff                      ', &
                     'dudt_hines                      ', &
                     'dudt_sso                        ', &
                     'dudt_cucall                     ', &
                     'dsddt_hdiff                     ', &
                     'dsvodt_hdiff                    '  &
                   /)
    n=6
    vwind_list%n=n
    vwind_list%var(1:n)=(/'dvdt_vdiff                      ', &
                          'dvdt_hines                      ', &
                          'dvdt_sso                        ', &
                          'dvdt_cucall                     ', &
                          'dsddt_hdiff                     ', &
                          'dsvodt_hdiff                    '  &
                        /)
    n=8
    temp_list%n=n
    temp_list%var(1:n)=(/'dtdt_vdiff                      ', &
                         'dstdt_hdiff                     ', &
                         'dtdt_rheat_sw                   ', &
                         'dtdt_rheat_lw                   ', &
                         'dtdt_hines                      ', &
                         'dtdt_sso                        ', &
                         'dtdt_cucall                     ', &
                         'dtdt_cloud                      '  &
                       /)
    n=3
    qhum_list%n=n
    qhum_list%var(1:n)=(/'dqdt_vdiff                      ', &
                         'dqdt_cucall                     ', &
                         'dqdt_cloud                      '  &
                       /)
    n=2
    xl_list%n=n
    xl_list%var(1:n)=(/'dxldt_vdiff                     ', &
                       'dxldt_cloud                     '  &
                     /)
    n=2
    xi_list%n=n
    xi_list%var(1:n)=(/'dxidt_vdiff                     ', &
                       'dxidt_cloud                     '  &
                     /)
    IF (p_parallel_io) THEN
      inml = open_nml ('namelist.echam')
      iunit = position_nml ('TDIAGCTL', inml, status=ierr)
      SELECT CASE (ierr)
      CASE (MISSING)
        CALL message('mo_diag_tendency_new, tdiag_init_stream', &
                    'no namelist tdiagctl in file namelist.echam, use default')
      CASE (POSITIONED)
        READ(iunit, tdiagctl)
      CASE DEFAULT
        CALL finish('mo_diag_tendency_new, tdiag_init_stream', &
                    'error in reading namelist tdiagctl')
      END SELECT
    END IF
    CALL p_bcast_event (puttdiag, p_io)
    CALL p_bcast (tdiagnam, p_io)
    !> Handle collective key words
    n_end=0
    DO i=1,n_tdiag_vars
      IF (n_end == n_tdiag_vars) EXIT  ! In this case, all variables are in list
      ichar=tolower(tdiagnam(i))
      SELECT CASE (TRIM(ichar))
!      SELECT CASE (ichar)
      CASE ('all')
        n_end=n_tdiag_vars
        c_vars(1:n_tdiag_vars)=c_tdiag_vars(1:n_tdiag_vars)
      CASE ('vdiff')
        DO j=1,vdiff_list%n
          IF (n_end == n_tdiag_vars) EXIT
          IF (st1_in_st2_idx(TRIM(vdiff_list%var(j)),c_vars(1:n_end))==0) THEN
            n_end=n_end+1
            c_vars(n_end)=vdiff_list%var(j)
          END IF
        END DO
      CASE ('hdiff')
        DO j=1,hdiff_list%n
          IF (n_end == n_tdiag_vars) EXIT
          IF (st1_in_st2_idx(TRIM(hdiff_list%var(j)),c_vars(1:n_end))==0) THEN
            n_end=n_end+1
            c_vars(n_end)=hdiff_list%var(j)
          END IF
        END DO
      CASE ('radheat')
        DO j=1,radheat_list%n
          IF (n_end == n_tdiag_vars) EXIT
          IF (st1_in_st2_idx(TRIM(radheat_list%var(j)),c_vars(1:n_end))==0) THEN
            n_end=n_end+1
            c_vars(n_end)=radheat_list%var(j)
          END IF
        END DO
      CASE ('gwspectrum')
        DO j=1,gwspectrum_list%n
          IF (n_end == n_tdiag_vars) EXIT
          IF (st1_in_st2_idx(TRIM(gwspectrum_list%var(j)),c_vars(1:n_end))==0) THEN
            n_end=n_end+1
            c_vars(n_end)=gwspectrum_list%var(j)
          END IF
        END DO
      CASE ('ssodrag')
        DO j=1,ssodrag_list%n
          IF (n_end == n_tdiag_vars) EXIT
          IF (st1_in_st2_idx(TRIM(ssodrag_list%var(j)),c_vars(1:n_end))==0) THEN
            n_end=n_end+1
            c_vars(n_end)=ssodrag_list%var(j)
          END IF
        END DO
      CASE ('cucall')
        DO j=1,cucall_list%n
          IF (n_end == n_tdiag_vars) EXIT
          IF (st1_in_st2_idx(TRIM(cucall_list%var(j)),c_vars(1:n_end))==0) THEN
            n_end=n_end+1
            c_vars(n_end)=cucall_list%var(j)
          END IF
        END DO
      CASE ('cloud')
        DO j=1,cloud_list%n
          IF (n_end == n_tdiag_vars) EXIT
          IF (st1_in_st2_idx(TRIM(cloud_list%var(j)),c_vars(1:n_end))==0) THEN
            n_end=n_end+1
            c_vars(n_end)=cloud_list%var(j)
          END IF
        END DO
      CASE ('uwind')
        DO j=1,uwind_list%n
          IF (n_end == n_tdiag_vars) EXIT
          IF (st1_in_st2_idx(TRIM(uwind_list%var(j)),c_vars(1:n_end))==0) THEN
            n_end=n_end+1
            c_vars(n_end)=uwind_list%var(j)
          END IF
        END DO
      CASE ('vwind')
        DO j=1,vwind_list%n
          IF (n_end == n_tdiag_vars) EXIT
          IF (st1_in_st2_idx(TRIM(vwind_list%var(j)),c_vars(1:n_end))==0) THEN
            n_end=n_end+1
            c_vars(n_end)=vwind_list%var(j)
          END IF
        END DO
      CASE ('temp')
        DO j=1,temp_list%n
          IF (n_end == n_tdiag_vars) EXIT
          IF (st1_in_st2_idx(TRIM(temp_list%var(j)),c_vars(1:n_end))==0) THEN
            n_end=n_end+1
            c_vars(n_end)=temp_list%var(j)
          END IF
        END DO
      CASE ('qhum')
        DO j=1,qhum_list%n
          IF (n_end == n_tdiag_vars) EXIT
          IF (st1_in_st2_idx(TRIM(qhum_list%var(j)),c_vars(1:n_end))==0) THEN
            n_end=n_end+1
            c_vars(n_end)=qhum_list%var(j)
          END IF
        END DO
      CASE ('xl')
        DO j=1,xl_list%n
          IF (n_end == n_tdiag_vars) EXIT
          IF (st1_in_st2_idx(TRIM(xl_list%var(j)),c_vars(1:n_end))==0) THEN
            n_end=n_end+1
            c_vars(n_end)=xl_list%var(j)
          END IF
        END DO
      CASE ('xi')
        DO j=1,xi_list%n
          IF (n_end == n_tdiag_vars) EXIT
          IF (st1_in_st2_idx(TRIM(xi_list%var(j)),c_vars(1:n_end))==0) THEN
            n_end=n_end+1
            c_vars(n_end)=xi_list%var(j)
          END IF
        END DO
      CASE ('dudt_vdiff','dvdt_vdiff','dtdt_vdiff','dqdt_vdiff','dxldt_vdiff', &
            'dxidt_vdiff','dsddt_hdiff','dsvodt_hdiff','dstdt_hdiff', &
            'dtdt_rheat_sw','dtdt_rheat_lw','dudt_hines', &
            'dvdt_hines','dtdt_hines','dudt_sso','dvdt_sso','dtdt_sso', &
            'dudt_cucall','dvdt_cucall','dtdt_cucall','dqdt_cucall', &
            'dtdt_cloud','dqdt_cloud','dxldt_cloud','dxidt_cloud')
        IF (n_end == n_tdiag_vars) EXIT
        IF (st1_in_st2_idx(TRIM(ichar),c_vars(1:n_end))==0) THEN
          n_end=n_end+1
          c_vars(n_end)=TRIM(ichar)
        END IF
      CASE ('end')
        EXIT
      CASE DEFAULT
        CALL finish ('mo_diag_tendency_new: init_tdiag', 'key word/variable '// &
                     TRIM(ichar)//' not allowed in tdiag list' )
      END SELECT
    END DO
    tdiagnam(1:n_end)=c_vars(1:n_end)
    tdiagnam(n_end+1:n_tdiag_vars)=''
    !> Open new stream
    CALL new_stream (tdiag,'tdiag',lrerun=.false.,interval=puttdiag)
    !> Add standard fields for post-processing:
    CALL add_stream_reference (tdiag, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (tdiag, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (tdiag, 'st'      ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (tdiag, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
    CALL add_stream_reference (tdiag, 'gboxarea','geoloc',lpost=.TRUE.)
    CALL set_stream_element_info(tdiag,'gboxarea',table=199,code=254)
    !> Add a 3d gp temperature 
    CALL add_stream_reference (tdiag, 'tm1'     ,'g1a'   ,lpost=.TRUE.)
    CALL set_stream_element_info(tdiag,'tm1','3d temperature at previous time step', &
                                 units='K',table=199,code=255)
    !> Set default stream properties
    CALL default_stream_setting (tdiag,lrerun=.false.,contnorest=.true., &
         laccu=.false.,lpost=.true.,table=199, code=AUTO)
    !> Add variables to tdiag stream if they are present in tdiagnam list
    !> dtudt, dvdt, dtdt for Hines and SSO drag
    lpost = st1_in_st2_proof( 'dudt_hines', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                          ,'dudt_hines'    &
                              ,tdiag_vars%dudt_hines          ,units='m/s/day' &
                              ,longname='u-wind tend. due to Hines g.w.'       &
                              ,code=13                                         )
    END IF
    lpost = st1_in_st2_proof( 'dvdt_hines', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                          ,'dvdt_hines'    &
                              ,tdiag_vars%dvdt_hines          ,units='m/s/day' &
                              ,longname='v-wind tend. due to Hines g.w.'       &
                              ,code=23                                         )
    END IF
    lpost = st1_in_st2_proof( 'dtdt_hines', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                          ,'dtdt_hines'  &
                              ,tdiag_vars%dtdt_hines          ,units='K/day' &
                              ,longname='temperature tend. due to Hines g.w.'&
                              ,code=3                                        )
    END IF
    lpost = st1_in_st2_proof( 'dudt_sso', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                          ,'dudt_sso'      &
                              ,tdiag_vars%dudt_sso            ,units='m/s/day' &
                              ,longname='u-wind tend. due to oro. waves'       &
                              ,code=14                                         )
    END IF
    lpost = st1_in_st2_proof( 'dvdt_sso', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                          ,'dvdt_sso'      &
                              ,tdiag_vars%dvdt_sso            ,units='m/s/day' &
                              ,longname='v-wind tend. due to oro. waves'       &
                              ,code=24                                         )
    END IF
    lpost = st1_in_st2_proof( 'dtdt_sso', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                          ,'dtdt_sso'    &
                              ,tdiag_vars%dtdt_sso            ,units='K/day' &
                              ,longname='temperature tend. due to oro. waves'&
                              ,code=4                                        )
    END IF
    !> solar radiation heating rate dtdt_rheat_sw
    !> thermal radiation dtdt_rheat_lw
    lpost = st1_in_st2_proof( 'dtdt_rheat_sw', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                      ,'dtdt_rheat_sw'  &
                              ,tdiag_vars%dtdt_rheat_sw   ,units='K/day'    &
                              ,longname='temperature tend. due to solar '   &
                              //'wavelength rad heating'                    &
                              ,code=62                                      )
    END IF
    lpost = st1_in_st2_proof( 'dtdt_rheat_lw', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                      ,'dtdt_rheat_lw' &
                              ,tdiag_vars%dtdt_rheat_lw   ,units='K/day'   &
                              ,longname='temperature tend. due to thermal '&
                              //'wavelength rad heating'                    &
                              ,code=72                                     )
    END IF
    !> tendency in u, v, T due to vdiff
    lpost = st1_in_st2_proof( 'dudt_vdiff', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                          ,'dudt_vdiff'    &
                              ,tdiag_vars%dudt_vdiff          ,units='m/s/day' &
                              ,longname='u-wind tend. due to vert. diff.'      &
                              ,code=11                                         )
    END IF
    lpost = st1_in_st2_proof( 'dvdt_vdiff', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                          ,'dvdt_vdiff'    &
                              ,tdiag_vars%dvdt_vdiff          ,units='m/s/day' &
                              ,code=21                                         )
    END IF
    lpost = st1_in_st2_proof( 'dtdt_vdiff', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                          ,'dtdt_vdiff'  &
                              ,tdiag_vars%dtdt_vdiff          ,units='K/day' &
                              ,longname='v-wind tend. due to vert. diff.'    &
                              ,code=1                                        )
    END IF
    lpost = st1_in_st2_proof( 'dqdt_vdiff', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                          ,'dqdt_vdiff'  &
                              ,tdiag_vars%dqdt_vdiff          ,units='1/day' &
                              ,longname='spec. hum. tend. due to vert. diff.'&
                              ,code=31                                       )
    END IF
    lpost = st1_in_st2_proof( 'dxldt_vdiff', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                          ,'dxldt_vdiff'  &
                              ,tdiag_vars%dxldt_vdiff         ,units='1/day'  &
                              ,longname='cloud water tend. due to vert. diff.'&
                              ,code=41                                        )
    END IF
    lpost = st1_in_st2_proof( 'dxidt_vdiff', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                          ,'dxidt_vdiff' &
                              ,tdiag_vars%dxidt_vdiff         ,units='1/day' &
                              ,longname='cloud ice tend. due to vert. diff.' &
                              ,code=51                                       )
    END IF
    lpost = st1_in_st2_proof( 'dsddt_hdiff', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                   ,'dsddt_hdiff'        &
                              ,tdiag_vars%dsddt_hdiff                        &
                              ,(/lc%nlev,2,lc%snsp/)                         &
                              ,(/lc%nlev,2,nsp/)                             & 
                              ,units='1/s/day'                               &
                              ,longname='sp. div tend. due to horiz. diff.'  &
                              ,code=97                        ,repr=SPECTRAL )
    END IF
    lpost = st1_in_st2_proof( 'dsvodt_hdiff', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                          ,'dsvodt_hdiff'&
                              ,tdiag_vars%dsvodt_hdiff                       &
                              ,(/lc%nlev,2,lc%snsp/)                         &
                              ,(/lc%nlev,2,nsp/)                             & 
                              ,units='1/s/day'                               &
                              ,longname='sp. vor tend. due to horiz. diff.'  &
                              ,code=87                        ,repr=SPECTRAL )
    END IF
    lpost = st1_in_st2_proof( 'dstdt_hdiff', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                          ,'dstdt_hdiff'  &
                              ,tdiag_vars%dstdt_hdiff                         &
                              ,(/lc%nlev,2,lc%snsp/)                          &
                              ,(/lc%nlev,2,nsp/)                              & 
                              ,units='K/day'                                  &
                              ,longname='sp. temp tend. due to horiz. diff.'  &
                              ,code=7                         ,repr=SPECTRAL  )
    END IF
    !> tendency in u, v, T due to cucall (cumulus clouds)
    lpost = st1_in_st2_proof( 'dudt_cucall', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                          ,'dudt_cucall'   &
                              ,tdiag_vars%dudt_cucall         ,units='m/s/day' &
                              ,longname='u-wind tend. due to convection'       &
                              ,code=15                                         )
    END IF
    lpost = st1_in_st2_proof( 'dvdt_cucall', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                          ,'dvdt_cucall'   &
                              ,tdiag_vars%dvdt_cucall         ,units='m/s/day' &
                              ,longname='v-wind tend. due to convection'       &
                              ,code=25                                         )
    END IF
    lpost = st1_in_st2_proof( 'dtdt_cucall', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                          ,'dtdt_cucall' &
                              ,tdiag_vars%dtdt_cucall         ,units='K/day' &
                              ,longname='temperature tend. due to convection'&
                              ,code=5                                        )
    END IF
    lpost = st1_in_st2_proof( 'dqdt_cucall', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                          ,'dqdt_cucall' &
                              ,tdiag_vars%dqdt_cucall         ,units='1/day' &
                              ,code=35                                       )
    END IF
    !> tendency in T due to cloud
    lpost = st1_in_st2_proof( 'dtdt_cloud', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                          ,'dtdt_cloud'  &
                              ,tdiag_vars%dtdt_cloud          ,units='K/day' &
                              ,longname='temp. tend. due to clouds'          &
                              ,code=6                                        )
    END IF
    lpost = st1_in_st2_proof( 'dqdt_cloud', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                          ,'dqdt_cloud'  &
                              ,tdiag_vars%dqdt_cloud          ,units='K/day' &
                              ,longname='spec. hum. tend. due to clouds'     &
                              ,code=36                                       )
    END IF
    lpost = st1_in_st2_proof( 'dxldt_cloud', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                          ,'dxldt_cloud' &
                              ,tdiag_vars%dxldt_cloud         ,units='K/day' &
                              ,longname='cloud water tend. due to clouds'    &
                              ,code=46                                       )
    END IF
    lpost = st1_in_st2_proof( 'dxidt_cloud', tdiagnam)
    IF (lpost) THEN
      CALL add_stream_element (tdiag                          ,'dxidt_cloud' &
                              ,tdiag_vars%dxidt_cloud         ,units='K/day' &
                              ,longname='cloud ice tend. due to clouds'      &
                              ,code=56                                       )
    END IF

  END SUBROUTINE init_tdiag
  SUBROUTINE set_tendency_gp2d(var_diag         ,var_tendency     ,kproma     &
                               ,kbdim           ,klev             ,mode       )
    INTEGER, INTENT(in)           :: kproma, kbdim, klev
    REAL(wp), INTENT(inout)       :: var_diag(kbdim,klev)     !> diagnosed tendency
    REAL(wp), INTENT(in)          :: var_tendency(kbdim,klev) !> tendency variable
    CHARACTER(LEN=3), INTENT(in)  :: mode
    SELECT CASE (mode)
    CASE ('set')
      var_diag(1:kproma,1:klev)=daylength * var_tendency(1:kproma,1:klev)
    CASE ('add')
      var_diag(1:kproma,1:klev)=var_diag(1:kproma,1:klev) + &
              daylength * var_tendency(1:kproma,1:klev)
    CASE ('sub')
      var_diag(1:kproma,1:klev)= - daylength * var_tendency(1:kproma,1:klev)
    CASE DEFAULT
      CALL finish('mo_diag_tendency_new, set_tendency_gp2d', &
                 'wrong mode, mode='//mode)
    END SELECT
  END SUBROUTINE set_tendency_gp2d
  SUBROUTINE set_tendency_sp3d(var_diag         ,var_tendency      &
                                   ,mode       )
    REAL(wp), INTENT(inout)       :: var_diag(lc%nlev,2,lc%snsp)     !> diagnosed tendency
    REAL(wp), INTENT(in)          :: var_tendency(lc%nlev,2,lc%snsp) !> tendency variable
    CHARACTER(LEN=3), INTENT(in)  :: mode
    SELECT CASE (mode)
    CASE ('add')
      var_diag=var_diag + &
              daylength * var_tendency
    CASE ('sub')
      var_diag= - daylength * var_tendency
    CASE DEFAULT
      CALL finish('mo_diag_tendency_new, set_tendency_sp3d', &
                 'wrong mode, mode='//mode)
    END SELECT
  END SUBROUTINE set_tendency_sp3d
END MODULE mo_diag_tendency_new
