!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_diag_amip2

  ! Calculates daily global averages (area-weighted) for AMIP2 table 4.
  ! Based on daystatz and daystat in ECHAM4f77.
  ! 
  ! L. Kornblueh, MPI October 2003, initial version
  ! L. Kornblueh, MPI, January 2004, corrections and clean runtime behavior

  USE mo_kind, ONLY: wp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_amip2_diag
  PUBLIC :: amip2_global_diag
  PUBLIC :: collect_amip2_diag

  ! Calculate/store daily  statistics

  ! Block AMIP2 statistics

  REAL(wp), ALLOCATABLE :: akez (:)    ! kinetic energy    (enek)       
  REAL(wp), ALLOCATABLE :: psuz (:)    ! surface pressure  (ps)         
  REAL(wp), ALLOCATABLE :: tmpz (:)    ! temperature       (ta)         
  REAL(wp), ALLOCATABLE :: amz  (:)    ! angular momentum  (moa)        
  REAL(wp), ALLOCATABLE :: tomz (:)    ! mountain torque                
  REAL(wp), ALLOCATABLE :: tovdz(:)    ! torque from *vdiff*            
  REAL(wp), ALLOCATABLE :: togwz(:)    ! torque from *gwdrag*           
  REAL(wp), ALLOCATABLE :: sstz (:)    ! sst over open sea              
  REAL(wp), ALLOCATABLE :: swgtz(:)

  ! temporary arrays to keep for tomz calculation with cyclic 
  ! boundary conditions, simplifies handling

  REAL(wp), ALLOCATABLE :: daphm1(:,:)
  REAL(wp), ALLOCATABLE :: dgeospm(:,:)

  ! Global AMIP2 statistics

  REAL(wp) :: ekac                    ! kinetic energy    (enek)       
  REAL(wp) :: psac                    ! surface pressure  (ps)         
  REAL(wp) :: tmac                    ! temperature       (ta)         
  REAL(wp) :: amac                    ! angular momentum  (moa)        
  REAL(wp) :: tomac                   ! mountain torque                
  REAL(wp) :: tovdac                  ! torque from *vdiff*            
  REAL(wp) :: togwac                  ! torque from *gwdrag*           
  REAL(wp) :: sstac                   ! sst over open sea              

  INTEGER :: amip2_output_unit

  CHARACTER(len=*), PARAMETER :: amip2_filename = 'AMIP2_diagnostics_table4' 

CONTAINS

  SUBROUTINE init_amip2_diag

    USE mo_control,       ONLY: ngl
    USE mo_decomposition, ONLY: ldc => local_decomposition

    LOGICAL, SAVE :: not_used = .TRUE.

    IF (not_used) THEN

      ALLOCATE (akez (ldc%ngpblks))
      ALLOCATE (psuz (ldc%ngpblks))
      ALLOCATE (tmpz (ldc%ngpblks))
      ALLOCATE (amz  (ldc%ngpblks))
      ALLOCATE (tomz (ngl))
      ALLOCATE (tovdz(ldc%ngpblks))
      ALLOCATE (togwz(ldc%ngpblks))
      ALLOCATE (sstz (ldc%ngpblks))
      ALLOCATE (swgtz(ldc%ngpblks))

      akez (:)  = 0.0_wp
      psuz (:)  = 0.0_wp
      tmpz (:)  = 0.0_wp
      amz (:)   = 0.0_wp
      tomz (:)  = 0.0_wp
      tovdz (:) = 0.0_wp
      togwz (:) = 0.0_wp
      sstz (:)  = 0.0_wp
      swgtz (:) = 0.0_wp

      ALLOCATE (daphm1(ldc%nproma,ldc%ngpblks))      
      ALLOCATE (dgeospm(ldc%nproma,ldc%ngpblks))      

      not_used = .FALSE.

    END IF

  END SUBROUTINE init_amip2_diag

  SUBROUTINE amip2_global_diag

    ! Calculates daily means of global averages of:
    !
    ! surface pressure
    ! kinetic energy     (per unit area)
    ! temperature
    ! angular momentum   (per unit area)
    ! torques            (per unit area)
    !
    ! daily_global_diag is called from scan1
    !
    ! U. Schlese, DKRZ, October 1998, initial version
    ! I. Kirchner, MPIM, November 2000, date/time control
    ! L. Kornblueh, MPIM, October 2003, packed in module and parallelized

    USE mo_io_units,        ONLY: nout
    USE mo_control,       ONLY: nlon, ngl
    USE mo_math_constants,ONLY: pi
    USE mo_physical_constants, ONLY: grav, earth_radius
    USE mo_gaussgrid,     ONLY: gl_budw
    USE mo_filename,      ONLY: find_next_free_unit
    USE mo_time_control,  ONLY: lstart, lresume, next_date, l_diagamip,     &
                                delta_time, l_putrerun, get_date_components, &
                                ndaylen
    USE mo_mpi,           ONLY: p_global_sum, p_parallel_io
    USE mo_transpose,     ONLY: gather_gp
    USE mo_decomposition, ONLY: gd => global_decomposition

    INTEGER  :: iy, id, im, hr, mn, se

    REAL(wp) :: zgps, zgtovd, zgtogw, zgtom, zgam, zgke, zgtmp, zgsst, zgswgt
    INTEGER  :: ioday, iomonth, ioyear
    REAL(wp) :: tosac, toac, amgac, dayl

    REAL(wp) :: zdtime, zg

    REAL(wp) :: zaphm1(nlon+1), zgeospm(nlon+1)

    REAL(wp), TARGET  :: gl_aphm1(nlon,ngl), gl_geospm(nlon,ngl)
    REAL(wp), POINTER :: gl(:,:)

    LOGICAL :: lexist

    INTEGER :: jg, jl

    zg = 1.0_wp/grav

    gl => gl_aphm1
    CALL gather_gp (gl, daphm1, gd)

    gl => gl_geospm
    CALL gather_gp (gl, dgeospm, gd)

    IF (p_parallel_io) THEN

      tomz(:) = 0.0_wp

      DO jg = 1, ngl
        zaphm1(1:nlon)  = gl_aphm1(1:nlon,jg)        
        zaphm1(nlon+1)  = zaphm1(1)
        zgeospm(1:nlon) = gl_geospm(1:nlon,jg)
        zgeospm(nlon+1) = zgeospm(1)

        DO jl = 1, nlon
          tomz(jg) = tomz(jg)-(zaphm1(jl)+zaphm1(jl+1))*0.5_wp &
               *(zgeospm(jl+1)-zgeospm(jl))
        END DO
      END DO

      tomz(:) = tomz(:)*gl_budw(:)*zg

      zgtom  = SUM(tomz)

    END IF

    ! 1. global summation

    zgps   = p_global_sum(psuz)
    zgtovd = p_global_sum(tovdz)
    zgtogw = p_global_sum(togwz)
    zgam   = p_global_sum(amz)
    zgke   = p_global_sum(akez)
    zgtmp  = p_global_sum(tmpz)
    zgsst  = p_global_sum(sstz)
    zgswgt = p_global_sum(swgtz)

    IF (p_parallel_io) THEN

      ! 1.1  blank block values for later use

      akez (:)  = 0.0_wp
      psuz (:)  = 0.0_wp
      tmpz (:)  = 0.0_wp
      amz (:)   = 0.0_wp
      tomz (:)  = 0.0_wp
      tovdz (:) = 0.0_wp
      togwz (:) = 0.0_wp
      sstz (:)  = 0.0_wp
      swgtz (:) = 0.0_wp

      zgtmp = zgtmp/zgps
      zgsst = zgsst/zgswgt

      ! 2. daily averages

      zdtime = delta_time

      IF (lstart) THEN

        ! 2.1  blank values at first timestep

        ekac   = 0.0_wp
        psac   = 0.0_wp
        tmac   = 0.0_wp
        amac   = 0.0_wp
        tomac  = 0.0_wp
        tovdac = 0.0_wp
        togwac = 0.0_wp
        sstac  = 0.0_wp

        amip2_output_unit = find_next_free_unit(77,99)
        OPEN (unit=amip2_output_unit, file=TRIM(amip2_filename))

        WRITE (nout,'(a)') 'AMIP2 diagnostics output file opened.'

      ELSE IF (lresume) THEN

        ! 2.2  read statistics from rerun file at beginning of run

        CALL rerun_read(ioday, iomonth, ioyear)

        WRITE (nout,'(a,a,i2.2,a,i2.2,a,i4.4)') &
             ' Global statistics rerun read', &
             ': day= ',ioday,' month= ',iomonth,' year= ',ioyear

        amip2_output_unit = find_next_free_unit(77,99)
        INQUIRE (file=TRIM(amip2_filename), exist=lexist)
        IF (lexist) THEN
          OPEN (unit=amip2_output_unit, file=TRIM(amip2_filename), &
               status='OLD', position='APPEND')
        ELSE
          OPEN (unit=amip2_output_unit, file=TRIM(amip2_filename))
        ENDIF

        WRITE (nout,'(a)') 'AMIP2 diagnostics output file opened.'

      ENDIF
#ifdef DEBUG
      ! 2.2.1 write instantanous values of statistics for AMIP2 output

      WRITE (nout,'(/,a)') 'AMIP2 global statistics:'
      WRITE (nout,'(a,e12.5)') 'Total kinetic energy   (enek) : ', zgke 
      WRITE (nout,'(a,e12.5)') 'Total angular momentum (moa)  : ', zgam
      WRITE (nout,'(a,e12.5)') 'Total surface torque   (torts): ', &
           zgtom+zgtovd+zgtogw
      WRITE (nout,'(a,f7.2)')  'Temperature            (ta)   : ', zgtmp
      WRITE (nout,'(a,f10.2)') 'Surface pressure       (ps)   : ', zgps
      WRITE (nout,'(a,f7.2)')  'SST over open ocean    (tso)  : ', zgsst
      WRITE (nout,'(/)')
#endif
      ! 2.3  accumulate statistics for AMIP2 output

      ekac   = ekac+zgke*zdtime          ! kinetic energy    (enek)
      psac   = psac+zgps*zdtime          ! surface pressure  (ps)
      tmac   = tmac+zgtmp*zdtime         ! temperature       (ta)
      amac   = amac+zgam*zdtime          ! angular momentum  (moa)
      tomac  = tomac+zgtom*zdtime        ! mountain torque
      tovdac = tovdac+zgtovd*zdtime      ! torque from *vdiff*
      togwac = togwac+zgtogw*zdtime      ! torque from *gwdrag*
      sstac  = sstac+zgsst*zdtime        ! sst over open sea

      ! 2.4 write to file at end of day and blank values
      !

      IF(l_diagamip) THEN

        CALL get_date_components(next_date,iy,im,id,hr,mn,se)

        dayl = REAL(ndaylen,wp)

        ekac   = ekac/dayl
        psac   = psac/dayl
        tmac   = tmac/dayl
        amac   = amac/dayl
        tomac  = tomac/dayl
        tovdac = tovdac/dayl
        togwac = togwac/dayl
        sstac  = sstac/dayl

        tosac = tovdac+togwac
        toac  = tomac+tosac
        amgac = amac*4.0_wp*pi*earth_radius**2

        WRITE(amip2_output_unit,'(2i3.2,i5.4,11e13.5)') &
             id, im, iy, ekac, psac, tmac, &
             amgac, amac, toac, tomac, tosac, tovdac, togwac, sstac

        WRITE (nout,'(a,i4,a,i2.2,a,i2.2,a)') &
             'AMIP2 diagnostics output written (',iy,'-',im,'-',id,').'

        ekac   = 0.0_wp
        psac   = 0.0_wp
        tmac   = 0.0_wp
        amac   = 0.0_wp
        tomac  = 0.0_wp
        tovdac = 0.0_wp
        togwac = 0.0_wp
        sstac  = 0.0_wp

      ENDIF

      ! 3. save statistics at end of run

      IF (l_putrerun) THEN

        CALL get_date_components(next_date,iy,im,id,hr,mn,se)

        CALL rerun_write(id, im, iy)

        WRITE (nout,'(a,a,i2.2,a,i2.2,a,i4.4)') &
             ' Global statistics rerun write', &
             ': day= ',id,' month= ',im,' year= ',iy

        CLOSE (amip2_output_unit)

      ENDIF

    ENDIF

  END SUBROUTINE amip2_global_diag

  SUBROUTINE collect_amip2_diag (kproma, kbdim, klev, klevp1, krow, & 
       ptm1,  pum1,    pvm1,   paphm1,  pgeospm,   &
       pustr, pustrgw, pustrm, pustrgwm,           &
       ptsm1m, pseaicem, ldland)

    ! Calculates daily means of block averages of:
    !
    ! surface pressure
    ! kinetic energy      (per unit area)
    ! temperature
    ! angular momentum    (per unit area)
    ! torques             (per unit area)
    ! sst over open ocean (excluding ice covered areas)
    !
    ! collect_amip2_diag is called from physc
    !
    ! U. Schlese, DKRZ, October 1998, initial version
    ! L. Kornblueh, MPIM, October 2003, packed in module and parallelized

    USE mo_control,      ONLY: nlon
    USE mo_geoloc,       ONLY: sqcst_2d, budw_2d
    USE mo_physical_constants,    ONLY: grav, earth_radius
    USE mo_time_control, ONLY: delta_time

    INTEGER, INTENT(in) :: kproma, kbdim, klev, klevp1, krow

    REAL(wp), INTENT(in) :: ptm1(kbdim,klev),      &
         pum1(kbdim,klev),     pvm1(kbdim,klev), &
         paphm1(kbdim,klevp1), pgeospm(kbdim),   &
         pustr(kbdim),         pustrgw(kbdim),   &
         pustrm(kbdim),        pustrgwm(kbdim),  &
         ptsm1m(kbdim),        pseaicem(kbdim)
    LOGICAL, INTENT(in) :: ldland(kbdim)

    REAL(wp) :: zdelp(kproma)

    REAL(wp) :: zdtime, zg

    INTEGER :: jk

    ! set up weighting factors

    zdtime = 1.0_wp/delta_time
    zg = 1.0_wp/grav

    DO jk = 1, klev

      zdelp(1:kproma)   = paphm1(1:kproma,jk+1)-paphm1(1:kproma,jk)      

      amz(krow)  = amz(krow)+SUM(zdelp(1:kproma)*pum1(1:kproma,jk) &
           *sqcst_2d(1:kproma,krow)*budw_2d(1:kproma,krow))*earth_radius*zg
      akez(krow) = akez(krow)+SUM(zdelp(1:kproma) &
           *(pum1(1:kproma,jk)**2+pvm1(1:kproma,jk)**2) &
           *budw_2d(1:kproma,krow))*0.5_wp*zg
      tmpz(krow) = tmpz(krow)+SUM(zdelp(1:kproma)*ptm1(1:kproma,jk) &
           *budw_2d(1:kproma,krow))

    END DO

    psuz(krow)  = SUM(paphm1(1:kproma,klevp1)*budw_2d(1:kproma,krow))
    tovdz(krow) = SUM(-(pustr(1:kproma)-pustrm(1:kproma)) &
         *sqcst_2d(1:kproma,krow)*budw_2d(1:kproma,krow))*zdtime*earth_radius
    togwz(krow) = SUM(-(pustrgw(1:kproma)-pustrgwm(1:kproma)) &
         *sqcst_2d(1:kproma,krow)*budw_2d(1:kproma,krow))*zdtime*earth_radius

    sstz(krow)  = SUM(ptsm1m(1:kproma)*budw_2d(1:kproma,krow), &
         MASK=(.NOT. ldland(1:kproma) .AND. pseaicem(1:kproma) < 0.5_wp))*REAL(nlon,wp)
    swgtz(krow) = SUM(budw_2d(1:kproma,krow), &
         MASK=(.NOT. ldland(1:kproma) .AND. pseaicem(1:kproma) < 0.5_wp))*REAL(nlon,wp)

    daphm1(1:kproma,krow) = paphm1(1:kproma,klevp1)
    dgeospm(1:kproma,krow) = pgeospm(1:kproma)

  END SUBROUTINE collect_amip2_diag

  SUBROUTINE rerun_write(iday, imon, iyear)

    USE mo_netcdf,        ONLY: FILE_INFO
    USE mo_io
    USE mo_filename,      ONLY: standard_rerun_file

    INTEGER, INTENT(in) :: iday, imon, iyear
    INTEGER             :: dims(1), fileID, date
    INTEGER             :: ekacID, psacID, tmacID, amacID, tomacID, tovdacID, togwacID, sstacID
    TYPE (FILE_INFO)  :: fileinfo

    date = iyear*10000 + imon*100 + iday
    fileinfo%opened    = .FALSE.

    CALL IO_OPEN(TRIM(standard_rerun_file)//'_diag_amip2', fileinfo, IO_WRITE)

    fileID = fileinfo%file_id

    CALL IO_PUT_ATT_INT (fileID, NF_GLOBAL, 'date', date)

    CALL IO_DEF_VAR(fileID, 'enek', NF_DOUBLE, 0, dims, ekacID)
    CALL IO_PUT_ATT_TEXT(fileID, ekacID, 'long_name', 'total kinetic energy')
    CALL IO_PUT_ATT_TEXT(fileID, ekacID, 'units', 'J/m**2')

    CALL IO_DEF_VAR(fileID, 'ps', NF_DOUBLE, 0, dims, psacID)
    CALL IO_PUT_ATT_TEXT(fileID, psacID, 'long_name', 'surface pressure')
    CALL IO_PUT_ATT_TEXT(fileID, psacID, 'units', 'Pa')

    CALL IO_DEF_VAR(fileID, 'ta', NF_DOUBLE, 0, dims, tmacID)
    CALL IO_PUT_ATT_TEXT(fileID, tmacID, 'long_name', 'temperature')
    CALL IO_PUT_ATT_TEXT(fileID, tmacID, 'units', 'K')

    CALL IO_DEF_VAR(fileID, 'moa', NF_DOUBLE, 0, dims, amacID)
    CALL IO_PUT_ATT_TEXT(fileID, amacID, 'long_name', 'total relative angular momentum')
    CALL IO_PUT_ATT_TEXT(fileID, amacID, 'units', 'kg/s')

    CALL IO_DEF_VAR(fileID, 'torts', NF_DOUBLE, 0, dims, tomacID)
    CALL IO_PUT_ATT_TEXT(fileID, tomacID, 'long_name', 'total surface torque')
    CALL IO_PUT_ATT_TEXT(fileID, tomacID, 'units', 'N/m')

    CALL IO_DEF_VAR(fileID, 'tovdac', NF_DOUBLE, 0, dims, tovdacID)
    CALL IO_PUT_ATT_TEXT(fileID, tovdacID, 'long_name', 'torque from *vdiff*')
    CALL IO_PUT_ATT_TEXT(fileID, tovdacID, 'units', 'N/m')

    CALL IO_DEF_VAR(fileID, 'togwac', NF_DOUBLE, 0, dims, togwacID)
    CALL IO_PUT_ATT_TEXT(fileID, togwacID, 'long_name', 'torque from *gwdrag*')
    CALL IO_PUT_ATT_TEXT(fileID, togwacID, 'units', 'N/m')

    CALL IO_DEF_VAR(fileID, 'tso', NF_DOUBLE, 0, dims, sstacID)
    CALL IO_PUT_ATT_TEXT(fileID, sstacID, 'long_name', 'SST over open sea')
    CALL IO_PUT_ATT_TEXT(fileID, sstacID, 'units', 'K')

    CALL IO_ENDDEF(fileID)

    CALL IO_PUT_VAR_DOUBLE1(fileID, ekacID, ekac)
    CALL IO_PUT_VAR_DOUBLE1(fileID, psacID, psac)
    CALL IO_PUT_VAR_DOUBLE1(fileID, tmacID, tmac)
    CALL IO_PUT_VAR_DOUBLE1(fileID, amacID, amac)
    CALL IO_PUT_VAR_DOUBLE1(fileID, tomacID, tomac)
    CALL IO_PUT_VAR_DOUBLE1(fileID, tovdacID, tovdac)
    CALL IO_PUT_VAR_DOUBLE1(fileID, togwacID, togwac)
    CALL IO_PUT_VAR_DOUBLE1(fileID, sstacID, sstac)

    CALL IO_CLOSE(fileinfo)

  END SUBROUTINE rerun_write

  SUBROUTINE rerun_read(iday, imon, iyear)

    USE mo_netcdf,        ONLY: FILE_INFO
    USE mo_io
    USE mo_filename,      ONLY: standard_rerun_file

    INTEGER, INTENT(out):: iday, imon, iyear
    INTEGER             :: fileID, date
    INTEGER             :: ekacID, psacID, tmacID, amacID, tomacID, tovdacID, togwacID, sstacID
    TYPE (FILE_INFO)  :: fileinfo

    fileinfo%opened = .FALSE.

    CALL IO_OPEN(TRIM(standard_rerun_file)//'_diag_amip2', fileinfo, IO_READ)

    fileID = fileinfo%file_id

    CALL IO_GET_ATT_INT(fileID, NF_GLOBAL, 'date', date)

    iyear = date / 10000
    imon = (date - iyear*10000) / 100
    iday =  date - iyear*10000 - imon*100

    CALL IO_INQ_VARID(fileID, 'enek', ekacID)
    CALL IO_INQ_VARID(fileID, 'ps', psacID)
    CALL IO_INQ_VARID(fileID, 'ta', tmacID)
    CALL IO_INQ_VARID(fileID, 'moa', amacID)
    CALL IO_INQ_VARID(fileID, 'torts', tomacID)
    CALL IO_INQ_VARID(fileID, 'tovdac', tovdacID)
    CALL IO_INQ_VARID(fileID, 'togwac', togwacID)
    CALL IO_INQ_VARID(fileID, 'tso', sstacID)

    CALL IO_GET_VAR_DOUBLE1(fileID, ekacID, ekac)
    CALL IO_GET_VAR_DOUBLE1(fileID, psacID, psac)
    CALL IO_GET_VAR_DOUBLE1(fileID, tmacID, tmac)
    CALL IO_GET_VAR_DOUBLE1(fileID, amacID, amac)
    CALL IO_GET_VAR_DOUBLE1(fileID, tomacID, tomac)
    CALL IO_GET_VAR_DOUBLE1(fileID, tovdacID, tovdac)
    CALL IO_GET_VAR_DOUBLE1(fileID, togwacID, togwac)
    CALL IO_GET_VAR_DOUBLE1(fileID, sstacID, sstac)

    CALL IO_CLOSE(fileinfo)

  END SUBROUTINE rerun_read

END MODULE mo_diag_amip2
