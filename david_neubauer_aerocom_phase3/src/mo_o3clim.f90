!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_o3clim

  !- Description:
  !
  !  This module provides a 3 dimensional monthly mean ozone climatology.
  !
  !  This module contains:
  !
  !  A) internal variables
  !  B) the subroutine su_o3clim to read the ozone file and
  !     initialize the ozone climatology
  !  C) the function o3clim to get a longitude height ozone section
  !  D) the subroutine cleanup_o3clim to deallocate memory
  !
  !- Author:
  !
  !  M.A. Giorgetta, MPI, May  2000
  !  U. Schulzweida, MPI, May  2002, blocking (nproma)
  !  A. Rhodin,      DWD, June 2002, new subroutine cleanup_o3clim
  !  U. Schulzweida, MPI, Jan  2004, read 3D monthly mean ozone climatology
  !  U. Schulzweida, MPI, May  2004, read annual cycle  monthly mean ozone
  !
  !=======================================================================

  USE mo_kind,       ONLY: dp
  USE mo_exception,  ONLY: finish, message, message_text, em_warn
  USE mo_interpo,    ONLY: wgt1_m,wgt2_m,nmw1_m,nmw2_m,wgt1,wgt2,nmw1cl,nmw2cl

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_o3clim_3, read_o3clim_4 
  PUBLIC :: pre_o3clim_3, pre_o3clim_4
  PUBLIC :: o3clim
  PUBLIC :: zozonec_x, cleanup_o3clim

  !=======================================================================
  !
  ! A)
  !
  ! unit number of ozone file
  INTEGER, PARAMETER :: nioz=21

  ! pressure level grid of ozone climatology in Pa
  INTEGER            :: noz          ! number of ozone full levels
  REAL(dp), ALLOCATABLE  :: pfoz(:)      ! full levels 
  REAL(dp), ALLOCATABLE  :: phoz(:)      ! half levels

  ! ozone climatology in mass mixing ratio in layers, latitudes and months
  REAL(dp), ALLOCATABLE  :: ozone3(:,:,:)  ! mass mixing ratio at full levels
  REAL(dp), ALLOCATABLE  :: ozone(:,:,:,:) ! mass mixing ratio at full levels
  REAL(dp), ALLOCATABLE, TARGET  :: zozonec_x(:,:,:)
  LOGICAL                :: lo3_set=.false.

  !=======================================================================

CONTAINS

  !=======================================================================
  !
  ! B)

  SUBROUTINE read_o3clim_3

    !- Description:
    !
    !  read_o3clim_3 provides a zonal mean and monthly mean
    !  climatology of ozone mass mixing ratio (g/g).
    !
    !  The NetCDF file accessed by unit nioz must contain ozone
    !  volume mixing ratios in ppmv at pressure levels.
    !
    !  The pressure levels are ordered from top to surface.
    !  At each level the data are given for ngl latitudes
    !  and all months.
    !
    !  The arrays are extended in time such that:
    !  month 0 corresponds to month 12 and month 13 to month 1
    !  in order to facilitate the interpolation in time.
    !
    !  After reading the data file the ozone concentration is
    !  converted to mass mixing ratio in g/g. 
    !
    !  Reading the NetCDF file is based on readozone in mo_midatm.f90
    !  (version ECHAM5.0.03) of U. Schulzweida
    !
    !- Author:
    !
    !  M.A. Giorgetta, MPI, May 2000
    !  L. Kornblueh, MPI, February 2011, fixes some bugs with respect to 
    !                       restart cycles an orthogonality in SW design
    !

    USE mo_control,       ONLY: ngl
    USE mo_mpi,           ONLY: p_pe, p_io, p_bcast
    USE mo_io,            ONLY: io_open_unit, io_close, io_read, &
         &                      io_var_id, io_file_id, ini_ozon
    USE mo_netCDF,        ONLY: io_inq_dimid, io_inq_dimlen,     &
         &                      io_inq_varid, io_get_var_double, &
         &                      io_get_vara_double
    USE mo_physical_constants, ONLY: amo3, amd
    USE mo_filename,      ONLY: NETCDF          

    ! ppmv2gg converts ozone from volume mixing ratio in ppmv
    ! to mass mixing ratio in g/g
    REAL(dp), PARAMETER :: ppmv2gg=1.e-6_dp*amo3/amd

    INTEGER               :: io_nlon  ! number of longitudes in NetCDF file
    INTEGER               :: io_ngl   ! number of latitudes in NetCDF file
    INTEGER, DIMENSION(4) :: io_start ! start index for NetCDF-read
    INTEGER, DIMENSION(4) :: io_count ! number of iterations for NetCDF-read

    INTEGER          :: jk, jkmax ! loop index

    CHARACTER(len=8) :: yformat

    ! Read ozone file
    ! ===============

    CALL message('','')
    WRITE(message_text,'(a,i2)') ' Read ozone climatology from unit ', nioz
    CALL message('',message_text)
    WRITE(message_text,*)' you chose zonal mean ozone boundary conditions.'// &
                           ' No transient boundary condition is possible'
    CALL message('read_o3clim_3',message_text,level=em_warn)

    IF (p_pe==p_io) THEN

       ini_ozon%format = NETCDF
       CALL io_open_unit (nioz, ini_ozon, io_read)
       io_file_id = ini_ozon%file_id

       ! Check resolution
       CALL io_inq_dimid  (io_file_id, 'lat', io_var_id)
       CALL io_inq_dimlen (io_file_id, io_var_id, io_ngl)
       CALL io_inq_dimid  (io_file_id, 'lon', io_var_id)
       CALL io_inq_dimlen (io_file_id, io_var_id, io_nlon)
       CALL io_inq_dimid  (io_file_id, 'level', io_var_id)
       CALL io_inq_dimlen (io_file_id, io_var_id, noz)

       IF (io_nlon /= 1.OR.io_ngl/=ngl) THEN
          WRITE(message_text,'(a,i0,i0)') &
               'read_o3clim_3: unexpected resolution ',io_nlon,io_ngl
          CALL message('',message_text)
          WRITE(message_text,'(a,i0)') &
               'expected number of latitudes = ',ngl
          CALL message('',message_text)
          WRITE(message_text,'(a,i0)') &
               'number of latitudes of ozone = ',io_ngl
          CALL message('',message_text)
          CALL finish ('read_o3clim_3','unexpected resolution')
       END IF

    END IF

    CALL p_bcast (noz, p_io)

    ! ALLOCATE memory
    IF ( .NOT. ALLOCATED(ozone3)) ALLOCATE(ozone3(noz,ngl,0:13))
    IF ( .NOT. ALLOCATED(pfoz))  ALLOCATE(pfoz(noz))
    IF ( .NOT. ALLOCATED(phoz))  ALLOCATE(phoz(noz+1))

    IF (p_pe == p_io) THEN

       ! read full pressure levels of ozone climatology
       CALL io_inq_varid      (io_file_id, 'level', io_var_id)
       CALL io_get_var_double (io_file_id, io_var_id, pfoz(:))

       WRITE(message_text,'(a,i0)') &
            'number of pressure levels of ozone climatology: ', noz
       CALL message('',message_text)
       CALL message('','pressure levels in [hPa]:')
       DO jk = 1, noz, 10
         jkmax = jk+10-1
         if (jkmax > noz) then
           jkmax = noz
         endif
         WRITE(yformat,'(a,i0,a)') '(', jkmax-jk+1, 'f8.2)'
         WRITE(message_text,TRIM(yformat)) pfoz(jk:jkmax)*0.01_dp
         CALL message('',message_text)
       ENDDO

       ! read ozone volume mixing ratio at full pressure levels in ppmv 
       CALL io_inq_varid (io_file_id, 'OZON', io_var_id)
       DO jk = 1, noz
          io_start(:) = (/ 1,   1, jk,  1 /)
          io_count(:) = (/ 1, ngl,  1, 12 /)
          ! for level jk: read ngl latitudes and 12 months
          CALL io_get_vara_double(io_file_id,io_var_id,io_start,io_count, &
               &                  ozone3(jk,:,1:12))
       END DO

       CALL io_close (ini_ozon)

    ENDIF

    CALL p_bcast (pfoz, p_io)
    CALL p_bcast (ozone3, p_io)


    ! Initialize arrays
    ! =================

    ! copy December to month 0
    ozone3(:,:,0)  = ozone3(:,:,12)

    ! copy January to month 13
    ozone3(:,:,13) = ozone3(:,:,1)

    ! convert from ppmv to g/g
    ozone3(:,:,:)  = ozone3(:,:,:)*ppmv2gg

    ! define half levels of ozone pressure grid
    ! upper boundary: ph =      0.Pa -> extrapolation of uppermost value
    ! lower boundary: ph = 125000.Pa -> extrapolation of lowermost value
    phoz(1)=0._dp
    phoz(2:noz)=(pfoz(1:noz-1)+pfoz(2:noz))/2._dp
    phoz(noz+1)=125000._dp

  END SUBROUTINE read_o3clim_3


  SUBROUTINE read_o3clim_4(rad_step_1, rad_step_2)

    !- Description:
    !
    !  read_o3clim_4 provides a monthly mean
    !  climatology of ozone mass mixing ratio (g/g).
    !
    !  The NetCDF file accessed by unit nioz must contain ozone
    !  volume mixing ratios in ppmv at pressure levels.
    !
    !  The pressure levels are ordered from top to surface.
    !  At each level the data are given for ngl latitudes
    !  and all months.
    !
    !  The arrays are extended in time such that:
    !  month 0 corresponds to month 12 and month 13 to month 1
    !  in order to facilitate the interpolation in time.
    !
    !  After reading the data file the ozone concentration is
    !  converted to mass mixing ratio in g/g. 
    !
    !  Reading the NetCDF file is based on readozone in mo_midatm.f90
    !  (version ECHAM5.0.03) of U. Schulzweida
    !
    !- Author:
    !
    !  M.A. Giorgetta, MPI, May 2000
    !  J.S. Rast, MPI, April 2010: modifications for new ozone clim. (cmip5)
    !  L. Kornblueh, MPI, February 2011, fixes some bugs with respect to 
    !                       restart cycles an orthogonality in SW design
    !  D. Zanchettin, MPI, July 2011, changes to read three-digits years
    !
    

    USE mo_control,         ONLY: ngl, nlon
    USE mo_mpi,             ONLY: p_pe, p_io, p_bcast
    USE mo_io,              ONLY: io_open, io_close, io_read, &
         &                        io_var_id
    USE mo_netCDF,          ONLY: io_inq_dimid, io_inq_dimlen,     &
         &                        io_inq_varid, io_get_var_double, &
         &                        io_get_vara_double, FILE_INFO
    USE mo_physical_constants, ONLY: amo3,amd
    USE mo_decomposition,   ONLY: lc => local_decomposition, &
                                  gl_dc => global_decomposition
    USE mo_transpose,       ONLY: scatter_gp
    USE mo_time_control,    ONLY: get_date_components
    USE mo_time_conversion, ONLY: time_days

    !INPUT PARMETES
    TYPE(time_days), INTENT(in), OPTIONAL :: rad_step_1, rad_step_2    

    !LOCAL VARIABLES
    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:,:,:)
    REAL(dp), POINTER :: gl_ozon(:,:,:,:)

    ! ppmv2gg converts ozone from volume mixing ratio in ppmv
    ! to mass mixing ratio in g/g
!    REAL(dp), PARAMETER :: ppmv2gg=1.e-6_dp*amo3/amd
    REAL(dp), PARAMETER :: vmr2mmr=amo3/amd

    INTEGER               :: io_nlon  ! number of longitudes in NetCDF file
    INTEGER               :: io_ngl   ! number of latitudes in NetCDF file
    INTEGER, DIMENSION(4) :: io_start ! start index for NetCDF-read
    INTEGER, DIMENSION(4) :: io_count ! number of iterations for NetCDF-read

    INTEGER               :: jk, jkmax ,i     ! loop index

    INTEGER :: IO_file_id0, IO_file_id1, IO_file_id2
    INTEGER :: icurrentyear, inextyear
    INTEGER       :: ihy0, ihy1, ihy2

    LOGICAL       :: lnewyear, lex0, lex1, lex2

    TYPE (FILE_INFO) :: fozon0, fozon1, fozon2

    CHARACTER(len=8) :: fn0, fn1, fn2
    CHARACTER(len=8) :: yformat

    ! Read ozone file
    ! ===============
    IF (PRESENT(rad_step_1) .AND. PRESENT(rad_step_2)) THEN
      CALL get_date_components (rad_step_1, year=icurrentyear)
      CALL get_date_components (rad_step_2, year=inextyear)
      lnewyear=icurrentyear/=inextyear
      
      IF (.NOT.lnewyear .AND. lo3_set) RETURN
    ENDIF

    IF (p_pe==p_io) THEN

      fozon0%opened = .FALSE.
      fozon1%opened = .FALSE.
      fozon2%opened = .FALSE.

      CALL set_years(ihy0, ihy1, ihy2)
      
      WRITE (fn0, '("ozon",i4.4)') ihy0
      WRITE (fn1, '("ozon",i4.4)') ihy1
      WRITE (fn2, '("ozon",i4.4)') ihy2

      CALL message('','')
      
      WRITE(message_text,*) 'fn0: ', TRIM(fn0),' fn1: ',TRIM(fn1), ' fn2: ',TRIM(fn2)
      CALL message('read_o3clim_4',message_text)

      INQUIRE (file=fn0, exist=lex0)
      INQUIRE (file=fn1, exist=lex1)
      INQUIRE (file=fn2, exist=lex2)

      IF ( .NOT. lex0 ) THEN
        WRITE (message_text,*) 'Could not open file <',fn0,'>'
        CALL message('',message_text)
        CALL finish ('read_o3clim_4', 'run terminated.')
      END IF

      IF ( .NOT. lex1 ) THEN
        WRITE (message_text,*) 'Could not open file <',fn1,'>'
        CALL message('',message_text)
        CALL finish ('read_o3clim_4', 'run terminated.')
      END IF

      IF ( .NOT. lex2 ) THEN
        WRITE (message_text,*) 'Could not open file <',fn2,'>'
        CALL message('',message_text)
        CALL finish ('read_o3clim_4', 'run terminated.')
      END IF
        
      CALL IO_open (fn0, fozon0, IO_READ)
      CALL IO_open (fn1, fozon1, IO_READ)
      CALL IO_open (fn2, fozon2, IO_READ)

      io_file_id0 = fozon0%file_id
      io_file_id1 = fozon1%file_id
      io_file_id2 = fozon2%file_id

      ! Check resolution
      CALL io_inq_dimid  (io_file_id1, 'lat', io_var_id)
      CALL io_inq_dimlen (io_file_id1, io_var_id, io_ngl)
      CALL io_inq_dimid  (io_file_id1, 'lon', io_var_id)
      CALL io_inq_dimlen (io_file_id1, io_var_id, io_nlon)
      CALL io_inq_dimid  (io_file_id1, 'plev', io_var_id)
      CALL io_inq_dimlen (io_file_id1, io_var_id, noz)

      IF (io_ngl/=ngl) THEN
        WRITE(message_text,'(a,i0,i0)') &
             'read_o3clim_4: unexpected resolution ',io_nlon,io_ngl
        CALL message('',message_text)
        WRITE(message_text,'(a,i0)') &
             'expected number of latitudes = ',ngl
        CALL message('',message_text)
        WRITE(message_text,'(a,i0)') &
             'number of latitudes of ozone = ',io_ngl
        CALL message('',message_text)
        CALL finish ('read_o3clim_4','unexpected resolution')
      END IF
       
      IF (io_nlon == 1) THEN
        CALL message('','Ozone data has zonal mean climatology')
      END IF
       
    END IF

    CALL p_bcast (noz, p_io)

    ! ALLOCATE memory
    IF ( .NOT. ALLOCATED(pfoz))  ALLOCATE(pfoz(noz))
    IF ( .NOT. ALLOCATED(phoz))  ALLOCATE(phoz(noz+1))

    !     Allocate memory for sst per PE

    IF ( .NOT. ALLOCATED(ozone)) ALLOCATE(ozone(lc%nproma,noz,lc%ngpblks,0:13))

    IF (p_pe == p_io) THEN

      ALLOCATE (zin(nlon,noz,ngl,0:13))

      ! read full pressure levels of ozone climatology
      CALL io_inq_varid      (io_file_id1, 'plev', io_var_id)
      CALL io_get_var_double (io_file_id1, io_var_id, pfoz(:))

       WRITE(message_text,'(a,i0)') &
            'number of pressure levels of ozone climatology: ', noz
       CALL message('',message_text)
       CALL message('','pressure levels in [hPa]:')
       DO jk = 1, noz, 10
          jkmax = jk+10-1
          IF (jkmax > noz) THEN
             jkmax = noz
         END IF
         WRITE(yformat,'(a,i0,a)') '(',jkmax-jk+1, 'f8.2)'
         WRITE(message_text,TRIM(yformat)) pfoz(jk:jkmax)*0.01_dp
         CALL message('',message_text)
       ENDDO

      ! read ozone volume mixing ratio at full pressure levels in ppmv 
      ! read ozone 12 month
      CALL io_inq_varid (io_file_id1, 'O3', io_var_id)
      DO jk = 1, noz
        io_start(:) = (/       1,   1, jk,  1 /)
        io_count(:) = (/ io_nlon, ngl,  1, 12 /)
        ! for level jk: read io_nlon longitudes, ngl latitudes and 12 months
        CALL io_get_vara_double(io_file_id1,io_var_id,io_start,io_count, &
             &                  zin(1:io_nlon,jk,:,1:12))
        IF (io_nlon == 1) THEN
          DO i = 2, nlon
            zin(i,jk,:,1:12) = zin(1,jk,:,1:12)
          END DO
        END IF
      END DO

      ! read ozone december of last year
      CALL io_inq_varid (io_file_id0, 'O3', io_var_id)
      DO jk = 1, noz
        io_start(:) = (/       1,   1, jk, 12 /)
        io_count(:) = (/ io_nlon, ngl,  1,  1 /)
        ! for level jk: read io_nlon longitudes, ngl latitudes and 1 months
        CALL io_get_vara_double(io_file_id0,io_var_id,io_start,io_count, &
             &                  zin(1:io_nlon,jk,:,0))
        IF (io_nlon == 1) THEN
          DO i = 2, nlon
            zin(i,jk,:,0) = zin(1,jk,:,0)
          END DO
        END IF
      END DO

      ! read ozone january of next year
      CALL io_inq_varid (io_file_id2, 'O3', io_var_id)
      DO jk = 1, noz
        io_start(:) = (/       1,   1, jk,  1 /)
        io_count(:) = (/ io_nlon, ngl,  1,  1 /)
        ! for level jk: read io_nlon longitudes, ngl latitudes and 1 months
        CALL io_get_vara_double(io_file_id2,io_var_id,io_start,io_count, &
             &                  zin(1:io_nlon,jk,:,13))
        IF (io_nlon == 1) THEN
          DO i = 2, nlon
            zin(i,jk,:,13) = zin(1,jk,:,13)
          END DO
        END IF
      END DO

      CALL io_close (fozon0)
      CALL io_close (fozon1)
      CALL io_close (fozon2)

!      ! convert from ppmv to g/g
!      zin(:,:,:,:) = MAX(zin(:,:,:,:)*ppmv2gg, 0._dp)
      ! convert from volume mixing ratio to mass mixing ratio
      zin(:,:,:,:) = MAX(zin(:,:,:,:)*vmr2mmr, 0._dp)

    ENDIF

    NULLIFY (gl_ozon)
    DO i = 0, 13
      IF (p_pe == p_io) gl_ozon => zin(:,:,:,i:i)
      CALL scatter_gp(gl_ozon, ozone(:,:,:,i:i), gl_dc)
    END DO

    IF (p_pe == p_io) THEN
      DEALLOCATE (zin)
    END IF

    CALL p_bcast (pfoz, p_io)

    ! define half levels of ozone pressure grid
    ! upper boundary: ph =      0.Pa -> extrapolation of uppermost value
    ! lower boundary: ph = 125000.Pa -> extrapolation of lowermost value
    phoz(1)=0._dp
    phoz(2:noz)=(pfoz(1:noz-1)+pfoz(2:noz))/2._dp
    phoz(noz+1)=125000._dp

    lo3_set=.true.

  END SUBROUTINE read_o3clim_4


  !=======================================================================
  !
  ! C)

  FUNCTION o3clim(krow,kproma,kbdim,klev,pph,ppf)

    !- Description:
    !
    !  The ozone zonal mean and/or monthly mean climatology is interpolated
    !  to the actual time of the integration and integrated from p=0 to
    !  the surface p=ps.
    !  The time interpolated profile of ozone is interpolated to the
    !  model full levels and integrated again fromp=0 to p=ps.
    !  Finally the ozone profile on the model levels is normalized such
    !  that the integrated amount on the model grid is identical to
    !  that on the grid of the climatology.
    !
    !- Author:
    !
    !  M.A. Giorgetta, MPI, May 2000

    ! INPUT
    ! -----

    INTEGER, INTENT(in)                      :: krow   ! local latitude index
    INTEGER, INTENT(in)                      :: kproma ! number of local longitudes
    INTEGER, INTENT(in)                      :: kbdim  ! first dimension of 2-d arrays
    INTEGER, INTENT(in)                      :: klev   ! number of levels
    REAL(dp), INTENT(in),DIMENSION(kbdim,klev)   :: ppf  ! full level pressure
    REAL(dp), INTENT(in),DIMENSION(kbdim,klev+1) :: pph  ! half level pressure

    ! OUTPUT
    ! ------
    REAL(dp), DIMENSION(kproma,klev)              :: o3clim ! ozone in g/g


    ! LOCAL
    ! -----

    ! pressure integrated ozone at half levels of ozone grid,
    ! and integral at surface
    REAL(dp), DIMENSION(kproma)                 :: zozintc

    ! time interpolated ozone at full levels of model grid,
    ! pressure integrated ozone at half levels of model,
    ! and integral at surface
    REAL(dp), DIMENSION(kproma,klev)            :: zozonem !kk
    REAL(dp), DIMENSION(kproma)                 :: zozintm

    REAL(dp) :: zdp1,zdp2 ! pressure weights in linear interpolation
    INTEGER  :: jl,jk,jkk ! loop indices
    INTEGER,DIMENSION(kproma)  :: jk1,jkn   ! first and last model level in interpolation

    INTEGER,DIMENSION(kproma)            :: kwork
    LOGICAL,DIMENSION(kproma)            :: kk_flag


!kk NEC version                                      Klaus Ketelsen
!kk   1. Remove EXIT from loop for vectorizing
!kk   2. DO jl=1,kproma as inner loop to get sufficiant vector length

!kk Now some loops have higher operation cout, but they are vectorizing
!kk The logic is the same as in the original version using EXIT in loops


    ! interpolate ozone profile to model grid
    ! ---------------------------------------
    ! set ozone concentration at levels above the uppermost level of
    ! the ozone climatology to the value in the uppermost level of 
    ! the ozone climatology

    jk1(:)     = 1
    kk_flag(:) = .TRUE.
    DO jk = 1,klev
       DO jl=1,kproma
          IF (ppf(jl,jk)<=pfoz(1) .AND. kk_flag(jl)) THEN 
             zozonem(jl,jk)=zozonec_x(jl,1,krow)
             jk1(jl)=jk+1
          ELSE
             kk_flag(jl) = .FALSE.
          END IF
       END DO
    END DO

    ! set ozone concentration at levels below the lowermost level of
    ! the ozone climatology to the value in the lowermost level of 
    ! the ozone climatology
    jkn(:)=klev
    kk_flag(:) = .TRUE.
    DO jk = klev,1,-1
       DO jl=1,kproma
          IF (ppf(jl,jk)>=pfoz(noz).AND. kk_flag(jl)) THEN
             zozonem(jl,jk)=zozonec_x(jl,noz,krow)
             jkn(jl)=jk-1
          ELSE
             kk_flag(jl) = .FALSE.
          END IF
       END DO
    ENDDO

    DO jk=1,klev
       kk_flag(:) = .TRUE.
       kwork(:)   = 1
       DO jkk = 1,noz
          DO jl=1,kproma
             IF(jk >= jk1(jl) .AND. jk <= jkn(jl))  THEN
                IF (ppf(jl,jk) <= pfoz(jkk) .AND. jkk >= kwork(jl) &
                                            .AND. kk_flag(jl)) THEN
                   kwork(jl)   = jkk
                   kk_flag(jl) = .FALSE.
                END IF
             END IF
          END DO
       END DO

       DO jl=1,kproma
          IF(jk >= jk1(jl) .AND. jk <= jkn(jl))  THEN
                jkk = kwork(jl)
                ! model level is in interval ]pfoz(jkk-1),pfoz(jkk)]
                ! -> make interpolation
                zdp1=pfoz(jkk)-ppf(jl,jk)
                zdp2=ppf(jl,jk)-pfoz(jkk-1)
                zozonem(jl,jk)=(zdp1*zozonec_x(jl,jkk-1,krow) &
                      +zdp2*zozonec_x(jl,jkk,krow))/ &
                     &      (zdp1+zdp2)
          END IF
       END DO
    END DO

       ! integrate ozone profile on grid of climatology
       ! from top to surface
       ! ----------------------------------------------
    zozintc=0._dp
    kk_flag(:) = .TRUE.
    jk1(:)     = 2
    DO jk=2,noz+1
          ! integrate layers of climatology above surface
       DO jl=1,kproma
          IF (phoz(jk)<=pph(jl,klev) .AND. kk_flag(jl) ) THEN
              zozintc(jl)=zozintc(jl)+ &
               &  zozonec_x(jl,jk-1,krow)*(phoz(jk)-phoz(jk-1))
              jk1(jl) = jk+1
          ELSE
              kk_flag(jl) = .FALSE.
          END IF
       END DO
    END DO
       ! integrate layer of climatology that is intersected
       ! by the surface from upper boundary to surface
    DO jl=1,kproma
       zozintc(jl)=zozintc(jl)+ &
            &  zozonec_x(jl,jk1(jl)-1,krow)*(pph(jl,klev)-phoz(jk1(jl)-1))
    END DO

       ! integrate ozone profile on grid of model
       ! from top to surface
       ! ----------------------------------------
    zozintm=0._dp
    DO jk=2,klev+1
       DO jl=1,kproma
         zozintm(jl)=zozintm(jl) + zozonem(jl,jk-1)*(pph(jl,jk)-pph(jl,jk-1))
       END DO
    END DO

       ! normalize interpolated ozone profile such that the
       ! ozone integral computed on the model grid is equal
       ! to that integrated on the grid of the climatology
       ! --------------------------------------------------
    DO jk=1,klev
       DO jl=1,kproma
         o3clim(jl,jk)=zozonem(jl,jk)/zozintm(jl) * zozintc(jl)
       END DO
    END DO


  END FUNCTION o3clim

  !=======================================================================

  SUBROUTINE pre_o3clim_4

    USE mo_decomposition,  ONLY: ldc => local_decomposition

    INTEGER  :: jl, jk, jrow
    INTEGER  :: ngpblks, nproma

    !  Executable statements

    ngpblks = ldc% ngpblks

    IF ( .NOT. ALLOCATED(zozonec_x)) THEN
       ALLOCATE(zozonec_x(ldc%nproma,noz,ldc%ngpblks))
       zozonec_x = 0._dp
    END IF

    DO jrow = 1, ngpblks

      IF ( jrow == ldc% ngpblks ) THEN
        nproma = ldc% npromz
      ELSE
        nproma = ldc% nproma
      END IF

      DO jk=1,noz
        DO jl=1,nproma
          zozonec_x(jl,jk,jrow) = wgt1_m*ozone(jl,jk,jrow,nmw1_m)+wgt2_m*ozone(jl,jk,jrow,nmw2_m)
        END DO
      END DO

    END DO

  END SUBROUTINE pre_o3clim_4
  !=======================================================================

  SUBROUTINE pre_o3clim_3(nmonth)

    USE mo_decomposition,  ONLY: ldc => local_decomposition
    USE mo_transpose,      ONLY: reorder

    ! LOCAL
    ! -----

    ! time interpolated ozone at full levels of ozone grid
    INTEGER, INTENT (IN) :: nmonth
    REAL(dp), DIMENSION(noz)                    :: zozone
    REAL(dp)    :: zzozonec_x(ldc% nglon, noz, ldc% nglat)

    INTEGER  :: jk, jrow, jglat
    INTEGER  :: nglat

    nglat = ldc% nglat      ! local number of latitudes

    IF ( .NOT. ALLOCATED(zozonec_x)) ALLOCATE(zozonec_x(ldc%nproma,noz,ldc%ngpblks))

    DO jrow = 1, nglat

      jglat = ldc% glat(jrow)

      ! Ozone at latitude jrow at ozone full levels
      ! -------------------------------------------
      IF ( nmonth == 0 ) THEN
        ! annual cycle switched on -> interpolate in time
        zozone(:) = wgt1*ozone3(:,jglat,nmw1cl)+wgt2*ozone3(:,jglat,nmw2cl)
      ELSE
        ! annual cycle switched off
        zozone(:) = ozone3(:,jglat,nmonth)
      ENDIF

      DO jk = 1,noz
        zzozonec_x(:,jk,jrow) = zozone(jk)
      END DO

    END DO

    CALL reorder(zozonec_x, zzozonec_x)

  END SUBROUTINE pre_o3clim_3
!------------------------------------------------------------------------------
  SUBROUTINE set_years(y1,y2,y3)

    USE mo_time_control,    ONLY: next_date, get_date_components

    INTEGER, INTENT(out) :: y1, y2, y3
    INTEGER :: yr, mo, dy, hr, mn, se

    CALL get_date_components(next_date, yr, mo, dy, hr, mn, se)

    y1 = yr - 1
    y2 = yr
    y3 = yr + 1

  END SUBROUTINE set_years
!------------------------------------------------------------------------------
  SUBROUTINE cleanup_o3clim
    !----------------------------
    ! deallocate module variables
    !----------------------------
    IF (ALLOCATED(zozonec_x)) DEALLOCATE (zozonec_x)
    IF (ALLOCATED(ozone3))    DEALLOCATE (ozone3)
    IF (ALLOCATED(ozone))     DEALLOCATE (ozone)
    IF (ALLOCATED(pfoz))      DEALLOCATE (pfoz)
    IF (ALLOCATED(phoz))      DEALLOCATE (phoz)

  END SUBROUTINE cleanup_o3clim
!------------------------------------------------------------------------------
END MODULE mo_o3clim

