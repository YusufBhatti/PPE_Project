!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_emi_volcano.f90
!!
!! \brief
!! handles emissions from continuously degassing and explosive volcanoes
!!
!! \author M. Schultz (FZ Juelich)
!!
!! \responsible_coder
!! M. Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# M. Schultz (FZ Juelich) - original code (2010-02-25)
!!
!! \limitations
!! input files have to be the reformatted ASCII files
!! (reformatted from original to include header lines)
!!
!! \details
!! Input data: AEROCOM climatology of SO2 emissions
!! Emissions are provided with lon and lat coordinates and a min and max altitude
!! for each cell. They are mapped onto the current ECHAM grid and stored as
!! a 3-dimensional boundary condition.
!! This code was originally developed by Philip Stier.
!! It has been almost entirely rewritten in order to allow for a commented
!! ascii file and to make use of the boundary condition concept.
!! Major conceptual change: emissions are now distributed onto the ECHAM grid once
!! at the model startup (formerly, the altitude->model level conversion was done
!! at every time step). Once the boundary condition scheme can store fields
!! on native altitude grids, this can be changed again, but I doubt that it makes
!! much difference.
!! The old code had a "feature" which has been corrected in this version: when two
!! lines in the ascii file would point to the same lon and lat coordinates, the 
!! minlevel and maxlevel fields used to be overwritten (so the first volcano would
!! emit at the levels of the second one). The new code can accumulate emission fluxes
!! for all volcanoes independently.
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
MODULE mo_hammoz_emi_volcano

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_emi_volcano
  PUBLIC :: read_emi_volcano
!!  PUBLIC :: read_ascii_file
!!  PUBLIC :: geoindex


  INTEGER            :: ibc_volcc, ibc_volce      ! Index to boundary conditions

CONTAINS

  SUBROUTINE init_emi_volcano(ktype)

  USE mo_exception,            ONLY: finish
  USE mo_boundary_condition,   ONLY: bc_find, bc_modify, BC_EVERYWHERE

  IMPLICIT NONE

  INTEGER, INTENT(in)          :: ktype  ! type of emissions to initialize:
                                         ! 1 = continuous degassing
                                         ! 2 = explosive

  CHARACTER(LEN=64)            :: bc_name
  INTEGER                      :: ibc     ! index to boundary condition

  IF (ktype < 1 .OR. ktype > 2) CALL finish('init_emi_volcano',      &
                                            'Invalid value of ktype! (Must be 1 or 2)')
  !-- locate corresponding boundary condition
  !   bc_define was called from init_emissions in mo_emi_interface
  IF (ktype == 1) THEN
    bc_name = 'VOLCC emissions of SO2'
  ELSE
    bc_name = 'VOLCE emissions of SO2'
  END IF
  CALL bc_find(bc_name, ibc)
  IF (ibc < 1) CALL finish('init_emi_volcano', 'Cannot find boundary condition for '//TRIM(bc_name))
  IF (ktype == 1) THEN
    ibc_volcc = ibc
  ELSE
    ibc_volce = ibc
  END IF
  !-- modify boundary condition to fit the specific needs of this module
  CALL bc_modify(ibc, bc_domain=BC_EVERYWHERE, bc_ndims=3, ef_actual_unit='kg m-2 s-1')

  END SUBROUTINE init_emi_volcano

  !---------------------------------------------------------------------------------------
  ! read the ascii file with the volcanic SO2 emission fluxes and distribute values on
  ! ECHAM 3D grid
  !!! NOTE: in future versions of HAMMOZ the values shall be stored on a native altitude 
  !!! grid and interpolated onto ECHAM grid at each time step (this corresponds to the ECHAM philosophy)
  !!! Need to decide on spacing of native grid (dz = 200 m ?, top = 12 km?)

  SUBROUTINE read_emi_volcano(kproma, kbdim, klev, ktype)

  USE mo_kind,                 ONLY: dp
  USE mo_exception,            ONLY: finish, message, message_text, em_info, em_error, em_warn
  USE mo_boundary_condition,   ONLY: bc_set
  USE mo_physical_constants,   ONLY: rgrav
  USE mo_memory_g3b,           ONLY: geosp
  USE mo_vphysc,               ONLY: vphysc
  USE mo_control,              ONLY: nlon, ngl, nlev
  USE mo_mpi,                  ONLY: p_parallel_io
  USE mo_transpose,            ONLY: scatter_gp, gather_gp
  USE mo_decomposition,        ONLY: gl_dc => global_decomposition, ldc => local_decomposition
  USE mo_gaussgrid,            ONLY: gridarea
  USE mo_submodel,             ONLY: emi_basepath
!!baustelle: workaround -- not all grheightm1 values defined when this routine is caleld with lfirst=true
USE mo_control,   ONLY: vct, nvclev

  IMPLICIT NONE

  INTEGER, INTENT(in)          :: kproma, kbdim, klev
  INTEGER, INTENT(in)          :: ktype  ! type of emissions to initialize:
                                         ! 1 = continuous degassing
                                         ! 2 = explosive

  CHARACTER(LEN=256)       :: cfilename, clabel
  CHARACTER(len=22)        :: formstr = '(2(1x,i4),3(1x,e10.4))'
  CHARACTER(len=30)        :: units
  INTEGER                  :: j, jk, jrow, jvolc, ibc, ilower, iupper, klon, klat
                              ! note: iupper and ilower refer to altitude, not level index!
  REAL(dp)                 :: zsumemi, zfac,                     &
                              zemivolc(kbdim, klev, ldc%ngpblks) 
  REAL(dp), TARGET         :: gl_emivolc(nlon, nlev, ngl),       &
                              gl_geosp(nlon, ngl),               &
                              gl_heights(nlon, nlev, ngl)
  REAL(dp), POINTER        :: gl2d(:,:), gl3d(:,:,:)
  REAL(dp), ALLOCATABLE    :: zlats(:), zlons(:), zemiss(:),     &    ! fields for reading
                              zheights_low(:), zheights_high(:)

  LOGICAL, ALLOCATABLE :: ll1(:) !SF #458 for replacing the former WHERE statement

  INTEGER, PARAMETER       :: maxlines=20000     ! max. number of lines in file
  INTEGER, PARAMETER       :: ncols   =5         ! number of columns in file
  INTEGER                  :: nrows              ! number of data values in file
  REAL(dp), ALLOCATABLE    :: zbuf(:,:)

!!baustelle!! - workaround
REAL(dp) :: as(nvclev), bs(nvclev), ph(nvclev), zf(klev)    ! vertical coordinate tables  (see mo_tpcore) [nvclev = klev+1, I hope]

 as(:) = vct(1:nvclev)
 bs(:) = vct(nvclev+1:2*nvclev)
 ph(:) = as(:)+bs(:)*1.e5_dp    ! nominal pressure values
 ph(1) = 1.e-19_dp
! convert to "height above surface"
 DO jk=1,klev
   zf(jk) = 8000._dp * (LOG(ph(klev+1)/1.e5_dp) - LOG(ph(jk)/1.e5_dp))
 END DO


  IF (ktype < 1 .OR. ktype > 2) CALL finish('init_emi_volcano',      &
                                            'Invalid value of ktype! (Must be 1 or 2)')
  !-- Initialisation
  gl_emivolc(:,:,:) = 0._dp

  !-- locate corresponding boundary condition
  IF (ktype == 1) THEN
    ibc = ibc_volcc
    cfilename = 'continuous_volcanos.dat'
    clabel = 'SO2 emissions from continuous volcanic degassing'
  ELSE
    ibc = ibc_volce
    cfilename = 'explosive_volcanos.dat'
    clabel = 'SO2 emissions from explosive volcanic eruptions'
  END IF
!>>SF the orig violates the symlink construction
  !cfilename = TRIM(emi_basepath)//'/'//TRIM(cfilename)
  cfilename = TRIM(cfilename)
!<<SF

  !-- gather global values of geopotential and gridbox height
  gl2d => gl_geosp
  CALL gather_gp (gl2d, geosp, gl_dc)
  gl3d => gl_heights
  CALL gather_gp (gl3d, vphysc%grheightm1, gl_dc)
!!baustelle!! -- workaround: overwrite gl_heights
DO jk=1,klev
gl_heights(:,jk,:) = zf(jk)
END DO 

  !-- read ascii file
  IF (p_parallel_io) THEN
    !-- calculate grid cell altitudes
    DO jk=1,klev
      gl_heights(:,jk,:) = gl_heights(:,jk,:) + gl_geosp(:,:)*rgrav
    END DO 
 
    CALL message('init_emi_volcano', 'Reading volcanic SO2 emissions from file '//TRIM(cfilename), &
                 level=em_info)

    ALLOCATE (zbuf(maxlines,ncols))
    CALL read_ascii_file(cfilename, ncols, maxlines, nrows, zbuf, formstr, 'init_emi_volcano', units)
    ! units should be 'kg y-1' for right calculations
    IF (TRIM(ADJUSTL(units)) .ne. 'kg y-1') THEN
      WRITE(message_text,'(a,a)') 'units of input file should be [kg y-1], units found are: ', TRIM(units)
      CALL finish('read_emi_volcano', message_text)
    ENDIF

    ! distribute zbuf into vectors
    ALLOCATE(zlons(nrows))
    ALLOCATE(zlats(nrows))
    ALLOCATE(zemiss(nrows))
    ALLOCATE(zheights_low(nrows))
    ALLOCATE(zheights_high(nrows))
    ALLOCATE(ll1(nrows)) !SF #458
    zlons(:)         = zbuf(1:nrows,1)
    zlats(:)         = zbuf(1:nrows,2)
    zemiss(:)        = zbuf(1:nrows,3)
    zheights_low(:)  = zbuf(1:nrows,4)
    zheights_high(:) = zbuf(1:nrows,5)
    ! change lons to 0..360:

    !>>SF #458: replaced the former WHERE statement by a logical array and MERGE
    ll1(:) = (zlons(:) < 0)
    zlons(:) = MERGE(zlons(:) + 360, zlons(:), ll1(:))
    !<<SF

    ! get total emissions
    zsumemi=0._dp
    DO j=1,nrows
      zsumemi = zsumemi + zemiss(j)
    END DO
    WRITE(message_text,'(a,e25.15,a)') 'Global total '//TRIM(clabel)//' = ',zsumemi,' kg(SO2) year-1'
    CALL message('read_emi_volcano', message_text, level=em_info)
    !-- distribute volcanic emissions onto ECHAM grid 
    DO jvolc=1, nrows
      !--- Find corresponding lat-lon indices:
      CALL geoindex(REAL(zlons(jvolc),dp),REAL(zlats(jvolc),dp),klon,klat)
      !--- Find model levels corresponding to min and max altitude of emissions
      ilower = -1
      DO jk = 1, klev
        IF (zheights_low(jvolc) <= gl_heights(klon,jk,klat)) ilower=jk+1
      END DO
      IF (ilower > klev) ilower = klev
      iupper = klev
      DO jk = klev, 1, -1
        IF (zheights_high(jvolc) >= gl_heights(klon,jk,klat)) iupper=jk
      END DO

      IF (ilower == -1) THEN
          WRITE(message_text,'(a,i0,e25.15,a,i0,i0,e25.15)') 'ilower=-1 for jvolc, zheights_low = ',jvolc,zheights_low(jvolc), &
                                                             'klon,klat,gl_heights=',klon,klat,gl_heights(klon,:,klat)
          CALL message('read_emi_volcano', message_text, level=em_warn)
      ENDIF

      IF (ilower < iupper) THEN
        WRITE(message_text,'(a,i0)') 'Error in emission levels! File '//TRIM(cfilename)//', line=',jvolc
        CALL message('read_emi_volcano', message_text, level=em_error)
      END IF
      !--- Enter portion of emission flux in each affected grid cell
      !    Unit conversion from [kg(SO2) grid_cell-1 year-1] to [kg(SO2) m-2 s-1]
      zfac = 1._dp/REAL(ilower-iupper+1, kind=dp)
      DO jk=iupper,ilower
        gl_emivolc(klon, jk, klat) = gl_emivolc(klon, jk, klat)             &
                                     + zfac*zemiss(jvolc)/(gridarea(klat)*86400._dp*365._dp)
      END DO
    END DO

  END IF      ! p_parallel_io

  !-- Scatter global emission field and set as boundary condition value
  gl3d => gl_emivolc
  CALL scatter_gp(gl3d, zemivolc, gl_dc)
  DO jrow = 1,ldc%ngpblks
    CALL bc_set(ibc, kproma, jrow, zemivolc(:,:,jrow))
  END DO

  !-- Deallocate temporal fields
  IF (p_parallel_io) THEN
     DEALLOCATE(zbuf,zheights_low,zheights_high,zemiss,zlats,zlons)
  END IF

  END SUBROUTINE read_emi_volcano

   
  SUBROUTINE read_ascii_file (filename, ncols, maxcount, ncount, buffer, formatstr, calling_routine, units)
  ! open and read a plain format ascii data file
  ! The file may contain comment lines (starting with '#')
  ! filename, number of columns and maximum number of lines must be provided,
  ! and a buffer array with sufficient capacity (dimension maxlines, ncols).
  ! The subroutine will return the number of lines that were actually read.
  ! Maximum number of columns is 32.
  !
  ! Note: the formatstr is not used at present as this led to an error on the IBM power 6

  USE mo_kind,          ONLY: dp
  USE mo_exception,     ONLY: finish
  USE mo_filename,      ONLY: find_next_free_unit
  USE mo_mpi,           ONLY: p_parallel_io

  IMPLICIT NONE

  ! Arguments
  CHARACTER(len=*), INTENT(in)    :: filename, formatstr, calling_routine
  INTEGER, INTENT(in)             :: ncols, maxcount
  INTEGER, INTENT(out)            :: ncount
  REAL(dp), INTENT(inout)         :: buffer(maxcount, ncols)
  CHARACTER(len=*), INTENT(inout) :: units

  ! Local variables
  INTEGER             :: iunit, ierr, ilines
  LOGICAL             :: lex
  INTEGER, PARAMETER  :: icolmax = 32
  REAL(dp)            :: zvalbuf(icolmax)
  CHARACTER(len=255)  :: linebuf

  buffer(:,:) = 0._dp
  ncount = 0           ! number of data lines
  ilines = 0           ! number of lines in file (including comments)
  linebuf = '                                                                                      '

  IF (p_parallel_io) THEN
    iunit = find_next_free_unit (20, 99)
    INQUIRE (file=filename,exist=lex)
    IF (.NOT. lex) THEN
      CALL finish( calling_routine,   &
            'Failed to open file '//TRIM(filename)//'. File does not exist or is not accessible.')
    END IF
    ! Open file
    OPEN (iunit,file=TRIM(filename),form='FORMATTED',iostat=ierr)
    IF (ierr /= 0) CALL finish('read_ascii_file', 'Could not open file '//TRIM(filename))
    REWIND(iunit)
    DO WHILE (ierr == 0)
      READ(iunit,FMT='(a)',IOSTAT=ierr) linebuf
      ilines = ilines + 1
      IF (ierr < 0) EXIT     ! eof reached
      IF (ierr > 0) THEN
        write(linebuf,'(i0)') ilines
        CALL finish( calling_routine, &
                'Failed to read line '//TRIM(linebuf)//' in file '//TRIM(filename) )
      END IF
      ! analyze ascii line
      IF (linebuf(1:1) == '#') CYCLE   ! comment line
      IF (linebuf(1:6) == 'units:') THEN
        units=TRIM(linebuf(7:))
        CYCLE
      ENDIF
!       READ(linebuf,FMT=formatstr,IOSTAT=ierr) zvalbuf(1:ncols)
      READ(linebuf,FMT=*,IOSTAT=ierr) zvalbuf(1:ncols)
      IF (ierr /= 0) THEN
        write(linebuf,'(i0)') ilines
        CALL finish( calling_routine, &
                     'Failed to read data in line '//TRIM(linebuf)//' in file '//TRIM(filename) )
      END IF
      ncount = ncount + 1
      IF (ncount > maxcount) THEN
        write(linebuf,'(i0)') maxcount
        CALL finish( calling_routine, &
                    'Error reading file '//TRIM(filename)//': number of lines exceeds '//linebuf )
      END IF
      buffer(ncount, 1:ncols) = zvalbuf(1:ncols)
    END DO

    CLOSE(iunit)
  END IF

  END SUBROUTINE read_ascii_file


  SUBROUTINE geoindex (plon, plat, klon, klat)

    !   *geoindex* calculates the corresponding lat-lon (klon,klat) index
    !              for given lat-lon [-90,90][0,360] coordinates (plon, plat)
    !
    !   Authors:
    !   ---------
    !   Philip Stier, MPI-MET            2002
    !
    !   Externals
    !   -----------
    !   none

    USE mo_kind,      ONLY: dp
    USE mo_exception, ONLY: finish
    USE mo_control,   ONLY: ngl, nlon
    USE mo_gaussgrid, ONLY: philat, philon

    IMPLICIT NONE

    REAL(dp),INTENT(IN)  :: plon, plat   ! Coordinates in degrees
    INTEGER, INTENT(OUT) :: klon, klat   ! Corresponding indices

    REAL(dp)             :: zxdif, zydif
    INTEGER              :: jlon, jlat


    IF((NINT(plon)<0).OR.(NINT(plon)>360).OR.(NINT(plat)<-90).OR.(NINT(plat)>90)) THEN
       CALL finish('geoindex:', 'Coordinates out of range')
    END IF

    klon=-999
    klat=-999

    zxdif=360._dp
    zydif=180._dp

    !--- 1) Search for closest longitude in [0,360}:

    DO jlon = 1, nlon
       IF(ABS(plon-philon(jlon)) < zxdif) THEN
          klon=jlon
          zxdif=ABS(plon-philon(jlon))
       END IF
    END DO

    !--- 2) Search for closest latitude [90,-90]:

    DO jlat = 1, ngl
       IF( ABS(plat-philat(jlat)) < zydif ) THEN
          klat=jlat
          zydif=ABS(plat-philat(jlat))
       END IF
    END DO


    IF(klon==-999 .OR. klat==-999) CALL finish('geoindex:', 'No index found')

  END SUBROUTINE geoindex

END MODULE mo_hammoz_emi_volcano

