!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_input_interface
IMPLICIT NONE

PUBLIC

INTEGER, PARAMETER :: INPUT_MODEL_CBALANCE = 1
INTEGER, PARAMETER :: INPUT_MODEL_JSBACH   = 2
INTEGER, PARAMETER :: INPUT_MODEL_ECHAM    = 3
INTEGER, PARAMETER :: INPUT_MODEL_ICON     = 4 ! Not yet in use

CHARACTER (LEN=32) :: dim_x_name, dim_y_name, dim_horiz_name

CONTAINS

  SUBROUTINE InputInitialize(model,debug,calendar,file_name,nprocx,nprocy) ! Present setup does not allow chunked reading
  USE mo_exception,    ONLY: finish, message
  USE mo_kind,         ONLY: dp
  USE mo_time_base,    only: get_calendar_type, JULIAN, CYL360, CYL365
  USE mo_time_control, ONLY: resume_date, dt_start, get_date_components, lresume, delta_time
  USE mo_mpi,          ONLY: p_pe, p_nprocs
  USE mo_control,      ONLY: nproca, nprocb, nprocio, nproma
  USE mo_filename,     ONLY: find_next_free_unit
  USE mo_io,           ONLY: IO_timestep
  USE mo_input_calendar, ONLY: INPUT_CALENDAR_GREGORIAN, INPUT_CALENDAR_DAY360, INPUT_CALENDAR_DAY365
  USE mo_input
  INTEGER,           INTENT(in) :: model
  INTEGER, OPTIONAL, INTENT(in) :: calendar, nprocx, nprocy
  CHARACTER (len=*), INTENT(in), OPTIONAL :: debug, file_name

  TYPE (input_file_list), POINTER :: IFile
  TYPE (input_dim_list),  POINTER :: dim_par, dim_lon, dim_lat
  LOGICAL, POINTER :: lmask(:), lAtPE(:,:)
  INTEGER :: mod_tim(6), ts, cal
  INTEGER :: nglobal, npts, nlarge
  INTEGER :: nx, ny, na, nblock, col_blocks, col_blocks_per_pe, n_hemisphere_lat
  INTEGER :: pe_row, pe_col, i, baseLo, baseHi, nlon, nlat, lons, lone, lonsz
  REAL(dp), POINTER :: slm(:,:)

    IF (nprocio /= 0) THEN 
      CALL message('WARNING','mo_input does presently not work together with IO_server. mo_input has been disabled!')
      RETURN
    ENDIF
 
    ! Determine correct time settings for mo_input
    SELECT CASE (get_calendar_type())
      CASE (JULIAN)
        cal = INPUT_CALENDAR_GREGORIAN
      CASE (CYL360)
        cal = INPUT_CALENDAR_DAY360
      CASE (CYL365)
        cal = INPUT_CALENDAR_DAY365
      CASE DEFAULT 
        CALL finish('InputInitialize','mo_input presently supports only Julian, Gregorian 360, 365 and 366 day calendars')
    END SELECT
    IF (PRESENT(calendar)) cal = calendar

    IF (lresume) then
      CALL get_date_components(resume_date,mod_tim(1),mod_tim(2),mod_tim(3),mod_tim(4),mod_tim(5),mod_tim(6))
    ELSE
      mod_tim = dt_start
    ENDIF
    IF (delta_time /= 0._dp) THEN
      ts = INT(delta_time+0.5)
    ELSE
      ts = IO_timestep
    ENDIF

    ! Intialize input module and read namelist specifications
    CALL InputInit(ts,mod_tim(1), mod_tim(2), mod_tim(3), mod_tim(4)*3600+ mod_tim(5)*60+ mod_tim(6), & 
                     dt_start(1),dt_start(2),dt_start(3),dt_start(4)*3600+dt_start(5)*60+dt_start(6), &
                   npe=p_nprocs,pe=p_pe,iope=0,sDbg=debug,calendar=cal)
    IF (model == INPUT_MODEL_CBALANCE) THEN
      CALL InputOptImport('namelist.cbalone',find_next_free_unit(30,100))
    ELSE
      CALL InputOptImport('namelist.echam,namelist.jsbach',find_next_free_unit(30,100))
    ENDIF

    ! Register needed standard dimensions
    CALL InputDimAdd('lon',alt_name='X',dim_ref=dim_lon,lCyclic=.TRUE.)
    CALL InputDimAdd('lat',alt_name='Y',dim_ref=dim_lat,ldec   =.TRUE.)
    CALL InputDimAdd('time',lRecSep=.TRUE.)
    CALL InputDimAdd('ntiles',alt_name='tiles') ! Special but commonly used JSBACH dimension

    ! Obtain dimension sizes and land-sea mask from initial file
    NULLIFY(slm)
    IF (PRESENT(file_name)) THEN
      IFile => InputFileAdd(TRIM(file_name))
    ELSE
      IFile => InputFileAdd('jsbach.nc')
    ENDIF

    CALL InputVarAdd('lon,lat',slm,'slm',IFile,dt_update=0)
    CALL InputGetData()
    CALL InputDimInq(ref=dim_lon,global=nlon)
    CALL InputDimInq(ref=dim_lat,global=nlat)
    
    ! Adjust setup if necessary
    IF (PRESENT(nprocx) .AND. PRESENT(nprocy)) THEN
      nx = nprocx
      ny = nprocy
    ELSE
      nx = nprocb
      ny = nproca
    ENDIF
    IF (nproma == 0) THEN
      na = nlon / nx ! In ECHAM nproma==0 is only corrected in the decomposition structures, not in the global variable
    ELSE
      na = nproma
    ENDIF

    ! Model specific settings and parallel setup
    dim_x_name = 'lon'
    dim_y_name = 'lat'
    dim_horiz_name = 'landpoints'
    SELECT CASE (model)
      CASE (INPUT_MODEL_CBALANCE)
        CALL InputDimPackedAdd('landpoints','lon,lat',lmask=slm>0.5_dp)
      CASE (INPUT_MODEL_JSBACH)
        CALL InputDimPackedAdd('landpoints','lon,lat',lmask=slm>0.5_dp,dim_ref=dim_par)
        nglobal = COUNT(slm>0.5_dp)
        npts    = nglobal/nx
        nlarge  = nglobal-npts*nx ! Missing points - this number of PEs should have a point more
        IF (p_pe < nlarge) THEN
          CALL InputDimLocalSet(lo=p_pe*(npts+1)+1,hi=(p_pe+1)*(npts+1),dim_ref=dim_par)
        ELSE
          CALL InputDimLocalSet(lo=nlarge+p_pe*npts+1,hi=nlarge+(p_pe+1)*npts,dim_ref=dim_par)
        ENDIF
      CASE (INPUT_MODEL_ECHAM)
        dim_x_name = 'nproma'
        dim_y_name = 'nblock'
        CALL InputDimAdd('nproma',Extent=na)
        pe_col           =     p_pe/ny
        pe_row           = MOD(p_pe,ny)
        n_hemisphere_lat = INT(CEILING(REAL(nlat,dp)/REAL(ny*2,dp)))
       ! lonsz            = INT(CEILING(REAL(nlon,dp)/REAL(nx  ,dp)))
        lonsz            = nlon/nx
        nlarge           = nlon-lonsz*nx
        IF ((nlon/nx/na)*nx*na == nlon) THEN ! Case is also covered by the general case below, but causes less overhead
          nblock = nlon*nlat/na 
          CALL InputDimAdd('nblock',Extent=nblock,dim_ref=dim_par)
          CALL InputDimEquivalence('nproma,nblock','lon,lat')
          col_blocks         = nlon / na
          col_blocks_per_pe  = nlon / na / nx
          baseLo = col_blocks_per_pe*pe_col + col_blocks * n_hemisphere_lat * pe_row
          baseHi = col_blocks_per_pe*pe_col - col_blocks * n_hemisphere_lat * pe_row + nblock
          ALLOCATE(lmask(nblock))
          lmask(:) = .FALSE.
          DO i=1,n_hemisphere_lat
            baseHi = baseHi - col_blocks
            lmask(baseLo+1:baseLo+col_blocks_per_pe) = .TRUE.
            lmask(baseHi+1:baseHi+col_blocks_per_pe) = .TRUE.
            baseLo = baseLo + col_blocks
          ENDDO
          CALL InputDimLocalSet(dim_ref=dim_par,lmask=lmask)
          DEALLOCATE(lmask)
          CALL InputDimPackedAdd('landpoints','nproma,nblock',lmask=slm>0.5_dp)
        ELSE
          baseLo =        n_hemisphere_lat * pe_row 
          baseHi = nlat - n_hemisphere_lat *(pe_row+1)
          ! lons   = lonsz * pe_col + 1
          ! lone   = MIN(lons+lonsz-1,nlon)
          IF (pe_col < nlarge) THEN
            lonsz = lonsz + 1
            lons  = pe_col*lonsz + 1
          ELSE
            lons = pe_col*lonsz + 1 + nlarge
          ENDIF
          lone = lons+lonsz-1
          ALLOCATE(lAtPE(nlon,nlat))
          lAtPE = .FALSE.
          lAtPE(lons:lone,baseLo+1:baseLo+n_hemisphere_lat) = .TRUE.
          lAtPE(lons:lone,baseHi+1:baseHi+n_hemisphere_lat) = .TRUE.
          CALL InputDimAdd('lon_iloc',Extent=lonsz) ! Just to get proper sizes of internal dimensions
          CALL InputDimLocalSetMulti('lon,lat',lmask=lAtPE)
          DEALLOCATE(lAtPE)
          CALL InputDimPackedAdd('landpoints','lon_iloc,lat_iloc',lmask=slm>0.5_dp)
          ! Preparation of reading of ECHAM grid-point variables.
          CALL InputDimAdd('nblock')
          CALL InputDimEquivalence('nproma,nblock','lon_iloc,lat_iloc')
        ENDIF
    END SELECT
    CALL InputDimSync()

    CALL message('InputInitialize','mo_input initialized')

  END SUBROUTINE InputInitialize

END MODULE mo_input_interface
