!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! mo_aero_volc_tab: use tabulated optical properties for volcanic aerosols
!!
!! @author Sebastian Rast, MPI Met, Hamburg
!!
!! $ID: n/a$
!!
!! @par Revision History
!! original source by J.S.Rast (2011-02-23), based on code 
!! by U. Schlese, M. Esch, S. Lorenz, and U. Niemeier (old mo_volc_data)
!!
MODULE mo_aero_volc_tab

  USE mo_kind,                ONLY: wp
  USE mo_decomposition,       ONLY: ldc=>local_decomposition, &
                                    global_decomposition
  USE mo_exception,           ONLY: finish,message
  USE mo_filename,            ONLY: find_next_free_unit
  USE mo_rrtm_params,         ONLY: nbndlw
  USE mo_radiation_parameters, &
                              ONLY: jpsw => nb_sw !number of solar radiation bands
  USE mo_mpi,                 ONLY: p_parallel_io, p_io, p_bcast
  USE mo_read_netcdf77,       ONLY: read_var_hs_nf77_3d  
  USE mo_transpose,           ONLY: scatter_gp
  USE mo_interpo,             ONLY: nmw1_m, nmw2_m, wgt1_m, wgt2_m
  USE mo_gaussgrid,           ONLY: philat
  USE mo_geoloc,              ONLY: ilat
  USE mo_control,             ONLY: nvclev, nlev, vct
  USE mo_time_control,        ONLY: radiation_date, get_date_components
  USE mo_time_conversion,     ONLY: day_in_year, year_len, day_len

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: su_aero_prop_ham, su_aero_prop_crow, & !< setup memory HAM/Crowley
            read_aero_volc_tables, & !< read general aerosol properties
            read_aero_prop_ham, read_aero_prop_crow, & !< read time dep. data
            add_aop_volc_ham, & !< add aerosol opt. prop. (HAM) to background
            add_aop_volc_crow, & !< add aerosol opt. prop. (Crowley) to backg. 
            cleanup_aero_volc_tab_ham, & !< clean up memory HAM aerosols
            cleanup_aero_volc_tab_crow !< clean up memory Crowley aerosols

  PUBLIC :: extrat    ! extinction ratio
  PUBLIC :: asym      ! asymmetry factor
  PUBLIC :: ssa       ! single scattering albedo

  CHARACTER(len=20),PARAMETER   :: fname='aero_volc_tables.dat'
  CHARACTER(len=9),PARAMETER    :: cfname_base_ham='aoddz_ham'
  CHARACTER(len=12),PARAMETER   :: cfname_base_crow='aodreff_crow'

  INTEGER, PARAMETER    :: nlat_crow=4, ntstep_crow=36
  REAL(wp), PARAMETER   :: dt_crow=1._wp/REAL(ntstep_crow,wp)
  REAL(wp), ALLOCATABLE :: extrat(:,:) !< extinction relative to ext at 550nm
  REAL(wp), ALLOCATABLE :: ssa   (:,:) !< single scattering albedo
  REAL(wp), ALLOCATABLE :: asym  (:,:) !< asymmetry factor
  REAL(wp), ALLOCATABLE :: aodh  (:,:,:,:) !< extinction at 550 nm (HAM)
  REAL(wp), ALLOCATABLE :: reffh (:,:,:,:) !< effective radius (HAM)
  REAL(wp), ALLOCATABLE :: fyear_crow(:)  !< fractional year (Crowley)
  REAL(wp), ALLOCATABLE :: aod_crow(:,:)  !< aod 550 nm (Crowley)
  REAL(wp), ALLOCATABLE :: reff_crow(:,:) !< effective radius (Crowley)
  REAL(wp), ALLOCATABLE :: pl_weights(:) !< weights for aod (Crowley)
  REAL(wp)              :: delta_reff, delta_reff_i, reff_min


  INTEGER               :: nrad !< number of different effective radii
  INTEGER               :: nlev_crow_min, nlev_crow_max !< min/max level
                           ! with aod weights > 0 for Crowley aerosols

  LOGICAL               :: laero_set_ham=.false.
  LOGICAL               :: laero_set_crow=.false.

CONTAINS

!! @par Description:
!! set up memory for optical aerosol parameters extracted from echam-HAM
!! simulations
!!
!! @par Revision History:
!! original source by J.S. Rast (2011-03-23)

SUBROUTINE su_aero_prop_ham
!
! !USES:
! !LOCAL VARIABLES

  INTEGER, PARAMETER             :: nmonths=12

! allocate memory for optical properties
  ALLOCATE(aodh(ldc%nproma,ldc%nlev,ldc%ngpblks,0:nmonths+1))
  ALLOCATE(reffh(ldc%nproma,ldc%nlev,ldc%ngpblks,0:nmonths+1))
END SUBROUTINE su_aero_prop_ham

!! @par Description:
!! set up memory for optical aerosol parameters given bye the data set
!! of T. Crowley
!!
!! @par Revision History:
!! original source by J.S. Rast (2011-04-06)

SUBROUTINE su_aero_prop_crow

  REAL(wp)             :: hyai(nvclev), hybi(nvclev)
  REAL(wp)             :: hyam(nlev), hybm(nlev)
  REAL(wp)             :: zp_jlev, zw_sum
  INTEGER              :: jlev
  CHARACTER(LEN=128)   :: cmessage

! allocate memory for optical properties (they will be read year by year)
  ALLOCATE(fyear_crow(ntstep_crow+1))
  ALLOCATE(aod_crow(ldc%nlat,ntstep_crow+1))
  ALLOCATE(reff_crow(ldc%nlat,ntstep_crow+1))
  ALLOCATE(pl_weights(nlev))
! define height profile for Crowley aerosols
! height levels in echam:
  hyai = vct(1:nvclev)
  hybi = vct(nvclev+1:2*nvclev)
  DO jlev=1,nlev
    hyam(jlev)=0.5_wp*(hyai(jlev)+hyai(jlev+1))
    hybm(jlev)=0.5_wp*(hybi(jlev)+hybi(jlev+1))
  END DO
! aerosol optical depth is distributed to echam levels with following weights
! [30hPa,40hPa[: 0.25; [40hPa,50hPa[: 0.5; [50hPa,60hPa[: 0.25.
! the distribution is based on a standard atmosphere with p_s=1013.25hPa
  pl_weights=0._wp
  nlev_crow_min=nlev
  nlev_crow_max=0
  DO jlev=1,nlev
    zp_jlev=hyam(jlev)+hybm(jlev)*101325._wp
    IF (zp_jlev>=3000._wp.AND.jlev<nlev_crow_min) nlev_crow_min=jlev
    IF (zp_jlev<6000._wp.AND.jlev>nlev_crow_max) nlev_crow_max=jlev
    IF (zp_jlev>=3000._wp.AND.zp_jlev<4000._wp) THEN
      pl_weights(jlev)=0.25_wp
    ELSE IF (zp_jlev>=4000._wp.AND.zp_jlev<5000._wp) THEN
      pl_weights(jlev)=0.5_wp
    ELSE IF (zp_jlev>=5000._wp.AND.zp_jlev<6000._wp) THEN
      pl_weights(jlev)=0.25_wp
    END IF
  END DO
! normalization
  zw_sum=SUM(pl_weights(nlev_crow_min:nlev_crow_max))
  IF (zw_sum<=0._wp) THEN
    CALL finish('su_aero_prop_crow(mo_aero_volc_tab)', 'there are no'// &
                ' echam levels suitable for volcanic aerosols')
  ELSE
    pl_weights(nlev_crow_min:nlev_crow_max)= &
              pl_weights(nlev_crow_min:nlev_crow_max)/zw_sum
  END IF
  IF (p_parallel_io) THEN
    CALL message('##########',' su_aero_prop_crow :##########')
    CALL message('','Model level dependent AOD weights for Crowley aerosols')
    CALL message('','model level    pressure/hPa      weight               ')
    DO jlev=1,nlev
      WRITE(cmessage,'(i4,5x,g13.3,3x,g13.3)')  &
           jlev,(hyam(jlev)+hybm(jlev)*101325._wp)*0.01_wp,pl_weights(jlev)
      CALL message('','|  '//TRIM(cmessage))
    END DO
  END IF

END SUBROUTINE su_aero_prop_crow

!! @par Description:
!! Reads lookup tables for extrapolation of aerosol optical properties given 
!! for 550 nm to other wavelengths
!!
!! @par Revision History:
!! original source by J.S. Rast (2011-03-23)

SUBROUTINE read_aero_volc_tables

! !LOCAL VARIABLES

  INTEGER, ALLOCATABLE  :: nbandmap(:)
  REAL(wp), ALLOCATABLE :: zextrat(:,:) !< extinction relative to ext at 550nm
  REAL(wp), ALLOCATABLE :: zssa   (:,:) !< single scattering albedo
  REAL(wp), ALLOCATABLE :: zasym  (:,:) !< asymmetry factor
  REAL(wp), ALLOCATABLE         :: reff(:)
  REAL(wp), ALLOCATABLE         :: zlambda(:), zl(:)
  REAL(wp)                      :: reff_dummy, eps
  INTEGER                       :: jptotal, iradunit, iread, irad, iw
  LOGICAL                       :: lex

  IF (p_parallel_io) THEN
    iradunit = find_next_free_unit (51,100)

    ! read data for lookup table from file "rad_table"
    INQUIRE (file=fname, exist=lex)
    IF (.NOT.lex) THEN
      CALL finish('read_aero_volc_tables(mo_aero_volc_tab.f90)', &
                  'file '//fname//' does not exist')
    END IF

    OPEN(UNIT=iradunit,FILE=fname,FORM='FORMATTED',STATUS='OLD',ACTION='READ')

    ! total number of bands
    jptotal=nbndlw+jpsw
    ALLOCATE(zlambda(jptotal),zl(jptotal))
    ALLOCATE(nbandmap(jptotal))
    IF (jptotal.NE.30) THEN
      CALL finish('read_aero_volc_tables(mo_aero_volc_tab)', &
         'expects different number of wavelength bands')
    END IF
    zlambda=(/0.23_wp,     0.30_wp,     0.40_wp,     0.55_wp,     0.70_wp, &
              1.00_wp,     1.27_wp,     1.46_wp,     1.78_wp,     2.05_wp, &
              2.32_wp,     2.79_wp,     3.47_wp,     8.00_wp,     3.60_wp, &
              4.00_wp,     4.20_wp,     4.60_wp,     5.20_wp,     6.15_wp, &
              7.00_wp,     7.85_wp,     8.85_wp,     9.75_wp,    11.20_wp, &
             13.20_wp,    15.10_wp,    18.00_wp,    30.00_wp,   100.00_wp /)
    ! verify wavelength of bands, read and count effective radii first
    READ(iradunit,*)
    READ(iradunit,*) (zl(iw),iw=1,jptotal)
    DO iw=1,jptotal
      IF (ABS(zlambda(iw)-zl(iw)).GT.1._wp) THEN
        CALL finish('read_aero_volc_tables(mo_aero_volc_tab)', &
        'file '//fname//' does not contain the correct wavelength bands'// &
        ' or the order is wrong')
      END IF
    END DO
    nbandmap=(/13,          12,           11,          10,          9, &
                8,           7,            6,           5,          4, &
                3,           2,            1,          14,    jpsw+16, &
          jpsw+15,     jpsw+14,      jpsw+13,     jpsw+12,    jpsw+11, &
          jpsw+10,      jpsw+9,       jpsw+8,      jpsw+7,     jpsw+6, &
           jpsw+5,      jpsw+4,       jpsw+3,      jpsw+2,     jpsw+1 /)
    DEALLOCATE(zlambda,zl)
    nrad=0
    DO
      READ(iradunit,*,iostat=iread) reff_dummy
      IF(iread/=0) EXIT
      nrad=nrad+1
      READ(iradunit,*,iostat=iread)
      IF(iread/=0) CALL finish('read_aero_volc_tables(mo_aero_volc_tab)', &
                  'file '//fname//' is corrupt')
      READ(iradunit,*,iostat=iread)
      IF(iread/=0) CALL finish('read_aero_volc_tables(mo_aero_volc_tab)', &
                  'file '//fname//' is corrupt')
      READ(iradunit,*,iostat=iread)
      IF(iread/=0) CALL finish('read_aero_volc_tables(mo_aero_volc_tab)', &
                  'file '//fname//' is corrupt')
    END DO
  END IF
  CALL p_bcast(nrad,p_io)
  CALL p_bcast(jptotal,p_io)
  ALLOCATE(extrat(jptotal,nrad))
  ALLOCATE(ssa(jptotal,nrad))
  ALLOCATE(asym(jptotal,nrad))
  IF (p_parallel_io) THEN
    ALLOCATE(reff(nrad))
    ALLOCATE(zextrat(jptotal,nrad))
    ALLOCATE(zssa(jptotal,nrad))
    ALLOCATE(zasym(jptotal,nrad))
    REWIND(iradunit)
    READ(iradunit,*)
    READ(iradunit,*) 
    DO irad=1,nrad
      READ(iradunit,*,iostat=iread) reff(irad)
      IF(iread/=0) CALL finish('read_aero_volc_tables(mo_aero_volc_tab)', &
                  'file '//fname//' is corrupt')
      READ(iradunit,*,iostat=iread) (zextrat(iw,irad),iw=1,jptotal)
      IF(iread/=0) CALL finish('read_aero_volc_tables(mo_aero_volc_tab)', &
                  'file '//fname//' is corrupt')
      READ(iradunit,*,iostat=iread) (zssa(iw,irad),iw=1,jptotal)
      IF(iread/=0) CALL finish('read_aero_volc_tables(mo_aero_volc_tab)', &
                  'file '//fname//' is corrupt')
      READ(iradunit,*,iostat=iread) (zasym(iw,irad),iw=1,jptotal)
      IF(iread/=0) CALL finish('read_aero_volc_tables(mo_aero_volc_tab)', &
                  'file '//fname//' is corrupt')
    END DO
    CLOSE(iradunit)
    !< reorder 
    DO iw=1,jptotal
      extrat(iw,1:nrad)=zextrat(nbandmap(iw),1:nrad)
! In Kinne's tables, co-singlescattering albedos are given
      ssa(iw,1:nrad)=1._wp-zssa(nbandmap(iw),1:nrad)
      asym(iw,1:nrad)=zasym(nbandmap(iw),1:nrad)
    END DO
    DEALLOCATE(zextrat,zssa,zasym,nbandmap)
    !< verify that reff is equidistant reff(1),reff(2),... 
    reff_min=reff(1)
    IF (nrad == 1) THEN
      delta_reff=1._wp
      delta_reff_i=1._wp
    ELSE
      delta_reff=reff(2)-reff(1)
      IF (delta_reff.LE.0._wp) THEN
        CALL finish('read_aero_volc_tables(mo_aero_volc_tab)', &
                    'table of effective radii not ascending in file '// &
                    fname)
      END IF 
      delta_reff_i=1._wp/delta_reff
      eps=1.e-7_wp*delta_reff
      DO irad=2,nrad-1
        IF (ABS(reff(irad+1)-reff(irad)-delta_reff).GT.eps) THEN
          CALL finish('read_aero_volc_tables(mo_aero_volc_tab)', &
                      'table of effective radii not equidistant in file '// &
                      fname)
        END IF
      END DO
    END IF ! (nrad == 1)
    !< Transform units from um to m
    delta_reff_i=delta_reff_i*1.e6_wp
    reff_min=reff_min*1.e-6_wp
    DEALLOCATE(reff)
  END IF ! (p_parallel_io)
  CALL p_bcast(delta_reff_i,p_io)
  CALL p_bcast(reff_min,p_io)
  CALL p_bcast(extrat,p_io)
  CALL p_bcast(ssa,p_io)
  CALL p_bcast(asym,p_io)

END SUBROUTINE read_aero_volc_tables
!! @par Description:
!! Read optical properties of aerosols provided by echam-HAM simulations.
!! Properties are given at 550 nm and must be interpolated using the lookup 
!! tables to the other wavelengths of echam.
!! 
!! @par Revision History:
!! original source by J.S. Rast (2011-03-23)

SUBROUTINE read_aero_prop_ham(rad_step_1, rad_step_2)

  USE mo_time_conversion,    ONLY: time_days
  
  !INPUT PARMETES
  TYPE(time_days), INTENT(in)   :: rad_step_1, rad_step_2

  !LOCAL VARIABLES
  INTEGER                       :: icurrentyear, inextyear
  INTEGER                       :: kyrm1, kyr, kyrp1
  LOGICAL                       :: lnewyear
  CHARACTER(len=5)              :: cyr
  CHARACTER(len=18)             :: cfname

  CALL get_date_components (rad_step_1, year=icurrentyear)
  CALL get_date_components (rad_step_2, year=inextyear)
  lnewyear=icurrentyear/=inextyear

  IF (.NOT.lnewyear .AND. laero_set_ham) return
  kyrm1=inextyear-1
  kyr=inextyear
  kyrp1=inextyear+1
  WRITE(cyr,*) kyrm1
  cfname=TRIM(cfname_base_ham)//'_'//TRIM(ADJUSTL(cyr))//'.nc'
  CALL aero_ham_read ( &
    ldc%nproma,           ldc%ngpblks,             ldc%nlev,            &
    aodh(:,:,:,0:0),      'AOD',                  reffh(:,:,:,0:0),    &
    'REFF',               12,                      12,                  &
    cfname                                                              )
  WRITE(cyr,*) kyr
  cfname=TRIM(cfname_base_ham)//'_'//TRIM(ADJUSTL(cyr))//'.nc'
  CALL aero_ham_read ( &
    ldc%nproma,           ldc%ngpblks,             ldc%nlev,            &
    aodh(:,:,:,1:12),     'AOD',                  reffh(:,:,:,1:12),   &
    'REFF',               1,                       12,                  &
    cfname                                                              )
  WRITE(cyr,*) kyrp1
  cfname=TRIM(cfname_base_ham)//'_'//TRIM(ADJUSTL(cyr))//'.nc'
  CALL aero_ham_read ( &
    ldc%nproma,           ldc%ngpblks,             ldc%nlev,            &
    aodh(:,:,:,13:13),    'AOD',                  reffh(:,:,:,13:13),  &
    'REFF',               1,                       1,                   &
    cfname                                                              )
  laero_set_ham=.true.
END SUBROUTINE read_aero_prop_ham

SUBROUTINE read_aero_prop_crow(rad_step_1, rad_step_2)

  USE mo_time_conversion,    ONLY: time_days
  
  !INPUT PARMETES
  TYPE(time_days), INTENT(in)   :: rad_step_1, rad_step_2

  !LOCAL VARIABLES
  REAL(wp), ALLOCATABLE         :: zaod(:,:), zreff(:,:)
  INTEGER                       :: icurrentyear, & !< current year
                                   inextyear, & !< year of next rad. time step
                                   iunit, iread, jd, ilat
  LOGICAL                       :: lnewyear, lex
  CHARACTER(len=5)              :: cyr
  CHARACTER(len=16)             :: cfname

  CALL get_date_components (rad_step_1, year=icurrentyear)
  CALL get_date_components (rad_step_2, year=inextyear)
  lnewyear=icurrentyear/=inextyear
  IF (.NOT.lnewyear .AND. laero_set_crow) return

  cfname=TRIM(cfname_base_crow)//'.dat'
  WRITE(cyr,*) inextyear

  IF (p_parallel_io) THEN
    ! allocate fields for original data on nlat_crow latitudes
    ALLOCATE(zaod(nlat_crow,ntstep_crow+1))
    ALLOCATE(zreff(nlat_crow,ntstep_crow+1))
    iunit = find_next_free_unit (51,100)

    ! read data for lookup table from file "rad_table"
    INQUIRE (file=cfname, exist=lex)
    IF (.NOT.lex) THEN
      CALL finish('read_aero_prop_crow(mo_aero_volc_tab.f90)', &
                  'file '//cfname//' does not exist')
    END IF

    OPEN(UNIT=iunit,FILE=cfname,FORM='FORMATTED',STATUS='OLD',ACTION='READ')

    ! skip header lines
    READ(iunit,*)
    READ(iunit,*)

    DO 
      READ(iunit,*,iostat=iread) fyear_crow(1), &
         zaod(1,1), zreff(1,1), &
         zaod(2,1), zreff(2,1), &
         zaod(3,1), zreff(3,1), &
         zaod(4,1), zreff(4,1)
      IF(iread>0) CALL finish('read_aero_prop_crow(mo_aero_volc_tab)', &
                  'file '//cfname//' is corrupt')
      IF(iread<0) CALL finish('read_aero_prop_crow(mo_aero_volc_tab)', &
                  'file '//cfname//' does not contain data for year '//cyr)
      IF (fyear_crow(1)>=inextyear) EXIT
    END DO
    DO jd=2,ntstep_crow+1
      READ(iunit,*,iostat=iread) fyear_crow(jd), &
         zaod(1,jd), zreff(1,jd), &
         zaod(2,jd), zreff(2,jd), &
         zaod(3,jd), zreff(3,jd), &
         zaod(4,jd), zreff(4,jd)
      IF(iread>0) CALL finish('read_aero_prop_crow(mo_aero_volc_tab)', &
                  'file '//cfname//' is corrupt')
      IF(iread<0) CALL finish('read_aero_prop_crow(mo_aero_volc_tab)', &
                  'file '//cfname//' does not contain sufficient data '// &
                  'for year '//cyr)
    END DO
    CLOSE(iunit)

    ! linear interpolation of nlat_crow zonal bands to echam latitudes
    ! - north and south of equator, the edges are 15 and 45 degrees
    ! - north and south of 45 degrees: constant values

    DO jd=1,ntstep_crow+1
      DO ilat=1,ldc%nlat
        IF(philat(ilat)>=45._wp) THEN
          aod_crow(ilat,jd)=zaod(1,jd)
          reff_crow(ilat,jd)=zreff(1,jd)
        ELSE IF(philat(ilat)>=15._wp .AND. philat(ilat)<45._wp) THEN
          aod_crow(ilat,jd)=((philat(ilat)-15._wp)*zaod(1,jd) +     &
                            (45._wp-philat(ilat))*zaod(2,jd))/30._wp
          reff_crow(ilat,jd)=((philat(ilat)-15._wp)*zreff(1,jd) +     &
                            (45._wp-philat(ilat))*zreff(2,jd))/30._wp
        ELSE IF(philat(ilat)>=-15._wp .AND. philat(ilat)<15._wp) THEN
          aod_crow(ilat,jd)=((philat(ilat)+15._wp)*zaod(2,jd) +     &
                            (15._wp-philat(ilat))*zaod(3,jd))/30._wp
          reff_crow(ilat,jd)=((philat(ilat)+15._wp)*zreff(2,jd) +     &
                            (15._wp-philat(ilat))*zreff(3,jd))/30._wp
        ELSE IF(philat(ilat)>=-45._wp .AND. philat(ilat)<-15._wp) THEN
          aod_crow(ilat,jd)=((philat(ilat)+45._wp)*zaod(3,jd) +     &
                            (-15._wp-philat(ilat))*zaod(4,jd))/30._wp
          reff_crow(ilat,jd)=((philat(ilat)+45._wp)*zreff(3,jd) +     &
                            (-15._wp-philat(ilat))*zreff(4,jd))/30._wp
        ELSE IF(philat(ilat)<-45._wp) THEN
          aod_crow(ilat,jd)=zaod(4,jd)
          reff_crow(ilat,jd)=zreff(4,jd)
        END IF
      END DO
    END DO

    ! set values without linear interpolation as alternativ method

!!$    DO jd=1,ntstep_crow+1
!!$      DO ilat=1,ldc%nlat
!!$        IF(philat(ilat)>=30._wp) THEN
!!$          aod_crow(ilat,jd)=zaod(1,jd)
!!$          reff_crow(ilat,jd)=zreff(1,jd)
!!$        ELSE IF(philat(ilat)>=0._wp .AND. philat(ilat)<30._wp) THEN
!!$          aod_crow(ilat,jd)=zaod(2,jd)
!!$          reff_crow(ilat,jd)=zreff(2,jd)
!!$        ELSE IF(philat(ilat)>=-30._wp .AND. philat(ilat)<0._wp) THEN
!!$          aod_crow(ilat,jd)=zaod(3,jd)
!!$          reff_crow(ilat,jd)=zreff(3,jd)
!!$        ELSE IF(philat(ilat)<-30._wp) THEN
!!$          aod_crow(ilat,jd)=zaod(4,jd)
!!$          reff_crow(ilat,jd)=zreff(4,jd)
!!$        END IF
!!$      END DO
!!$    END DO
! transform effective radii into from um into m
   reff_crow=1.e-6_wp*reff_crow
   DEALLOCATE(zaod,zreff)

   END IF

   CALL p_bcast(fyear_crow,p_io)
   CALL p_bcast(aod_crow,p_io)
   CALL p_bcast(reff_crow,p_io)
   laero_set_crow=.true.
END SUBROUTINE read_aero_prop_crow  

!! @par Description:
!! Add optical properties of volcanic aerosols provided from echam-HAM 
!! simulations to the optical properties of backgroud (e.g. Kinne) aerosols
!!
!! @par Revision History:
!! original source by J.S. Rast (2011-03-24)

SUBROUTINE add_aop_volc_ham( &
           & kproma           ,kbdim                 ,klev             ,&
           & krow             ,nb_lw                 ,nb_sw            ,&
           & paer_tau_lw_vr   ,paer_tau_sw_vr        ,paer_piz_sw_vr   ,&
           & paer_cg_sw_vr                                              )

! !INPUT PARAMETERS:
  INTEGER, INTENT(in)    :: kproma, kbdim, klev, krow
  INTEGER, INTENT(in)    :: nb_lw, nb_sw !< number of wavelengths for IR, VIS
  
! !OUTPUT PARAMETERS:
  REAL(wp),INTENT(inout),DIMENSION(kbdim,klev,nb_lw):: &
   paer_tau_lw_vr        !<aerosol optical depth (far IR)
  REAL(wp),INTENT(inout),DIMENSION(kbdim,klev,nb_sw):: &
   paer_tau_sw_vr,   &   !<aerosol optical depth (solar), sum_i(tau_i)
   paer_piz_sw_vr,   &   !<weighted sum of single scattering albedos, 
                         !<sum_i(tau_i*omega_i)
   paer_cg_sw_vr         !<weighted sum of asymmetry factors, 
                         !<sum_i(tau_i*omega_i*g_i)

! !LOCAL VARIABLES
  INTEGER               :: iw,jk,jki,jl,isum
  REAL(wp)              :: zaod(kbdim,klev), zreff(kbdim,klev)
  REAL(wp)              :: zaod_sw(kbdim,klev,nb_sw)
  REAL(wp)              :: zaod_lw(kbdim,klev,nb_lw)
  REAL(wp)              :: zssa_sw(kbdim,klev,nb_sw)
  REAL(wp)              :: zssa_lw(kbdim,klev,nb_lw)
  REAL(wp)              :: zasy_sw(kbdim,klev,nb_sw)
  INTEGER               :: ireff(kbdim,klev),icount(kbdim,klev)

! 1. Time interpolation of extinction and effective radius
  DO jk=1,klev
    DO jl=1,kproma
      zaod(jl,jk)=wgt1_m*aodh(jl,jk,krow,nmw1_m)+ &
                  wgt2_m*aodh(jl,jk,krow,nmw2_m)
      zreff(jl,jk)=wgt1_m*reffh(jl,jk,krow,nmw1_m)+ &
                   wgt2_m*reffh(jl,jk,krow,nmw2_m)
    END DO
  END DO
! 2. calculation of wavelength dependent aod
  ireff(1:kproma,1:klev)=INT((zreff(1:kproma,1:klev)-reff_min)*delta_reff_i)+1
  icount=0
  WHERE (ireff(1:kproma,1:klev).LT.1)
    ireff(1:kproma,1:klev)=1
    icount(1:kproma,1:klev)=1
  END WHERE
  isum=SUM(icount)
  IF(isum.GT.0) THEN
    CALL message('add_aop_volc_ham(mo_aero_volc_tab.f90)', &
                 'effective radius smaller than minium eff. radius in table')
  END IF
  icount=0
  WHERE (ireff(1:kproma,1:klev).GT.nrad)
    ireff(1:kproma,1:klev)=nrad
    icount(1:kproma,1:klev)=1
  END WHERE
  isum=SUM(icount)
  IF(isum.GT.0) THEN
    CALL message('add_aop_volc_ham(mo_aero_volc_tab.f90)', &
                 'effective radius larger than maximum eff. radius in table')
  END IF
  DO iw=1,nb_sw
    DO jk=1,klev
      DO jl=1,kproma
        zaod_sw(jl,jk,iw)=zaod(jl,jk)*extrat(iw,ireff(jl,jk))
        zssa_sw(jl,jk,iw)=ssa(iw,ireff(jl,jk))
        zasy_sw(jl,jk,iw)=asym(iw,ireff(jl,jk))
      END DO
    END DO
  END DO

  DO jk=1,klev
    jki=klev-jk+1
    WHERE (zaod_sw(1:kproma,jki,1:nb_sw)>0._wp) 
      paer_cg_sw_vr(1:kproma,jk,1:nb_sw)=paer_tau_sw_vr(1:kproma,jk,1:nb_sw)*&
       paer_piz_sw_vr(1:kproma,jk,1:nb_sw)*paer_cg_sw_vr(1:kproma,jk,1:nb_sw)+&
       zaod_sw(1:kproma,jki,1:nb_sw)*zssa_sw(1:kproma,jki,1:nb_sw)*&
       zasy_sw(1:kproma,jki,1:nb_sw)
      paer_piz_sw_vr(1:kproma,jk,1:nb_sw)=paer_tau_sw_vr(1:kproma,jk,1:nb_sw)*&
       paer_piz_sw_vr(1:kproma,jk,1:nb_sw)+&
       zaod_sw(1:kproma,jki,1:nb_sw)*zssa_sw(1:kproma,jki,1:nb_sw)
      paer_tau_sw_vr(1:kproma,jk,1:nb_sw)=paer_tau_sw_vr(1:kproma,jk,1:nb_sw)+&
       zaod_sw(1:kproma,jki,1:nb_sw)
      paer_piz_sw_vr(1:kproma,jk,1:nb_sw)=paer_piz_sw_vr(1:kproma,jk,1:nb_sw)/&
       paer_tau_sw_vr(1:kproma,jk,1:nb_sw)
      paer_cg_sw_vr(1:kproma,jk,1:nb_sw)=paer_cg_sw_vr(1:kproma,jk,1:nb_sw)/&
       (paer_tau_sw_vr(1:kproma,jk,1:nb_sw)* &
       paer_piz_sw_vr(1:kproma,jk,1:nb_sw))
    END WHERE
  END DO
! 4. Optical properties at thermal wavelengths
  DO iw=1,nb_lw
    DO jk=1,klev
      DO jl=1,kproma
        zaod_lw(jl,jk,iw)=zaod(jl,jk)*extrat(iw+nb_sw,ireff(jl,jk))
        zssa_lw(jl,jk,iw)=ssa(iw+nb_sw,ireff(jl,jk))
      END DO
    END DO
  END DO

  DO jk=1,klev
    jki=klev-jk+1
    paer_tau_lw_vr(1:kproma,jk,1:nb_lw)=paer_tau_lw_vr(1:kproma,jk,1:nb_lw)+ &
          zaod_lw(1:kproma,jki,1:nb_lw)*(1._wp-zssa_lw(1:kproma,jki,1:nb_lw))
  END DO      
END SUBROUTINE add_aop_volc_ham
!! @par Description:
!! Add optical properties of volcanic aerosols provided by T. Crowley
!!
!! @par Revision History:
!! original source by J.S. Rast (2011-04-08)

SUBROUTINE add_aop_volc_crow( &
           & kproma           ,kbdim                 ,klev             ,&
           & krow             ,nb_lw                 ,nb_sw            ,&
           & paer_tau_lw_vr   ,paer_tau_sw_vr        ,paer_piz_sw_vr   ,&
           & paer_cg_sw_vr                                              )

! !INPUT PARAMETERS:
  INTEGER, INTENT(in)    :: kproma, kbdim, klev, krow
  INTEGER, INTENT(in)    :: nb_lw, nb_sw !< number of wavelengths for IR, VIS
  
! !OUTPUT PARAMETERS:
  REAL(wp),INTENT(inout),DIMENSION(kbdim,klev,nb_lw):: &
   paer_tau_lw_vr        !<aerosol optical depth (far IR)
  REAL(wp),INTENT(inout),DIMENSION(kbdim,klev,nb_sw):: &
   paer_tau_sw_vr,   &   !<aerosol optical depth (solar), sum_i(tau_i)
   paer_piz_sw_vr,   &   !<weighted sum of single scattering albedos, 
                         !<sum_i(tau_i*omega_i)
   paer_cg_sw_vr         !<weighted sum of asymmetry factors, 
                         !<sum_i(tau_i*omega_i*g_i)

! !LOCAL VARIABLES
  INTEGER               :: ireff(kbdim)
  INTEGER               :: iyear, ihour, imint, isec, idayyr, &
                           jl, jk, jki, iw
  INTEGER               :: itw_1, itw_2
  REAL(wp)              :: zaod(kbdim),zreff(kbdim)
  REAL(wp)              :: rday_len, rdaysyr, fyear
  REAL(wp)              :: zaod_sw(kbdim,klev,nb_sw)
  REAL(wp)              :: zaod_lw(kbdim,klev,nb_lw)
  REAL(wp)              :: zssa_sw(kbdim,klev,nb_sw)
  REAL(wp)              :: zssa_lw(kbdim,klev,nb_lw)
  REAL(wp)              :: zasy_sw(kbdim,klev,nb_sw)
  REAL(wp)              :: rtw_1, rtw_2


! Calculate weights for time interpolation

  CALL get_date_components(radiation_date, year=iyear, hour=ihour, &
                           minute=imint, second=isec)
  idayyr=day_in_year(radiation_date)-1
  rdaysyr=1._wp/REAL(year_len(radiation_date),wp)
  rday_len=1._wp/REAL(day_len(),wp)
  isec=3600*ihour+60*imint+isec
  fyear=REAL(iyear,wp)+(REAL(idayyr+rday_len*REAL(isec,wp),wp))*rdaysyr
  itw_1=FLOOR((fyear-fyear_crow(1))/dt_crow)+1
  itw_2=itw_1+1
  rtw_1=(fyear_crow(itw_2)-fyear)/dt_crow
  rtw_2=1._wp-rtw_1

  DO jl=1,kproma  
    zreff(jl)=reff_crow(ilat(jl,krow),itw_1)*rtw_1+ &
             reff_crow(ilat(jl,krow),itw_2)*rtw_2
    zaod(jl)=aod_crow(ilat(jl,krow),itw_1)*rtw_1+ &
             aod_crow(ilat(jl,krow),itw_2)*rtw_2
  END DO
  ireff(1:kproma)=INT((zreff(1:kproma)-reff_min)*delta_reff_i)+1
! Optical properties at solar wave lengths
  DO iw=1,nb_sw
    DO jk=1,klev
      DO jl=1,kproma
        zaod_sw(jl,jk,iw)=zaod(jl)*pl_weights(jk)*extrat(iw,ireff(jl))
        zssa_sw(jl,jk,iw)=ssa(iw,ireff(jl))
        zasy_sw(jl,jk,iw)=asym(iw,ireff(jl))
      END DO
    END DO
  END DO

  DO jk=1,klev
    jki=klev-jk+1
    WHERE (zaod_sw(1:kproma,jki,1:nb_sw)>0._wp) 
      paer_cg_sw_vr(1:kproma,jk,1:nb_sw)=paer_tau_sw_vr(1:kproma,jk,1:nb_sw)*&
      paer_piz_sw_vr(1:kproma,jk,1:nb_sw)*paer_cg_sw_vr(1:kproma,jk,1:nb_sw)+&
        zaod_sw(1:kproma,jki,1:nb_sw)*zssa_sw(1:kproma,jki,1:nb_sw)*&
        zasy_sw(1:kproma,jki,1:nb_sw)
      paer_piz_sw_vr(1:kproma,jk,1:nb_sw)=paer_tau_sw_vr(1:kproma,jk,1:nb_sw)*&
      paer_piz_sw_vr(1:kproma,jk,1:nb_sw)+&
        zaod_sw(1:kproma,jki,1:nb_sw)*zssa_sw(1:kproma,jki,1:nb_sw)
      paer_tau_sw_vr(1:kproma,jk,1:nb_sw)=paer_tau_sw_vr(1:kproma,jk,1:nb_sw)+&
        zaod_sw(1:kproma,jki,1:nb_sw)
      paer_piz_sw_vr(1:kproma,jk,1:nb_sw)=paer_piz_sw_vr(1:kproma,jk,1:nb_sw)/&
        paer_tau_sw_vr(1:kproma,jk,1:nb_sw)
      paer_cg_sw_vr(1:kproma,jk,1:nb_sw)=paer_cg_sw_vr(1:kproma,jk,1:nb_sw)/&
        (paer_tau_sw_vr(1:kproma,jk,1:nb_sw)* &
         paer_piz_sw_vr(1:kproma,jk,1:nb_sw))
    END WHERE
  END DO
! Optical properties at thermal wavelengths
  DO iw=1,nb_lw
    DO jk=1,klev
      DO jl=1,kproma
        zaod_lw(jl,jk,iw)=zaod(jl)*pl_weights(jk)*extrat(iw+nb_sw,ireff(jl))
        zssa_lw(jl,jk,iw)=ssa(iw+nb_sw,ireff(jl))
      END DO
    END DO
  END DO

  DO jk=1,klev
    jki=klev-jk+1
    paer_tau_lw_vr(1:kproma,jk,1:nb_lw)=paer_tau_lw_vr(1:kproma,jk,1:nb_lw)+ &
          zaod_lw(1:kproma,jki,1:nb_lw)*(1._wp-zssa_lw(1:kproma,jki,1:nb_lw))
  END DO      
END SUBROUTINE add_aop_volc_crow

!! @par Description:
!! Deallocate fields for internal rerun
!! 
!! @par Revision History:
!! original source by J.S. Rast (2011-04-04)
SUBROUTINE cleanup_aero_volc_tab_ham

  IF(ALLOCATED(aodh)) DEALLOCATE(aodh)
  IF(ALLOCATED(reffh)) DEALLOCATE(reffh)
  IF(ALLOCATED(extrat)) DEALLOCATE(extrat)
  IF(ALLOCATED(ssa)) DEALLOCATE(ssa)
  IF(ALLOCATED(asym)) DEALLOCATE(asym)

  laero_set_ham=.false.

END SUBROUTINE cleanup_aero_volc_tab_ham

!! @par Description:
!! Deallocate fields for internal rerun
!! 
!! @par Revision History:
!! original source by J.S. Rast (2011-05-30)
SUBROUTINE cleanup_aero_volc_tab_crow

  IF(ALLOCATED(fyear_crow))  DEALLOCATE(fyear_crow)
  IF(ALLOCATED(aod_crow))    DEALLOCATE(aod_crow)
  IF(ALLOCATED(reff_crow))   DEALLOCATE(reff_crow)
  IF(ALLOCATED(pl_weights))  DEALLOCATE(pl_weights)
  IF(ALLOCATED(extrat))      DEALLOCATE(extrat)
  IF(ALLOCATED(ssa))         DEALLOCATE(ssa)
  IF(ALLOCATED(asym))        DEALLOCATE(asym)

  laero_set_crow=.false.

END SUBROUTINE cleanup_aero_volc_tab_crow


SUBROUTINE aero_ham_read (&
  kbdim,           ngpblks,            klev,           &
  ext,             cext,               reff,           &
  creff,           imnthb,             imnthe,         &
  cfname                                               )

  !INPUT/OUTPUT PARAMETERS
  INTEGER, INTENT(in) :: kbdim,  & !> maximum block length
                         ngpblks,& !> number of blocks
                         klev,   & !> number of model levels
                         imnthb, & !> first month to read
                         imnthe    !> last month to read
  CHARACTER(len=*), INTENT(in)    :: cext, &!> name of extincetion var in file
                                     creff,&!> name of effective radius var
                                     cfname !> name of netcdf file
  REAL(wp), INTENT(out) :: ext(kbdim,klev,ngpblks,imnthb:imnthe),&!< extinction
                           reff(kbdim,klev,ngpblks,imnthb:imnthe) !< eff. rad.
  !LOCAL VARIABLES

  LOGICAL                :: lex
  INTEGER                :: j,ierr,ilon
  REAL(wp), POINTER          :: field_3dm(:,:,:)
  REAL(wp), POINTER          :: field_3d(:,:,:)

  IF (p_parallel_io) THEN
    INQUIRE (file=TRIM(cfname), exist=lex)
    IF (.NOT. lex) THEN
      CALL finish('read_aero_prop_ham(mo_aero_volc_tab)', &
                  'file '//TRIM(cfname)//' does not exist')
    END IF
    ALLOCATE(field_3dm(1,ldc%nlev,ldc%nlat))
    ALLOCATE(field_3d(ldc%nlon,ldc%nlev,ldc%nlat))
  END IF
  DO j=imnthb,imnthe
    IF (p_parallel_io) THEN
      CALL read_var_hs_nf77_3d(TRIM(cfname),       'lon',      'mlev',     &
                               'lat',              'time',     j,          &
                               TRIM(cext),         field_3dm,  ierr        )
      DO ilon=1,ldc%nlon
        field_3d(ilon,1:ldc%nlev,1:ldc%nlat)=field_3dm(1,1:ldc%nlev,1:ldc%nlat)
      END DO
    END IF
    CALL scatter_gp(field_3d,ext(:,:,:,j),global_decomposition)
    IF (p_parallel_io) THEN
      CALL read_var_hs_nf77_3d(TRIM(cfname),       'lon',      'mlev',     &
                               'lat',              'time',     j,          &
                               TRIM(creff),        field_3dm,  ierr        )
      DO ilon=1,ldc%nlon
        field_3d(ilon,1:ldc%nlev,1:ldc%nlat)=field_3dm(1,1:ldc%nlev,1:ldc%nlat)
      END DO
    END IF
    CALL scatter_gp(field_3d,reff(:,:,:,j),global_decomposition)
  END DO
  IF (p_parallel_io) THEN
    DEALLOCATE(field_3d,field_3dm)
  END IF
END SUBROUTINE aero_ham_read

END MODULE mo_aero_volc_tab
