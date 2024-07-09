!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_volc_data

  !
  !  read and initialize data for volcanic aerosols
  !
  !   U.Schlese, MPI Hamburg, June-2007
  !   M.Esch,    MPI Hamburg, June-2007
  !   S.Lorenz,  MPI Hamburg, Nov-2007, Feb-2008
  !

  USE mo_kind,         ONLY: dp
  USE mo_filename,     ONLY: find_next_free_unit
  USE mo_control,      ONLY: ngl

  IMPLICIT NONE

  SAVE

  PRIVATE

  ! public variables 
  !
  ! quantities for 50 effective radii (2-100nm) (every 2nm)
  !            and 22 wavelength intervals (SW:1-6, LW:7-22)

  PUBLIC :: jprad1    ! lookup tables lower bound for eff. radii
  PUBLIC :: jprad2    ! lookup tables upper bound for eff. radii
  PUBLIC :: jpwave    ! number of wavelength intervals
  PUBLIC :: extrat    ! extinction ratio
  PUBLIC :: asym      ! asymmetry factor
  PUBLIC :: ssa       ! single scattering albedo

  ! volcanic aerosol data from T.Crowley

  PUBLIC :: jpz       ! number of zonal bands (4: 90-30N,30-0N,0-30S,30-90S)
  PUBLIC :: jpd       ! number of aerosol data sets (108, 36 per year)

  ! on zonal model grid:

  PUBLIC :: aod, reff         ! interpolated zonal forcing data

  ! public subroutines

  PUBLIC :: init_volc_tables  ! initialize lookup tables
  PUBLIC :: read_volc_data    ! read volcanic aerosol data
  PUBLIC :: interp_volc       ! time interpolate volcanic aerosol data

  INTEGER, PARAMETER :: jprad1 =  1   ! reff/2
  INTEGER, PARAMETER :: jprad2 =  50  ! reff/2
  INTEGER, PARAMETER :: jpwave =  22  ! wavenumbers (SW+LW)
  INTEGER, PARAMETER :: jpz    =   4  ! zonal bands
  INTEGER, PARAMETER :: jpd    =  36  ! data records of one year
  

  REAL(dp) :: extrat(jprad1:jprad2,jpwave)
  REAL(dp) :: ssa   (jprad1:jprad2,jpwave)
  REAL(dp) :: asym  (jprad1:jprad2,jpwave)

  REAL(dp), ALLOCATABLE  :: aod(:,:), reff(:,:)  !  (ngl,0:37)
! REAL(dp) :: aod(ngl,0:jpd+1), reff(ngl,0:jpd+1)
 
!-----------------------------------------------------------------------

CONTAINS

  SUBROUTINE init_volc_tables
  ! 
  ! called from setphys <- initialize <- control <- master
  ! intialize lookup tables for aerosol radiation parameters
  !

    INTEGER jr, jw, iradunit
    LOGICAL lex

      iradunit = find_next_free_unit (51,100)
 
    ! read data for lookup table from file "rad_table"
    INQUIRE (file='rad_table', exist=lex)
    OPEN(UNIT=iradunit,FILE='rad_table',FORM='FORMATTED',STATUS='OLD',ACTION='READ')

    READ(iradunit,*) ! skip header
    DO jw=1,jpwave
      READ (iradunit,*) (extrat(jr,jw),jr=1,50)
    END DO

    READ(iradunit,*)
    READ(iradunit,*)
    READ(iradunit,*)
    READ(iradunit,*)
  

    DO jw=1,jpwave
      READ (iradunit,*) (ssa(jr,jw),jr=1,50)
    END DO 

    READ(iradunit,*)
    READ(iradunit,*)
    READ(iradunit,*)
    READ(iradunit,*)


    DO jw=1,jpwave
      READ (iradunit,*) (asym(jr,jw),jr=1,50)
    END DO

    CLOSE(iradunit)

  END SUBROUTINE init_volc_tables

  SUBROUTINE read_volc_data

   !  reads  the Crowley aerosol data 

   !  called from setphys
   ! must be read when initial forecast time is already known !!


    USE mo_time_control,    ONLY: next_date, get_date_components
    USE mo_gaussgrid,       ONLY: philat

    INTEGER   :: iyr, jd, ivolcunit, jl
    REAL(dp)  :: zaod  (jpz,0:jpd+1)       ! local aerosol optical depth
    REAL(dp)  :: zreff (jpz,0:jpd+1)       ! local effective radius
    REAL(dp)  :: yearvolc  (0:jpd+1)       ! year of the volcano

    CALL get_date_components(next_date, iyr)

    ivolcunit = find_next_free_unit (51,100)

    OPEN(UNIT=ivolcunit,FILE='volc_data',FORM='FORMATTED',STATUS='OLD',ACTION='READ')

!   search for first record of the year before current forecast year
!   
!    skip headers

     READ(ivolcunit,*)
     READ(ivolcunit,*)

  1  READ(ivolcunit,*) yearvolc(1), zaod(1,1), zreff(1,1), zaod(2,1), zreff(2,1), &
                     zaod(3,1), zreff(3,1), zaod(4,1), zreff(4,1)

    IF(INT(yearvolc(1)).NE.iyr-1) GO TO 1

!   read one year plus one record before and after resp.
!   
!   go to the last record of the year  before current year

    DO jd=2,jpd-1
      READ(ivolcunit,*)
    END DO

!   read data

    DO jd=0,jpd+1
      READ(ivolcunit,*) yearvolc(jd), zaod(1,jd), zreff(1,jd), zaod(2,jd), zreff(2,jd), &
                      zaod(3,jd), zreff(3,jd), zaod(4,jd), zreff(4,jd)
    END DO

    WRITE(*,*) ' Year of first record of volc_data: ', yearvolc(1)

!   linear interpolation between jpz=4 zonal bands to latitudes:
!    - north and south of equator edges are: 15 and 45 degrees
!    - north and south of 45 deg: constant values

    DO jd=0,jpd+1
      DO jl=1,ngl
        IF(philat(jl).GE.45._dp) THEN
          aod(jl,jd)=zaod(1,jd)
          reff(jl,jd)=zreff(1,jd)
        ELSE IF(philat(jl).GE.15._dp .AND. philat(jl).LT.45._dp) THEN
          aod (jl,jd)=( ( philat(jl)-15._dp)* zaod(1,jd) +              &
                        ( 45._dp-philat(jl))* zaod(2,jd) )/30._dp
          reff(jl,jd)=( ( philat(jl)-15._dp)*zreff(1,jd) +              &
                        ( 45._dp-philat(jl))*zreff(2,jd) )/30._dp
        ELSE IF(philat(jl).GE.-15._dp .AND. philat(jl).LT.15._dp) THEN
          aod (jl,jd)=( ( philat(jl)+15._dp)* zaod(2,jd) +              &
                        ( 15._dp-philat(jl))* zaod(3,jd) )/30._dp
          reff(jl,jd)=( ( philat(jl)+15._dp)*zreff(2,jd) +              &
                        ( 15._dp-philat(jl))*zreff(3,jd) )/30._dp
        ELSE IF(philat(jl).GE.-45._dp .AND. philat(jl).LT.-15._dp) THEN
          aod (jl,jd)=( ( philat(jl)+45._dp)* zaod(3,jd) +              &
                        (-15._dp-philat(jl))* zaod(4,jd) )/30._dp
          reff(jl,jd)=( ( philat(jl)+45._dp)*zreff(3,jd) +              &
                        (-15._dp-philat(jl))*zreff(4,jd) )/30._dp
        ELSE IF(philat(jl).LT.-45._dp) THEN
          aod(jl,jd)=zaod(4,jd)
          reff(jl,jd)=zreff(4,jd)
        END IF
      END DO
    END DO
    WRITE(*,*) ' Interpolated AOD data of last (jpd=36) record in %:'
    WRITE(*,'(8f12.5)') aod(:,jpd)*100._dp
      
  END SUBROUTINE read_volc_data
 
  SUBROUTINE interp_volc(kproma, kbdim, klev, krow, pswext,pswssa,pswasy,pext)

   !  called from rad_int
   ! interpolate in time and distribute to 2-D field
   !

    USE mo_time_control,    ONLY: next_date
    USE mo_time_conversion, ONLY: day_in_year
    USE mo_geoloc,          ONLY: ilat
    USE mo_memory_g1a,      ONLY: tm1, qm1, xlm1, xim1
    USE mo_memory_g3b,      ONLY: geosp
    USE mo_scan_buffer,     ONLY: alnpr, alpha
    USE mo_physical_constants,       ONLY: vtmpc1

    ! INPUT
    ! -----
 
    ! local dimensions
    INTEGER, INTENT(in)                      :: kproma ! number of local longitudes
    INTEGER, INTENT(in)                      :: kbdim  ! first dimension of 2-d arrays
    INTEGER, INTENT(in)                      :: klev   ! 

    ! gauss grid description
    INTEGER, INTENT(in)                      :: krow   ! sequential index

    ! output variables for volcanic aerosols
    REAL(dp), INTENT(out), DIMENSION(kbdim,klev,6) :: pswext   ! extinction SW
    REAL(dp), INTENT(out), DIMENSION(kbdim,klev,6) :: pswasy   ! asymmetry factor SW
    REAL(dp), INTENT(out), DIMENSION(kbdim,klev,6) :: pswssa   ! single scattering albedo SW
    REAL(dp), INTENT(out), DIMENSION(kbdim,klev,16) :: pext    ! extinction for longwave

!   local variables
    REAL(dp) :: zext_z(kbdim,jpwave), zasym_z(kbdim,jpwave), zssa_z(kbdim,jpwave)
    REAL(dp) :: zwgt1, zwgt2
    REAL(dp) :: ztv(kbdim,klev), zgeop(kbdim,klev), zalpha(kbdim,klev)
    REAL(dp) :: zww(klev), ztotal(kbdim), zgeopf(kbdim,klev)
    REAL(dp) :: zwthick(kbdim,klev)

    INTEGER   :: idayv       ! day in year
    INTEGER   :: iivolc      ! index in volcanic data
    INTEGER   :: jl, jk, jw, ireff, iizon, iday10 
    INTEGER   :: itop, ibot  ! top/bottom layers for aerosol distribution

!
!   Vertical distribution profile for volcanic aerosol forcing
!   Claudia Timmreck, Stephan Lorenz Feb 21, 2008
!

!   decide vertical distribution profile
!     pin10 - Nov 2007 - top 5 levels 1:2:4:2:1
!     pin9  - Feb 2008 - uppermost level not used, levels 2-4: 1:2:1
!     levels 2-4 are effected - caution, fits to 19 levels only

    itop=2 !!!!!!!
    ibot=4 !!!!!!!

    zww(2)=1.00_dp
    zww(3)=2.00_dp
    zww(4)=1.00_dp

!   assumes 36 steps of forcing data per year
!   inexact at end of year (365/366 days, not 360)

    idayv=day_in_year(next_date)
    iivolc=(idayv-1)/10+1
    IF(iivolc.EQ.37) iivolc=36

    iday10=MOD(idayv-1,10)
    zwgt2=REAL(iday10,dp)*0.1_dp
    zwgt1=1.0_dp-zwgt2
 
  DO jw = 1,jpwave
    DO jl = 1, kproma
      iizon=ilat(jl,krow)

! divide reff  by 2 to project the radii to 50 indices 
      ireff = NINT((zwgt1*reff(iizon,iivolc)+zwgt2*reff(iizon,iivolc+1))*100._dp*0.5_dp)
! catch 
      ireff = MIN(ireff,jprad2)
      ireff = MAX(ireff,jprad1)
!
      zext_z(jl,jw)=(zwgt1*aod(iizon,iivolc)+zwgt2*aod(iizon,iivolc+1))*extrat(ireff,jw)
      zasym_z(jl,jw)=asym(ireff,jw)
      zssa_z(jl,jw)=ssa(ireff,jw)
    END DO
  END DO

! put data into levels itop to ibot, weight with layer thickness in z

! virtual temperature:
    ztv(1:kproma,:) = tm1(1:kproma,:,krow)*(1._dp+vtmpc1*qm1(1:kproma,:,krow) &
                      -(xlm1(1:kproma,:,krow)+xim1(1:kproma,:,krow)))

! geopotential at half levels: 
    zalpha(:,:)=0.0_dp
    CALL geopot(zgeop,ztv,alnpr(:,:,krow),zalpha,geosp(:,krow),kbdim,kproma)

! geopotential at full levels (needed for uppermost level):

    CALL geopot(zgeopf,ztv,alnpr(:,:,krow),alpha(:,:,krow),geosp(:,krow),kbdim,kproma)

! total thickness of layers containing aerosol (weighted with profile)

   ztotal(:)=0._dp

   DO jk = itop,ibot
     DO jl = 1,kproma
       ztotal(jl) = ztotal(jl)+(zgeop(jl,jk-1)-zgeop(jl,jk))*zww(jk)
     END DO
   END DO

  DO jk = itop,ibot 
    DO jl = 1,kproma
      zwthick(jl,jk) = (zgeop(jl,jk-1)-zgeop(jl,jk))*zww(jk)/ztotal(jl)
    END DO
  END DO

! copy to 3d array to be used in the radiation subroutines
! !!! invert order of longwave data with respect to wavelength index  

  DO jw=1,6
    DO jk = itop,ibot  
      DO jl = 1, kproma
        pswext(jl,jk,jw) = zext_z(jl,jw)*zwthick(jl,jk)
        pswssa(jl,jk,jw) = zssa_z(jl,jw)
        pswasy(jl,jk,jw) = zasym_z(jl,jw)
      END DO
    END DO
  END DO

  DO jw=1,16             
    DO jk = itop,ibot 
      DO jl = 1, kproma
        pext(jl,jk,16+1-jw) = zext_z(jl,jw+6)*zwthick(jl,jk)
      END DO
    END DO
  END DO

  END SUBROUTINE interp_volc

END MODULE mo_volc_data
