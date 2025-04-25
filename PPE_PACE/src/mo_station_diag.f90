!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_station_diag
!----------------------------------------------------
! PURPOSE:
!
!   CFMIP5, cfSites in K. Taylor's table
!
!   module for writing high frequency point and profile output 
!   at individual sites
!  
! Strategy: write one output file for each site
!    
!           instantaneous counterparts of accumulated point variables
!           must be stored in additional fields
!
!           initialize pointers to 2-D and 3-D stream elements
!           to avoid calls to get_stream_element
! 
!           where 3-D variables are not available in output streams, 
!           collect profiles from 2-D slices in physc.f90  
!
!
! Entry points:
!     - station_diag_nml is called from initialize.f90 (reads namelist)
!     - station_diag is called from scan1.f90 
!     - collect_station_diag is called from physc.f90
!
!
! Other:
!   init_station_diag (called once from station_diag) calls:
!            - init_svars  
!                  (initialize info on variables)                
!            - init_pointers
!                   (initialize pointers to stream elements)        
!            - station_find 
!                  (determine PE and coordinates (row, proma) of stations)
!            - open_station_file (open output file, define variables, etc.)
!            - write_lon_lat ( write coordinates to output file)               
!                    
!---------------------------------------------------------
   USE mo_linked_list,   ONLY: t_stream
   USE mo_kind,          ONLY: wp

   IMPLICIT NONE

   PRIVATE

! do cfSites output?
   LOGICAL :: lostation = .FALSE. 

!-------------------------------------------------
! free IO/unit ?
   INTEGER, PARAMETER :: iun=2123
!
!-------------------------------------------------
! count number of calls, reset on restart
   INTEGER :: ncall 
!
!-------------------------------------------------
! number of stations to process (from pointlocation.txt)
   INTEGER, PARAMETER :: nstation = 124
!
!-------------------------------------------------
!  number of point variables
   INTEGER, PARAMETER :: nvar_0d = 32
!
!-------------------------------------------------
!  maximum number of profile variables from streams 16+ten
   INTEGER, PARAMETER :: nvar_1d_max = 30
!  
!  actual number of profile variables 
!  ( lower than nvar_1d_max only if not all streams/variables
!    are available)
   INTEGER :: nvar_1d
!  number of profile variables minus number of tendency variables 
    INTEGER :: nvar_1d_s

!
!-------------------------------------------------
!  number of additional profile variables from collect_station_diag
   INTEGER, PARAMETER :: nvar_1dp = 2 
!
!-------------------------------------------------
! has this module been initialized?
   LOGICAL :: module_is_initialized = .FALSE.

!-------------------------------------------------

  TYPE file_type 
    CHARACTER(LEN=150) :: name
   ! file handle
    INTEGER :: ncid
   !axis ids
    INTEGER :: time_dimid
    INTEGER :: mlev_dimid, ilev_dimid
   ! special variable ids
    INTEGER :: lon_varid, lat_varid 
    INTEGER :: time_varid 
  END TYPE file_type 

  TYPE station_type
   INTEGER id
   REAL(wp) lon, lat
   CHARACTER(LEN=150) :: name
   INTEGER :: np, nr  ! grid box indices 
   LOGICAL :: on_this_pe
   ! output file
  TYPE(file_type) :: ofile
  END TYPE station_type 

  TYPE var_type2d
    CHARACTER(LEN=150) :: iname, oname
    CHARACTER(LEN=150)  :: unit
    INTEGER :: varid
   ! pointer to variable
    REAL(wp), DIMENSION(:,:), POINTER :: ptr2d
  END TYPE var_type2d


  TYPE var_type3d
    CHARACTER(LEN=150) :: iname, oname
    CHARACTER(LEN=150)  :: unit
    INTEGER :: varid
   ! pointer to output stream
    TYPE (t_stream)  ,POINTER   :: stream 
   ! pointer to variable
    REAL(wp), DIMENSION(:,:,:), POINTER :: ptr3d
   ! is this an ilev (mlev+1) variable?
    LOGICAL :: lnlevp1
  END TYPE var_type3d


 ! additional 3-d variables (ones that are not part of a stream)
  INTEGER, PARAMETER :: p_wap   =  2
  INTEGER, PARAMETER :: p_zg    =  1


! nlev, station, variable
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: profiles
  CHARACTER(LEN=150), DIMENSION(nvar_1dp) :: oname1d 
  CHARACTER(LEN=150), DIMENSION(nvar_1dp) :: unit1d 
  INTEGER,  DIMENSION(nvar_1dp) :: varid1d
  

! information on stations
  TYPE(station_type), DIMENSION(nstation) :: st

! information on 2-d variables  
  TYPE(var_type2d),  DIMENSION(nvar_0d ) :: v2d

! information on 3-d variables  
  TYPE(var_type3d),  DIMENSION(nvar_1d_max ) :: v3d


  PUBLIC :: lostation, station_diag, station_diag_nml, & 
            collect_station_diag, init_station_diag, cleanup_station_diag

!------------------------------------------------

CONTAINS

 SUBROUTINE station_diag
! called  each timestep 

  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_time_control,       ONLY: next_date, get_date_components

  IMPLICIT NONE

!local vars
  INTEGER :: year, month, day, hour, minute, second
  REAL(wp) :: yyyymmdd

    IF (.NOT. module_is_initialized ) THEN 
     CALL finish('station_diag: ', 'not initialized')
    END IF

    CALL get_date_components(next_date, year, month, day, hour, minute, second)

    IF ( MOD(minute,30) .EQ. 0 ) THEN
       
      ! day in format %Y%m%d.%f
       yyyymmdd = REAL(ABS(year),wp)*10000._wp+ REAL(month,wp)*100._wp+ REAL(day,wp)  &
                + (REAL(hour,wp)*3600._wp+REAL(minute,wp)*60._wp+REAL(second,wp))/86400._wp
                
             IF (year < 0) THEN
               yyyymmdd = -yyyymmdd
             END IF  

             write(message_text,'(A14,3(I2.2,A))') ' station diag ', &
                         hour,":",minute,":",second
             call message('', message_text )
 

        ncall = ncall + 1

        CALL station_write( yyyymmdd )

    END IF
 END SUBROUTINE station_diag

!--------------------------------------------------
 SUBROUTINE station_write( yyyymmdd )

  USE mo_control,          ONLY : nlev, nlevp1
  USE mo_time_control,     ONLY : delta_time
  USE mo_netcdf
  USE mo_geoloc,           ONLY : sqcst_2d
  USE mo_scan_buffer,      ONLY: vervel

  IMPLICIT NONE
      
  REAL(wp), INTENT(IN) :: yyyymmdd

!local

  REAL(wp) :: v1dt, odelta_time
  REAL(wp), DIMENSION(nlev)   :: v2dt
  REAL(wp), DIMENSION(nlevp1) :: v2dtp1
  INTEGER :: n,i, time_varid
  INTEGER :: ncid, varid
  INTEGER, DIMENSION(1) :: start1, count1
  INTEGER, DIMENSION(2) :: start2, count2, count2p1

    start1 = (/ ncall /)
    count1 = (/ 1 /)

    start2 = (/ 1, ncall /)
    count2 = (/ nlev, 1 /)
    count2p1 = (/ nlevp1, 1 /)

    odelta_time=1._wp/delta_time

    DO n=1, nstation      
       

      IF ( st(n)%on_this_pe ) THEN
          
        ncid =  st(n)%ofile%ncid

      ! fill in time variable
          time_varid=st(n)%ofile%time_varid 
          CALL nf_check(NF_PUT_VARA_DOUBLE(ncid, time_varid, start1, count1, yyyymmdd))


      ! write point variable to file
         DO i = 1, nvar_0d 
              
            v1dt = v2d(i)%ptr2d(st(n)%np,st(n)%nr) 
            varid = v2d(i)%varid

            CALL nf_check(NF_PUT_VARA_DOUBLE( ncid, varid, start1, count1, v1dt  )) 

         END DO


      ! write  profile variable to file
   
         DO i = 1, nvar_1d

           IF ( v3d(i)%lnlevp1 ) THEN

             v2dtp1 =  v3d(i)%ptr3d(st(n)%np,:,st(n)%nr)        
             varid = v3d(i)%varid

             CALL nf_check(NF_PUT_VARA_DOUBLE( ncid, varid, start2, count2p1, v2dtp1  ))

           ELSE

             IF (v3d(i)%iname == 'um1' .OR. v3d(i)%iname == 'vm1') THEN
                v2dt =  v3d(i)%ptr3d(st(n)%np,:,st(n)%nr) / sqcst_2d(st(n)%np,st(n)%nr)
             ELSE
                v2dt =  v3d(i)%ptr3d(st(n)%np,:,st(n)%nr)
             END IF
             varid = v3d(i)%varid

             CALL nf_check(NF_PUT_VARA_DOUBLE( ncid, varid, start2, count2, v2dt  ))
 
           END IF  

         END DO



      ! write additional profile variable to file
   
         DO i = 1, nvar_1dp-1

             v2dt = profiles(:,n,i)        
             varid = varid1d(i)

            CALL nf_check(NF_PUT_VARA_DOUBLE( ncid, varid, start2, count2, v2dt  )) 

         END DO


            v2dt = vervel(st(n)%np,:,st(n)%nr)
            varid = varid1d(p_wap)

            CALL nf_check(NF_PUT_VARA_DOUBLE( ncid, varid, start2, count2, v2dt  )) 

      END IF

    END DO


 END SUBROUTINE station_write

!-----------------------------------------------------------------

 SUBROUTINE station_diag_nml
! read namelist only, init_station_diag (below) is called from station_diag
  USE mo_mpi,         ONLY: p_parallel, p_parallel_io, p_bcast, p_io
  USE mo_namelist,    ONLY: position_nml, POSITIONED, open_nml

  IMPLICIT NONE 

      INTEGER :: ierr, inml, iunit

      NAMELIST /stationctl/ lostation

     
     IF (p_parallel_io) THEN
         inml = open_nml ('namelist.echam')
         iunit = position_nml ('STATIONCTL', inml, status=ierr)
         SELECT CASE (ierr)
         CASE (POSITIONED)
            READ (iunit, stationctl)          
         END SELECT
      END IF


     IF (p_parallel) THEN
         CALL p_bcast (lostation, p_io)
     END IF


 END SUBROUTINE station_diag_nml
!------------------------------------------------------------------------

 SUBROUTINE init_station_diag

   USE mo_decomposition, ONLY: ldc => local_decomposition

 ! initialize information on output variables
  CALL init_svars 
 ! set pointers to variables
  CALL init_pointers
 ! determine PE and coordinates (row, proma) of stations  
  CALL station_find
 ! open output file for each PE 
  CALL open_station_file
 ! write lon, lat arrays
  CALL write_lon_lat
 ! allocate
  ALLOCATE(profiles(ldc%nlev, nstation, nvar_1dp-1)) 
 ! save present value for de-accumulating variables 
 !  (non-zero after restarts) 
   ncall=0
   module_is_initialized = .TRUE.


 END SUBROUTINE  init_station_diag

!----------------------------------------

 SUBROUTINE cleanup_station_diag

  CALL close_station_file

  DEALLOCATE(profiles)

 END SUBROUTINE cleanup_station_diag 

!-------------------------------------
 SUBROUTINE init_svars

  USE mo_exception,          ONLY : message, message_text, finish
  USE mo_memory_gl,          ONLY : gl
  USE mo_memory_g3b,         ONLY : g3b
  USE mo_memory_g1a,         ONLY : g1a
  USE mo_memory_g2a,         ONLY : g2a
  USE mo_memory_cfdiag,      ONLY : cfdiag, locfdiag
  USE mo_diag_tendency_new,  ONLY : tdiag_vars
  USE mo_control,            ONLY : ltdiag

  IMPLICIT NONE
! initialize information on output variables
  INTEGER ::  i, j

   i=1
  v2d(i)%iname='temp2'
  v2d(i)%oname='tas'

   i=i+1 
  v2d(i)%iname='tsurf_na'
  v2d(i)%oname='ts'

   i=i+1
  v2d(i)%iname='aps'
  v2d(i)%oname='ps'

   i=i+1
  v2d(i)%iname='u10'
  v2d(i)%oname='uas'
 
   i=i+1
  v2d(i)%iname='v10'
  v2d(i)%oname='vas'

   i=i+1
  v2d(i)%iname='wind10_na'
  v2d(i)%oname='sfcWind'

   i=i+1 
  v2d(i)%iname='aprl_na'
  v2d(i)%oname='aprl'

   i=i+1 
  v2d(i)%iname='aprc_na'
  v2d(i)%oname='prc'

   i=i+1 
  v2d(i)%iname='aprs_na'
  v2d(i)%oname='prsn'

   i=i+1
  v2d(i)%iname='evap_na'
  v2d(i)%oname='evspsbl'
 
   i=i+1
  v2d(i)%iname='ustr_na'
  v2d(i)%oname='tauu'
 
   i=i+1
  v2d(i)%iname='vstr_na'
  v2d(i)%oname='tauv'

   i=i+1
  v2d(i)%iname='ahfl_na'
  v2d(i)%oname='mhfls'

   i=i+1
  v2d(i)%iname='ahfs_na'
  v2d(i)%oname='mhfss'
 
   i=i+1
  v2d(i)%iname='trads_na'
  v2d(i)%oname='trads'

   i=i+1
  v2d(i)%iname='tradsu_na'
  v2d(i)%oname='mrlus'

   i=i+1
  v2d(i)%iname='srads_na'
  v2d(i)%oname='srads'

   i=i+1
  v2d(i)%iname='sradsu_na'
  v2d(i)%oname='mrsus'

   i=i+1 !185
  v2d(i)%iname='srafs_na'
  v2d(i)%oname='srafs'

   i=i+1 !186
  v2d(i)%iname='trafs_na'
  v2d(i)%oname='trafs'
 
   i=i+1 !184
  v2d(i)%iname='srad0d_na'
  v2d(i)%oname='rsdt'
 
   i=i+1 !203
  v2d(i)%iname='srad0u_na'
  v2d(i)%oname='mrsut'
 
   i=i+1 !179
  v2d(i)%iname='trad0_na'
  v2d(i)%oname='mrlut'
 
   i=i+1 !188
  v2d(i)%iname='traf0_na'
  v2d(i)%oname='mrlutcs'

   i=i+1 !187
  v2d(i)%iname='sraf0_na'
  v2d(i)%oname='sraf0'

   i=i+1 !230
  v2d(i)%iname='qvi_na'
  v2d(i)%oname='prw'

   i=i+1 !164
  v2d(i)%iname='aclcov_na'
  v2d(i)%oname='aclcov'

   i=i+1 !231
  v2d(i)%iname='xlvi_na'
  v2d(i)%oname='clwvi'

   i=i+1 !150
  v2d(i)%iname='xivi_na'
  v2d(i)%oname='clivi'

   i=i+1 !217
  v2d(i)%iname='topmax'
  v2d(i)%oname='cct'

   i=i+1
  v2d(i)%iname='slm'
  v2d(i)%oname='slm'

   i=i+1
  v2d(i)%iname='geosp'
  v2d(i)%oname='geosp'


   IF ( i .GT. nvar_0d ) THEN
     CALL finish('mo_station_diag: ', 'increase nvar_0d')
   END IF 

   DO i=1,nvar_0d 
     v2d(i)%unit=unit_info(v2d(i)%iname,g3b) 
   END DO


! profile variables
      i=1
    v3d(i)%iname="aclc"
    v3d(i)%oname="cl"
    v3d(i)%stream => g3b
   
      i=i+1
    v3d(i)%iname="xl"
    v3d(i)%oname="clw"
    v3d(i)%stream => gl

      i=i+1
    v3d(i)%iname="xi"
    v3d(i)%oname="cli"
    v3d(i)%stream => gl
  
        i=i+1
    v3d(i)%iname="q"
    v3d(i)%oname="hus"
    v3d(i)%stream => gl

      i=i+1
    v3d(i)%iname="tm1"
    v3d(i)%oname="ta"
    v3d(i)%stream => g1a
 
      i=i+1
    v3d(i)%iname="um1"
    v3d(i)%oname="ua"
    v3d(i)%stream => g2a

      i=i+1
    v3d(i)%iname="vm1"
    v3d(i)%oname="va"
    v3d(i)%stream => g2a
 
      i=i+1
    v3d(i)%iname="relhum"
    v3d(i)%oname="hur"
    v3d(i)%stream => g3b


   IF ( locfdiag ) THEN

      i=i+1
    v3d(i)%iname="irlu"
    v3d(i)%oname="rlu"
    v3d(i)%stream => cfdiag
 
      i=i+1
    v3d(i)%iname="irsu"
    v3d(i)%oname="rsu"
    v3d(i)%stream => cfdiag

      i=i+1
    v3d(i)%iname="irld"
    v3d(i)%oname="rld"
    v3d(i)%stream => cfdiag

      i=i+1
    v3d(i)%iname="irsd"
    v3d(i)%oname="rsd"
    v3d(i)%stream => cfdiag

      i=i+1
    v3d(i)%iname="irlucs"
    v3d(i)%oname="rlucs"
    v3d(i)%stream => cfdiag
   
      i=i+1
    v3d(i)%iname="irsucs"
    v3d(i)%oname="rsucs"
    v3d(i)%stream => cfdiag

      i=i+1
    v3d(i)%iname="irldcs"
    v3d(i)%oname="rldcs"
    v3d(i)%stream => cfdiag

      i=i+1
    v3d(i)%iname="irsdcs"
    v3d(i)%oname="rsdcs"
    v3d(i)%stream => cfdiag

      i=i+1
    v3d(i)%iname="imc"
    v3d(i)%oname="mc"
    v3d(i)%stream => cfdiag

   ELSE
   
     WRITE(message_text,*) 'radiation flux profile dignostics turned off'
     CALL message('mo_station_diag ', message_text )
   
   END IF

    nvar_1d_s = i

   DO j = 1, nvar_1d_s
     v3d(j)%unit   =  unit_info(v3d(j)%iname,v3d(j)%stream)
   END DO



   IF ( ltdiag ) THEN
     IF (ASSOCIATED(tdiag_vars%dqdt_vdiff)) THEN
        i=i+1
       v3d(i)%iname="dqdt_vdiff"
       v3d(i)%oname="dqdt_vdiff"
       v3d(i)%ptr3d => tdiag_vars%dqdt_vdiff
       v3d(i)%unit = '1/day'
     END IF
     IF (ASSOCIATED(tdiag_vars%dqdt_cucall)) THEN
        i=i+1
       v3d(i)%iname="dqdt_cucall"
       v3d(i)%oname="dqdt_cucall"
       v3d(i)%ptr3d => tdiag_vars%dqdt_cucall
       v3d(i)%unit = '1/day'
     END IF
     IF (ASSOCIATED(tdiag_vars%dqdt_cloud)) THEN
        i=i+1
       v3d(i)%iname="dqdt_cloud"
       v3d(i)%oname="dqdt_cloud"
       v3d(i)%ptr3d => tdiag_vars%dqdt_cloud
       v3d(i)%unit = '1/day'
     END IF
!!$     IF (ASSOCIATED(tdiag_vars%dqdt_tpcore)) THEN
!!$        i=i+1
!!$       v3d(i)%iname="dqdt_tpcore"
!!$       v3d(i)%oname="dqdt_tpcore"
!!$       v3d(i)%ptr3d => tdiag_vars%dqdt_tpcore
!!$       v3d(i)%unit = '1/day'
!!$     END IF
     IF (ASSOCIATED(tdiag_vars%dtdt_vdiff)) THEN
        i=i+1
       v3d(i)%iname="dtdt_vdiff"
       v3d(i)%oname="dtdt_vdiff"
       v3d(i)%ptr3d => tdiag_vars%dtdt_vdiff
       v3d(i)%unit = 'K/day'
     END IF
     IF (ASSOCIATED(tdiag_vars%dtdt_cucall)) THEN
        i=i+1
       v3d(i)%iname="dtdt_cucall"
       v3d(i)%oname="dtdt_cucall"
       v3d(i)%ptr3d => tdiag_vars%dtdt_cucall
       v3d(i)%unit = 'K/day'
     END IF
     IF (ASSOCIATED(tdiag_vars%dtdt_cloud)) THEN
        i=i+1
       v3d(i)%iname="dtdt_cloud"
       v3d(i)%oname="dtdt_cloud"
       v3d(i)%ptr3d => tdiag_vars%dtdt_cloud
       v3d(i)%unit = 'K/day'
     END IF
     IF (ASSOCIATED(tdiag_vars%dtdt_rheat_sw)) THEN
        i=i+1
       v3d(i)%iname="dtdt_rheat_sw"
       v3d(i)%oname="dtdt_rheat_sw"
       v3d(i)%ptr3d => tdiag_vars%dtdt_rheat_sw
       v3d(i)%unit = 'K/day'
     END IF
     IF (ASSOCIATED(tdiag_vars%dtdt_rheat_lw)) THEN
        i=i+1
       v3d(i)%iname="dtdt_rheat_lw"
       v3d(i)%oname="dtdt_rheat_lw"
       v3d(i)%ptr3d => tdiag_vars%dtdt_rheat_lw
       v3d(i)%unit = 'K/day'
     END IF
     IF (ASSOCIATED(tdiag_vars%dtdt_hines)) THEN
        i=i+1
       v3d(i)%iname="dtdt_hines"
       v3d(i)%oname="dtdt_hines"
       v3d(i)%ptr3d => tdiag_vars%dtdt_hines
       v3d(i)%unit = 'K/day'
     END IF
     IF (ASSOCIATED(tdiag_vars%dtdt_sso)) THEN
        i=i+1
       v3d(i)%iname="dtdt_sso"
       v3d(i)%oname="dtdt_sso"
       v3d(i)%ptr3d => tdiag_vars%dtdt_sso
       v3d(i)%unit = 'K/day'
     END IF

   END IF

   nvar_1d = i




   IF ( nvar_1d .GT. nvar_1d_max ) THEN
     CALL finish('mo_station_diag: ', 'increase nvar_1d_max')
   END IF 


! additional profile variables
 
    oname1d(p_wap)    = "wap"
    unit1d(p_wap)     = "Pa s-1"

    oname1d(p_zg)     = "zg"
    unit1d(p_zg)      = "m"
 


 END SUBROUTINE init_svars
!---------------------------------------------------
 SUBROUTINE init_pointers

   USE mo_memory_g3b
   USE mo_memory_base,      ONLY: get_stream_element
   USE mo_memory_base,      ONLY: get_stream_element_info,  memory_info
   USE mo_exception,        ONLY : message_text, finish
   USE mo_control,          ONLY: nlev, nlevp1

! local

   INTEGER :: i
   CHARACTER(LEN=150) :: cname
   TYPE(t_stream), POINTER :: stream
   TYPE (memory_info) :: info

! set pointers 
  DO i = 1, nvar_0d 
     cname=TRIM(v2d(i)%iname)
     CALL get_stream_element (g3b, cname, v2d(i)%ptr2d)  
  END DO
 
  DO i = 1, nvar_1d_s 
     cname=TRIM(v3d(i)%iname)
     stream => v3d(i)%stream
     CALL get_stream_element (stream, cname, v3d(i)%ptr3d)

     CALL  get_stream_element_info( stream, TRIM(cname), info)
      
     !  set lnlevp1 
     IF (info%klev .EQ. nlev ) THEN
       v3d(i)%lnlevp1 = .FALSE.
     ELSE IF (info%klev .EQ. nlevp1 ) THEN
       v3d(i)%lnlevp1 = .TRUE.
     ELSE
         WRITE(message_text,*) TRIM(v3d(i)%iname), ' vertical grid not supported'
         CALL finish('mo_station_diag: ', message_text)
     END IF

  END DO

 END SUBROUTINE init_pointers

!---------------------------------------------------

 SUBROUTINE  collect_station_diag(kbdim, klev, krow, pgeom1 ) 
!!, pvervel )
!, pgeom1 )
   USE mo_exception,          ONLY: finish
   

   INTEGER, INTENT(IN) :: kbdim, klev, krow

   REAL(wp), DIMENSION( kbdim,klev), INTENT(IN ) ::   pgeom1

!, &
!                                                      pvervel
!, &
!                                                      pgeom1


! local
   INTEGER :: n, jk

    IF (.NOT. module_is_initialized ) THEN 
     CALL finish('collect_station_diag: ', 'not initialized')
    END IF

    DO n=1, nstation      

      IF ( st(n)%on_this_pe ) THEN
        IF ( st(n)%nr .EQ. krow ) THEN     
          DO jk = 1,klev
            profiles( jk, n, p_zg ) = pgeom1 (st(n)%np, jk) 
          END DO
        END IF
      END IF

    END DO
    

 END SUBROUTINE  collect_station_diag




!---------------------------------------------------


 SUBROUTINE station_find

! called during initialization
! determine PE and coordinates (row, proma) of stations  

  USE mo_exception,          ONLY: message, message_text

  USE mo_geoloc,        ONLY: philon_2d, philat_2d
  USE mo_decomposition, ONLY: ldc=>local_decomposition
  USE mo_mpi,         ONLY: p_parallel, p_parallel_io, p_bcast, &
                              p_io, p_pe, p_nprocs

!local
  REAL(wp) :: lon, lat
  INTEGER :: id, n, jl, nproma, jrow, np, ngpblks
  REAL(wp) :: dl, dlt
  CHARACTER(LEN=150) :: sname

  REAL(wp), ALLOCATABLE, DIMENSION(:, :) :: recvbuf 
  REAL(wp), DIMENSION(nstation ) :: sendbuf

  INTEGER, DIMENSION(nstation ) :: onpe
  INTEGER, DIMENSION(nstation ) :: myjl, myjrow

  ALLOCATE(recvbuf(nstation, p_nprocs))

! initialize
   onpe(:) = -999
   st(:)%np = -999
   st(:)%nr = -999
   recvbuf(:,:) = -999._wp 
   ngpblks = ldc % ngpblks
   myjl(:) = -999
   myjrow(:) = -999

! read in station locations
    IF (p_parallel_io) THEN
     OPEN(iun,FILE="pointlocations.txt")
    END IF

   DO n=1,nstation
     IF (p_parallel_io) THEN
           !READ(iun, *) st(n)%id,  st(n)%lon,  st(n)%lat, st(n)%name
           READ(iun, *) id,  lon,  lat, sname
           IF ( lon .LT. 0._wp ) THEN 
               lon = 360._wp +  lon
           END IF
             write(message_text,'(I4.4,2F7.2,2X,A30)') id,  lon,  lat, sname
             call message('', message_text )
     END IF
 
     IF (p_parallel) THEN
         CALL p_bcast (id, p_io)
         CALL p_bcast (lon, p_io)
         CALL p_bcast (lat, p_io)
         CALL p_bcast (sname, p_io)
     END IF

     st(n)%id=id 
     st(n)%lon=lon 
     st(n)%lat=lat 
     st(n)%name=sname
   END DO

    IF (p_parallel_io) THEN
     CLOSE(iun)
    END IF


! determine distance from closest grid point on each PE
   DO n=1,nstation
      dl=999999._wp
     DO jrow = 1, ngpblks        
       IF ( jrow == ngpblks ) THEN
         nproma = ldc% npromz
       ELSE
         nproma = ldc% nproma
       END IF

       DO jl=1,nproma      
           dlt = (philat_2d(jl,jrow) - st(n)%lat)**2 &
             + ( philon_2d(jl,jrow) - st(n)%lon)**2 
          IF ( dlt .LT.  dl ) THEN
             dl = dlt
             myjl(n) = jl
             myjrow(n) = jrow
          END IF     
       END DO
     END DO

       sendbuf(n) =  dl
      
    END DO

! communicate all the distances to the IO PE
   CALL p_gather_real_1d2d (sendbuf, recvbuf, p_io)
 
! let IO PE decide which PE contains closest grid point
    IF (p_parallel_io) THEN 
      DO n=1,nstation
          dl=999999._wp
        DO np=1,p_nprocs
          IF ( recvbuf(n,np) .GT. 0._wp .AND. recvbuf(n,np) .LE. dl ) THEN
            dl=recvbuf(n,np)
            onpe(n)=np-1
          END IF 
         END DO
      END DO
    END IF


! let other PEs know the decisions
     IF (p_parallel) THEN
         CALL p_bcast (onpe, p_io)
     END IF


! acknowledge decisions
    DO n=1,nstation
      st(n)%on_this_pe = .false.
      IF ( onpe(n) .EQ. p_pe ) THEN
        st(n)%on_this_pe = .true.
      END IF
    END DO


! on PEs which have stations:
! save the indices of the closest grid points 

    DO n=1,nstation
      IF ( st(n)%on_this_pe ) THEN
             st(n)%np=myjl(n)
             st(n)%nr=myjrow(n)
      END IF
    END DO


  DEALLOCATE(recvbuf)


!!$! for debug:
!!$    do n=1,nstation
!!$     if ( st(n)%on_this_pe ) then
!!$      write(message_text,*) "YYYY station nr, pe ", n, p_pe, st(n)%on_this_pe
!!$      call message('', message_text, all_print=.true.  )   
!!$      write(message_text,*) "YYYY station   lon lat  ",  st(n)%lon, st(n)%lat
!!$      call message('', message_text, all_print=.true.  )   
!!$      write(message_text,*) "YYYY grid point lon lat ",  philon_2d(st(n)%np,st(n)%nr), philat_2d(st(n)%np,st(n)%nr)
!!$    call message('', message_text, all_print=.true.  )   
!!$      write(message_text,*) "YYYY grid point (local pe) ",  st(n)%np, st(n)%nr
!!$      call message('', message_text, all_print=.true. )   
!!$      write(message_text,*) "YYYY -----------"
!!$      call message('', message_text, all_print=.true. )
!!$     end if
!!$    enddo
!!$    
!!$ CALL finish('mo_station_diag: ', 'TEST')

 END SUBROUTINE station_find


!---------------------------------------------------------------------
 SUBROUTINE open_station_file

   USE mo_exception,          ONLY: message, message_text
   USE mo_netcdf
   USE mo_mpi,                ONLY: p_pe 
   USE mo_control,            ONLY: nlev, nlevp1
   USE mo_time_control,       ONLY: next_date, timelabel_type,  str_date
   USE mo_filename,           ONLY: out_expname

   IMPLICIT NONE

!local

  INTEGER :: n, i, ncid, varid, lat_varid, lon_varid, time_varid
  INTEGER :: mlev_dimid, ilev_dimid, time_dimid
  
  INTEGER, DIMENSION(1) :: dimids1d 
  INTEGER, DIMENSION(2) :: dimids2d 
  CHARACTER(LEN=150) :: dumname
  CHARACTER(len=19) :: time_label



    time_label = TRIM(str_date(timelabel_type,next_date))

    station_loop: DO n=1,nstation

     pe_if: IF ( st(n)%on_this_pe ) THEN

    ! determine filename
       WRITE(dumname,'(A8,I4.4,A1)') "_cfSites_", n
         dumname =TRIM(out_expname)//"_"//TRIM(time_label)//TRIM(dumname)//'.nc'
        
     ! open output file
        CALL nf_check( NF_CREATE(TRIM(dumname), NF_CLOBBER, ncid))
        st(n)%ofile%name =  dumname

        WRITE(message_text,*) ' Open cfSites file ',  TRIM(dumname), ncid
        CALL message('', message_text, all_print=.true. )


     ! declare dimensions
      
       CALL nf_check( NF_DEF_DIM(ncid, 'mlev', nlev, mlev_dimid) )
       CALL nf_check( NF_DEF_DIM(ncid, 'ilev', nlevp1, ilev_dimid) )
       CALL nf_check( NF_DEF_DIM(ncid, 'time', NF_UNLIMITED, time_dimid) )

       dimids1d=(/ time_dimid /)

      !declare special variables

       CALL nf_check( NF_DEF_VAR (ncid, 'lon', NF_DOUBLE, 0, 0, lon_varid ) )
       CALL nf_check( NF_DEF_VAR (ncid, 'lat', NF_DOUBLE, 0, 0, lat_varid ) )
       CALL nf_check( NF_DEF_VAR (ncid, 'time', NF_DOUBLE, 1, dimids1d, time_varid ) )

       CALL nf_check (NF_PUT_ATT_TEXT(ncid,lon_varid,'long_name',9, 'longitude'))
       CALL nf_check (NF_PUT_ATT_TEXT(ncid,lon_varid,'units', 12, 'degrees_east'))
       CALL nf_check (NF_PUT_ATT_TEXT(ncid,lon_varid,'standard_name',9,'longitude'))
       CALL nf_check (NF_PUT_ATT_TEXT(ncid,lat_varid,'long_name',8, 'latitude'))
       CALL nf_check (NF_PUT_ATT_TEXT(ncid,lat_varid,'units',13, 'degrees_north'))
       CALL nf_check (NF_PUT_ATT_TEXT(ncid,lat_varid,'standard_name',8, 'latitude'))

       CALL nf_check (NF_PUT_ATT_TEXT(ncid,time_varid,'units', 16,'day as %Y%m%d.%f'))
       CALL nf_check (NF_PUT_ATT_TEXT(ncid,time_varid,'calendar',19,'proleptic_gregorian'))


    ! declare point variables
       DO i=1,nvar_0d

          dumname =  v2d(i)%oname
          CALL nf_check( NF_DEF_VAR(ncid, TRIM(dumname), NF_DOUBLE, 1, &
                     dimids1d,  varid) )
          v2d(i)%varid=varid
  
          CALL put_unit( ncid, varid, v2d(i)%unit , v2d(i)%oname )

       END DO

    ! declare profile variables
       DO i=1,nvar_1d
         
         IF ( v3d(i)%lnlevp1 ) THEN
          dimids2d=(/ ilev_dimid, time_dimid /)
         ELSE  
          dimids2d=(/ mlev_dimid, time_dimid /)
         END IF

 
         dumname =  v3d(i)%oname
         CALL nf_check( NF_DEF_VAR(ncid, TRIM(dumname), NF_DOUBLE, 2, &
                    dimids2d,  varid) )
         v3d(i)%varid=varid

         CALL put_unit( ncid, varid, v3d(i)%unit , v3d(i)%oname )


       END DO

    ! declare additional profile variables

          dimids2d=(/ mlev_dimid, time_dimid /)

       DO i=1,nvar_1dp

          dumname =  oname1d(i)
          CALL nf_check( NF_DEF_VAR(ncid, TRIM(dumname), NF_DOUBLE, 2, &
                     dimids2d,  varid) )
          varid1d(i)=varid

          CALL put_unit( ncid, varid, unit1d(i) , oname1d(i))

       END DO



       CALL nf_check( NF_ENDDEF(ncid) )

      st(n)%ofile%ncid = ncid
      st(n)%ofile%lat_varid = lat_varid
      st(n)%ofile%lon_varid = lon_varid
      st(n)%ofile%time_varid = time_varid
      st(n)%ofile%mlev_dimid = mlev_dimid
      st(n)%ofile%ilev_dimid = ilev_dimid
      st(n)%ofile%time_dimid = time_dimid


        WRITE(message_text,*) 'mo_station_diag: opened file on PE ', st(n)%ofile%name,  p_pe
        CALL message('', message_text, all_print=.true. )

     END IF pe_if 

   END DO station_loop

 END SUBROUTINE open_station_file

!---------------------------------------------------------------------
 SUBROUTINE write_lon_lat

   USE mo_geoloc,        ONLY: philon_2d, philat_2d
   USE mo_netcdf

   IMPLICIT NONE

!local
   REAL(wp) :: olon, olat
   INTEGER :: n

    ! fill lon, lat arrays
      DO n=1, nstation
         IF ( st(n)%on_this_pe ) THEN

           olon = philon_2d(st(n)%np,st(n)%nr)
           olat = philat_2d(st(n)%np,st(n)%nr)

          CALL nf_check(NF_PUT_VAR_DOUBLE( st(n)%ofile%ncid,  st(n)%ofile%lon_varid, olon  ))
          CALL nf_check(NF_PUT_VAR_DOUBLE( st(n)%ofile%ncid,  st(n)%ofile%lat_varid, olat  ))

        END IF
     END DO

    
 END SUBROUTINE write_lon_lat

!---------------------------------------------------------------------
 SUBROUTINE close_station_file

   USE mo_exception,          ONLY: message, message_text
   USE mo_netcdf
   USE mo_mpi,                ONLY: p_pe 

   IMPLICIT NONE

   INTEGER :: n
  
   
   DO n=1,nstation
     IF ( st(n)%on_this_pe ) THEN 
      CALL nf_check( NF_CLOSE(st(n)%ofile%ncid))
      
       WRITE(message_text,*) 'mo_station_diag: closed file on PE ', st(n)%ofile%name,  p_pe
       CALL message('', message_text, all_print=.true. )

     END IF
   END DO 

 END SUBROUTINE close_station_file

!---------------------------------------------------------------------
FUNCTION unit_info(iname,stream)  RESULT (unit)

  USE mo_memory_base,      ONLY: get_stream_element_info,  memory_info

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN ) :: iname
  TYPE (t_stream)  ,POINTER   :: stream 

  CHARACTER(LEN=150) :: unit

!local
  TYPE (memory_info) :: info

  CALL get_stream_element_info  (stream, TRIM(iname), info)
  unit=info%units
 
  IF ( LEN(TRIM(unit )) == 0 ) THEN
   unit="1"
  END IF

END FUNCTION unit_info
!---------------------------------------------------------------------
SUBROUTINE put_unit ( ncid, varid, unit, oname )

  USE mo_netcdf 
  USE mo_exception,          ONLY: message_text,finish

  INTEGER, INTENT (IN ) :: ncid, varid
  CHARACTER(LEN=150), INTENT(IN ) :: unit, oname

!local
  INTEGER :: lent

          lent=LEN(TRIM(unit))
          IF ( lent == 0 ) THEN
           WRITE(message_text,*)  "UNIT NOT SPECIFIED ", oname
           CALL finish('mo_station_diag: ', message_text)
          END IF
          CALL nf_check (NF_PUT_ATT_TEXT(ncid,varid,'units', lent, TRIM(unit)))

END SUBROUTINE put_unit
!---------------------------------------------------------------------
! temporary :use from module_..
 SUBROUTINE p_gather_real_1d2d (sendbuf, recvbuf, p_dest, comm)

    USE mo_mpi,         ONLY:p_all_comm, p_real_dp

    IMPLICIT NONE

    REAL(wp),          INTENT(inout) :: sendbuf(:), recvbuf(:,:)
    INTEGER,           INTENT(in) :: p_dest
    INTEGER, OPTIONAL, INTENT(in) :: comm

    INTEGER :: p_error

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

     CALL MPI_GATHER(sendbuf, SIZE(sendbuf), p_real_dp, &
                     recvbuf, SIZE(sendbuf), p_real_dp, &
                     p_dest, p_comm, p_error)
#else
     recvbuf(:,LBOUND(recvbuf,2)) = sendbuf(:)
#endif
 END SUBROUTINE p_gather_real_1d2d



END MODULE mo_station_diag
