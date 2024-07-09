!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_cosp_offline
! TODO: 
!     stop when restart does not coincide with output file freq.
!     split output files
!     check that all allocated arrays are deallocated
!     allow arbitrary restart times (tm1 variables), also in mo_station_diag
!----------------------------------------------------
! PURPOSE:
!
!  write out COSP input in orbital curtain format
!  for running the COSP simulators in OFFLINE mode
!
!  see  cf3hr in K. Taylor's CMIP table
! 
! Strategy:   
!           determine indices of points along A-Train orbit 
!             - one month at a time for each timestep
!               curtains are centered at each timestep
!             - curtains can span several PEs
!               each full curtain consists of segments (i.e. part
!                of a curtain that is stored on a particular PE).
!
!           for accumulated stream variables calculate instanteneous value 
!           using value from previous timestep 
!
!           initialize pointers to 2-D and 3-D stream elements
!           to avoid calls to get_stream_element
! 
!           where 3-D variables are not available in existing streams, 
!           write them to a new stream called cospoffl
!
! Entry points:
!     - cosp_offline_nml is called from initialize.f90 (reads namelist)
!     - cosp_offline is called from stepon.f90 (outside openMP loop)
!     - construct_stream_cospoffl, destruct_stream_cospoffl are called from mo_memory_streams.f90
!     - collect data in mo_radiation.f90, physc.f90, cloud.f90, cuflx.f90 
!     - control.f90: cosp_offline_init (needs geoloc, etc.)
!
! Other:
!
!  -cosp_offline 
!        calls:
!          - init_cosp_offline (initialize curtains anew for each month, set pointers to variables, see below)
!          - cosp_offl_write (sample  model data along curtain, write to file, see below)   
!          - cosp_offline_finalize (close file, deallocate vars)
!
!     -init_cosp_offline (called once every model month from cosp_offline) 
!        calls:
!            - init_svars  
!                  (initialize info on variables)                
!            - init_pointers
!                   (initialize pointers to stream elements)        
!            - read_track  
!                   ( read in satellite orbit data)
!            - get_indidices
!                  (determine PE and coordinates (row, proma) of points along satellite track)
!            - open_output_file (open output file, define variables, etc.)
!
!      
!
!       -cosp_offl_write:
!          calls:  
!                 - collect_cu1d_wrap (which calls collect_cu1d)
!                 - collect_cu2d_wrap (which calls collect_cu2d)
!                 ( assemble entire curtains from segments, cu1d: points, cu2d: profiles)  
!                   ind_sort_cu1d: sort points according to orbit time 
!                                  since initially segments that are inserted into the 
!                                  final curtains are ordered based on pe id.)
!                                        
!------------------------------------------------------------------------------------------------


   USE mo_linked_list,   ONLY: t_stream
   USE mo_kind,          ONLY: dp
   USE mo_linked_list,   ONLY: t_stream, NETCDF

   IMPLICIT NONE

   INCLUDE 'netcdf.inc'

   PRIVATE

! do cf3hr output?
   LOGICAL :: locospoffl = .FALSE. 
! output frequ for 2-d offline testing only
   INTEGER :: offl2dout = -1
!-------------------------------------------------
! free IO/unit ?
   INTEGER, PARAMETER :: iun=2124
!-------------------------------------------------
! orbital locations and times for one month sampled at model time steps
!  first dimension: location; second dimension: time step   
   REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: slon, slat, stime
!
!-------------------------------------------------
! count number of calls, reset on restart
   INTEGER :: ncall = 1 
!
!-------------------------------------------------
!  number of point variables
   INTEGER, PARAMETER :: nvar_0d = 6
!
!-------------------------------------------------
!  maximum number of profile variables from streams
!   (not so nice, try linked list next time.) 
   INTEGER, PARAMETER :: nvar_1d_max = 20
!  
!  actual number of profile variables 
!  ( lower than nvar_1d_max only if not all streams/variables
!    are available)
   INTEGER :: nvar_1d
!
!-------------------------------------------------
! has this module been initialized?
   INTEGER ::  cur_month = -999
!-------------------------------------------------

  TYPE file_type 
    CHARACTER(LEN=150) :: name
   ! file handle
    INTEGER :: ncid
   !axis ids
    INTEGER :: point_dimid
    INTEGER :: lev_dimid
    INTEGER :: hydro_dimid
   ! special variable ids
    INTEGER :: year_varid, month_varid, day_varid 
    INTEGER :: hour_varid, minute_varid, second_varid 
    INTEGER :: lon_varid, lat_varid, tim_varid
   ! maximum number of points in an entire curtain
   ! for the whole month
    INTEGER  :: gmax 
  END TYPE file_type 

  TYPE curtain_type
   ! number of points in the entire curtain
   ! for the interval [t-dt/2,t+dt/2] for this particular time step
   !  note that each curtain can span several PEs 
   INTEGER :: gmax
   ! number of points in  the curtain segment 
   ! that is found on the local PE
   INTEGER :: pmax 
   ! longest curtain segment on any PE
   INTEGER :: pamax 
   ! starting point for inserting curtain segment into the entire
   ! curtain (for this particular time step)
   INTEGER :: start
   !
   LOGICAL, ALLOCATABLE, DIMENSION(: ) :: on_this_pe
   ! co-ordinates
   INTEGER, ALLOCATABLE, DIMENSION(: ) :: myjl, myjrow
   ! exact time in orbit 
   REAL(dp), ALLOCATABLE, DIMENSION(: ) :: mytime
  END TYPE curtain_type

  TYPE var_type2d
    CHARACTER(LEN=150) :: iname, oname
    CHARACTER(LEN=150)  :: unit
    INTEGER :: varid
    LOGICAL :: laccu  
   ! pointer to output stream
    TYPE (t_stream)  ,POINTER   :: stream 
   ! pointer to variable
    REAL(dp), DIMENSION(:,:), POINTER :: ptr2d
   ! value of accumulated 2-d variable from previous call
   ! curtain t+1 !!
    REAL(dp), DIMENSION(:), ALLOCATABLE :: tm1
  END TYPE var_type2d

  TYPE var_type3d
    CHARACTER(LEN=150) :: iname, oname
    CHARACTER(LEN=150)  :: unit
    INTEGER :: varid
    LOGICAL :: laccu
   ! pointer to output stream
    TYPE (t_stream)  ,POINTER   :: stream 
   ! pointer to variable
    REAL(dp), DIMENSION(:,:,:), POINTER :: ptr3d
   ! value of accumulated 2-d variable from previous call
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: tm1
   ! is this an ilev (mlev+1) variable?
    LOGICAL :: lnlevp1
  END TYPE var_type3d

! output file
  TYPE(file_type) :: ofile

! information on 2-d variables  
  TYPE(var_type2d),  DIMENSION(nvar_0d ) :: v2d

! information on 3-d variables  
  TYPE(var_type3d),  DIMENSION(nvar_1d_max ) :: v3d

! curtain  
  TYPE(curtain_type), POINTER, DIMENSION(:) :: cu

! stream to collect extra vars
  TYPE (t_stream), PUBLIC, POINTER :: cospoffl

  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cospoffl_lsrain ! Flux large scale cloud rain [kg m^-2 s^-1]
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cospoffl_lssnow ! Flux large scale cloud snow [kg m^-2 s^-1]
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cospoffl_ccrain ! Flux convective cloud rain [kg m^-2 s^-1]
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cospoffl_ccsnow ! Flux convective cloud snow [kg m^-2 s^-1]

  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cospoffl_reffl  ! Liquid water droplet effective radius [um]
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cospoffl_reffi  ! Ice crystal effective radius [um]

  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cospoffl_geom1  ! geopotential height [m3 s-2]
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cospoffl_geohm1 ! geopotential height at interfaces  [m3 s-2]

  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cospoffl_p      ! pressure [Pa]
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cospoffl_ph     ! pressure at interfaces  [Pa]

  LOGICAL :: module_is_initialized


  PUBLIC :: locospoffl, cosp_offline, cosp_offline_nml, & 
            construct_stream_cospoffl, destruct_stream_cospoffl, &
            init_cosp_offline

!------------------------------------------------

CONTAINS

 SUBROUTINE cosp_offline
! called  each timestep 

  USE mo_exception,          ONLY: message, message_text
  USE mo_time_control,       ONLY: l_putrerun, lbreak, lstop, &
                                   next_date, get_date_components

  IMPLICIT NONE

!local vars
  INTEGER :: year, month, day, hour, minute, second


  INTEGER :: n,i, nr, np
 
  TYPE(curtain_type), POINTER :: cur, next 


    CALL get_date_components(next_date, year, month, day, hour, minute, second)

  ! find points along track, initialize variables, open output file
    IF ( month .NE. cur_month .OR. .NOT. module_is_initialized ) THEN 
       ncall = 1
     IF ( module_is_initialized ) THEN
       CALL cosp_offline_finalize 
     END IF
     ! each full month
     CALL init_cosp_offline
     cur_month = month
    END IF

     cur => cu(ncall)
   

             write(message_text,'(A14,3(I2.2,A))') ' cosp offl ', &
                         hour,":",minute,":",second
             call message('', message_text )
 

       CALL cosp_offl_write(  cur, year, month, day, hour, minute, second )



   ! for accumulated 2-d variables, save value at present time step
   !  but for the next curtain
    IF ( .NOT. ( lstop .OR. (l_putrerun .AND. lbreak) )) THEN

         next => cu(ncall+1)

     DO n=1, next%pmax
       IF ( next%on_this_pe(n) ) THEN 
          np = next%myjl(n)
          nr = next%myjrow(n)
         DO i=1,nvar_0d
           IF (  v2d(i)%laccu ) THEN
             v2d(i)%tm1(n)=v2d(i)%ptr2d(np,nr)
           END IF
         END DO
         DO i=1,nvar_1d
          IF (  v3d(i)%laccu ) THEN
             v3d(i)%tm1(:,n)=v3d(i)%ptr3d(np,:,nr)
          END IF
        END DO
      END IF
    END DO
   END IF


       ncall = ncall + 1

  ! close output file
    IF ((l_putrerun .AND. lbreak) .OR. lstop) THEN
      CALL cosp_offline_finalize
    END IF

write(message_text,* )  'end  '
call message('mo_cosp_offline ', message_text )

 END SUBROUTINE cosp_offline

!--------------------------------------------------
 SUBROUTINE cosp_offl_write(  cur, year, month, day, hour, minute, second )

  USE mo_control,          ONLY: nlev
  USE mo_mpi,              ONLY: p_pe, p_parallel_io,  &
                                  p_nprocs, p_io

  IMPLICIT NONE

  INTEGER, INTENT (IN ) :: year, month, day, hour, minute, second
  TYPE(curtain_type), INTENT (IN ), POINTER :: cur

!local

  INTEGER :: ivar

  REAL(dp), ALLOCATABLE, DIMENSION(:)   :: cu1d
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: cu2d
  REAL(dp), ALLOCATABLE, DIMENSION(:)   :: etime
  INTEGER,  ALLOCATABLE, DIMENSION(:)   :: new_ind

  ! see who is sending part of the curtain (this does not include
  !   the IO PE, since it does not need to send to itself)
  !  isen=.TRUE.   :  Will send something.  
  !  isen=.FALSE.  :  Nope. Got nothing to send.
  LOGICAL :: isen
  ! see who has part of the curtain (like isen, but might include IO PE)
  LOGICAL :: igot
  !  igot=.TRUE.    : I saw a satellite! 
  !  igot=.FALSE.   : Nope. No A-train here. 

  INTEGER :: isen_pe, is, ie
  INTEGER, DIMENSION( p_nprocs ) :: sending_PEs, istart, iend 
 
   ALLOCATE(cu1d(ofile%gmax))
   ALLOCATE(cu2d(nlev, ofile%gmax))
   ALLOCATE(etime(ofile%gmax))
   ALLOCATE(new_ind(ofile%gmax))

   etime(:) = -999._dp
   new_ind(:) = -999

   CALL nc_date_write(  cur, year, month, day, hour, minute, second )

! see who has part of the curtain, who is going to send something 
!   and where it will go 
     isen=.FALSE.
     igot=.FALSE.
     isen_pe = -999
     sending_PEs(:) = -999
     is=-999
     ie=-999
    IF (  cur%pmax .GT. 0  ) THEN
      is=cur%start
      ie=is+cur%pmax-1
      igot=.TRUE.
      IF ( .NOT. p_parallel_io ) THEN
        isen=.TRUE.
        isen_pe = p_pe
      END IF
    END IF
    CALL p_gather_int_0d1d (isen_pe, sending_PEs, p_io)
    CALL p_gather_int_0d1d (is, istart, p_io)
    CALL p_gather_int_0d1d (ie, iend, p_io)


! exact orbit time

      CALL collect_cu1d ( etime, cur, ofile%gmax, igot, isen, &
                      sending_PEs, istart, iend, opt='tim' )      

      CALL set_ind_sort_cu( etime, new_ind, ofile%gmax )

      CALL ind_sort_cu1d(etime,  new_ind, ofile%gmax )

      CALL nc_write_cu1d( ofile%tim_varid, etime, ofile%gmax, cur%gmax )


! lon, lat
      CALL collect_cu1d ( cu1d, cur, ofile%gmax, igot, isen, &
                      sending_PEs, istart, iend, opt='lon' )

      CALL ind_sort_cu1d(cu1d,  new_ind, ofile%gmax )

      CALL nc_write_cu1d( ofile%lon_varid, cu1d, ofile%gmax, cur%gmax )


      CALL collect_cu1d ( cu1d, cur, ofile%gmax, igot, isen, &
                      sending_PEs, istart, iend, opt='lat' )
       
      CALL ind_sort_cu1d(cu1d,  new_ind, ofile%gmax )

      CALL nc_write_cu1d( ofile%lat_varid, cu1d, ofile%gmax, cur%gmax )
 

! point variables
    DO ivar=1, nvar_0d
     
     CALL collect_cu1d_wrap( ivar, cu1d, cur, ofile%gmax, igot, isen, &
                      sending_PEs, istart, iend )

     CALL ind_sort_cu1d(cu1d,  new_ind, ofile%gmax )

     CALL nc_write_cu1d( v2d(ivar)%varid, cu1d, ofile%gmax, cur%gmax )

   END DO

! profile variables

    DO ivar=1, nvar_1d
      
     CALL collect_cu2d_wrap( ivar, cu2d, cur, nlev, ofile%gmax, igot, isen, &
                      sending_PEs, istart, iend )

     CALL ind_sort_cu2d(cu2d,  new_ind, nlev, ofile%gmax )

     CALL nc_write_cu2d( v3d(ivar)%varid, cu2d, nlev, ofile%gmax, cur%gmax )

    END DO

    
   DEALLOCATE(cu1d, cu2d)
   DEALLOCATE(etime)   
   DEALLOCATE(new_ind)

 END SUBROUTINE cosp_offl_write

!-----------------------------------------------------------------
 SUBROUTINE collect_cu1d_wrap( ivar, cu1d, cur, gmax, igot, isen, &
                         sending_PEs, istart, iend )


   USE mo_mpi,              ONLY: p_nprocs
   USE mo_time_control,     ONLY: delta_time

   IMPLICIT NONE

   TYPE(curtain_type), INTENT (IN ), POINTER :: cur

   INTEGER, INTENT(IN) :: ivar, gmax
   LOGICAL, INTENT(IN) :: igot, isen

   INTEGER, DIMENSION( p_nprocs ), INTENT(IN) :: sending_PEs, &
                                                 istart,      &
                                                 iend

   REAL(dp), DIMENSION(gmax), INTENT(INOUT ) :: cu1d

! local

    INTEGER :: i
    REAL(dp) :: odelta_time

    REAL(dp), DIMENSION(gmax)  :: cu1d_tm1



       cu1d(:) = -999._dp
       cu1d_tm1(:) = -999._dp

     CALL collect_cu1d( cu1d, cur, gmax, igot, isen, &
                      sending_PEs, istart, iend, ivar=ivar )


    IF (  v2d(ivar)%laccu ) THEN
      odelta_time=1._dp/delta_time
      CALL collect_cu1d( cu1d_tm1, cur, gmax, igot, isen, &
                    sending_PEs, istart, iend, ivar=ivar, opt='tm1' )


     DO i=1, gmax
       IF ( cu1d(i) .NE.  -999._dp .AND. cu1d_tm1(i) .NE.  -999._dp ) THEN 
        cu1d(i)=odelta_time*(cu1d(i)-cu1d_tm1(i))
       ELSE
        cu1d(i)=-999._dp
       END IF
     END DO
    END IF


 END SUBROUTINE collect_cu1d_wrap

!------------------------------------------------------------------------------------

 SUBROUTINE collect_cu1d( cu1d, cur, gmax, igot, isen, &
                         sending_PEs, istart, iend, ivar, opt )
! collect points along curtain
!  default: v2d(ivar), but also options for lat, lon, tm1

  USE mo_mpi,              ONLY: p_pe, p_parallel_io, p_global_comm, &
                                  p_nprocs, p_io,  p_send, p_recv
  USE mo_decomposition,    ONLY: gdc => global_decomposition
  USE mo_geoloc,           ONLY: philon_2d, philat_2d


  IMPLICIT NONE  
  
   TYPE(curtain_type), INTENT (IN ), POINTER :: cur

   INTEGER, INTENT(IN) :: gmax
   LOGICAL, INTENT(IN) :: igot, isen

   INTEGER, DIMENSION( p_nprocs ), INTENT(IN) :: sending_PEs, &
                                                 istart,      &
                                                 iend

   REAL(dp), DIMENSION(gmax), INTENT(INOUT ) :: cu1d

   ! get v2d(ivar)
   INTEGER, OPTIONAL, INTENT(IN) :: ivar

   ! special options: curtain from last time step, lon, lat
   CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: opt

   

!local
   REAL(dp), ALLOCATABLE, DIMENSION(:) :: sendbuf1d, recvbuf1d
  
   INTEGER :: i, np ,nr, ie, is, je, pe

   CHARACTER(LEN=3) :: lopt
   LOGICAL :: got_opt, got_ivar 

   got_opt =  PRESENT( opt )
   got_ivar = PRESENT( ivar )   

   IF ( got_opt ) THEN
      lopt=opt
   ELSE
      lopt='var'
   END IF

   CALL collect_arg_check(got_opt, got_ivar, lopt, 'collect_cu1d')

   ALLOCATE(sendbuf1d(cur%pamax))
   ALLOCATE(recvbuf1d(cur%pamax))  

         sendbuf1d(:) = -999._dp
         recvbuf1d(:) = -999._dp        

! fill send buffer for all PEs which have part of the curtain
!  this includes pe_io
    IF (  igot ) THEN
      SELECT CASE (TRIM(lopt))
      CASE('lon')
        DO i=1,  cur%pmax  
          np = cur%myjl(i)
          nr = cur%myjrow(i)
          IF ( philon_2d(np,nr) .LE. 180._dp ) THEN
            sendbuf1d(i) = philon_2d(np,nr) 
          ELSE
           sendbuf1d(i) = philon_2d(np,nr)  - 360._dp
          END IF
        END DO 
      CASE('lat')
        DO i=1,  cur%pmax  
          np = cur%myjl(i)
          nr = cur%myjrow(i)
          sendbuf1d(i) = philat_2d(np,nr) 
        END DO 
      CASE('tm1')
        DO i=1,  cur%pmax  
          sendbuf1d(i) = v2d(ivar)%tm1(i)
        END DO
      CASE('tim')
        DO i=1,  cur%pmax
          sendbuf1d(i) = cur%mytime(i) 
        END DO 
      CASE ('var')
        DO i=1,  cur%pmax  
          np = cur%myjl(i)
          nr = cur%myjrow(i)
          sendbuf1d(i) = v2d(ivar)%ptr2d(np,nr) 
        END DO
      END SELECT
    END IF 


! send
    IF ( isen ) THEN
       CALL p_send(sendbuf1d, p_io, 700+p_pe, comm=p_global_comm) 
    END IF


! receive, assemble, and write to output file
   
    IF ( p_parallel_io ) THEN
      DO i = 1, p_nprocs
           recvbuf1d(:) = -999._dp
           pe    = gdc(i)%pe
           IF ( ANY(sending_PEs .EQ. pe)  ) THEN

            CALL p_recv(recvbuf1d, pe, 700+pe, comm=p_global_comm)

           END IF
           IF ( pe .EQ. p_io ) THEN
              recvbuf1d(:) = sendbuf1d(:)
           END IF

           is = istart(i)
           ie = iend(i)
           je = ie - is +1 ! pmax

           IF ( is .GT. 0 ) THEN
             cu1d(is:ie) = recvbuf1d(1:je)
           END IF
       END DO   

     END IF

!!$      if ( p_parallel_io ) then
!!$       write(655,*) "---------- ",i
!!$       write(655,*) " is, ie, je ", is, ie, je 
!!$       write(655,*) recvbuf1d(:)
!!$       write(655,*) "======================================"
!!$       write(655,*) cu1d(:)
!!$       write(655,*) "======================================"
!!$     end if

 END SUBROUTINE collect_cu1d

!-----------------------------------------------------------------
 SUBROUTINE collect_cu2d_wrap( ivar, cu2d, cur, nlev, gmax, igot, isen, &
                         sending_PEs, istart, iend )


   USE mo_mpi,              ONLY: p_nprocs
   USE mo_time_control,     ONLY: delta_time

   IMPLICIT NONE

   TYPE(curtain_type), INTENT (IN ), POINTER :: cur

   INTEGER, INTENT(IN) :: ivar, nlev, gmax
   LOGICAL, INTENT(IN) :: igot, isen

   INTEGER, DIMENSION( p_nprocs ), INTENT(IN) :: sending_PEs, &
                                                 istart,      &
                                                 iend

   REAL(dp), DIMENSION(nlev, gmax), INTENT(INOUT ) :: cu2d

! local

    INTEGER :: i, n
    REAL(dp) :: odelta_time

    REAL(dp), DIMENSION(nlev, gmax)  :: cu2d_tm1



       cu2d(:,:) = -999._dp
       cu2d_tm1(:,:) = -999._dp

     CALL collect_cu2d( cu2d, cur, nlev, gmax, igot, isen, &
                      sending_PEs, istart, iend, ivar=ivar )


    IF (  v3d(ivar)%laccu ) THEN
      odelta_time=1._dp/delta_time
      CALL collect_cu2d( cu2d_tm1, cur, nlev, gmax, igot, isen, &
                    sending_PEs, istart, iend, ivar=ivar, opt='tm1' )


     DO i=1, gmax
       DO n=1,nlev
         IF ( cu2d(n,i) .NE.  -999._dp .AND. cu2d_tm1(n,i) .NE.  -999._dp ) THEN 
          cu2d(n,i)=odelta_time*(cu2d(n,i)-cu2d_tm1(n,i))
         ELSE
          cu2d(n,i)=-999._dp
         END IF
       END DO 
     END DO
    END IF


 END SUBROUTINE collect_cu2d_wrap



!-----------------------------------------------------------------

 SUBROUTINE collect_cu2d( cu2d, cur, nlev, gmax, igot, isen, &
                         sending_PEs, istart, iend, ivar, opt )
! collect points along curtain
!  default: v2d(ivar), but also options for lat, lon, tm1

  USE mo_mpi,              ONLY: p_pe, p_parallel_io, p_global_comm, &
                                  p_nprocs, p_io,  p_send, p_recv
  USE mo_decomposition,    ONLY: gdc => global_decomposition


  IMPLICIT NONE  
  
   TYPE(curtain_type), INTENT (IN ), POINTER :: cur

   INTEGER, INTENT(IN) :: nlev, gmax
   LOGICAL, INTENT(IN) :: igot, isen

   INTEGER, DIMENSION( p_nprocs ), INTENT(IN) :: sending_PEs, &
                                                 istart,      &
                                                 iend

   REAL(dp), DIMENSION(nlev, gmax), INTENT(INOUT ) :: cu2d

   ! get v3d(ivar)
   INTEGER, OPTIONAL, INTENT(IN) :: ivar

   ! special options: curtain from last time step
   CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: opt

   

!local
   REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: sendbuf2d, recvbuf2d
  
   INTEGER :: i, in,  n, np ,nr, ie, is, je, pe

   CHARACTER(LEN=3) :: lopt
   
   LOGICAL :: got_opt, got_ivar
  
      
   got_opt = PRESENT( opt )
   got_ivar = PRESENT( ivar )

   IF ( got_opt ) THEN
      lopt=opt
   ELSE
      lopt='var'
   END IF

   CALL collect_arg_check(got_opt, got_ivar, lopt, 'collect_cu2d')
  

   ALLOCATE(sendbuf2d(nlev,cur%pamax))
   ALLOCATE(recvbuf2d(nlev,cur%pamax))  

         sendbuf2d(:,:) = -999._dp
         recvbuf2d(:,:) = -999._dp        

! fill send buffer for all PEs which have part of the curtain
!  this includes pe_io
!  also remove uppermost level from variables dimensioned klevp1 
    IF (  igot ) THEN
      SELECT CASE (TRIM(lopt))
      CASE('tm1')
        IF ( v3d(ivar)%lnlevp1 ) THEN
          DO n=1,  cur%pmax  
            DO i=2,nlev+1
             sendbuf2d(i-1,n) = v3d(ivar)%tm1(i,n)
            END DO
          END DO
        ELSE
          DO n=1,  cur%pmax  
            DO i=1,nlev
             sendbuf2d(i,n) = v3d(ivar)%tm1(i,n)
            END DO
          END DO
        END IF 
      CASE ('var')
        IF ( v3d(ivar)%lnlevp1 ) THEN
          DO n=1,  cur%pmax
              np = cur%myjl(n)
              nr = cur%myjrow(n)
            DO i=2,nlev+1
              sendbuf2d(i-1,n) = v3d(ivar)%ptr3d(np,i,nr)
            END DO 
          END DO
         ELSE
           DO n=1,  cur%pmax
               np = cur%myjl(n)
               nr = cur%myjrow(n)
             DO i=1,nlev
               sendbuf2d(i,n) = v3d(ivar)%ptr3d(np,i,nr)
             END DO 
           END DO
         END IF
      END SELECT
    END IF 


! send
    IF ( isen ) THEN
       CALL p_send(sendbuf2d, p_io, 1400+p_pe, comm=p_global_comm) 
    END IF


! receive, assemble, and write to output file
   
    IF ( p_parallel_io ) THEN
      DO i = 1, p_nprocs
           recvbuf2d(:,:) = -999._dp
           pe    = gdc(i)%pe
           IF ( ANY(sending_PEs .EQ. pe)  ) THEN

            CALL p_recv(recvbuf2d, pe, 1400+pe, comm=p_global_comm)

           END IF
           IF ( pe .EQ. p_io ) THEN
              recvbuf2d(:,:) = sendbuf2d(:,:)
           END IF

           is = istart(i)
           ie = iend(i)
           je = ie - is +1 ! pmax

           IF ( is .GT. 0 ) THEN
             DO in=1,nlev
               cu2d(in,is:ie) = recvbuf2d(in,1:je)
             END DO
           END IF
       END DO   

     END IF

!!$      if ( p_parallel_io ) then
!!$       write(655,*) "---------- ",i
!!$       write(655,*) " is, ie, je ", is, ie, je 
!!$       write(655,*) recvbuf1d(:)
!!$       write(655,*) "======================================"
!!$       write(655,*) cu1d(:)
!!$       write(655,*) "======================================"
!!$     end if

 END SUBROUTINE collect_cu2d

!---------------------------------------------------------------------
 SUBROUTINE cosp_offline_nml
! read namelist only, init_cosp_offline (below) is called from cosp_offline
  USE mo_mpi,         ONLY: p_parallel, p_parallel_io, p_bcast, p_io
  USE mo_namelist,    ONLY: position_nml, POSITIONED, open_nml

  IMPLICIT NONE 

      INTEGER :: ierr, inml, iunit

      NAMELIST /cospofflctl/ locospoffl, offl2dout

     IF (p_parallel_io) THEN
         inml = open_nml ('namelist.echam')
         iunit = position_nml ('COSPOFFLCTL', inml, status=ierr)
         SELECT CASE (ierr)
         CASE (POSITIONED)
            READ (iunit, cospofflctl)          
         END SELECT
      END IF


     IF (p_parallel) THEN
         CALL p_bcast (locospoffl, p_io)
         CALL p_bcast (offl2dout, p_io)
     END IF


 END SUBROUTINE cosp_offline_nml
!------------------------------------------------------------------------

 SUBROUTINE init_cosp_offline

   USE mo_cosp_simulator, ONLY: locosp, Lisccp_sim
   USE mo_exception,      ONLY: message_text, finish
   USE mo_time_control,   ONLY: next_date, get_date_components

   IMPLICIT NONE

   INTEGER :: year, month, day, hour, minute, second

!local

    IF ( .NOT. ( locosp .AND. Lisccp_sim  ) ) THEN
        WRITE(message_text,*) 'mo_cosp_offline needs locosp and lisccp_sim = true'
        CALL finish('mo_cosp_offline ', message_text)
    END IF 

  CALL get_date_components(next_date, year, month, day, hour, minute, second)

 ! read and sample satellite track 
  CALL read_track(year, month, day, hour, minute, second)
 ! determine PE and coordinates (row, proma) of points along track  
  CALL get_indices
 ! initialize information on output variables
  CALL init_svars 
 ! set pointers to variables
  CALL init_pointers
 ! open output file 
  CALL open_output_file

  module_is_initialized = .TRUE.
  cur_month = month

 END SUBROUTINE  init_cosp_offline

!-----------------------------------------------------------------------------------------

 SUBROUTINE read_track(year, month, day, hour, minute, second)
  ! PURPOSE:
  !  determine satellite orbit coordinates (slon, slat) and times (stime) 
  !  during each model time step 
  !
  !  use satellite curtain centered around present time (see point 2 below)
  !
  ! SEQUENCE:
  !  (1) read time info for one month from orbital data file 
  !    and convert to seconds since last day of previous month (dsf)
  !  (2) calculate sampling times (i.e. model time steps)  for the entire 
  !      next month, again in seconds since last day of the previous month  
  !  (3) find times closest to model time step in orbital data (indarr)
  !  (4) set up slice (curtain) for each timestep 
  !  (5) determine satellite orbit coordinates (slon, slat)  and times (stime)
  !      during model time steps (2d-arrays!) 
  !
  

   USE mo_netcdf,           ONLY: nf_check
   USE mo_time_control,     ONLY: delta_time
   USE mo_mpi,              ONLY: p_parallel, p_parallel_io, p_bcast, p_io
   USE mo_exception,        ONLY: message, message_text, finish

   IMPLICIT NONE

   INTEGER, INTENT (IN ) :: year, month, day, hour, minute, second

!local
   ! model 
   REAL(dp), ALLOCATABLE, DIMENSION(:) ::  dsm !  seconds since last day of previous month model 

   ! orbit
   INTEGER, ALLOCATABLE, DIMENSION(:) :: fyear0, fmonth0
   INTEGER, ALLOCATABLE, DIMENSION(:) :: fyear, fmonth, fday, fhour, fmin
   REAL, ALLOCATABLE, DIMENSION(:) :: fsec
   REAL, ALLOCATABLE, DIMENSION(:) :: flon, flat

   REAL, ALLOCATABLE, DIMENSION(:) :: dsf ! seconds since last day of previous month orbit

   ! indices of orbit times closest to model time steps
   INTEGER, ALLOCATABLE, DIMENSION(:) :: indarr 
   
   
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: indarr2d
   
   INTEGER :: nmax1, nmax, nlmax 
   INTEGER :: mdays
   INTEGER :: ncid, varid, nloctot, nlength, dimid, i, j
   INTEGER :: tfirst, tlen, mfirst, mlast
   INTEGER :: off, offs, offe
   REAL(dp) :: d, dn
 
   INTEGER, DIMENSION(1) :: starta, counta
  
  WRITE(message_text,* )  'read A-TRAIN orbit data'
  CALL message('mo_cosp_offline ', message_text )

  p_parallel_io_if_1:  IF (p_parallel_io) THEN


  !-----------------------------------------------------------------------------
  ! (1) get orbital time data for one month and convert to
  !     seconds since last day of previous month
  !-----------------------------------------------------------------------------
      CALL nf_check(nf_open('cloudsat_orbit_2008.nc', nf_nowrite, ncid))

      CALL nf_check(nf_inq_dimid(ncid, 'location', dimid))
      CALL nf_check(nf_inq_dimlen(ncid, dimid, nloctot))

      ALLOCATE( fmonth0( nloctot ), fyear0( nloctot ) )

       CALL nf_check(nf_inq_varid(ncid, 'year', varid))
       CALL nf_check(nf_get_var_int(ncid, varid, fyear0))
       CALL nf_check(nf_inq_varid(ncid, 'month', varid))
       CALL nf_check(nf_get_var_int(ncid, varid, fmonth0))


        tlen = 0
        DO i=1,nloctot
          IF ( year .EQ. fyear0(i) .AND. month .EQ. fmonth0(i) ) THEN
            IF (   tlen .EQ. 0 ) THEN
               tfirst = i
            END IF
            tlen = tlen + 1
          END IF
        END DO
  
        ! need offset since we want curtains to be centered around time step
        off = INT(0.01_dp + delta_time/86400._dp * REAL(tlen,dp)/28._dp) 

        offs = MIN(off, tfirst +1)
        offe = MIN(2*off, nloctot - tlen)

        starta(1) = tfirst - offs
        nlength  = tlen + offe
   
        ! indices of first/last of month in fday, fhour, fmin, fsec (see below) 
        mfirst = offs + 1
        mlast =  tlen + offs 

       IF ( tlen .LT. 1 ) THEN
        WRITE(message_text,*) 'no points in orbital file for this month'
        CALL finish('mo_cosp_offline ', message_text)
       END IF

       ALLOCATE( fyear(nlength) )
       ALLOCATE( fmonth(nlength) )
       ALLOCATE( fday(nlength) )
       ALLOCATE( fhour(nlength) )
       ALLOCATE( fmin(nlength) )
       ALLOCATE( fsec(nlength) )

       ALLOCATE( dsf(nlength) )
 

       counta(1) = nlength
       CALL nf_check(nf_inq_varid(ncid, 'year', varid))
       CALL nf_check(nf_get_vara_int(ncid, varid, starta, counta, fyear))
       CALL nf_check(nf_inq_varid(ncid, 'month', varid))
       CALL nf_check(nf_get_vara_int(ncid, varid, starta, counta, fmonth))
       CALL nf_check(nf_inq_varid(ncid, 'day', varid))
       CALL nf_check(nf_get_vara_int(ncid, varid, starta, counta, fday))
       CALL nf_check(nf_inq_varid(ncid, 'hour', varid))
       CALL nf_check(nf_get_vara_int(ncid, varid, starta, counta, fhour))
       CALL nf_check(nf_inq_varid(ncid, 'minute', varid))
       CALL nf_check(nf_get_vara_int(ncid, varid, starta, counta, fmin))
       CALL nf_check(nf_inq_varid(ncid, 'second', varid))
       CALL nf_check(nf_get_vara_real(ncid, varid, starta, counta, fsec))


       ! seconds since last day of previous month for orbit 
        DO i = 1, mfirst-1
          dsf(i) =   3600._dp * REAL(fhour(i), dp) &
                     + 60._dp * REAL(fmin(i), dp) +  REAL(fsec(i),dp)
        END DO

        DO i = mfirst, mlast 
          dsf(i) =  86400._dp * REAL(fday(i), dp) &
                    + 3600._dp * REAL(fhour(i), dp) &
                    + 60._dp * REAL(fmin(i), dp) + REAL(fsec(i),dp)
        END DO

        DO i=mlast+1, nlength
          dsf(i) =  86400._dp * REAL(fday(mlast)+1, dp) &
                    + 3600._dp * REAL(fhour(i), dp) &
                    + 60._dp * REAL(fmin(i), dp) + REAL(fsec(i),dp)
        END DO
 
    
      ! number of days in the month (from orbital file)
      mdays=MAXVAL(fday(mfirst:mlast))

 !----------------------------------------------------------------------------------
 ! (2) calculate model sampling times in seconds since last day of previous month 
 !     for one month based on model time step
 !--------------------------------------------------------------------------------

      ! number of model time steps remaining in this month 
       nmax1=(mdays-day+1)*INT(86400._dp/delta_time) + 2
       ALLOCATE(dsm(nmax1))
       ALLOCATE(indarr(nmax1))

     ! seconds since last day of previous month for model time steps
     ! offset by half time step, so that sampling will be centered
     ! at present step

      dsm(1) =  86400._dp*day + 3600._dp*hour + 60._dp*minute + REAL(second,dp) - 0.5_dp * delta_time
      DO i=2, nmax1  
         dsm(i) =  dsm(i-1) + delta_time
      END DO

  !---------------------------------------------------------------------------------
  ! (3) find indices of times closest to model time steps in orbital data (indarr)
  !---------------------------------------------------------------------------------
  
        i=1
       DO j=1,  nmax1 
           d  = 9999999._dp      
           dn = ABS( dsm(j) - dsf(i) )
          DO  WHILE ( dn .LE. d .AND. i .LT. nlength )
             d = dn
             i=i+1
             dn =  ABS( dsm(j) - dsf(i) )     
          END DO
              i=MAX(i-1,1)
           indarr(j) = i
       END DO


  !---------------------------------------------------------------------------------
  ! (4) set up slice (curtain) for each timestep (indarr2d)
  !---------------------------------------------------------------------------------

  ! maximum number of orbit locations in one time step
       
       nlmax=0
     DO i = 1,  nmax1-1
       nlmax=MAX(indarr(i+1)-indarr(i),nlmax)
     END DO

     ALLOCATE(indarr2d(nlmax,nmax1-1))
   
     indarr2d(:,:) = -999

     DO i=1, nmax1 - 1
         j=1
       DO WHILE ( j .LE.  indarr(i+1)-indarr(i) )  
         indarr2d(j,i) =  indarr(i)+j-1
         j=j+1
       END DO
     END DO

      nmax=nmax1-1 


! debug++

!!$open(1112,file="debu.out")
!!$write(1112,*) "delta_time", delta_time
!!$write(1112,*) "day, hour, minute -- delta_time ", day, hour, minute, "--", delta_time 
!!$do i=1, nmax-1 !!!!!
!!$write(1112,*) "-------------", i
!!$    do j=1,nlmax
!!$     !write(1112,*)  indarr2d(j,i), fday(indarr(i)), fhour(indarr(i)), fmin(indarr(i)), fsec(indarr(i))
!!$     if (  indarr2d(j,i) .ne. -999 ) then
!!$     write(1112,*)  j, indarr2d(j,i), ' -- ',fday(indarr2d(j,i)), fhour(indarr2d(j,i)), fmin(indarr2d(j,i)), fsec(indarr2d(j,i))
!!$     else
!!$      write(1112,*)  indarr2d(j,i)
!!$     end if 
!!$    end do
!!$end do
!!$close(1112)
! debug --  
     
   END IF p_parallel_io_if_1

 !------------------------------------------------------------------------------
 !  (5) determine satellite orbit coordinates (slon, slat) and times (stime) 
 !      at model time steps 
 !-----------------------------------------------------------------------------
     IF (p_parallel) THEN
       CALL p_bcast(nmax, p_io)
       CALL p_bcast(nlmax, p_io)
     END IF

 
       ALLOCATE(slon(nlmax,nmax))
       ALLOCATE(slat(nlmax,nmax))
       ALLOCATE(stime(nlmax,nmax))
    

       slon(:,:) = -999._dp
       slat(:,:) = -999._dp
       stime(:,:) = -999._dp

  p_parallel_io_if_2:  IF (p_parallel_io) THEN    
    ! read lon/lat from orbit file

       ALLOCATE( flon(nlength) )
       ALLOCATE( flat(nlength) )

       CALL nf_check(nf_inq_varid(ncid, 'lon', varid))
       CALL nf_check(nf_get_vara_real(ncid, varid, starta, counta, flon))
       CALL nf_check(nf_inq_varid(ncid, 'lat', varid))
       CALL nf_check(nf_get_vara_real(ncid, varid, starta, counta, flat))



      DO i=1,nmax 
         DO j=1,nlmax
           IF ( indarr2d(j,i) .NE. -999 ) THEN
             slon(j,i)  = REAL(flon(indarr2d(j,i)), dp)
            IF ( slon(j,i) .LT. 0._dp ) THEN 
              slon(j,i) = 360._dp +  slon(j,i)
            END IF
             slat(j,i)  = REAL(flat(indarr2d(j,i)), dp)
             ! day in format %Y%m%d.%f
             stime(j,i)     = REAL(ABS(fyear(indarr2d(j,i))),dp)*10000._dp+ REAL(fmonth(indarr2d(j,i)),dp)*100._dp &
                           + REAL(fday(indarr2d(j,i)),dp) + (REAL(fhour(indarr2d(j,i)),dp)*3600._dp  &
                           + REAL(fmin(indarr2d(j,i)),dp)*60._dp + REAL(fsec(indarr2d(j,i)),dp))/86400._dp
           END IF
         END DO
       END DO
     
        DEALLOCATE( indarr )
        DEALLOCATE( indarr2d )
        DEALLOCATE ( flon, flat )
        DEALLOCATE(dsm)  
        DEALLOCATE(fyear0, fmonth0) 
        DEALLOCATE(fyear, fmonth) 
        DEALLOCATE(fday, fhour, fmin, fsec, dsf)
 
   END IF p_parallel_io_if_2


 ! inform others
   CALL p_bcast(slon, p_io)   
   CALL p_bcast(slat, p_io)   
   CALL p_bcast(stime, p_io)   


 END SUBROUTINE read_track
!----------------------------------------

 SUBROUTINE cosp_offline_finalize

    USE mo_mpi,              ONLY: p_parallel_io
    USE mo_netcdf,           ONLY: nf_check
    USE mo_exception,        ONLY: message, message_text

    IMPLICIT NONE

   IF ( p_parallel_io ) THEN

    CALL nf_check( NF_CLOSE(ofile%ncid))
 
    WRITE(message_text,*) 'mo_cosp_offline: closed file', ofile%name
    CALL message('', message_text)

  END IF

 END SUBROUTINE cosp_offline_finalize 

!-------------------------------------
 SUBROUTINE init_svars

  USE mo_exception,          ONLY: finish
  USE mo_memory_gl,          ONLY: gl
  USE mo_memory_g3b,         ONLY: g3b
  USE mo_memory_g1a,         ONLY: g1a
  USE mo_cosp_simulator,     ONLY: scosp


  IMPLICIT NONE

! initialize information on output variables
  INTEGER ::  i, j

   i=1  
  v2d(i)%iname='aps'
  v2d(i)%oname='psfc'
  v2d(i)%stream => g3b

   i=i+1
  v2d(i)%iname='u10'
  v2d(i)%oname='u_wind'
  v2d(i)%stream => g3b

   i=i+1
  v2d(i)%iname='v10'
  v2d(i)%oname='v_wind'
  v2d(i)%stream => g3b

   i=i+1
  v2d(i)%iname='slm'
  v2d(i)%oname='landmask'
  v2d(i)%stream => g3b

    i=i+1
  v2d(i)%iname='geosp'
  v2d(i)%oname='geosp'
  v2d(i)%stream => g3b

   i=i+1
  v2d(i)%iname='cosp_sunlit'
  v2d(i)%oname='sunlit'
  v2d(i)%stream => scosp


   IF ( i .GT. nvar_0d ) THEN
     CALL finish('mo_cosp_offline: ', 'increase nvar_0d')
   END IF 


   IF ( .NOT. module_is_initialized ) THEN
     DO i=1,nvar_0d 
       ALLOCATE( v2d(i)%tm1( ofile%gmax ))
     END DO
   END IF

   DO i=1,nvar_0d 
     v2d(i)%tm1(:)=0._dp
     v2d(i)%laccu=accu_info(v2d(i)%iname,g3b) 
     v2d(i)%unit=unit_info(v2d(i)%iname,g3b) 
   END DO


! profile variables
      i=1
    v3d(i)%iname="aclc"
    v3d(i)%oname="tca"
    v3d(i)%stream => g3b
   
   ! vmr? mmr? check later, not needed for cfmip
      i=i+1
    v3d(i)%iname="ao3"
    v3d(i)%oname="mr_ozone"
    v3d(i)%stream => g3b

      i=i+1
    v3d(i)%iname="xl"
    v3d(i)%oname="mr_lsliq"
    v3d(i)%stream => gl

      i=i+1
    v3d(i)%iname="xi"
    v3d(i)%oname="mr_lsice"
    v3d(i)%stream => gl
  
        i=i+1
    v3d(i)%iname="q"
    v3d(i)%oname="qv"
    v3d(i)%stream => gl

      i=i+1
    v3d(i)%iname="tm1"
    v3d(i)%oname="T_abs"
    v3d(i)%stream => g1a
 
      i=i+1
    v3d(i)%iname="relhum"
    v3d(i)%oname="rh"
    v3d(i)%stream => g3b

     i=i+1
    v3d(i)%iname="cisccp_tau3d"
    v3d(i)%oname="dtau_s"
    v3d(i)%stream => scosp

    i=i+1
    v3d(i)%iname="cisccp_emi3d"
    v3d(i)%oname="dem_s"
    v3d(i)%stream => scosp


     i=i+1
    v3d(i)%iname="reffl"
    v3d(i)%oname="Reffl"
    v3d(i)%stream => scosp

     i=i+1
    v3d(i)%iname="reffi"
    v3d(i)%oname="Reffi"
    v3d(i)%stream => scosp

     i=i+1
    v3d(i)%iname="fl_lsrain"
    v3d(i)%oname="fl_lsrain"
    v3d(i)%stream => cospoffl

     i=i+1
    v3d(i)%iname="fl_lssnow"
    v3d(i)%oname="fl_lssnow"
    v3d(i)%stream => cospoffl

     i=i+1
    v3d(i)%iname="fl_ccrain"
    v3d(i)%oname="fl_ccrain"
    v3d(i)%stream => cospoffl

     i=i+1
    v3d(i)%iname="fl_ccsnow"
    v3d(i)%oname="fl_ccsnow"
    v3d(i)%stream => cospoffl

     i=i+1
    v3d(i)%iname="geom1"
    v3d(i)%oname="geom"
    v3d(i)%stream => cospoffl

     i=i+1
    v3d(i)%iname="geohm1"
    v3d(i)%oname="geohm"
    v3d(i)%stream => cospoffl

     i=i+1
    v3d(i)%iname="pfull"
    v3d(i)%oname="pfull"
    v3d(i)%stream => cospoffl

     i=i+1
    v3d(i)%iname="phalf"
    v3d(i)%oname="phalf"
    v3d(i)%stream => cospoffl

    
    nvar_1d = i

   DO j = 1, nvar_1d
     v3d(j)%laccu  =  accu_info(v3d(j)%iname,v3d(j)%stream)
     v3d(j)%unit   =  unit_info(v3d(j)%iname,v3d(j)%stream)
   END DO


   IF ( nvar_1d .GT. nvar_1d_max ) THEN
     CALL finish('mo_cosp_offline: ', 'increase nvar_1d_max')
   END IF 


 END SUBROUTINE init_svars
!---------------------------------------------------
 SUBROUTINE init_pointers

   USE mo_memory_g3b
   USE mo_memory_base,      ONLY: get_stream_element
   USE mo_memory_base,      ONLY: get_stream_element_info,  memory_info
   USE mo_exception,        ONLY: message_text, finish
   USE mo_control,          ONLY: nlev, nlevp1


   IMPLICIT NONE
! local

   INTEGER :: i
   CHARACTER(LEN=150) :: cname
   TYPE(t_stream), POINTER :: stream
   TYPE (memory_info) :: info

! set pointers 
  DO i = 1, nvar_0d 
     cname=TRIM(v2d(i)%iname)
     stream => v2d(i)%stream
     CALL get_stream_element (stream, cname, v2d(i)%ptr2d)  
  END DO
 
  DO i = 1, nvar_1d 
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
         WRITE(message_text,*) TRIM(v3d(i)%iname), &
           ': var not found or vert grid not supported ', info%klev
         CALL finish('mo_cosp_offline: ', message_text)
     END IF

     ! allocate v3d tm1
     IF  ( v3d(i)%laccu ) THEN
       IF (info%klev .GT. 0 ) THEN
         IF ( .NOT.  module_is_initialized  ) THEN 
           ALLOCATE ( v3d(i)%tm1(  info%klev ,ofile%gmax ) )
         END IF  
       ELSE
         WRITE(message_text,*) TRIM(v3d(i)%iname), ' info%klev not assigned'
         CALL finish('mo_cosp_offline: ', message_text)
       END IF
     END IF
  END DO

 END SUBROUTINE init_pointers


!---------------------------------------------------


 SUBROUTINE get_indices

! called during initialization
! determine PE and coordinates (row, proma) along curtain 

  USE mo_exception,     ONLY: message, message_text

  USE mo_geoloc,        ONLY: philon_2d, philat_2d
  USE mo_decomposition, ONLY: ldc=>local_decomposition, &
                              gdc => global_decomposition
  USE mo_mpi,           ONLY: p_parallel, p_parallel_io, p_bcast, &
                              p_io, p_pe, p_nprocs, p_send, p_recv, &
                              p_global_comm 


  IMPLICIT NONE
!local

 ! the following include multiple locations in the same grid box
  INTEGER, ALLOCATABLE, DIMENSION(: ) :: onpe0
  INTEGER, ALLOCATABLE, DIMENSION(:,: ) :: onpe1          
  INTEGER, ALLOCATABLE, DIMENSION(:,: ) :: myjl1, myjrow1

 ! keep these, the others are duplicates
  LOGICAL, ALLOCATABLE, DIMENSION(:,: ) :: lkeep
  
 ! number of points in curtain segment (after removing duplicates)
  INTEGER,  ALLOCATABLE, DIMENSION(: )   :: nlkeep1 

 ! number of points in entire curtain (after removing duplicates) 
  INTEGER,  ALLOCATABLE, DIMENSION(: )   :: ngkeep  

 ! number of time steps from orbit file
  INTEGER :: nmax

 ! number of points from orbit file
  INTEGER :: nlmax

 ! number of points in the largest curtain in the whole month
  INTEGER :: ngkeepmax

 ! number of points in the largest curtain segment
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nlkeepmax

 ! start index for assembling curtain segments into curtain
  INTEGER,  ALLOCATABLE, DIMENSION(: )  :: nlstart 

 ! number of points that is send "nl send" 
  INTEGER :: nlsend

  
  ! other variables
  REAL(dp), ALLOCATABLE, DIMENSION(:, :) :: recvbuf 
  REAL(dp), ALLOCATABLE, DIMENSION(: )   :: sendbuf
  INTEGER,  ALLOCATABLE, DIMENSION(: )   :: irecvbuf 

  INTEGER :: i, n, jl, nproma, jrow, np, ngpblks, j
  INTEGER :: nn, jl1, jr1, nlk, pe, isu, nlkeepmax1, nlstart1
  REAL(dp) :: dl, dlt

  nlmax = SIZE(slat, 1)
  nmax  = SIZE(slat, 2)

  ALLOCATE( onpe0(nlmax) )
  ALLOCATE ( onpe1(nlmax, nmax) )
  ALLOCATE ( myjl1(nlmax, nmax) )
  ALLOCATE ( myjrow1(nlmax, nmax ) )
  ALLOCATE ( lkeep(nlmax, nmax) )

  ALLOCATE(nlkeep1(nmax))
  ALLOCATE(ngkeep(nmax))
  ALLOCATE(nlkeepmax(nmax))
  ALLOCATE(nlstart(nmax))

  ALLOCATE(recvbuf(nlmax, p_nprocs))
  ALLOCATE(sendbuf(nlmax))
  ALLOCATE(irecvbuf( p_nprocs ))

 
  ! initialize
   ngkeep(:)=-999
   ngkeepmax=-999
   nlstart(:) = -999
   myjl1(:,:) = -999
   myjrow1(:,:) = -999
   lkeep(:,:) = .FALSE.
   onpe1(:,:) = -999    

  WRITE(message_text,* )  'calculate satellite curtains for ', nmax, ' steps'
  CALL message('', message_text )


  step_loop1: DO i=1, nmax

   recvbuf(:,:) = -999._dp
   sendbuf(:) = -999._dp
   onpe0(:)  = -999
   irecvbuf(:) = 0 
   ngpblks = ldc % ngpblks

    IF ( MOD(i,20) .EQ. 0 ) THEN
      WRITE(message_text,* ) '     step ', i, ' of ',   nmax
      CALL message('  ', message_text )
    END IF

! determine distance from closest grid point on each PE
   DO n=1,nlmax
      dl=999999._dp
     DO jrow = 1, ngpblks        
       IF ( jrow == ngpblks ) THEN
         nproma = ldc% npromz
       ELSE
         nproma = ldc% nproma
       END IF

       DO jl=1,nproma      
           dlt = (philat_2d(jl,jrow) - slat(n,i))**2 &
             + ( philon_2d(jl,jrow) - slon(n,i))**2 
          IF ( dlt .LT.  dl ) THEN
             dl = dlt
             myjl1(n,i) = jl
             myjrow1(n,i) = jrow
          END IF     
       END DO
     END DO

       sendbuf(n) =  dl
      
    END DO
 
! communicate all the distances to the IO PE
   CALL p_gather_real_1d2d (sendbuf, recvbuf, p_io)
 
! let IO PE decide which PE contains closest grid point
    IF (p_parallel_io) THEN 
      DO n=1, nlmax
          dl=999999._dp
        DO np=1,p_nprocs
          IF ( recvbuf(n,np) .GT. 0._dp .AND. recvbuf(n,np) .LE. dl .AND.  &
                     recvbuf(n,np) .LT.  999999._dp ) THEN
            dl=recvbuf(n,np)
            onpe0(n)=np-1
          END IF 
         END DO
      END DO
    END IF

  

! let other PEs know the decisions
     IF (p_parallel) THEN
         CALL p_bcast (onpe0, p_io)
     END IF

! remove multiple points in the same grid box
      nn = 0
      jl1 = -999
      jr1 = -999
   DO n=1, nlmax
     IF ( onpe0(n) .EQ. p_pe ) THEN


       IF (  myjl1(n,i) .NE. jl1 .OR.  myjrow1(n,i) .NE. jr1 ) THEN
        IF ( slon (n,i) .NE. -999 ) THEN
          lkeep(n,i) = .TRUE.
          nn=nn+1
        END IF
       END IF
         jl1 =  myjl1(n,i)
         jr1 = myjrow1(n,i)
      END IF 


   END DO


! save number of points for this pe and this time step 
!  (length of the curtain segment on the local pe after removing 
!   duplicates)

    nlkeep1(i) = nn


! sum up the number of points from all the individual segments 
    CALL p_gather_int_0d1d (nn, irecvbuf, p_io)
     isu =  SUM(irecvbuf)
     ngkeepmax = MAX(isu,ngkeepmax)
     nlkeepmax1 = MAXVAL(irecvbuf)      

!  ... and let everybody know 
     IF (p_parallel) THEN
         CALL p_bcast (isu, p_io)
         CALL p_bcast (ngkeepmax, p_io)
         CALL p_bcast (nlkeepmax1, p_io)
     END IF

     ngkeep(i) = isu
     nlkeepmax(i)  = nlkeepmax1

! later, the curtain segments will need to be assembled into one entire curtain
!  so, here we determine a starting index where the curtain
!     segment will be inserted into the entire curtain

   IF ( p_parallel_io ) THEN
     ! first calculate nlstart for all the other PEs 
     !  and let them know the result

        nlstart1=nlkeep1(i) ! reserve first part for io PE (if necessary)
      DO j=1,  p_nprocs
          pe  = gdc(j)%pe
         IF ( pe .NE. p_io ) THEN
           IF ( irecvbuf(j) .GT. 0 ) THEN 
              nlsend = nlstart1 + 1
              nlstart1 =  nlstart1 + irecvbuf(j)
           ELSE
             nlsend = -999
           END IF

         CALL p_send(nlsend, pe, 500+pe,  comm=p_global_comm)

         END IF
      END DO
   ELSE

         CALL p_recv(nlstart1, p_io, 500+p_pe,  comm=p_global_comm )  

   END IF

    ! now also set nlstart for p_io
    IF ( p_parallel_io ) THEN
      IF ( nlkeep1(i) .GT. 0 ) THEN 
        nlstart1=1
      ELSE
        nlstart1 = -999
      END IF
    END IF

    nlstart(i) =  nlstart1   
 
    DO n=1, nlmax
     onpe1(n,i) = onpe0(n)
    END DO

   END DO step_loop1


  DEALLOCATE ( recvbuf )
  DEALLOCATE ( irecvbuf )
  DEALLOCATE ( sendbuf )


   ofile%gmax = ngkeepmax 


! save final coordinates after removing multiple points in the same gridbox

    IF (  module_is_initialized ) THEN
      DEALLOCATE ( cu )
    END IF


    ALLOCATE( cu(nmax))

    DO i=1,nmax
     
      nlk=MAX(nlkeep1(i),1)

      ALLOCATE ( cu(i)%myjl(nlk) )
      ALLOCATE ( cu(i)%myjrow(nlk ) )
      ALLOCATE ( cu(i)%mytime(nlk) )
      ALLOCATE ( cu(i)%on_this_pe(nlk ) )

      cu(i)%on_this_pe(:) = .FALSE.
      cu(i)%myjl(:) = -999
      cu(i)%myjrow(:) = -999
      cu(i)%mytime(:) = -999._dp

      cu(i)%pmax = nlkeep1(i)
      cu(i)%gmax = ngkeep(i)
      cu(i)%pamax =  nlkeepmax(i)
      cu(i)%start = nlstart(i)

    END DO


  DO i=1, nmax 
      nn=1
    DO n=1, nlmax
      IF ( onpe1(n,i) .EQ. p_pe ) THEN
        IF ( lkeep(n,i) ) THEN
          cu(i)%on_this_pe(nn) = .TRUE.
          cu(i)%myjl(nn)    = myjl1(n,i)
          cu(i)%myjrow(nn)  = myjrow1(n,i)
          cu(i)%mytime(nn) = stime(n,i)
          nn=nn+1
        END IF 
      END IF
    END DO
  END DO 


!!$!debug ++
!!$nn=1500 + p_pe
!!$ write(nn,*)  "-----------------------------------------------"
!!$  do i=1, 2!!!nmax
!!$     write(nn,*)  " i, cu(i)%on_this_pe(1) ",   i, cu(i)%on_this_pe(1),  nlkeep1(i)
!!$   do n=1, nlkeep1(i)
!!$     if (cu(i)%on_this_pe(n) ) then
!!$      write(nn,*) i, p_pe, " -- ", cu(i)%myjl(n), cu(i)%myjrow(n), "  -- ",  &
!!$         philon_2d(cu(i)%myjl(n), cu(i)%myjrow(n)), philat_2d(cu(i)%myjl(n), cu(i)%myjrow(n))
!!$
!!$     end if 
!!$   end do
!!$  end do  



  DEALLOCATE ( lkeep, nlkeep1, ngkeep, nlstart, nlkeepmax )
  DEALLOCATE ( myjl1, myjrow1, onpe1 ) 
  DEALLOCATE ( slon, slat, stime )



 END SUBROUTINE get_indices


!---------------------------------------------------------------------
 SUBROUTINE open_output_file

   USE mo_exception,          ONLY: message, message_text
   USE mo_netcdf,             ONLY: nf_check
   USE mo_mpi,                ONLY: p_parallel_io
   USE mo_control,            ONLY: nlev
   USE mo_time_control,       ONLY: next_date, timelabel_type,  str_date
   USE mo_filename,           ONLY: out_expname

   IMPLICIT NONE

!local

  INTEGER :: i, ncid, varid, year_varid, month_varid, day_varid 
  INTEGER :: hour_varid, minute_varid, second_varid, tim_varid
  INTEGER :: lev_dimid, point_dimid, hydro_dimid
  
  INTEGER, DIMENSION(1) :: dimids1d 
  INTEGER, DIMENSION(2) :: dimids2d 
  CHARACTER(LEN=150) :: dumname
  CHARACTER(len=19) :: time_label

  pe_io_if: IF ( p_parallel_io ) THEN    

    time_label = TRIM(str_date(timelabel_type,next_date))

    ! determine filename
         dumname =TRIM(out_expname)//'_'//TRIM(time_label)//'_cfOff.nc'
        
     ! open output file
        CALL nf_check( NF_CREATE(TRIM(dumname), NF_CLOBBER, ncid))
        ofile%name =  dumname

        WRITE(message_text,*) ' Open cfOff file ',  TRIM(dumname), ncid
        CALL message('', message_text )


     ! declare dimensions
       CALL nf_check( NF_DEF_DIM(ncid, 'point', NF_UNLIMITED, point_dimid) )
       CALL nf_check( NF_DEF_DIM(ncid, 'level', nlev, lev_dimid) )
       CALL nf_check( NF_DEF_DIM(ncid, 'hydro', 2, hydro_dimid) )

       
       ofile% point_dimid = point_dimid
       ofile% lev_dimid   = lev_dimid
       ofile% hydro_dimid = hydro_dimid


       dimids1d=(/ point_dimid /)


      !declare time variables

       CALL nf_check( NF_DEF_VAR (ncid, 'year'  , NF_INT, 1, dimids1d, year_varid ) )
       CALL nf_check( NF_DEF_VAR (ncid, 'month' , NF_INT, 1, dimids1d, month_varid ) )
       CALL nf_check( NF_DEF_VAR (ncid, 'day'   , NF_INT, 1, dimids1d, day_varid ) )
       CALL nf_check( NF_DEF_VAR (ncid, 'hour',   NF_INT, 1, dimids1d, hour_varid ) )
       CALL nf_check( NF_DEF_VAR (ncid, 'minute', NF_INT, 1, dimids1d, minute_varid ) )
       CALL nf_check( NF_DEF_VAR (ncid, 'second', NF_REAL, 1, dimids1d, second_varid ) )
       CALL nf_check( NF_DEF_VAR (ncid, 'orbtime', NF_REAL, 1, dimids1d, tim_varid ) )

       ofile% year_varid    =  year_varid
       ofile% month_varid   =  month_varid
       ofile% day_varid     =  day_varid 
       ofile% hour_varid    =  hour_varid
       ofile% minute_varid  =  minute_varid
       ofile% second_varid  =  second_varid
       ofile% tim_varid  =  tim_varid
      
       CALL nf_check (NF_PUT_ATT_TEXT(ncid, year_varid   ,'long_name', 4, 'year'  ))
       CALL nf_check (NF_PUT_ATT_TEXT(ncid, month_varid  ,'long_name', 5, 'month' ))
       CALL nf_check (NF_PUT_ATT_TEXT(ncid, day_varid    ,'long_name', 3, 'day'   ))
       CALL nf_check (NF_PUT_ATT_TEXT(ncid, hour_varid   ,'long_name', 4, 'hour'  ))
       CALL nf_check (NF_PUT_ATT_TEXT(ncid, minute_varid ,'long_name', 6,'minute' ))
       CALL nf_check (NF_PUT_ATT_TEXT(ncid, second_varid ,'long_name', 6,'second' ))
       CALL nf_check (NF_PUT_ATT_TEXT(ncid, tim_varid    ,'long_name', 30,'corresponding exact orbit time' ))

    ! declare lon and lat variables
          dumname='lon'
       CALL nf_check( NF_DEF_VAR(ncid, TRIM(dumname), NF_REAL, 1, &
                  dimids1d,  varid) )
       ofile% lon_varid = varid

      
          dumname='lat'
       CALL nf_check( NF_DEF_VAR(ncid, TRIM(dumname), NF_REAL, 1, &
                  dimids1d,  varid) )
       ofile% lat_varid = varid

    ! declare point variables
       DO i=1,nvar_0d

          dumname =  v2d(i)%oname


          CALL nf_check( NF_DEF_VAR(ncid, TRIM(dumname), NF_REAL, 1, &
                     dimids1d,  varid) )
          v2d(i)%varid=varid
  
          CALL put_unit( ncid, varid, v2d(i)%unit , v2d(i)%oname )

       END DO


    ! declare profile variables
       DO i=1,nvar_1d
         
          dimids2d=(/ lev_dimid, point_dimid /)

         dumname =  v3d(i)%oname


         CALL nf_check( NF_DEF_VAR(ncid, TRIM(dumname), NF_REAL, 2, &
                    dimids2d,  varid) )
         v3d(i)%varid=varid

         CALL put_unit( ncid, varid, v3d(i)%unit , v3d(i)%oname )


       END DO

       CALL nf_check( NF_ENDDEF(ncid) )

       ofile%ncid=ncid

        WRITE(message_text,*) 'mo_cosp_offline: opened file  ', ofile%name
        CALL message('', message_text )

   END IF pe_io_if

 END SUBROUTINE open_output_file
!------------------------------------------------------------------
  SUBROUTINE nc_date_write( cur, year, month, day, hour, minute, second )

    USE mo_mpi,              ONLY: p_parallel_io
    USE mo_netcdf,           ONLY: nf_check

   IMPLICIT NONE
 
    INTEGER, INTENT (IN ) :: year, month, day, hour, minute, second
    TYPE(curtain_type), INTENT (IN ), POINTER :: cur


! local
     INTEGER, ALLOCATABLE, DIMENSION(:) :: v1di
     REAL, ALLOCATABLE, DIMENSION(:) :: v1dr ! not dp
     INTEGER:: ncid, varid, start1, gmax, gcmax 

     INTEGER, DIMENSION(1) :: start1d, count1d

   IF ( p_parallel_io ) THEN    
      gmax = ofile%gmax
      gcmax = cur%gmax

      ALLOCATE(v1di(gmax))
      ALLOCATE(v1dr(gmax))

       ncid =  ofile%ncid
       start1 = (ncall-1) * gmax + 1
       start1d = (/ start1 /)
       count1d = (/ gmax /)

!!$     write(333,*)  "-------------------------------------------"
!!$     write(333,*) "hour, min, sec ", hour, minute, second
!!$     write(333,*)  start1, gmax


       varid=  ofile% year_varid
       v1di = curfilli1(year,gmax,gcmax) 
      CALL nf_check(NF_PUT_VARA_INT(ncid, varid, start1d, count1d, v1di ))

       varid=  ofile% month_varid
       v1di =  curfilli1(month,gmax,gcmax) 

     CALL nf_check(NF_PUT_VARA_INT(ncid, varid, start1d, count1d, v1di ))

       varid=  ofile% day_varid
       v1di = curfilli1(day,gmax,gcmax) 
     CALL nf_check(NF_PUT_VARA_INT(ncid, varid, start1d, count1d, v1di ))

       varid=  ofile% hour_varid
       v1di = curfilli1(hour,gmax,gcmax) 
!!$ write(333,*) v1di
     CALL nf_check(NF_PUT_VARA_INT(ncid, varid, start1d, count1d, v1di ))

       varid=  ofile% minute_varid
       v1di = curfilli1(minute,gmax,gcmax) 
     CALL nf_check(NF_PUT_VARA_INT(ncid, varid, start1d, count1d, v1di ))

        varid=  ofile% second_varid
        v1dr=curfillr1(second,gmax,gcmax)
     CALL nf_check(NF_PUT_VARA_REAL(ncid, varid, start1d, count1d, v1dr ))


     DEALLOCATE(v1di, v1dr)
   END IF


  
  END SUBROUTINE nc_date_write
!------------------------------------------------------------------
FUNCTION curfilli1(in, gmax, gcmax)  RESULT (v1)

  IMPLICIT NONE

  INTEGER, INTENT(IN ) :: in, gmax, gcmax

  INTEGER, DIMENSION(gmax) :: v1

  v1(1:gcmax) = in
  IF ( gcmax+1 .LE. gmax ) THEN
   v1( gcmax+1:gmax ) = -999 
  END IF

END FUNCTION curfilli1
!------------------------------------------------------------------
FUNCTION curfillr1(in, gmax, gcmax)  RESULT (v1)

  IMPLICIT NONE

  INTEGER,  INTENT(IN ) :: in
  INTEGER, INTENT(IN ) :: gmax, gcmax

  REAL, DIMENSION(gmax) :: v1


  v1(1:gcmax) = REAL(in)
  IF ( gcmax+1 .LE. gmax ) THEN
   v1( gcmax+1:gmax ) = -999. 
  END IF

END FUNCTION curfillr1
!------------------------------------------------------------------
 SUBROUTINE nc_write_cu1d(  varid, cu1d, gmax, gcmax )

   USE mo_mpi,              ONLY: p_parallel_io
   USE mo_netcdf,           ONLY: nf_check


   IMPLICIT NONE

   INTEGER, INTENT(IN ) ::   varid, gmax, gcmax 

   REAL(dp), DIMENSION( gmax ) :: cu1d
   REAL, DIMENSION( gmax ) :: cu1o


! local
   INTEGER :: ncid,  start1
   INTEGER, DIMENSION(1) :: start1d, count1d


   IF ( p_parallel_io ) THEN 
     ncid =  ofile%ncid

      cu1o(1:gcmax)=REAL(cu1d(1:gcmax)) 
       IF ( gcmax+1 .LE. gmax ) THEN  
        cu1o(gcmax+1:gmax)=-999. 
       END IF

       start1 = (ncall-1) * gmax + 1
       start1d = (/ start1 /)
       count1d = (/ gmax /)

      CALL nf_check(NF_PUT_VARA_REAL( ncid, varid, start1d, count1d, cu1o  )) 

   END IF


 END SUBROUTINE nc_write_cu1d
!------------------------------------------------------------------
 SUBROUTINE nc_write_cu2d(  varid, cu2d, nlev, gmax, gcmax )

   USE mo_mpi,              ONLY: p_parallel_io
   USE mo_netcdf,           ONLY: nf_check


   IMPLICIT NONE

   INTEGER, INTENT(IN ) ::   varid, nlev, gmax, gcmax 

   REAL(dp), DIMENSION( nlev, gmax ) :: cu2d
   REAL, DIMENSION( nlev, gmax ) :: cu2o


! local
   INTEGER :: i, ncid,  start2
   INTEGER, DIMENSION(2) :: start2d, count2d


   IF ( p_parallel_io ) THEN 
     ncid =  ofile%ncid
     
      DO i=1, nlev
        cu2o(i,1:gcmax)=REAL(cu2d(i,1:gcmax))
      END DO 
       IF ( gcmax+1 .LE. gmax ) THEN
         DO i=1, nlev   
           cu2o(i, gcmax+1:gmax)=-999.
         END DO 
       END IF

       start2 = (ncall-1) * gmax + 1
       start2d = (/ 1, start2 /)
       count2d = (/ nlev, gmax /)

      CALL nf_check(NF_PUT_VARA_REAL( ncid, varid, start2d, count2d, cu2o  )) 

   END IF


 END SUBROUTINE nc_write_cu2d
!---------------------------------------------------------------------
  SUBROUTINE construct_stream_cospoffl

    USE mo_memory_base,   ONLY: new_stream, add_stream_element, &
                                default_stream_setting 
    USE mo_control,       ONLY: nlev
    USE mo_time_event,    ONLY: io_time_event

    IMPLICIT NONE

    LOGICAL :: lpost

     IF  (offl2dout .LE. 0) THEN
      CALL new_stream (cospoffl,'cospoffl',filetype=NETCDF, interval=io_time_event(25,'months','last',0))
      lpost=.FALSE.
     ELSE ! testing 2-D offline
      CALL new_stream (cospoffl,'cospoffl',filetype=NETCDF, interval=io_time_event(offl2dout,'minutes','last',0))
      lpost=.TRUE. 
     END IF



    CALL default_stream_setting (cospoffl,  units     = '',        &
                                         lrerun    = .FALSE. ,     &
                                         lpost     = lpost  ,     &
                                         contnorest = .TRUE.,      &
                                         laccu     = .FALSE.      )


    CALL add_stream_element (cospoffl, 'fl_lsrain', cospoffl_lsrain,  &
                             longname='Large-scale rain flux', units='kg m-2 s-1')
    CALL add_stream_element (cospoffl, 'fl_lssnow', cospoffl_lssnow,  &
                             longname='Large-scale snow flux', units='kg m-2 s-1')
    CALL add_stream_element (cospoffl, 'fl_ccrain', cospoffl_ccrain,  &
                             longname='Convective rain flux', units='kg m-2 s-1')
    CALL add_stream_element (cospoffl, 'fl_ccsnow', cospoffl_ccsnow,  &
                             longname='Convective snow flux', units='kg m-2 s-1')

    CALL add_stream_element (cospoffl, 'geom1', cospoffl_geom1, &
         longname = 'geopotential height', &
         units = 'm3 s-2')
    CALL add_stream_element (cospoffl, 'geohm1', cospoffl_geohm1, &
         klev=nlev+1,                                                 &
         longname = 'geopotential height at interfaces', &
         units = 'm3 s-2')
    CALL add_stream_element (cospoffl, 'pfull', cospoffl_p, &
         longname = 'pressure', &
         units = 'Pa')
    CALL add_stream_element (cospoffl, 'phalf', cospoffl_ph, &
         klev=nlev+1,                                                 &
         longname = 'pressure at interfaces', &
         units = 'Pa')


  END  SUBROUTINE construct_stream_cospoffl
!---------------------------------------------------------------------
  SUBROUTINE destruct_stream_cospoffl

    USE mo_memory_base,   ONLY: delete_stream

    IMPLICIT NONE

    CALL delete_stream (cospoffl)

  END  SUBROUTINE destruct_stream_cospoffl

!---------------------------------------------------------------------
FUNCTION accu_info(iname,stream)  RESULT (laccu)

  USE mo_memory_base,      ONLY: get_stream_element_info,  memory_info
  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN ) :: iname
  TYPE (t_stream)  ,POINTER   :: stream 

  LOGICAL  :: laccu

!local
  TYPE (memory_info) :: info

  CALL get_stream_element_info  (stream, TRIM(iname), info)
  laccu=info%laccu

END FUNCTION accu_info
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
           CALL finish('mo_cosp_offline: ', message_text)
          END IF
          CALL nf_check (NF_PUT_ATT_TEXT(ncid,varid,'units', lent, TRIM(unit)))

END SUBROUTINE put_unit
!---------------------------------------------------------------------
SUBROUTINE collect_arg_check(got_opt, got_ivar, lopt, rname)

 USE mo_exception,          ONLY: message_text, finish
 
 IMPLICIT NONE

 LOGICAL, INTENT(IN) :: got_opt, got_ivar
 CHARACTER(LEN=3), INTENT(IN)  :: lopt 
 CHARACTER(LEN=*), INTENT(IN) :: rname
 
   IF (.NOT. ( got_opt .OR. got_ivar ) ) THEN
       WRITE(message_text,*) rname, ': neither lopt nor ivar are set'
       CALL finish('mo_cosp_offline ', message_text )
   END IF
   IF ( got_opt .AND. got_ivar ) THEN
     IF ( lopt .EQ. 'lon' .OR.  lopt .EQ. 'lat' ) THEN
       WRITE(message_text,*) rname, ': ivar and opt=lat or opt=lon'
       CALL finish('mo_cosp_offline ', message_text )
     END IF 
   END IF
   IF ( lopt .EQ. 'var' .OR. lopt .EQ. 'tm1' ) THEN
     IF ( .NOT.  got_ivar  ) THEN
       WRITE(message_text,*) rname, ': need to set ivar'
       CALL finish('mo_cosp_offline ', message_text )
     END IF
   END IF

END SUBROUTINE collect_arg_check
!---------------------------------------------------------------------
SUBROUTINE set_ind_sort_cu( etime, ind, gmax)
  IMPLICIT NONE
  INTEGER, INTENT( IN ) :: gmax
  REAL(dp), DIMENSION(gmax),  INTENT( IN ) :: etime
  INTEGER,  DIMENSION(gmax),  INTENT( INOUT ) :: ind

  REAL(dp), DIMENSION(gmax) :: ttime
  INTEGER :: i, ii(1)

  ttime=etime
  WHERE ( ttime .EQ. -999._dp ) ttime = 1.E20_dp
  ind(:)=-999 
   
   i=1
  DO WHILE ( MINVAL(ttime) .LT. 1.E20_dp)
   ii=MINLOC(ttime)
   ind(i)=ii(1)
   i=i+1
   ttime(ii(1)) = 1.E20_dp
  END DO

END SUBROUTINE set_ind_sort_cu
!---------------------------------------------------------------------
SUBROUTINE ind_sort_cu1d( cu1d , ind, gmax)
  IMPLICIT NONE
  INTEGER, INTENT( IN ) :: gmax
  REAL(dp), DIMENSION(gmax),  INTENT( INOUT ) :: cu1d
  INTEGER,  DIMENSION(gmax),  INTENT( IN ) :: ind

  REAL(dp), DIMENSION(gmax) :: cu1dt
  INTEGER :: i

  cu1dt(:) = cu1d(:)
  DO i=1, gmax
    IF ( ind(i) .NE. -999 ) THEN
     cu1d(i) = cu1dt(ind(i))
    ELSE
      cu1d(i) = -999._dp 
    END IF
  END DO

END SUBROUTINE ind_sort_cu1d
!---------------------------------------------------------------------

SUBROUTINE ind_sort_cu2d( cu2d , ind, nlev, gmax)
  IMPLICIT NONE
  INTEGER, INTENT( IN ) :: gmax, nlev
  REAL(dp), DIMENSION(nlev, gmax),  INTENT( INOUT ) :: cu2d
  INTEGER,  DIMENSION(gmax),  INTENT( IN ) :: ind

  REAL(dp), DIMENSION(nlev,gmax) :: cu2dt
  INTEGER :: i

  cu2dt(:,:) = cu2d(:,:)
  DO i=1, gmax
    IF ( ind(i) .NE. -999 ) THEN
     cu2d(:,i) = cu2dt(:,ind(i))
    ELSE
      cu2d(:,i) = -999._dp 
    END IF
  END DO

END SUBROUTINE ind_sort_cu2d

!---------------------------------------------------------------------
! temporary :use from module_..
 SUBROUTINE p_gather_real_1d2d (sendbuf, recvbuf, p_dest, comm)

    USE mo_mpi,         ONLY:p_all_comm, p_real_dp

    IMPLICIT NONE

    REAL(dp),          INTENT(inout) :: sendbuf(:), recvbuf(:,:)
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

!--------------------------------------------------------------------

 SUBROUTINE p_gather_int_0d1d (sendbuf, recvbuf, p_dest, comm)

    USE mo_mpi,         ONLY:p_all_comm, p_int

    IMPLICIT NONE

    INTEGER,           INTENT(inout) :: sendbuf, recvbuf(:)
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

     CALL MPI_GATHER(sendbuf,  1, p_int, &
                     recvbuf,  1, p_int, &
                     p_dest, p_comm, p_error)
#else
     recvbuf(:) = sendbuf
#endif
 END SUBROUTINE p_gather_int_0d1d



END MODULE mo_cosp_offline
