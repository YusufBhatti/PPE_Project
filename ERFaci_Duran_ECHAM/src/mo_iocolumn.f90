!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_iocolumn
  !> Purpose is to do all IO/Memory of Single Column Model,
  !> Additional output for Single column, 
  !> Generic interfaces for Assimilation/nudging of any model field
  !> in grid space
  !>
  !> Authors:
  !>
  !> A. Rhodin,      MPI, November  1999, original source
  !> I. Kirchner,    MPI, December  2000, time control update
  !> A. Rhodin,      DWD, March     2003, nproma blocking
  !> L. Kornblueh,   MPI, August    2009, bug fix with respect to pole filter
  !> S. K. Cheedela, MPI, January   2010, rewritten
  !>
  !>code needs some final polishing
  USE mo_netcdf
  USE mo_kind,          ONLY: wp
  USE mo_time_control,  ONLY:delta_time
  USE mo_mpi,           ONLY: p_pe, p_io
  USE mo_exception,     ONLY: finish, message
  USE mo_decomposition, ONLY: ldc => local_decomposition
  USE mo_linked_list,   ONLY: t_stream
  USE mo_memory_base,   ONLY: delete_stream, add =>add_stream_element, &
                              default_stream_setting, HYBRID_H, HYBRID, AUTO

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: get_ext_state
  PUBLIC :: get_ext_tend
  PUBLIC :: init_forcing
  PUBLIC :: read_initial3
  PUBLIC :: read_initial2
  PUBLIC :: get_col_location
  PUBLIC :: get_col_vct
  PUBLIC :: construct_scm
  PUBLIC :: destruct_scm

  INTERFACE get_ext_state
    MODULE PROCEDURE get_ext_state3
    MODULE PROCEDURE get_ext_state2
  END INTERFACE  

  INTERFACE get_ext_tend
    MODULE PROCEDURE get_ext_tend3
    MODULE PROCEDURE get_ext_tend2
  END INTERFACE 
  
  CHARACTER(len= 32),PUBLIC    :: forcingfile=""
  REAL(wp),ALLOCATABLE,PUBLIC  :: ftimes(:)
  REAL(wp),ALLOCATABLE         :: forcingtimes(:)
  REAL(wp)                     :: ftime_len 
  REAL(wp),PUBLIC              :: fstart_time
  INTEGER                      :: nftime
  INTEGER                      :: nfid      !netcdf file id of forcing data
  INTEGER                      :: nvarid
  
  !additional variables for output
  !can be shifted to a new module mo_memory_scm
  TYPE(t_stream),POINTER,PUBLIC :: scm 

  REAL(wp),POINTER,PUBLIC       :: scm_u(:,:,:)
  REAL(wp),POINTER,PUBLIC       :: scm_v(:,:,:)
  REAL(wp),POINTER,PUBLIC       :: scm_omega(:,:,:)
  REAL(wp),POINTER,PUBLIC       :: scm_div(:,:,:)
  REAL(wp),POINTER,PUBLIC       :: scm_t(:,:,:)
  REAL(wp),POINTER,PUBLIC       :: scm_aps(:,:)

  REAL(wp),POINTER,PUBLIC       :: ddt_u_phy(:,:,:)
  REAL(wp),POINTER,PUBLIC       :: ddt_v_phy(:,:,:)
  REAL(wp),POINTER,PUBLIC       :: ddt_t_phy(:,:,:)
  REAL(wp),POINTER,PUBLIC       :: ddt_q_phy(:,:,:)
  REAL(wp),POINTER,PUBLIC       :: ddt_xl_phy(:,:,:)
  REAL(wp),POINTER,PUBLIC       :: ddt_xi_phy(:,:,:)

  REAL(wp),POINTER,PUBLIC       :: ddt_u_adv(:,:,:)
  REAL(wp),POINTER,PUBLIC       :: ddt_v_adv(:,:,:)
  REAL(wp),POINTER,PUBLIC       :: ddt_t_adv(:,:,:)
  REAL(wp),POINTER,PUBLIC       :: ddt_q_adv(:,:,:)
  REAL(wp),POINTER,PUBLIC       :: ddt_xl_adv(:,:,:)
  REAL(wp),POINTER,PUBLIC       :: ddt_xi_adv(:,:,:)
  REAL(wp),POINTER,PUBLIC       :: pressure(:,:,:) 

CONTAINS
  
  SUBROUTINE init_forcing
  INTEGER :: ntimeid
  INTEGER :: st(1),ct(1)

    IF(forcingfile/="") THEN
      IF(p_pe == p_io) THEN  
        
        !>get dimensions of time, allocate
        !>find ftimes, finit_time from attribute and ftime_len
        CALL nf_check(nf_open(forcingfile, NF_NOWRITE, nfid),forcingfile)
        
        CALL nf_check(nf_inq_dimid(nfid,"time",ntimeid),forcingfile)
        CALL nf_check(nf_inq_dimlen(nfid,ntimeid,nftime),forcingfile)
        
        ALLOCATE(ftimes(nftime))   
        ALLOCATE(forcingtimes(nftime))       
        st= (/1/)
        ct= (/nftime/)
        CALL nf_check(nf_inq_varid(nfid,"time",ntimeid),forcingfile)
        CALL nf_check(nf_get_vara_double(nfid,ntimeid,st,ct,ftimes),forcingfile)
        
        !>find fstart_time from attribute
        CALL io_get_att_double(nfid,ntimeid,'start',fstart_time) 

        ftimes = ftimes/86400 !converted to fraction of days
        ftime_len = ftimes(nftime)-ftimes(1)
        
      END IF  
    END IF
  END SUBROUTINE init_forcing

  SUBROUTINE get_ext_state3(varname,state3inout,relaxnt,interpn,optn)
  ! 
  !>To get external state, interpolate and relax based on options 
  !
  CHARACTER (len=*),INTENT(in) :: varname
  INTEGER,INTENT(in)           :: relaxnt
  INTEGER,INTENT(in)           :: interpn
  REAL(wp),INTENT(inout)       :: state3inout(:,:,:)
  INTEGER,OPTIONAL,INTENT(in)  :: optn        !future option for any data assimilation
  REAL(wp),POINTER             :: readstate3(:,:,:)
  INTEGER                      :: type=1 

  type = 1

  ALLOCATE(readstate3(ldc%nglon,size(state3inout,2),ldc%nglat)) 
                                              ! must be same dim as varinout

  !idea is if full gcm read 3d field and if column model read 1d..
  
  CALL read_ext_var3(varname,readstate3,interpn) 

  !Execute relaxation
 
  !gather varinout if necessary may be can avoid this part of gather
  IF(PRESENT(optn)) type=optn
  
  SELECT CASE (type)
    
    CASE(1)
      IF(relaxnt<=delta_time) THEN 
        state3inout = readstate3
      ELSE
        state3inout = state3inout + (readstate3-state3inout)*delta_time/relaxnt
      ENDIF

    CASE(2)
     !reserved for any new data assimilation algorithm now cgils relax

  END SELECT

  DEALLOCATE(readstate3)
  
  END SUBROUTINE get_ext_state3 
   
  SUBROUTINE get_ext_state2(varname,state2inout,relaxnt,interpn,optn)
  ! 
  !To get external state, interpolate and relax based on options 
  !
  CHARACTER (len=*),INTENT(in) :: varname
  INTEGER,INTENT(in)           ::relaxnt
  INTEGER,INTENT(in)           ::interpn
  REAL(wp),INTENT(inout)       :: state2inout(:,:)
  INTEGER,OPTIONAL,INTENT(in)  ::optn          ! future option for any data assimilation
  REAL(wp),POINTER             ::readstate2(:,:)
  INTEGER                      ::type = 1


  ALLOCATE(readstate2(ldc%nglon,ldc%nglat))    ! must be same dim as varinout

  !idea is if full gcm read 3d field and if column model read 1d..
  
  CALL read_ext_var2(varname,readstate2,interpn) 

  !Execute relaxation
 
  !gather varinout if necessary may be can avoid this part of gather
  IF(PRESENT(optn)) type=optn
  
  SELECT CASE (type)
    
    CASE(1)
      IF(relaxnt<=delta_time) THEN 
        state2inout = readstate2
      ELSE
        state2inout = state2inout + (readstate2-state2inout)*delta_time/relaxnt 
      ENDIF

    CASE(2)
     !reserved for any new data assimilation algorithm
  END SELECT

  DEALLOCATE(readstate2)
  
  END SUBROUTINE get_ext_state2 

  SUBROUTINE get_ext_tend3(tendname,tend3inout,interpn,tendoptn)
  CHARACTER (len=*),INTENT(in) ::tendname
  REAL(wp),INTENT(inout)       ::tend3inout(:,:,:)
  INTEGER,INTENT(in)           ::interpn
  INTEGER,INTENT(in)           ::tendoptn
  
  REAL(wp),ALLOCATABLE         ::readtend3(:,:,:)

  ALLOCATE(readtend3(ldc%nglon,ldc%nlev, ldc%nglat)) !same dim as tend3inout 

  CALL read_ext_var3(tendname,readtend3,interpn) 

  IF(tendoptn == 0) THEN
    tend3inout = readtend3
  ELSE
    tend3inout = readtend3 + tend3inout 
  ENDIF
  
  DEALLOCATE(readtend3)

  END SUBROUTINE get_ext_tend3

  SUBROUTINE get_ext_tend2(tendname,tend2inout,interpn,tendoptn)
  CHARACTER (len=*),INTENT(in) ::tendname
  REAL(wp),INTENT(inout)       ::tend2inout(:,:)
  INTEGER,INTENT(in)           ::interpn
  INTEGER,INTENT(in)           ::tendoptn
  
  REAL(wp),ALLOCATABLE         ::readtend2(:,:)

  ALLOCATE(readtend2(ldc%nglon,ldc%nglat)) !same dim as tend3inout 

  CALL read_ext_var2(tendname,readtend2,interpn) 

  IF(tendoptn == 0) THEN
    tend2inout = readtend2
  ELSE
    tend2inout = readtend2 + tend2inout 
  ENDIF
  
  DEALLOCATE(readtend2)

  END SUBROUTINE get_ext_tend2
  
  SUBROUTINE read_ext_var3(varname,var3out,interpn)    
  ! Reads "name" from file , scaters data if decomposition 
  ! scatter and gather can be done here or  
  ! currently only done for column model... 
  USE mo_time_control,    ONLY:current_date, get_date_components
  USE mo_time_base,       ONLY:julian_date,Set_JulianDay

  CHARACTER (len=*) ,INTENT(in) ::varname
  INTEGER,INTENT(in)            ::interpn
  REAL(wp),INTENT(inout)        ::var3out(:,:,:) 
  !local
  REAL(wp),ALLOCATABLE          ::readvar2(:,:)

  REAL(wp),POINTER              ::intvar3(:,:,:)
  REAL(wp)                      ::ifac1,ifac2              !interpolation factors
  REAL(wp)                      ::mtime
  INTEGER,ALLOCATABLE           ::st(:),ct(:)
  INTEGER                       ::i,nlevels
  INTEGER                       :: yr,mon,day,hr,min,sec
  TYPE(julian_date)             :: dayout

  IF(p_pe == p_io) THEN
   
    CALL get_date_components(current_date,yr,mon,day,hr,min,sec)
    sec = hr*3600+min*60+sec   ! necessary as below subroutine has only sec in input
    CALL Set_JulianDay(yr,mon,day,sec,dayout)
    mtime =  dayout%day+dayout%fraction    !mtime = 3800._wp + fstart_time

    IF (interpn==1) THEN !cyclic case
      !modify mtime for comparision with ftime,
      mtime = mod(mtime,ftime_len) !issue w.r to 0.5 in julianday
      forcingtimes = ftimes
    ELSE
      forcingtimes = fstart_time + ftimes 
    ENDIF
       
    nlevels = size(var3out,2)
    IF(ldc%col_1d) THEN !> in column mode
      ALLOCATE(intvar3(1,nlevels,1))
      ALLOCATE(st(2))
      ALLOCATE(ct(2)) 
      IF(mtime > forcingtimes(1) .AND. mtime < forcingtimes(nftime) ) THEN 
        DO i=1,nftime-1
          IF (mtime > forcingtimes(i) .AND. mtime <= forcingtimes(i+1) ) EXIT            
        END DO
        ALLOCATE(readvar2(nlevels,2))
        readvar2 = -999.0_wp
        st =(/1,i/) 
        ct =(/nlevels,2/) 
        CALL nf_check(nf_inq_varid(nfid,varname,nvarid))
        CALL nf_check(nf_get_vara_double(nfid,nvarid,st,ct,readvar2))
        IF(ANY(readvar2==-999.0_wp)) &
          & CALL finish('read_ext_var3','ERROR in reading, check dimensions')
        CALL get_ifactors(forcingtimes(i),forcingtimes(i+1),mtime,ifac1,ifac2)
        intvar3(1,:,1) = ifac1*readvar2(:,1) + ifac2*readvar2(:,2)
      ELSEIF (mtime <= forcingtimes(1)) THEN 
        st =(/1,1/) 
        ct =(/nlevels,1/) 
        ALLOCATE(readvar2(nlevels,1))
        readvar2 = -999.0_wp
        CALL nf_check(nf_inq_varid(nfid,varname,nvarid))
        CALL nf_check(nf_get_vara_double(nfid,nvarid,st,ct,readvar2))
        IF(ANY(readvar2==-999.0_wp)) &
          & CALL finish('read_ext_var2','ERROR in reading, check dimensions')
        intvar3(1,:,1) = readvar2(:,1)
      ELSE !IF (mtime > forcingtimes(nftime))
        st =(/1,nftime/) 
        ct =(/nlevels,1/) 
        ALLOCATE(readvar2(nlevels,1))
        readvar2 = -999.0_wp
        CALL nf_check(nf_inq_varid(nfid,varname,nvarid))
        CALL nf_check(nf_get_vara_double(nfid,nvarid,st,ct,readvar2))
        IF(ANY(readvar2==-999.0_wp)) &
          & CALL finish('read_ext_var2','ERROR in reading, check dimensions')
        intvar3(1,:,1) = readvar2(:,1)
      ENDIF
   
      var3out = intvar3

      DEALLOCATE(readvar2)
      DEALLOCATE(intvar3)

    ELSE
     !full model 3d casse still needs to be fixed
     !CALL scatter_gp(intvar3,outvar3,gdc) !still needs to be fixed
     !DEALLOCATE(readvar2)
     !NULLIFY(intvar3)
    END IF
       
  END IF
 
 END SUBROUTINE read_ext_var3

 SUBROUTINE read_ext_var2(varname,var2out,interpn)    
  ! Reads "name" from file , scaters data if decomposition 
  ! scatter and gather can be done here or  
  ! currently only done for column model... 
  USE mo_time_control,    ONLY:current_date, get_date_components
  USE mo_time_base,       ONLY:julian_date,Set_JulianDay

  CHARACTER (len=*) ,INTENT(in) ::varname
  INTEGER,INTENT(in)      ::interpn
  REAL(wp),INTENT(out)    ::var2out(:,:) 
  !local
  REAL(wp),ALLOCATABLE    ::readvar1(:)

  REAL(wp),POINTER        ::intvar2(:,:)
  REAL(wp)                ::ifac1,ifac2              !interpolation factors
  REAL(wp)                ::mtime
  INTEGER,ALLOCATABLE     ::st(:),ct(:)
  INTEGER                 ::i
  INTEGER                 :: yr,mon,day,hr,min,sec
  TYPE(julian_date)       :: dayout

  IF(p_pe == p_io) THEN
   
    CALL get_date_components(current_date,yr,mon,day,hr,min,sec)
    sec = hr*3600+min*60+sec   ! necessary as below subroutine has only sec in input
    CALL Set_JulianDay(yr,mon,day,sec,dayout)
    mtime =  dayout%day+dayout%fraction    !mtime = 3800._wp + fstart_time

    IF (interpn==1) THEN !cyclic case
      !modify mtime for comparision with ftime,
      mtime = mod(mtime,ftime_len) !issue w.r to 0.5 in julianday
      forcingtimes = ftimes
    ELSE
      forcingtimes = fstart_time + ftimes 
    ENDIF
       
    IF(ldc%col_1d) THEN !> in column mode
      ALLOCATE(intvar2(1,1))
      ALLOCATE(st(1))
      ALLOCATE(ct(1)) 
      IF(mtime > forcingtimes(1) .AND. mtime < forcingtimes(nftime)) THEN 
        DO i=1,nftime-1
          IF (mtime > forcingtimes(i) .AND. mtime <= forcingtimes(i+1) ) EXIT            
        END DO
        ALLOCATE(readvar1(2))
        readvar1 = -999.0_wp
        st =(/i/) 
        ct =(/2/) 
        CALL nf_check(nf_inq_varid(nfid,varname,nvarid))
        CALL nf_check(nf_get_vara_double(nfid,nvarid,st,ct,readvar1))
        IF(ANY(readvar1==-999.0_wp)) &
          & CALL finish('read_ext_var2','ERROR in reading, check dimensions')
        CALL get_ifactors(forcingtimes(i),forcingtimes(i+1),mtime,ifac1,ifac2)
        intvar2(1,1) = ifac1*readvar1(1) + ifac2*readvar1(2)
      ELSEIF (mtime <= forcingtimes(1)) THEN 
        st =(/1/) 
        ct =(/1/) 
        ALLOCATE(readvar1(1))
        readvar1 = -999.0_wp
        CALL nf_check(nf_inq_varid(nfid,varname,nvarid))
        CALL nf_check(nf_get_vara_double(nfid,nvarid,st,ct,readvar1))
        IF(ANY(readvar1==-999.0_wp)) &
          & CALL finish('read_ext_var2','ERROR in reading, check dimensions')
        intvar2(1,1) = readvar1(1)
      ELSE !IF (mtime > forcingtimes(nftime))
        st =(/nftime/) 
        ct =(/1/) 
        ALLOCATE(readvar1(1))
        readvar1 = -999.0_wp
        CALL nf_check(nf_inq_varid(nfid,varname,nvarid))
        CALL nf_check(nf_get_vara_double(nfid,nvarid,st,ct,readvar1))
        IF(ANY(readvar1==-999.0_wp)) &
          & CALL finish('read_ext_var2','ERROR in reading, check dimensions')
        intvar2(1,1) = readvar1(1)
      ENDIF
    
      var2out = intvar2

      DEALLOCATE(readvar1)
      DEALLOCATE(intvar2)

    ELSE
     !full model 3d casse still needs to be fixed
     !CALL scatter_gp(intvar2,outvar2,gdc) !still needs to be fixed
     !DEALLOCATE(readvar3)
     !NULLIFY(intvar2)
    END IF
       
  END IF
  
  END SUBROUTINE read_ext_var2  

  SUBROUTINE read_initial2(varname,var2out,readoptn)    
  ! Reads initial data for single column model from scm forcing file 
  !.error will be based on if readoption is essential or optional 
  CHARACTER (len=*) ,INTENT(in)         ::varname
  REAL(wp),INTENT(inout)                ::var2out(:,:) 
  CHARACTER(len=*),OPTIONAL,INTENT(in)  :: readoptn ! reading option for error message 
                                                    ! can be eighter essential or optional  
  !local
  REAL(wp),ALLOCATABLE    ::readvar1(:)
  INTEGER                 ::st(1), ct(1)
  IF(p_pe == p_io) THEN
  st =(/1/) 
  ct =(/1/) 
  ALLOCATE(readvar1(1))
  readvar1 = -999.0_wp
  IF ( readoptn == "essential" ) THEN 
    CALL nf_check(nf_inq_varid(nfid,varname,nvarid),varname)
    CALL nf_check(nf_get_vara_double(nfid,nvarid,st,ct,readvar1),varname)
    IF(ANY(readvar1==-999.0_wp)) CALL finish('ERROR in reading essential inital data',varname)
    var2out(1,1) = readvar1(1)
  ELSE!optional case
    !no read error checks done
    !CALL io_inq_varid(nfid,varname,nvarid)
    !CALL io_get_vara_double(nfid,nvarid,st,ct,readvar1)
    IF(ANY(readvar1==-999.0_wp)) THEN
      CALL message('ERROR in reading optional inital data',varname)
    ELSE
      var2out(1,1) = readvar1(1)
    ENDIF
  ENDIF
  ENDIF
      
  ! additional else can be added for 3D case.
  END SUBROUTINE read_initial2

  SUBROUTINE read_initial3(varname,var3out,readoptn)    
  ! Reads initial data for single column model from scm forcing file 
  ! error will be based on if readoption is essential or optional 
  CHARACTER (len=*) ,INTENT(in)         ::varname
  REAL(wp),INTENT(inout)                ::var3out(:,:,:) 
  CHARACTER(len=*),OPTIONAL,INTENT(in)  :: readoptn ! reading option for error message 
                                                    ! can be eighter essential or optional  
  !local
  REAL(wp),ALLOCATABLE                  ::readvar2(:,:)
  INTEGER                               ::st(2), ct(2)
  
  IF (p_pe == p_io) THEN
    st =(/1,1/) 
    ct =(/ldc%nlev,1/) 
    ALLOCATE(readvar2(ldc%nlev,1))
    readvar2 = -999.0_wp
    IF ( readoptn == "essential" ) THEN
      CALL nf_check(nf_inq_varid(nfid,varname,nvarid),varname)
      CALL nf_check(nf_get_vara_double(nfid,nvarid,st,ct,readvar2),varname)
      IF(ANY(readvar2==-999.0_wp)) CALL finish('ERROR in reading essential inital data',varname)
      var3out(1,:,1) = readvar2(:,1)
    ELSE!optional case
      !no read error checks done
      !CALL nf_check(nf_inq_varid(nfid,varname,nvarid))
      !CALL nf_check(nf_get_vara_double(nfid,nvarid,st,ct,readvar2))
      IF(ANY(readvar2==-999.0_wp)) THEN
        CALL message('ERROR in reading optional inital data',varname)
      ELSE
        var3out(1,:,1) = readvar2(:,1)
      ENDIF
    ENDIF
  ENDIF
      
  ! additional else can be added for 3D case.
  END SUBROUTINE read_initial3

  SUBROUTINE get_ifactors(t1,t2,t,ifac1,ifac2)
  ! Computes linear interpolation factors a,b for interpolating 
  ! forcing at t1,t2 to model time t--- f(t)= a*f(t1)+b*f(t2) 
  REAL(wp),INTENT(IN)   :: t1,t2  !time in julian days or any standard days since
  REAL(wp),INTENT(IN)   :: t 
  REAL(wp),INTENT(OUT)  :: ifac1,ifac2
  
  ifac1 = (t2-t)/(t2-t1)
  ifac2 = (t-t1)/(t2-t1)
  
  !check for factors
   
  END SUBROUTINE get_ifactors

  SUBROUTINE get_col_location(dlat,dlon)
  REAL(wp), INTENT(out) :: dlat(1),dlon(1) !latitude and longitude in degrees
  CALL nf_check(nf_inq_varid(nfid, "lon", nvarid))
  CALL nf_check(nf_get_var_double(nfid,nvarid,dlon))
  CALL nf_check(nf_inq_varid(nfid, "lat", nvarid))
  CALL nf_check(nf_get_var_double(nfid,nvarid,dlat))
  ! Read dlat dlon from forcing file, check if single column
  ! find index of location on current horizontal grid
  END SUBROUTINE get_col_location

  SUBROUTINE get_col_vct(nlev, nlevp1, nvclev, vct)
  ! not sure if P-bcast additional complication required for working in 
  !  full model configuration!!
  ! no start and count used for reading can cause problems
  INTEGER,INTENT(inout) :: nlev,nlevp1
  INTEGER,INTENT(inout) :: nvclev
  REAL(wp),POINTER      :: vct(:)
  REAL(wp),ALLOCATABLE  :: vct_a(:),vct_b(:)
  IF (forcingfile/="") THEN
    CALL nf_check(nf_inq_dimid(nfid,"nlev",nvarid))
    CALL nf_check(nf_inq_dimlen(nfid,nvarid,nlev)) 
    ALLOCATE(vct_a(nlev+1))
    ALLOCATE(vct_b(nlev+1))
    CALL nf_check(nf_inq_varid(nfid,"vct_a",nvarid))
    CALL nf_check(nf_get_var_double(nfid,nvarid,vct_a))
    CALL nf_check(nf_inq_varid(nfid,"vct_b",nvarid))
    CALL nf_check(nf_get_var_double(nfid,nvarid,vct_b))
  
    !IF (p_parallel) CALL p_bcast (nvclev, p_io)
    WRITE (*,*) 'Model hybrid levels changed to ',nlev,' from forcingfile',forcingfile
    !>Warning to remind setting the timestep size, if vertical resolution changed
    IF(nlev+1>nvclev) THEN
      CALL message('Warning','Set the timestep size using namelist option delta_time')
      CALL message('','else could lead to instability')
    END IF
    
    nvclev   = nlev+1
    nlevp1   = nvclev
    
    NULLIFY(vct)
    ALLOCATE(vct(2*nvclev))

    vct(1:nvclev) = vct_a
    vct(nvclev+1:2*nvclev) = vct_b

    DEALLOCATE(vct_a)
    DEALLOCATE(vct_b)

  ENDIF
  !IF (p_parallel) CALL p_bcast (vct, p_io)
  
  END SUBROUTINE get_col_vct

  SUBROUTINE construct_scm
  ! construct output streams of ECHAM SCM 

  CALL default_stream_setting( scm,lpost = .TRUE., code= AUTO, leveltype= HYBRID) ! code? table? othervars LEveltype

  !add additional variables to output stream
  CALL add (scm,'u',         scm_u,     longname='Zonal Velocity',                             units='m/s')
  CALL add (scm,'v',         scm_v,     longname='Meridional Velocity',                        units='m/s')
  CALL add (scm,'omega',     scm_omega, longname='Vertical Pressure Velocity',                 units='Pa/s')
  CALL add (scm,'div',       scm_div,   longname='Horizontal divergence',                      units='1/s')
  CALL add (scm,'t',         scm_t,     longname='Temperature',                                units='K')
   
  CALL add (scm,'ddt_u_adv', ddt_u_adv, longname='LS Advective Zonal Wind Tendency',           units='m/s')
  CALL add (scm,'ddt_v_adv', ddt_v_adv, longname='LS Advective Meridional Wind Tendency',      units='m/s')
  CALL add (scm,'ddt_t_adv', ddt_t_adv, longname='LS Advective Temparature Tendency',          units='K/s')
  CALL add (scm,'ddt_q_adv', ddt_q_adv, longname='LS Advective Water Vapor Specific Humidity', units='kg/kg')
  CALL add (scm,'ddt_xl_adv',ddt_xl_adv,longname='LS Advective Liquid Water Specific Humidity',units='kg/kg')
  CALL add (scm,'ddt_xi_adv',ddt_xi_adv,longname='Ls Advective Ice Specific Humidity',         units='kg/kg')
   
  CALL add (scm,'ddt_u_phy', ddt_u_phy, longname='Physics Zonal Wind Tendency',                units='m/s2')
  CALL add (scm,'ddt_v_phy', ddt_v_phy, longname='Physics Meridional Wind Tendency',           units='m/s2')
  CALL add (scm,'ddt_t_phy', ddt_t_phy, longname='Physics Temparature Tendency',               units='K/s')
  CALL add (scm,'ddt_q_phy', ddt_q_phy, longname='Physics Water Vapor Specific Humidity',      units='kg/kg/s')
  CALL add (scm,'ddt_xl_phy',ddt_xl_phy,longname='Physics Liquid Water Specific Humidity',     units='kg/kg/s')
  CALL add (scm,'ddt_xi_phy',ddt_xi_phy,longname='Physics Ice Specific Humidity',              units='kg/kg/s')
  
  CALL add (scm,'aps',       scm_aps,   longname='surface pressure',           lpost= .false., units='m/s') 
  CALL add (scm,'pressure',  pressure,  longname='Pressure',               leveltype=HYBRID_H, units='Pa')

  END SUBROUTINE construct_scm

  SUBROUTINE destruct_scm

  CALL delete_stream(scm)

  END SUBROUTINE destruct_scm

  END MODULE mo_iocolumn
