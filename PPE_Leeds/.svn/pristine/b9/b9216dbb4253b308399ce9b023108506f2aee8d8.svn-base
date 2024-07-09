!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
module mo_cbalone_io

  USE mo_kind,            ONLY: dp
  USE mo_cbalone_memory,  ONLY: check_err
  USE mo_jsbach,          ONLY: debug, debug_Cconservation

  implicit none

  public :: getDailyData
  public :: writeSingleTimeStep

  private

  TYPE netcdfInfo_type
     character(len=64) :: varName  = ''
     integer           :: varType  = -1   !! E.g. NF_FLOAT or NF_DOUBLE
     integer           :: dimType  = -1   !! BOX_TYPE or  TILES_TYPE (see below)  
     integer           :: varId    = -1   !! Is set when calling netcdf definition routine
     real(dp),pointer  :: array_box(:)     => NULL() !! For access to the data (only if NF_FLOAT  and BOX_TYPE)
     real(dp),pointer  :: array_tiles(:,:) => NULL() !! For access to box data (only if NF_FLOAT  and TILES_TYPE)
  end TYPE netcdfInfo_type

  !! For the definition of Netcdf dimensions: 
  integer,parameter :: BOX_TYPE = 3          !! lat x lon x time
  integer,parameter :: TILES_TYPE = 4        !! lat x lon x tiles x time

  !! Module variables

  integer :: numberOfOutputVariables
  TYPE(netcdfInfo_type),pointer :: netcdfInfo(:)

contains
  
  ! --- getDailyData()  -----------------------------------------------------------------------------------
  !
  ! This routine reads in a time series of daily forcing data, called LAI_yDayMean, NPP_yDayMean, topSoilTemp_yDayMean
  ! Runoff_yDayMean, and alpha_yDayMean. The routine automatically detects whether data are packed or in LonLat format.
  !
  SUBROUTINE getDailyData(filename, grid, run_dynveg, run_disturbance,with_nitrogen, with_yasso, &
                          cbalance_diag, nbalance_offline, dynveg_clim, ntiles, nday)

    USE mo_cbalone_memory, ONLY: cbal_offline_type, Nbal_offline_type, dynveg_clim_type
    USE mo_jsbach_grid,    ONLY: grid_type

    include 'netcdf.inc'
    character(len=*),intent(in)            :: filename
    type(grid_type),       INTENT(in)      :: grid
    LOGICAL,               INTENT(in)      :: run_dynveg
    LOGICAL,               INTENT(in)      :: run_disturbance
    LOGICAL,               INTENT(in)      :: with_nitrogen
    LOGICAL,               INTENT(in)      :: with_yasso   
    TYPE(cbal_offline_type), INTENT(inout) :: cbalance_diag
    TYPE(Nbal_offline_type), INTENT(inout) :: nbalance_offline
    TYPE(dynveg_clim_type),INTENT(inout)   :: dynveg_clim 
    integer,intent(in)                     :: ntiles
    integer,intent(out)                    :: nday               !! number of days in the considered month
    ! local variables

    logical :: lai_found, npp_found, temp_found, alpha_found,runoff_found
    logical :: yasso_precip_found, yasso_temp_found 
    logical :: ave_npp5_found, bio_exist_found, rel_hum_found
    logical :: max_wind10_found , prev_day_max_wind10_found
    logical :: prev_day_mean_wind10_found
    logical :: prev_day_temp_max_found, prev_day_temp_min_found, prev_day_precip_mean_found !, dew_point_temp_found
    logical :: prev_day_mean_vol_moist_found

    real,     allocatable :: hlp2D_sp(:,:),hlp3D_sp(:,:,:),hlp4D_sp(:,:,:,:)
    real(dp), allocatable :: hlp2D(:,:),hlp3D(:,:,:),hlp4D(:,:,:,:)

    ! error status return
    integer  ::  iret
    ! netCDF id
    integer  ::  ncid
    ! dimensions
    integer  ::  ndim,nvar,nattr,nunlimdim
    integer,allocatable,dimension(:)  ::  ndimlen
    character(len=50),allocatable,dimension(:)  ::  fdimnam
    ! variables
    integer ::  nin(100)
    integer,allocatable ,dimension(:)  ::  nvartyp
    integer,allocatable,dimension(:)   ::  nvardim
    integer,allocatable,dimension(:,:) ::  nvardimid
    integer,allocatable,dimension(:)   ::  nvaratt
    character(len=128),allocatable,dimension(:)  ::  fvarnam
    ! indices
    integer  ::  i, ii, itile, day

    ! preparations

    lai_found                     = .false.
    npp_found                     = .false.
    temp_found                    = .false.
    alpha_found                   = .false.
    yasso_temp_found              = .false.
    yasso_precip_found            = .false.
    runoff_found                  = .false.
    ave_npp5_found                = .false.
    bio_exist_found               = .false.
    prev_day_max_wind10_found     = .false.
    max_wind10_found              = .false.
    rel_hum_found                 = (.not. run_dynveg) .and. (.not. run_disturbance)
    prev_day_mean_wind10_found    = (.not. run_dynveg) .and. (.not. run_disturbance)
    prev_day_temp_max_found       = (.not. run_dynveg) .and. (.not. run_disturbance)
    prev_day_temp_min_found       = (.not. run_dynveg) .and. (.not. run_disturbance)
!    dew_point_temp_found          = (.not. run_dynveg) .and. (.not. run_disturbance)
    prev_day_precip_mean_found    = (.not. run_dynveg) .and. (.not. run_disturbance)
    prev_day_mean_vol_moist_found = (.not. run_dynveg) .and. (.not. run_disturbance)
    if (run_disturbance) then
       rel_hum_found                 = .not. ASSOCIATED(dynveg_clim%rel_hum_air)
       prev_day_mean_wind10_found    = .not. associated(dynveg_clim%prev_day_mean_wind10)
       prev_day_temp_max_found       =  .not. associated(dynveg_clim%prev_day_temp_max)
       prev_day_temp_min_found       =  .not. associated(dynveg_clim%prev_day_temp_min)
!       dew_point_temp_found          = .not. associated(dynveg_clim%dew_point_temp)
       prev_day_precip_mean_found    = .not. associated(dynveg_clim%prev_day_precip_mean)
       prev_day_mean_vol_moist_found = .not. associated(dynveg_clim%prev_day_mean_vol_moist)
    endif

    ! open the input-file

    if(debug) write (*,*) 'getDailyData(): read from: ',filename
    iret = nf_open(filename,nf_nowrite,ncid)
    call check_err(iret,"Error in getDailyData() when opening file "//trim(filename))

    ! check what is in the input-file

    iret = nf_inq(ncid,ndim,nvar,nattr,nunlimdim)
    call check_err(iret,"Error in getDailyData() when inquiring contents of "//trim(filename))
    if (debug) write (*,*) 'getDailyData(): ndim,nvar,nattr,nunlimdim',ndim,nvar,nattr,nunlimdim

    ! get the dimension name and length

    allocate (ndimlen(ndim))
    allocate (fdimnam(ndim))
    do i=1,ndim
       iret = nf_inq_dim(ncid,i,fdimnam(i),ndimlen(i))
       call check_err(iret,"Error in getDailyData() when inquiring dimension of "//trim(fdimnam(i))//" in "//trim(filename))
       if (debug) write (*,*) 'getDailyData(): i,fdimnam,ndimlen',i,trim(fdimnam(i)),ndimlen(i)
    end do

    ! set ndate (number of output intervals of the monthly file) and noutput_per_day

    nday = 0
    do i=1,ndim
       if (fdimnam(i) == 'time') nday = ndimlen(i)
    end do
    if (debug) write (*,*) 'getDailyData(): nday: ',nday
    if (nday == 0) then
       write(*,*) 'getDailyData(): variable time not found in '//trim(filename)
       stop
    end if

    ! get variable names, types and shapes

    allocate (fvarnam(nvar))
    allocate (nvartyp(nvar))
    allocate (nvardim(nvar))
    allocate (nvardimid(nvar,100))
    allocate (nvaratt(nvar))
    do i=1,nvar
       iret = nf_inq_var(ncid,i,fvarnam(i),nvartyp(i),nvardim(i),nin,nvaratt(i))
       call check_err(iret)
       do ii=1,nvardim(i)
          nvardimid(i,ii)=nin(ii)
       end do
       if (debug) write (*,*) 'getDailyData(): i,fvarnam,nvartyp,nvardim,nvardimid,nvaratt',&
            i,trim(fvarnam(i)),nvartyp(i),nvardim(i),nvardimid(i,1),nvaratt(i)
    end do

    ! get data
    
    do i=1,nvar
       select case (trim(fvarnam(i)))
       case('LAI_yDayMean','NPP_yDayMean','alpha_yDayMean','ave_npp5','bio_exist')
          select case (nvardim(i))
          case(3) !! The array has 3 dimensions => We find packed data
             if (ndimlen(nvardimid(i,1)) /= grid%nland .or. ndimlen(nvardimid(i,2)) /= ntiles) then
                write(*,*) "getDailyData(): Dimension error for variable "//trim(fvarnam(i))
                write(*,*) "getDailyData(): nland=",grid%nland," but ndimlen(nvardimid(i,1)=",ndimlen(nvardimid(i,1))
                write(*,*) "getDailyData(): ntiles=",ntiles," but ndimlen(nvardimid(i,2)=",ndimlen(nvardimid(i,2))
                stop
             end if
             if (ndimlen(nvardimid(i,3)) /= nday) then
                write(*,*) "getDailyData(): Number of time steps wrong for input variable "//trim(fvarnam(i))
                stop
             end if

             allocate (hlp3D(grid%nland,ntiles,nday))
             select case (nvartyp(i))
             case(NF_REAL)
                allocate (hlp3D_sp(grid%nland,ntiles,nday))
                iret = nf_get_var_real(ncid,i,hlp3D_sp)
                call check_err(iret,"Error(1) in getDailyData() reading variable "//trim(fvarnam(i))//" from "//trim(filename))
                hlp3D = REAL(hlp3D_sp)
                deallocate(hlp3D_sp)
             case(NF_DOUBLE)
                iret = nf_get_var_double(ncid,i,hlp3D)
                call check_err(iret,"Error(1) in getDailyData() reading variable "//trim(fvarnam(i))//" from "//trim(filename))
             case default
                write(*,*) "ERROR: unsupported type of variable ", trim(fvarnam(i)), " in file ", trim(filename)
                stop
             end select
          case(4)!! The array has 4 dimensions => We find latlon data
             if(ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat &
                  .or. ndimlen(nvardimid(i,3)) /= ntiles) then
                write(*,*) "getDailyData(): Dimension error for variable "//trim(fvarnam(i))//" in file "//trim(filename)//":"
                write(*,*) "getDailyData(): nlon=",grid%nlon," but found nndimlen(nvardimid(i,1))=",ndimlen(nvardimid(i,1))
                write(*,*) "getDailyData(): nlat=",grid%nlat," but found nndimlen(nvardimid(i,2))=",ndimlen(nvardimid(i,2))
                write(*,*) "getDailyData(): ntiles=",ntiles," but ndimlen(nvardimid(i,3)=",ndimlen(nvardimid(i,3))
                stop
             end if
             if (ndimlen(nvardimid(i,4)) /= nday) then
                write(*,*) "getDailyData(): Number of time steps wrong for input variable "//&
                     trim(fvarnam(i))//" in file "//trim(filename)//":"
                stop
             end if
             allocate (hlp4D(grid%nlon,grid%nlat,ntiles,nday))
             select case (nvartyp(i))
             case(NF_REAL)
                allocate (hlp4D_sp(grid%nlon,grid%nlat,ntiles,nday))
                iret = nf_get_var_real(ncid,i,hlp4D_sp)
                call check_err(iret,"Error(2)in getDailyData() reading variable "//trim(fvarnam(i))//" from "//trim(filename))
                hlp4D=real(hlp4d_sp)
                deallocate(hlp4D_sp)
             case(NF_DOUBLE)
                iret = nf_get_var_double(ncid,i,hlp4D)
                call check_err(iret,"Error(2)in getDailyData() reading variable "//trim(fvarnam(i))//" from "//trim(filename))
             case default
                write(*,*) "ERROR: unsupported type of variable ", trim(fvarnam(i)), " in file ", trim(filename)
                stop
             end select

             ! pack the array
             allocate (hlp3D(grid%nland,ntiles,nday))
             do itile=1,ntiles
                do day=1,nday
                   hlp3D(:,itile,day) = pack(hlp4D(:,:,itile,day),mask=grid%mask)
                end do
             end do
             deallocate(hlp4D)
          case default
             write(*,*) "getDailyData(): ERROR: Variable "//trim(fvarnam(i))//&
                          " from file "//trim(filename)//" is neither of type 'packed' nor 'LonLat'."
             stop
          end select
       case('topSoilTemp_yDayMean','Runoff_yDayMean','rel_hum_air','rel_hum_air_yDay','prev_day_max_wind10','max_wind10', &
            'WindSpeed_yDayMean', 'prev_day_mean_wind10',&
            'pseudo_temp_yDay','pseudo_precip_yDay','pseudo_temp2_30d','pseudo_precip_30d','max_2m_temp','min_2m_temp',&
            'prev_day_temp_max','prev_day_temp_min','prev_day_precip_mean','precip','prev_day_mean_vol_moist')
             ! 'dew_point_temp'
          select case (nvardim(i))
          case(2) !! The array has 2 dimensions => We find packed data
             if (ndimlen(nvardimid(i,1)) /= grid%nland) then
                write(*,*) "getDailyData(): Dimension error for variable "//trim(fvarnam(i))//" in file "//trim(filename)//":"
                write(*,*) "getDailyData(): nland=",grid%nland," but ndimlen(nvardimid(i,1)=",ndimlen(nvardimid(i,1))
                stop
             end if
             if (ndimlen(nvardimid(i,2)) /= nday) then
                write(*,*) "getDailyData(): Number of time steps wrong for input variable "//trim(fvarnam(i))//&
                     " in file "//trim(filename)//":"
                write(*,*) "getDailyData(): nday=",nday," but ndimlen(nvardimid(i,2))=",ndimlen(nvardimid(i,2))
                stop
             end if
             allocate (hlp2D(grid%nland,nday))
             select case (nvartyp(i))
             case(NF_REAL)
               allocate (hlp2D_sp(grid%nland,nday))
               iret = nf_get_var_real(ncid,i,hlp2D_sp)
               call check_err(iret,"Error(3) in getDailyData() reading variable "//trim(fvarnam(i))//" from "//trim(filename))
               hlp2D=real(hlp2d_sp)
               deallocate(hlp2D_sp)
             case(NF_DOUBLE)
               iret = nf_get_var_double(ncid,i,hlp2D)
               call check_err(iret,"Error(3) in getDailyData() reading variable "//trim(fvarnam(i))//" from "//trim(filename))
             case default
                write(*,*) "ERROR: unsupported type of variable ", trim(fvarnam(i)), " in file ", trim(filename)
                stop
             end select

          case(3)!! The array has 3 dimensions => We find latlon data
             if(ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) then
                write(*,*) "getDailyData(): Dimension error for variable "//trim(fvarnam(i))//" in file "//trim(filename)//":"
                write(*,*) "getDailyData(): nlon=",grid%nlon," but found nndimlen(nvardimid(i,1))=",ndimlen(nvardimid(i,1))
                write(*,*) "getDailyData(): nlat=",grid%nlat," but found nndimlen(nvardimid(i,2))=",ndimlen(nvardimid(i,2))
                stop
             end if
             if (ndimlen(nvardimid(i,3)) /= nday) then
                write(*,*) "getDailyData(): Number of time steps wrong for input variable "//&
                     trim(fvarnam(i))//" in file "//trim(filename)//":"
                write(*,*) "getDailyData(): nday=",nday," but ndimlen(nvardimid(i,3))=",ndimlen(nvardimid(i,3))
                stop
             end if
             allocate (hlp3D(grid%nlon,grid%nlat,nday))
             select case (nvartyp(i))
             case(NF_REAL)
                allocate (hlp3D_sp(grid%nlon,grid%nlat,nday))
                iret = nf_get_var_real(ncid,i,hlp3D_sp)
                call check_err(iret,"Error(4) in getDailyData() when reading variable "//trim(fvarnam(i))//" from "//trim(filename))
                hlp3D=real(hlp3D_sp)
                deallocate(hlp3D_sp)
             case(NF_DOUBLE)
                iret = nf_get_var_double(ncid,i,hlp3D)
                call check_err(iret,"Error(4) in getDailyData() when reading variable "//trim(fvarnam(i))//" from "//trim(filename))
             case default
                write(*,*) "ERROR: unsupported type of variable ", trim(fvarnam(i)), " in file ", trim(filename)
                stop
             end select

             ! pack the array
             allocate (hlp2D(grid%nland,nday))
             do day=1,nday
                hlp2D(:,day) = pack(hlp3D(:,:,day),mask=grid%mask)
             end do
             deallocate(hlp3D)
          case(4)!! The array has 4 dimensions => We find latlon data,with only one level
             if(ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) then
                write(*,*) "getDailyData(): Dimension error for variable "//trim(fvarnam(i))//" in file "//trim(filename)//":"
                write(*,*) "getDailyData(): nlon=",grid%nlon," but found nndimlen(nvardimid(i,1))=",ndimlen(nvardimid(i,1))
                write(*,*) "getDailyData(): nlat=",grid%nlat," but found nndimlen(nvardimid(i,2))=",ndimlen(nvardimid(i,2))
                stop
             end if
             if (ndimlen(nvardimid(i,4)) /= nday) then
                write(*,*) "getDailyData(): Number of time steps wrong for input variable "//&
                     trim(fvarnam(i))//" in file "//trim(filename)//":"
                write(*,*) "getDailyData(): nday=",nday," but ndimlen(nvardimid(i,3))=",ndimlen(nvardimid(i,4))
                stop
             end if
             if (ndimlen(nvardimid(i,3)) /= 1) then
                write(*,*) "getDailyData(): Number of tiles wrong for input variable "//&
                     trim(fvarnam(i))//" in file "//trim(filename)//":"
                write(*,*) "getDailyData(): should be 1 but is ",ndimlen(nvardimid(i,3))
                stop
             end if

             allocate (hlp4D(grid%nlon,grid%nlat,1,nday))
             select case (nvartyp(i))
             case(NF_REAL)
                allocate (hlp4D_sp(grid%nlon,grid%nlat,1,nday))
                iret = nf_get_var_real(ncid,i,hlp4D_sp)
                call check_err(iret,"Error(9) in getDailyData() when reading variable "//trim(fvarnam(i))//" from "//trim(filename))
                hlp4D=real(hlp4D_sp)
                deallocate(hlp4D_sp)
             case(NF_DOUBLE)
                iret = nf_get_var_double(ncid,i,hlp4D)
                call check_err(iret,"Error(5) in getDailyData() when reading variable "//trim(fvarnam(i))//" from "//trim(filename))
                write(*,*) "grid%nlon,grid%nlat,nday=",grid%nlon,grid%nlat,nday
             case default
                write(*,*) "ERROR(abcdef): unsupported type of variable ", trim(fvarnam(i)), " in file ", trim(filename)
                stop
             end select

             ! pack the array
             allocate (hlp2D(grid%nland,nday))
             do day=1,nday
                hlp2D(:,day) = pack(hlp4D(:,:,1,day),mask=grid%mask)
             end do
             deallocate(hlp4D)
          case default
             write(*,*) "getDailyData(): ERROR: Variable "//trim(fvarnam(i))//&
                  " from file "//trim(filename)//" is neither of type 'packed' nor 'LonLat'."
             stop
          end select
       end select

       !! copy data to correct variable

       select case(fvarnam(i))
       case('LAI_yDayMean')
          cbalance_diag%LAI_yDayMean(:,:,1:nday) = hlp3D(:,:,1:nday)
          deallocate(hlp3D)
          lai_found=.true.
       case('NPP_yDayMean')
          cbalance_diag%NPP_yDayMean(:,:,1:nday) = hlp3D(:,:,1:nday)
          deallocate(hlp3D)
          npp_found=.true.
       case('topSoilTemp_yDayMean')
          do itile=1,ntiles
             cbalance_diag%topSoilTemp_yDayMean(:,itile,1:nday) = hlp2D(:,1:nday)
          end do
          deallocate(hlp2D)
          temp_found=.true.
       case('alpha_yDayMean')
          cbalance_diag%alpha_yDayMean(:,:,1:nday) = hlp3D(:,:,1:nday)
          deallocate(hlp3D)
          alpha_found=.true.
       case('pseudo_temp_yDay','pseudo_temp2_30d')       ! Meteorological data for Yasso
          if (with_yasso) then 
             do itile=1,ntiles
                cbalance_diag%pseudo_temp_yDay(:,itile,1:nday) = hlp2D(:,1:nday)
             end do
          end if
          deallocate(hlp2D)
          yasso_temp_found=.true.
       case('pseudo_precip_yDay','pseudo_precip_30d')      ! Meteorological data for Yasso
          if (with_yasso) then 
             do itile=1,ntiles
                cbalance_diag%pseudo_precip_yDay(:,itile,1:nday) = hlp2D(:,1:nday)
             end do
          end if
          deallocate(hlp2D)
          yasso_precip_found=.true.
       case('Runoff_yDayMean')
          if(with_nitrogen) then
             do itile=1,ntiles
                nbalance_offline%Runoff_yDayMean(:,itile,1:nday) = hlp2D(:,1:nday)
             end do
          endif
          deallocate(hlp2D)
          runoff_found=.true.
       case('ave_npp5')
          IF (run_dynveg) dynveg_clim%ave_npp5(:,:,1:nday) = hlp3D(:,:,1:nday)
          deallocate(hlp3D)
          ave_npp5_found=.true.
       case('bio_exist')
          IF (run_dynveg) dynveg_clim%bio_exist(:,:,1:nday) = hlp3D(:,:,1:nday)
          deallocate(hlp3D)
          bio_exist_found=.true.
       case('rel_hum_air_yDay','rel_hum_air')
          if (.not. rel_hum_found) then
             if (associated(dynveg_clim%rel_hum_air)) dynveg_clim%rel_hum_air(:,1:nday) = hlp2D(:,1:nday)
          endif
          deallocate(hlp2D)
          rel_hum_found=.true.
       case('prev_day_max_wind10')
          IF (run_dynveg .or. run_disturbance) dynveg_clim%prev_day_max_wind10(:,1:nday) = hlp2D(:,1:nday)
          deallocate(hlp2D)
          prev_day_max_wind10_found=.true.
       case('max_wind10')
          IF (run_dynveg .or. run_disturbance) dynveg_clim%max_wind10(:,1:nday) = hlp2D(:,1:nday)
          deallocate(hlp2D)
          max_wind10_found=.true.
       case('WindSpeed_yDayMean')
          if (.not. prev_day_mean_wind10_found) then 
             if (associated(dynveg_clim%prev_day_mean_wind10)) dynveg_clim%prev_day_mean_wind10(:,1:nday) = hlp2D(:,1:nday)
          endif
          deallocate(hlp2D)
          prev_day_mean_wind10_found = .true.
       case('prev_day_mean_wind10')
          if (.not. prev_day_mean_wind10_found) then 
             if (associated(dynveg_clim%prev_day_mean_wind10)) dynveg_clim%prev_day_mean_wind10(:,1:nday) = hlp2D(:,1:nday)
          endif
          deallocate(hlp2D)
          prev_day_mean_wind10_found = .true.
       case('max_2m_temp')
          if (.not. prev_day_temp_max_found) then 
             if (associated(dynveg_clim%prev_day_temp_max)) dynveg_clim%prev_day_temp_max(:,1:nday)= hlp2D(:,1:nday)
          endif
          deallocate(hlp2D)
          prev_day_temp_max_found=.true. 
       case('prev_day_temp_max')
          if (.not. prev_day_temp_max_found) then 
             if (associated(dynveg_clim%prev_day_temp_max)) dynveg_clim%prev_day_temp_max(:,1:nday)= hlp2D(:,1:nday)
          endif
          deallocate(hlp2D)
          prev_day_temp_max_found=.true. 
       case('min_2m_temp')
          if (.not. prev_day_temp_min_found) then 
             if (associated(dynveg_clim%prev_day_temp_min)) dynveg_clim%prev_day_temp_min(:,1:nday)= hlp2D(:,1:nday)
          endif
          deallocate(hlp2D)
          prev_day_temp_min_found=.true. 
       case('prev_day_temp_min')
          if (.not. prev_day_temp_min_found) then 
             if (associated(dynveg_clim%prev_day_temp_min)) dynveg_clim%prev_day_temp_min(:,1:nday)= hlp2D(:,1:nday)
          endif
          deallocate(hlp2D)
          prev_day_temp_min_found=.true. 
!       case('dew_point_temp')
!          if (.not. dew_point_temp_found) then 
!             if (associated(dynveg_clim%dew_point_temp)) dynveg_clim%dew_point_temp(:,1:nday)= hlp2D(:,1:nday)
!          endif
!          deallocate(hlp2D)
!          dew_point_temp_found=.true.
       case('prev_day_precip_mean','precip')
          if (.not. prev_day_precip_mean_found) then
             if (associated(dynveg_clim%prev_day_precip_mean)) dynveg_clim%prev_day_precip_mean(:,1:nday)=hlp2D(:,1:nday)
          endif
          deallocate(hlp2D)
          prev_day_precip_mean_found= .true.
       case('prev_day_mean_vol_moist')
          if (.not. prev_day_mean_vol_moist_found) then
             if (associated(dynveg_clim%prev_day_mean_vol_moist)) dynveg_clim%prev_day_mean_vol_moist(:,1:nday)=hlp2D(:,1:nday)
          endif
          deallocate(hlp2D)
          prev_day_mean_vol_moist_found= .true.
       end select

    end do

    ! close the file

    iret = nf_close(NCID)
    call check_err(iret,"Error(4) in getDailyData() when closing "//trim(filename))

    ! check whether all data have been found

    if(.not. LAI_found) then
       write(*,*) "ERROR: Did not find variable LAI_yDayMean in "//trim(filename)//"."
       stop
    end if
    if(.not. NPP_found) then
       write(*,*) "ERROR: Did not find variable NPP_yDayMean in "//trim(filename)//"."
       stop
    end if
    if(.not. temp_found) then
       if (with_yasso) then
          ! with yasso only needed for technical reasons
          cbalance_diag%topSoilTemp_yDayMean(:,:,:) = 277._dp
       else
          write(*,*) "ERROR: Did not find variable topSoilTemp_yDayMean in "//trim(filename)//"."
          stop
       end if
    end if
    if(.not. alpha_found) then
       if (with_yasso) then
          ! with yasso only needed for technical reasons
          cbalance_diag%alpha_yDayMean(:,:,:) = 1._dp
       else
          write(*,*) "ERROR: Did not find variable alpha_yDayMean in "//trim(filename)//"."
          stop
       end if
    end if
    if(.not. runoff_found .and. with_nitrogen) then
       write(*,*) "ERROR: Did not find variable Runoff_yDayMean in "//trim(filename)//"."
       stop
    end if
    if(.not. yasso_precip_found .and. with_yasso) then
       write(*,*) "ERROR: Did not find variable pseudo_precip_yDay  in "//trim(filename)//"."
       stop
    end if
    if(.not. yasso_temp_found .and. with_yasso) then
       write(*,*) "ERROR: Did not find variable pseudo_temp_yDay  in "//trim(filename)//"."
       stop
    end if
    IF (run_dynveg) THEN
       if(.not. ave_npp5_found) then
          write(*,*) "ERROR: Did not find variable ave_npp5 in "//trim(filename)//"."
          stop
       end if
       if(.not. bio_exist_found) then
         dynveg_clim%bio_exist = 1._dp
         write(*,*) "WARNING: Did not find variable bio_exist in "//trim(filename)//"."
         write(*,*) "WARNING: Assuming bio_exist = 1"
!StW        stop
       end if
    endif
    if (run_disturbance) then
       if(.not. rel_hum_found) then
          write(*,*) "ERROR: Did not find variable rel_hum_air_yDay in "//trim(filename)//"."
          stop
       end if
       if(.not. prev_day_max_wind10_found) then
          write(*,*) "ERROR: Did not find variable prev_day_max_wind10 in "//trim(filename)//"."
          stop
       end if
       if(.not. max_wind10_found) then
          write(*,*) "ERROR: Did not find variable max_wind10 in "//trim(filename)//"."
          stop
       end if
       if (.not. prev_day_mean_wind10_found) then
          write(*,*) "ERROR: Did not find variable WindSpeed_yDayMean in "//trim(filename)//"."
          stop
       end if
       if (.not.  prev_day_temp_max_found) then
          write(*,*) "ERROR: Did not find variable max_2m_temp/prev_day_temp_max in "//trim(filename)//"."
          stop
       end if
       if (.not.  prev_day_temp_min_found) then
          write(*,*) "ERROR: Did not find variable min_2m_temp/prev_day_temp_min in "//trim(filename)//"."
          stop
       end if
!       if (.not. dew_point_temp_found) then
!          write(*,*) "Warning: Did not find variable dew_point_temp in "//trim(filename)//" using Tmin-4 instead."
!       end if
       if (.not. prev_day_precip_mean_found) then
          write(*,*) "ERROR: Did not find variable precip/prev_day_precip_mean in "//trim(filename)//"."
          stop
       end if

    endif

    ! deallocate local arrays

    deallocate(fdimnam,ndimlen)
    deallocate(fvarnam,nvartyp,nvardim,nvardimid,nvaratt)
  end subroutine getDailyData

  ! --- initNetcdfInfo() -------------------------------------------------------------------------------------------

  SUBROUTINE initNetcdfInfo(grid, cbalance, cbalance_diag, nbalance, nbalance_offline, dynveg, surface, &
                            landcover_change, landuse_transitions, with_nitrogen, with_yasso, run_dynveg, &
                            run_disturbance, dynveg_feedback, use_external_landcover_maps, with_landuse_transitions, lcc_scheme)

    USE mo_cbalone_memory,         ONLY: cbal_offline_type, Nbal_offline_type
    USE mo_cbal_bethy,             ONLY: cbalance_type
    USE mo_cbal_landcover_change,  ONLY: landcover_change_type, &
                                         landuse_transitions_type
    USE mo_cbal_bethy,             ONLY: nbalance_type
    USE mo_land_surface,           ONLY: land_surface_type
    USE mo_dynveg,                 ONLY: dynveg_type
    USE mo_jsbach_grid,            ONLY: grid_type
    USE mo_disturbance,            ONLY: disturbance, dist_opts
    USE mo_disturbance_thonicke,   ONLY: fire_TH_diag, state_FTH

    TYPE(cbal_offline_type),       INTENT(inout) :: cbalance_diag
    TYPE(cbalance_type),           INTENT(inout) :: cbalance
    TYPE(grid_type),               INTENT(in) :: grid
    TYPE(nbalance_type),           INTENT(in) :: nbalance
    TYPE(nbal_offline_type),       INTENT(in) :: nbalance_offline
    TYPE(dynveg_type),             INTENT(in) :: dynveg
    TYPE(land_surface_type),       INTENT(in) :: surface
    TYPE(landcover_change_type),   INTENT(in) :: landcover_change
    TYPE(landuse_transitions_type),INTENT(in) :: landuse_transitions
    LOGICAL,                       INTENT(in) :: with_nitrogen
    LOGICAL,                       INTENT(in) :: with_yasso   
    LOGICAL,                       INTENT(in) :: run_dynveg
    LOGICAL,                       INTENT(in) :: run_disturbance
    LOGICAL,                       INTENT(in) :: dynveg_feedback
    LOGICAL,                       INTENT(in) :: use_external_landcover_maps                  
    LOGICAL,                       INTENT(in) :: with_landuse_transitions
    INTEGER,                       INTENT(in) :: lcc_scheme


    INCLUDE 'netcdf.inc'

    INTEGER :: n



    numberOfOutputVariables = 30
    IF (with_nitrogen)                                             numberOfOutputVariables = numberOfOutputVariables + 52
    
    IF (with_yasso)                                                numberOfOutputVariables = numberOfOutputVariables + 26
    IF (use_external_landcover_maps .OR. with_landuse_transitions) THEN
       SELECT CASE (lcc_scheme)
          CASE (1)
                                                                   numberOfOutputVariables = numberOfOutputVariables +  3
            IF (with_nitrogen)                                     numberOfOutputVariables = numberOfOutputVariables +  3
          CASE (2)
            IF (with_landuse_transitions) THEN
                                                                   numberOfOutputVariables = numberOfOutputVariables + 20 + 12
            ELSE
                                                                   numberOfOutputVariables = numberOfOutputVariables + 20
            ENDIF
       END SELECT
    ENDIF
    IF (with_landuse_transitions)                                  numberOfOutputVariables = numberOfOutputVariables + 23
    IF (dynveg_feedback .OR. with_landuse_transitions)             numberOfOutputVariables = numberOfOutputVariables +  1
    IF (use_external_landcover_maps .AND. with_nitrogen)           numberOfOutputVariables = numberOfOutputVariables +  4
    IF (run_dynveg)                                                numberOfOutputVariables = numberOfOutputVariables +  7
    IF (run_disturbance) THEN
                                                                   numberOfOutputVariables = numberOfOutputVariables +  8
       IF (dist_opts%fire_algorithm == 3) THEN
                                                                   numberOfOutputVariables = numberOfOutputVariables +  7
          IF (dist_opts%ldiag)                                     numberOfOutputVariables = numberOfOutputVariables + 14
          IF (ASSOCIATED(fire_TH_diag%avg_FRP_gridcell))           numberOfOutputVariables = numberOfOutputVariables +  1
          IF (ASSOCIATED(fire_TH_diag%avg_numfires_gridcell))      numberOfOutputVariables = numberOfOutputVariables +  1
       ENDIF
    ENDIF
    IF (debug_Cconservation) THEN
                                                                   numberOfOutputVariables = numberOfOutputVariables +  1
       IF (with_landuse_transitions)                               numberOfOutputVariables = numberOfOutputVariables + 14
       IF (use_external_landcover_maps .OR. with_landuse_transitions) &
                                                                   numberOfOutputVariables = numberOfOutputVariables +  1
       IF (run_dynveg)                                             numberOfOutputVariables = numberOfOutputVariables +  2
    END IF

    allocate(netcdfInfo(1:numberOfOutputVariables))
    n=0

    n=n+1
    netcdfInfo(n)%varName = "Cpool_green"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_tiles => cbalance%Cpool_green

    n=n+1
    netcdfInfo(n)%varName = "Cpool_reserve"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_tiles => cbalance%Cpool_reserve

    n=n+1
    netcdfInfo(n)%varName = "Cpool_woods"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_tiles => cbalance%Cpool_woods

    IF (.NOT. with_yasso) THEN
       n=n+1
       netcdfInfo(n)%varName = "Cpool_litter_green_ag"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%Cpool_litter_green_ag

       n=n+1
       netcdfInfo(n)%varName = "Cpool_litter_wood_ag"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%Cpool_litter_wood_ag

       n=n+1
       netcdfInfo(n)%varName = "Cpool_litter_wood_bg"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%Cpool_litter_wood_bg

       n=n+1
       netcdfInfo(n)%varName = "Cpool_litter_green_bg"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%Cpool_litter_green_bg

       n=n+1
       netcdfInfo(n)%varName = "Cpool_slow"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%Cpool_slow

    ELSE
       n=n+1
       netcdfInfo(n)%varName = "YCpool_acid_ag1"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%YCpool_acid_ag1
 
       n=n+1
       netcdfInfo(n)%varName = "YCpool_water_ag1"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%YCpool_water_ag1

       n=n+1
       netcdfInfo(n)%varName = "YCpool_ethanol_ag1"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%YCpool_ethanol_ag1
 
       n=n+1
       netcdfInfo(n)%varName = "YCpool_nonsoluble_ag1"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%YCpool_nonsoluble_ag1
    
       n=n+1
       netcdfInfo(n)%varName = "YCpool_acid_bg1"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%YCpool_acid_bg1
 
       n=n+1
       netcdfInfo(n)%varName = "YCpool_water_bg1"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%YCpool_water_bg1
 
       n=n+1
       netcdfInfo(n)%varName = "YCpool_ethanol_bg1"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%YCpool_ethanol_bg1
 
       n=n+1
       netcdfInfo(n)%varName = "YCpool_nonsoluble_bg1"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%YCpool_nonsoluble_bg1
    
       n=n+1
       netcdfInfo(n)%varName = "YCpool_humus_1"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%YCpool_humus_1

       n=n+1
       netcdfInfo(n)%varName = "avg_YCpool_acid_ag1"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_YCpool_acid_ag1

       n=n+1
       netcdfInfo(n)%varName = "avg_YCpool_water_ag1"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_YCpool_water_ag1

       n=n+1
       netcdfInfo(n)%varName = "avg_YCpool_ethanol_ag1"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_YCpool_ethanol_ag1

       n=n+1
       netcdfInfo(n)%varName = "avg_YCpool_nonsoluble_ag1"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_YCpool_nonsoluble_ag1
    
       n=n+1
       netcdfInfo(n)%varName = "avg_YCpool_acid_bg1"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_YCpool_acid_bg1

       n=n+1
       netcdfInfo(n)%varName = "avg_YCpool_water_bg1"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_YCpool_water_bg1

       n=n+1
       netcdfInfo(n)%varName = "avg_YCpool_ethanol_bg1"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_YCpool_ethanol_bg1

       n=n+1
       netcdfInfo(n)%varName = "avg_YCpool_nonsoluble_bg1"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_YCpool_nonsoluble_bg1
   
       n=n+1
       netcdfInfo(n)%varName = "avg_YCpool_humus_1"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_YCpool_humus_1


       ! WOODY POOLS (SIZE CLASS 2) 
       n=n+1
       netcdfInfo(n)%varName = "YCpool_acid_ag2"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%YCpool_acid_ag2
 
       n=n+1
       netcdfInfo(n)%varName = "YCpool_water_ag2"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%YCpool_water_ag2

       n=n+1
       netcdfInfo(n)%varName = "YCpool_ethanol_ag2"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%YCpool_ethanol_ag2
 
       n=n+1
       netcdfInfo(n)%varName = "YCpool_nonsoluble_ag2"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%YCpool_nonsoluble_ag2
    
       n=n+1
       netcdfInfo(n)%varName = "YCpool_acid_bg2"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%YCpool_acid_bg2
 
       n=n+1
       netcdfInfo(n)%varName = "YCpool_water_bg2"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%YCpool_water_bg2
 
       n=n+1
       netcdfInfo(n)%varName = "YCpool_ethanol_bg2"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%YCpool_ethanol_bg2
 
       n=n+1
       netcdfInfo(n)%varName = "YCpool_nonsoluble_bg2"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%YCpool_nonsoluble_bg2
    
       n=n+1
       netcdfInfo(n)%varName = "YCpool_humus_2"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance%YCpool_humus_2

       n=n+1
       netcdfInfo(n)%varName = "avg_YCpool_acid_ag2"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_YCpool_acid_ag2

       n=n+1
       netcdfInfo(n)%varName = "avg_YCpool_water_ag2"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_YCpool_water_ag2

       n=n+1
       netcdfInfo(n)%varName = "avg_YCpool_ethanol_ag2"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_YCpool_ethanol_ag2

       n=n+1
       netcdfInfo(n)%varName = "avg_YCpool_nonsoluble_ag2"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_YCpool_nonsoluble_ag2
    
       n=n+1
       netcdfInfo(n)%varName = "avg_YCpool_acid_bg2"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_YCpool_acid_bg2

       n=n+1
       netcdfInfo(n)%varName = "avg_YCpool_water_bg2"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_YCpool_water_bg2

       n=n+1
       netcdfInfo(n)%varName = "avg_YCpool_ethanol_bg2"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_YCpool_ethanol_bg2

       n=n+1
       netcdfInfo(n)%varName = "avg_YCpool_nonsoluble_bg2"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_YCpool_nonsoluble_bg2
    
       n=n+1
       netcdfInfo(n)%varName = "avg_YCpool_humus_2"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_YCpool_humus_2

    END IF

    n=n+1
    netcdfInfo(n)%varName = "LAI_previousDayMean"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_tiles => cbalance_diag%LAI_previousDayMean

    n=n+1
    netcdfInfo(n)%varName = "box_Cpools_total"
    netcdfInfo(n)%dimType = BOX_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_box => cbalance_diag%box_Cpools_total

    n=n+1
    netcdfInfo(n)%varName = "box_test_Ccons"
    netcdfInfo(n)%dimType = BOX_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_box => cbalance_diag%box_test_Ccons

    IF (debug_Cconservation) THEN
       n=n+1
       netcdfInfo(n)%varName = "test_Cconserv"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles =>  cbalance_diag%testCconserv
    END IF

    n=n+1
    netcdfInfo(n)%varName = "avg_Cpool_green"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_tiles => cbalance_diag%avg_Cpool_green

    n=n+1
    netcdfInfo(n)%varName = "avg_Cpool_reserve"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE 
    netcdfInfo(n)%array_tiles => cbalance_diag%avg_Cpool_reserve

    n=n+1
    netcdfInfo(n)%varName = "avg_Cpool_woods"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_tiles => cbalance_diag%avg_Cpool_woods 

    IF (.NOT. with_yasso) THEN
       n=n+1
       netcdfInfo(n)%varName = "avg_Cpool_litter_green_ag"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_Cpool_litter_green_ag

       n=n+1
       netcdfInfo(n)%varName = "avg_Cpool_litter_wood_ag"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_Cpool_litter_wood_ag 

       n=n+1
       netcdfInfo(n)%varName = "avg_Cpool_litter_wood_bg"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_Cpool_litter_wood_bg 

       n=n+1
       netcdfInfo(n)%varName = "avg_Cpool_litter_green_bg"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_Cpool_litter_green_bg

       n=n+1
       netcdfInfo(n)%varName = "avg_Cpool_slow"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => cbalance_diag%avg_Cpool_slow
    END IF

    n=n+1
    netcdfInfo(n)%varName = "avg_soil_respiration"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_tiles =>  cbalance_diag%avg_soil_respiration

    n=n+1
    netcdfInfo(n)%varName = "avg_NPP_yDayMean"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_tiles =>  cbalance_diag%avg_NPP_yDayMean

    n=n+1
    netcdfInfo(n)%varName = "avg_NPP_act_yDayMean"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_tiles =>  cbalance_diag%avg_NPP_act_yDayMean

    n=n+1
    netcdfInfo(n)%varName = "avg_NPP_flux_correction"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_tiles =>  cbalance_diag%avg_NPP_flux_correction

    n=n+1
    netcdfInfo(n)%varName = "avg_excess_NPP"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_tiles =>  cbalance_diag%avg_excess_NPP

    n=n+1
    netcdfInfo(n)%varName = "avg_root_exudates"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_tiles =>  cbalance_diag%avg_root_exudates

    n=n+1
    netcdfInfo(n)%varName = "avg_Cflux_herbivory"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_tiles =>  cbalance_diag%avg_Cflux_herbivory

    n=n+1
    netcdfInfo(n)%varName = "avg_Cflux_herbivory_LG"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_tiles =>  cbalance_diag%avg_Cflux_herbivory_LG

    n=n+1
    netcdfInfo(n)%varName = "avg_Cflux_herbivory_2_atm"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_tiles =>  cbalance_diag%avg_Cflux_herbivory_2_atm

    n=n+1
    netcdfInfo(n)%varName = "avg_box_NEP"
    netcdfInfo(n)%dimType = BOX_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_box =>  cbalance_diag%avg_box_NEP

    n=n+1
    netcdfInfo(n)%varName = "landindex"
    netcdfInfo(n)%dimType = BOX_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_box =>  grid%landindex

    IF (with_nitrogen) THEN
       n=n+1
       netcdfInfo(n)%varName = "Npool_green"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance%Npool_green

       n=n+1
       netcdfInfo(n)%varName = "Npool_mobile"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance%Npool_mobile

       n=n+1
       netcdfInfo(n)%varName = "Npool_woods"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance%Npool_woods

       n=n+1
       netcdfInfo(n)%varName = "Npool_litter_green_ag"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance%Npool_litter_green_ag

       n=n+1
       netcdfInfo(n)%varName = "Npool_litter_wood_ag"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance%Npool_litter_wood_ag

       n=n+1
       netcdfInfo(n)%varName = "Npool_litter_wood_bg"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance%Npool_litter_wood_bg

       n=n+1
       netcdfInfo(n)%varName = "Npool_litter_green_bg"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance%Npool_litter_green_bg
  
       n=n+1
       netcdfInfo(n)%varName = "Npool_slow"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance%Npool_slow

       n=n+1
       netcdfInfo(n)%varName = "avg_Npool_green"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_Npool_green

       n=n+1
       netcdfInfo(n)%varName = "avg_Npool_mobile"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE 
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_Npool_mobile

       n=n+1
       netcdfInfo(n)%varName = "avg_Npool_woods"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_Npool_woods 

       n=n+1
       netcdfInfo(n)%varName = "avg_Npool_litter_green_ag"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_Npool_litter_green_ag

       n=n+1
       netcdfInfo(n)%varName = "avg_Npool_litter_wood_ag"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_Npool_litter_wood_ag 

       n=n+1
       netcdfInfo(n)%varName = "avg_Npool_litter_wood_bg"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_Npool_litter_wood_bg 

       n=n+1
       netcdfInfo(n)%varName = "avg_Npool_litter_green_bg"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_Npool_litter_green_bg

       n=n+1
       netcdfInfo(n)%varName = "avg_Npool_slow"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_Npool_slow

       n=n+1
       netcdfInfo(n)%varName = "SMINN_pool"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance%SMINN_pool

       n=n+1
       netcdfInfo(n)%varName = "avg_SMINN_pool"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_SMINN_pool

       n=n+1
       netcdfInfo(n)%varName = "avg_SMINN_herbivory"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_SMINN_herbivory

       n=n+1
       netcdfInfo(n)%varName = "redFact_Nlimitation"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance%redFact_Nlimitation

       n=n+1
       netcdfInfo(n)%varName = "avg_minNflux_litter_green_ag"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_minNflux_litter_green_ag

       n=n+1
       netcdfInfo(n)%varName = "avg_minNflux_litter_wood_ag"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_minNflux_litter_wood_ag

       n=n+1
       netcdfInfo(n)%varName = "avg_minNflux_litter_wood_bg"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_minNflux_litter_wood_bg

       n=n+1
       netcdfInfo(n)%varName = "avg_minNflux_litter_green_bg"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_minNflux_litter_green_bg

       n=n+1
       netcdfInfo(n)%varName = "avg_minNflux_slow"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_minNflux_slow

       n=n+1
       netcdfInfo(n)%varName = "avg_Nplant_demand"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_Nplant_demand

       n=n+1
       netcdfInfo(n)%varName = "avg_Nsoil_demand"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_Nsoil_demand

       n=n+1
       netcdfInfo(n)%varName = "avg_Ntotal_demand"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_Ntotal_demand

       n=n+1
       netcdfInfo(n)%varName = "ndep_to_sminn"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance%ndep_to_sminn

       n=n+1
       netcdfInfo(n)%varName = "ndep_forc"
       netcdfInfo(n)%dimType = BOX_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_box => nbalance%ndep_forc

       n=n+1
       netcdfInfo(n)%varName = "nfert_forc_2d"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance%nfert_forc_2d

       n=n+1
       netcdfInfo(n)%varName = "nfert_forc"
       netcdfInfo(n)%dimType = BOX_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_box => nbalance%nfert_forc

       n=n+1
       netcdfInfo(n)%varName = "avg_nfert_to_sminn"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_nfert_to_sminn

       n=n+1
       netcdfInfo(n)%varName = "box_Npools_total"
       netcdfInfo(n)%dimType = BOX_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_box => nbalance_offline%box_Npools_total

       n=n+1
       netcdfInfo(n)%varName = "avg_soil_respiration_pot"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles =>  cbalance_diag%avg_soil_respiration_pot
    
       n=n+1
       netcdfInfo(n)%varName = "box_test_Ncons"
       netcdfInfo(n)%dimType = BOX_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_box => nbalance_offline%box_test_Ncons

!!$    n=n+1
!!$    netcdfInfo(n)%varName = "nfix_to_sminn"
!!$    netcdfInfo(n)%dimType = TILES_TYPE
!!$    netcdfInfo(n)%varType = NF_DOUBLE
!!$    netcdfInfo(n)%array_tiles => nbalance%nfix_to_sminn
!!$
!!$    n=n+1
!!$    netcdfInfo(n)%varName = "nfix_to_sminn_npp"
!!$    netcdfInfo(n)%dimType = TILES_TYPE
!!$    netcdfInfo(n)%varType = NF_DOUBLE
!!$    netcdfInfo(n)%array_tiles => nbalance%nfix_to_sminn_npp
!!$
!!$    n=n+1
!!$    netcdfInfo(n)%varName = "nfix_to_sminn_et"
!!$    netcdfInfo(n)%dimType = TILES_TYPE
!!$    netcdfInfo(n)%varType = NF_DOUBLE
!!$    netcdfInfo(n)%array_tiles => nbalance%nfix_to_sminn_et

!!$    n=n+1
!!$    netcdfInfo(n)%varName = "sminn_to_denit"
!!$    netcdfInfo(n)%dimType = TILES_TYPE
!!$    netcdfInfo(n)%varType = NF_DOUBLE
!!$    netcdfInfo(n)%array_tiles => nbalance%sminn_to_denit
!!$
!!$    n=n+1
!!$    netcdfInfo(n)%varName = "sminn_leach"
!!$    netcdfInfo(n)%dimType = TILES_TYPE
!!$    netcdfInfo(n)%varType = NF_DOUBLE
!!$    netcdfInfo(n)%array_tiles => nbalance%sminn_leach

       n=n+1
       netcdfInfo(n)%varName = "avg_nfix_to_sminn"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_nfix_to_sminn

       n=n+1
       netcdfInfo(n)%varName = "avg_sminn_to_denit"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_sminn_to_denit

       n=n+1
       netcdfInfo(n)%varName = "avg_sminn_leach"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_sminn_leach

       n=n+1
       netcdfInfo(n)%varName = "avg_N2O_emissions_ExternalN"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_N2O_emissions_ExternalN
    
       n=n+1
       netcdfInfo(n)%varName = "avg_N2O_emissions_mineraliz"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_N2O_emissions_mineraliz

       n=n+1
       netcdfInfo(n)%varName = "avg_N2O_emissions_nfert"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_N2O_emissions_nfert

       n=n+1
       netcdfInfo(n)%varName = "avg_N2O_emissions_grazing"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_N2O_emissions_grazing
       
       n=n+1
       netcdfInfo(n)%varName = "avg_NetEcosyst_N_flux"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%avg_NetEcosyst_N_flux

       n=n+1
       netcdfInfo(n)%varName = "sum_N_pools"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance%sum_N_pools

!!$    n=n+1
!!$    netcdfInfo(n)%varName = "NPP_act_yDayMean"
!!$    netcdfInfo(n)%dimType = TILES_TYPE
!!$    netcdfInfo(n)%varType = NF_FLOAT
!!$    netcdfInfo(n)%array_tiles =>  cbalance%NPP_act_yDayMean

       n=n+1
       netcdfInfo(n)%varName = "test_Nconserv"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles =>  nbalance_offline%testNconserv

       n=n+1
       netcdfInfo(n)%varName = "test_Ngreen"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%test_Ngreen

       n=n+1
       netcdfInfo(n)%varName = "test_Nwoods"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%test_Nwoods

       n=n+1
       netcdfInfo(n)%varName = "test_Nlitter_green_ag"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%test_Nlitter_green_ag

       n=n+1
       netcdfInfo(n)%varName = "test_Nlitter_green_bg"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%test_Nlitter_green_bg

       n=n+1
       netcdfInfo(n)%varName = "test_Nlitter_wood_ag"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%test_Nlitter_wood_ag

       n=n+1
       netcdfInfo(n)%varName = "test_Nlitter_wood_bg"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%test_Nlitter_wood_bg

       n=n+1
       netcdfInfo(n)%varName = "test_Nslow"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => nbalance_offline%test_Nslow
    END IF

    IF (use_external_landcover_maps .or. with_landuse_transitions) THEN

       IF (debug_Cconservation) THEN
          n=n+1
          netcdfInfo(n)%varName = "LCC_testCconserv"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box => landcover_change%LCC_testCconserv
       END IF

       IF (lcc_scheme==1) THEN

          n=n+1
          netcdfInfo(n)%varName = "LCC_flux_box_C2atmos"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box =>  landcover_change%LCC_flux_box_C2atmos

          n=n+1
          netcdfInfo(n)%varName = "LCC_flux_box_C2litterGreenPools"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box =>  landcover_change%LCC_flux_box_C2litterGreenPools

          n=n+1
          netcdfInfo(n)%varName = "LCC_flux_box_C2litterWoodPool"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box =>  landcover_change%LCC_flux_box_C2litterWoodPool
          
       ENDIF
    
       IF (lcc_scheme==2) THEN    ! LCC after Grand Slam protocol
          
          ! Fluxes to atmosphere

          n=n+1
          netcdfInfo(n)%varName = "LCC_flux_C2atmos"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          IF (with_landuse_transitions) THEN
            netcdfInfo(n)%array_box => landuse_transitions%C2atmos_LUtrans
          ELSE
            netcdfInfo(n)%array_box => landcover_change%C2atmos
          ENDIF

          n=n+1
          netcdfInfo(n)%varName = "LCC_flux_box_C2atmos"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box =>  landcover_change%LCC_flux_box_C2atmos

          n=n+1
          netcdfInfo(n)%varName = "LCC_flux_box_C_onSite_2_atmos"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box =>  landcover_change%boxC_onSite_2_atmos      

          n=n+1
          netcdfInfo(n)%varName = "LCC_flux_box_C_paper_2_atmos"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box =>  landcover_change%boxC_paper_2_atmos

          n=n+1
          netcdfInfo(n)%varName = "LCC_flux_box_C_construction_2_atmos"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box =>  landcover_change%boxC_construction_2_atmos

          n=n+1
          netcdfInfo(n)%varName = "LCC_flux_C_onSite_2_atmos"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box =>  landcover_change%C_onSite_2_atmos

          n=n+1
          netcdfInfo(n)%varName = "LCC_flux_C_paper_2_atmos"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box =>  landcover_change%C_paper_2_atmos

          n=n+1
          netcdfInfo(n)%varName = "LCC_flux_C_construction_2_atmos"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box =>  landcover_change%C_construction_2_atmos

          ! LCC induced transfer to anthro pools
          n=n+1
          netcdfInfo(n)%varName = "LCC_flux_box_C_2_onSite"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box =>  landcover_change%boxC_2_onSite 
          
          n=n+1
          netcdfInfo(n)%varName = "LCC_flux_box_C_2_paper"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box =>  landcover_change%boxC_2_paper 

          n=n+1
          netcdfInfo(n)%varName = "LCC_flux_box_C_2_construction"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box =>  landcover_change%boxC_2_construction 

          n=n+1
          netcdfInfo(n)%varName = "LCC_flux_C_2_onSite"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box =>  landcover_change%C_2_onSite 

          n=n+1
          netcdfInfo(n)%varName = "LCC_flux_C_2_paper"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box =>  landcover_change%C_2_paper 

          n=n+1
          netcdfInfo(n)%varName = "LCC_flux_C_2_construction"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box =>  landcover_change%C_2_construction 

          ! Pools

          n=n+1
          netcdfInfo(n)%varName = "LCC_box_C_onSite"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box => cbalance_diag%avg_Cpool_onSite

          n=n+1
          netcdfInfo(n)%varName = "LCC_box_C_paper"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box => cbalance_diag%avg_Cpool_paper

          n=n+1
          netcdfInfo(n)%varName = "LCC_box_C_construction"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box => cbalance_diag%avg_Cpool_construction

          n=n+1
          netcdfInfo(n)%varName = "LCC_C_onSite"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box => cbalance%Cpool_onSite

          n=n+1
          netcdfInfo(n)%varName = "LCC_C_paper"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box => cbalance%Cpool_paper

          n=n+1
          netcdfInfo(n)%varName = "LCC_C_construction"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box => cbalance%Cpool_construction

          if (with_landuse_transitions) then  
             n=n+1
             netcdfInfo(n)%varName = "LCC_flux_box_C_paper_harvest_2_atmos"
             netcdfInfo(n)%dimType = BOX_TYPE
             netcdfInfo(n)%varType = NF_DOUBLE
             netcdfInfo(n)%array_box =>  landcover_change%boxC_paper_harv_2_atmos

             n=n+1
             netcdfInfo(n)%varName = "LCC_flux_box_C_construction_harvest_2_atmos"
             netcdfInfo(n)%dimType = BOX_TYPE
             netcdfInfo(n)%varType = NF_DOUBLE
             netcdfInfo(n)%array_box =>  landcover_change%boxC_construction_harv_2_atmos

             n=n+1
             netcdfInfo(n)%varName = "LCC_flux_C_paper_harvest_2_atmos"
             netcdfInfo(n)%dimType = BOX_TYPE
             netcdfInfo(n)%varType = NF_DOUBLE
             netcdfInfo(n)%array_box =>  landcover_change%C_paper_harv_2_atmos

             n=n+1
             netcdfInfo(n)%varName = "LCC_flux_C_construction_harvest_2_atmos"
             netcdfInfo(n)%dimType = BOX_TYPE
             netcdfInfo(n)%varType = NF_DOUBLE
             netcdfInfo(n)%array_box =>  landcover_change%C_construction_harv_2_atmos

             n=n+1
             netcdfInfo(n)%varName = "LCC_flux_box_C_2_paper_harvest"
             netcdfInfo(n)%dimType = BOX_TYPE
             netcdfInfo(n)%varType = NF_DOUBLE
             netcdfInfo(n)%array_box =>  landcover_change%boxC_2_paper_harv

             n=n+1
             netcdfInfo(n)%varName = "LCC_flux_box_C_2_construction_harvest"
             netcdfInfo(n)%dimType = BOX_TYPE
             netcdfInfo(n)%varType = NF_DOUBLE
             netcdfInfo(n)%array_box =>  landcover_change%boxC_2_construction_harv 

             n=n+1
             netcdfInfo(n)%varName = "LCC_flux_C_2_paper_harvest"
             netcdfInfo(n)%dimType = BOX_TYPE
             netcdfInfo(n)%varType = NF_DOUBLE
             netcdfInfo(n)%array_box =>  landcover_change%C_2_paper_harv 

             n=n+1
             netcdfInfo(n)%varName = "LCC_flux_C_2_construction_harvest"
             netcdfInfo(n)%dimType = BOX_TYPE
             netcdfInfo(n)%varType = NF_DOUBLE
             netcdfInfo(n)%array_box =>  landcover_change%C_2_construction_harv 

             n=n+1
             netcdfInfo(n)%varName = "LCC_box_C_paper_harvest"
             netcdfInfo(n)%dimType = BOX_TYPE
             netcdfInfo(n)%varType = NF_DOUBLE
             netcdfInfo(n)%array_box => cbalance_diag%avg_Cpool_paper_harvest

             n=n+1
             netcdfInfo(n)%varName = "LCC_box_C_construction_harvest"
             netcdfInfo(n)%dimType = BOX_TYPE
             netcdfInfo(n)%varType = NF_DOUBLE
             netcdfInfo(n)%array_box => cbalance_diag%avg_Cpool_construction_harvest

             n=n+1
             netcdfInfo(n)%varName = "LCC_C_paper_harvest"
             netcdfInfo(n)%dimType = BOX_TYPE
             netcdfInfo(n)%varType = NF_DOUBLE
             netcdfInfo(n)%array_box => cbalance%Cpool_paper_harvest

             n=n+1
             netcdfInfo(n)%varName = "LCC_C_construction_harvest"
             netcdfInfo(n)%dimType = BOX_TYPE
             netcdfInfo(n)%varType = NF_DOUBLE
             netcdfInfo(n)%array_box => cbalance%Cpool_construction_harvest

          ENDIF

       ENDIF

    ENDIF

    IF (with_landuse_transitions) THEN

        n=n+1
        netcdfInfo(n)%varName = "TransMtrx_CROP_2_PAST"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%TransMtrx_CROP_2_PAST
 
        n=n+1
        netcdfInfo(n)%varName = "TransMtrx_PAST_2_CROP"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%TransMtrx_PAST_2_CROP 
 
        n=n+1
        netcdfInfo(n)%varName = "TransMtrx_NATL_2_PAST"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%TransMtrx_NATL_2_PAST 
 
        n=n+1
        netcdfInfo(n)%varName = "TransMtrx_PAST_2_NATL"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%TransMtrx_PAST_2_NATL 
 
        n=n+1
        netcdfInfo(n)%varName = "TransMtrx_CROP_2_NATL"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%TransMtrx_CROP_2_NATL 
 
        n=n+1
        netcdfInfo(n)%varName = "TransMtrx_NATL_2_CROP"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%TransMtrx_NATL_2_CROP 
 
        n=n+1
        netcdfInfo(n)%varName = "TransMtrx_FRST_2_PAST"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%TransMtrx_FRST_2_PAST 
 
        n=n+1
        netcdfInfo(n)%varName = "TransMtrx_PAST_2_FRST"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%TransMtrx_PAST_2_FRST 
 
        n=n+1
        netcdfInfo(n)%varName = "TransMtrx_GRAS_2_PAST"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%TransMtrx_GRAS_2_PAST 
 
        n=n+1
        netcdfInfo(n)%varName = "TransMtrx_PAST_2_GRAS"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%TransMtrx_PAST_2_GRAS 
 
        n=n+1
        netcdfInfo(n)%varName = "TransMtrx_FRST_2_CROP"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%TransMtrx_FRST_2_CROP 
 
        n=n+1
        netcdfInfo(n)%varName = "TransMtrx_CROP_2_FRST"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%TransMtrx_CROP_2_FRST 
 
        n=n+1
        netcdfInfo(n)%varName = "TransMtrx_GRAS_2_CROP"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%TransMtrx_GRAS_2_CROP 
 
        n=n+1
        netcdfInfo(n)%varName = "TransMtrx_CROP_2_GRAS"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%TransMtrx_CROP_2_GRAS 
 
        n=n+1
        netcdfInfo(n)%varName = "Grass_coverFract_lastYear"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%Grass_coverFract_lastYear
        
        n=n+1
        netcdfInfo(n)%varName = "NatWood_coverFract_lastYear"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%NatWood_coverFract_lastYear
        
        n=n+1
        netcdfInfo(n)%varName = "Pasture_coverFract_lastYear"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%Pasture_coverFract_lastYear
 
        n=n+1
        netcdfInfo(n)%varName = "Crop_coverFract_lastYear"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%Crop_coverFract_lastYear
        
        IF (debug_Cconservation) THEN
           n=n+1
           netcdfInfo(n)%varName = "Test_NATL_2_PAST"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => landuse_transitions%Test_NATL_2_PAST
 
           n=n+1
           netcdfInfo(n)%varName = "Test_PAST_2_NATL"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => landuse_transitions%Test_PAST_2_NATL
 
           n=n+1
           netcdfInfo(n)%varName = "Test_CROP_2_NATL"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => landuse_transitions%Test_CROP_2_NATL
 
           n=n+1
           netcdfInfo(n)%varName = "Test_NATL_2_CROP"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => landuse_transitions%Test_NATL_2_CROP
 
           n=n+1
           netcdfInfo(n)%varName = "Test_CROP_2_PAST"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => landuse_transitions%Test_CROP_2_PAST
 
           n=n+1
           netcdfInfo(n)%varName = "Test_PAST_2_CROP"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => landuse_transitions%Test_PAST_2_CROP
 
           n=n+1
           netcdfInfo(n)%varName = "Test_FRST_2_PAST"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => landuse_transitions%Test_FRST_2_PAST
 
           n=n+1
           netcdfInfo(n)%varName = "Test_PAST_2_FRST"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => landuse_transitions%Test_PAST_2_FRST
 
           n=n+1
           netcdfInfo(n)%varName = "Test_GRAS_2_PAST"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => landuse_transitions%Test_GRAS_2_PAST
 
           n=n+1
           netcdfInfo(n)%varName = "Test_PAST_2_GRAS"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => landuse_transitions%Test_PAST_2_GRAS
 
           n=n+1
           netcdfInfo(n)%varName = "Test_FRST_2_CROP"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => landuse_transitions%Test_FRST_2_CROP
 
           n=n+1
           netcdfInfo(n)%varName = "Test_CROP_2_FRST"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => landuse_transitions%Test_CROP_2_FRST
 
           n=n+1
           netcdfInfo(n)%varName = "Test_GRAS_2_CROP"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => landuse_transitions%Test_GRAS_2_CROP
 
           n=n+1
           netcdfInfo(n)%varName = "Test_CROP_2_GRAS"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => landuse_transitions%Test_CROP_2_GRAS
        END IF
 
        n=n+1
        netcdfInfo(n)%varName = "CROP_2_NATL_ignored"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%CROP_2_NATL_ignored
 
        n=n+1
        netcdfInfo(n)%varName = "PAST_2_NATL_ignored"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%PAST_2_NATL_ignored 
 
        n=n+1
        netcdfInfo(n)%varName = "Box_harvest"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%Box_harvest 
 
        n=n+1
        netcdfInfo(n)%varName = "Box_flux_harvest"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%Box_flux_harvest 
 
        n=n+1
        netcdfInfo(n)%varName = "Box_flux_harvest_2atmos"
        netcdfInfo(n)%dimType = BOX_TYPE
        netcdfInfo(n)%varType = NF_DOUBLE
        netcdfInfo(n)%array_box => landuse_transitions%Box_flux_harvest_2atmos

    END IF

    IF (use_external_landcover_maps) then

       if (lcc_scheme == 1) THEN 
          if (with_nitrogen) THEN
             n=n+1
             netcdfInfo(n)%varName = "LCC_flux_box_N2atmos"
             netcdfInfo(n)%dimType = BOX_TYPE
             netcdfInfo(n)%varType = NF_DOUBLE
             netcdfInfo(n)%array_box =>  landcover_change%LCC_flux_box_N2atmos

             n=n+1
             netcdfInfo(n)%varName = "LCC_flux_box_N2litterGreenPools"
             netcdfInfo(n)%dimType = BOX_TYPE
             netcdfInfo(n)%varType = NF_DOUBLE
             netcdfInfo(n)%array_box =>  landcover_change%LCC_flux_box_N2litterGreenPools

             n=n+1
             netcdfInfo(n)%varName = "LCC_flux_box_N2litterWoodPool"
             netcdfInfo(n)%dimType = BOX_TYPE
             netcdfInfo(n)%varType = NF_DOUBLE
             netcdfInfo(n)%array_box =>  landcover_change%LCC_flux_box_N2litterWoodPool

             n=n+1
             netcdfInfo(n)%varName = "LCC_flux_box_N2SMINNpool"
             netcdfInfo(n)%dimType = BOX_TYPE
             netcdfInfo(n)%varType = NF_DOUBLE
             netcdfInfo(n)%array_box =>  landcover_change%LCC_flux_box_N2SMINNpool

          END IF
       endif

    end if

    IF (dynveg_feedback .OR. with_landuse_transitions) THEN
       n=n+1
       netcdfInfo(n)%varName = "cover_fract"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => surface%cover_fract
    END IF

    IF (run_dynveg) THEN
       n=n+1
       netcdfInfo(n)%varName = "cover_fract_pot"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => surface%cover_fract_pot

       n=n+1
       netcdfInfo(n)%varName = "act_fpc"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles =>  dynveg%act_fpc

       n=n+1
       netcdfInfo(n)%varName = "pot_fpc"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles =>  dynveg%pot_fpc

       n=n+1
       netcdfInfo(n)%varName = "bare_fpc"
       netcdfInfo(n)%dimType = BOX_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_box =>  dynveg%bare_fpc

       n=n+1
       netcdfInfo(n)%varName = "desert_fpc"
       netcdfInfo(n)%dimType = BOX_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_box =>  dynveg%desert_fpc

       n=n+1
       netcdfInfo(n)%varName = "max_green_bio"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles =>  dynveg%max_green_bio

       n=n+1
       netcdfInfo(n)%varName = "sum_green_bio_memory"
       netcdfInfo(n)%dimType = BOX_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_box =>  dynveg%sum_green_bio_memory

       IF (debug_Cconservation) THEN
          n=n+1
          netcdfInfo(n)%varName = "dynveg_testCconserv_1"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box => dynveg%dynveg_testCconserv_1

          n=n+1
          netcdfInfo(n)%varName = "dynveg_testCconserv_2"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_DOUBLE
          netcdfInfo(n)%array_box => dynveg%dynveg_testCconserv_2
       END IF

     END IF

     IF (run_disturbance) THEN

       n=n+1
       netcdfInfo(n)%varName = "fuel"
       netcdfInfo(n)%dimType = BOX_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_box => disturbance%fuel

       n=n+1
       netcdfInfo(n)%varName = "burned_frac"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles =>  disturbance%burned_frac

       n=n+1
       netcdfInfo(n)%varName = "damaged_frac"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles =>  disturbance%damaged_frac

       n=n+1
       netcdfInfo(n)%varName = "carbon_2_GreenLitterPools"
       netcdfInfo(n)%dimType = BOX_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_box => disturbance%carbon_2_GreenLitterPools

       n=n+1
       netcdfInfo(n)%varName = "carbon_2_WoodLitterPools"
       netcdfInfo(n)%dimType = BOX_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_box => disturbance%carbon_2_WoodLitterPools

!       n=n+1
!       netcdfInfo(n)%varName = "carbon_2_atmos"
!       netcdfInfo(n)%dimType = BOX_TYPE
!       netcdfInfo(n)%varType = NF_DOUBLE
!       netcdfInfo(n)%array_box => disturbance%carbon_2_atmos

       n=n+1
       netcdfInfo(n)%varName = "box_burned_frac_avg"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => disturbance%box_burned_frac_avg

       n=n+1
       netcdfInfo(n)%varName = "box_damaged_frac_avg"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_tiles => disturbance%box_damaged_frac_avg

       n=n+1
       netcdfInfo(n)%varName = "box_fire_CO2_flux_2_atm"
       netcdfInfo(n)%dimType = BOX_TYPE
       netcdfInfo(n)%varType = NF_DOUBLE
       netcdfInfo(n)%array_box =>  disturbance%box_CO2_flux_2_atmos

       IF (dist_opts%fire_algorithm == 3) THEN
         n=n+1
         netcdfInfo(n)%varName = "box_burned_frac_diag_avg"
         netcdfInfo(n)%dimType = TILES_TYPE
         netcdfInfo(n)%varType = NF_DOUBLE
         netcdfInfo(n)%array_tiles => disturbance%box_burned_frac_diag_avg

         n=n+1
         netcdfInfo(n)%varName = "frac_1hr_wood"
         netcdfInfo(n)%dimType = TILES_TYPE
         netcdfInfo(n)%varType = NF_DOUBLE
         netcdfInfo(n)%array_tiles => state_FTH%frac_1hr_wood

         n=n+1
         netcdfInfo(n)%varName = "frac_10hr_wood"
         netcdfInfo(n)%dimType = TILES_TYPE
         netcdfInfo(n)%varType = NF_DOUBLE
         netcdfInfo(n)%array_tiles => state_FTH%frac_10hr_wood

         n=n+1
         netcdfInfo(n)%varName = "frac_100hr_wood"
         netcdfInfo(n)%dimType = TILES_TYPE
         netcdfInfo(n)%varType = NF_DOUBLE
         netcdfInfo(n)%array_tiles => state_FTH%frac_100hr_wood

         n=n+1
         netcdfInfo(n)%varName = "frac_1000hr_wood"
         netcdfInfo(n)%dimType = TILES_TYPE
         netcdfInfo(n)%varType = NF_DOUBLE
         netcdfInfo(n)%array_tiles => state_FTH%frac_1000hr_wood

         n=n+1
         netcdfInfo(n)%varName = "frac_litter_wood_new"
         netcdfInfo(n)%dimType = TILES_TYPE
         netcdfInfo(n)%varType = NF_DOUBLE
         netcdfInfo(n)%array_tiles => cbalance%frac_litter_wood_new
         
         n=n+1
         netcdfInfo(n)%varName = "NI_acc"
         netcdfInfo(n)%dimType = BOX_TYPE
         netcdfInfo(n)%varType = NF_DOUBLE
         netcdfInfo(n)%array_box => state_FTH%NI_acc
         
         if (dist_opts%ldiag) then
           n=n+1
           netcdfInfo(n)%varName = "avg_FDI"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => fire_TH_diag%avg_FDI

           n=n+1
           netcdfInfo(n)%varName = "avg_ROS"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box =>fire_TH_diag%avg_ROS

           n=n+1
           netcdfInfo(n)%varName = "avg_fuel_1hr_total"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => fire_TH_diag%avg_fuel_1hr_total

           n=n+1
           netcdfInfo(n)%varName = "avg_fuel_10hr_total"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => fire_TH_diag%avg_fuel_10hr_total 

           n=n+1
           netcdfInfo(n)%varName = "avg_fuel_100hr_total"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => fire_TH_diag%avg_fuel_100hr_total

           n=n+1
           netcdfInfo(n)%varName = "avg_HumanIgnition"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => fire_TH_diag%avg_HumanIgnition

           n=n+1
           netcdfInfo(n)%varName = "avg_LightningIgnition"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => fire_TH_diag%avg_LightningIgnition

           n=n+1
           netcdfInfo(n)%varName = "avg_FireDuration"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => fire_TH_diag%avg_FireDuration

           n=n+1
           netcdfInfo(n)%varName = "avg_numfire"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => fire_TH_diag%avg_numfire

           n=n+1
           netcdfInfo(n)%varName = "avg_fire_intensity"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => fire_TH_diag%avg_fire_intensity

           n=n+1
           netcdfInfo(n)%varName = "avg_postfire_mortality"
           netcdfInfo(n)%dimType = TILES_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_tiles => fire_TH_diag%avg_postfire_mortality
           
           n=n+1
           netcdfInfo(n)%varName = "IR_ROS"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => fire_TH_diag%IR_ROS

           n=n+1
           netcdfInfo(n)%varName = "avg_vegetation_height"
           netcdfInfo(n)%dimType = TILES_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_tiles => fire_TH_diag%avg_vegetation_height

           n=n+1
           netcdfInfo(n)%varName = "avg_carbon_2_atmos"
           netcdfInfo(n)%dimType = TILES_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_tiles => fire_TH_diag%avg_carbon_2_atmos

         ENDIF

         IF (ASSOCIATED(fire_TH_diag%avg_FRP_gridcell)) THEN
           n=n+1
           netcdfInfo(n)%varName = "avg_FRP_gridcell"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => fire_TH_diag%avg_FRP_gridcell
         ENDIF

         IF (ASSOCIATED(fire_TH_diag%avg_numfires_gridcell)) THEN
           n=n+1
           netcdfInfo(n)%varName = "avg_numfires_gridcell"
           netcdfInfo(n)%dimType = BOX_TYPE
           netcdfInfo(n)%varType = NF_DOUBLE
           netcdfInfo(n)%array_box => fire_TH_diag%avg_numfires_gridcell
         ENDIF

       ENDIF

    END IF

    if(n /= numberOfOutputVariables) stop 'initNetcdfInfo(): programming error!!'

  end subroutine initNetcdfInfo

  ! --- writeSingleTimeStep() --------------------------------------------------------------------------------------------

  SUBROUTINE writeSingleTimeStep(filename, grid, surface, cbalance, cbalance_diag, nbalance, nbalance_offline, dynveg, &
                                 landcover_change, landuse_transitions, with_nitrogen, with_yasso, run_dynveg,&
                                 run_disturbance, &
                                 dynveg_feedback, use_external_landcover_maps, with_landuse_transitions, lcc_scheme, &
                                 run_year_first, timeStep, restart, finish)

    USE mo_cbalone_memory,        ONLY: cbal_offline_type, Nbal_offline_type
    USE mo_cbal_landcover_change, ONLY: landcover_change_type, &
                                        landuse_transitions_type
    USE mo_cbal_bethy,            ONLY: nbalance_type, cbalance_type
    USE mo_land_surface,          ONLY: land_surface_type
    USE mo_jsbach_grid,           ONLY: grid_type
    USE mo_dynveg,                ONLY: dynveg_type

    INCLUDE 'netcdf.inc'
    CHARACTER(len=*),              INTENT(in) :: filename
    TYPE(grid_type),               INTENT(in) :: grid
    TYPE(land_surface_type),       INTENT(in) :: surface
    TYPE(cbalance_type),           INTENT(inout) :: cbalance
    TYPE(cbal_offline_type),       INTENT(inout) :: cbalance_diag
    TYPE(nbalance_type),           INTENT(in) :: nbalance
    TYPE(Nbal_offline_type),       INTENT(in) :: nbalance_offline
    TYPE(dynveg_type),             INTENT(in) :: dynveg
    TYPE(landcover_change_type),   INTENT(in) :: landcover_change
    TYPE(landuse_transitions_type),INTENT(in) :: landuse_transitions
    LOGICAL,  INTENT(in)                  :: with_nitrogen
    LOGICAL,  INTENT(in)                  :: with_yasso   
    LOGICAL,  INTENT(in)                  :: run_dynveg
    LOGICAL,  INTENT(in)                  :: run_disturbance
    LOGICAL,  INTENT(in)                  :: dynveg_feedback
    LOGICAL,  INTENT(in)                  :: use_external_landcover_maps
    LOGICAL,  INTENT(in)                  :: with_landuse_transitions
    INTEGER,  INTENT(in)                  :: lcc_scheme
    INTEGER,  INTENT(in)                  :: run_year_first
    REAL(dp), INTENT(in)                  :: timeStep !! indicator for time
    LOGICAL,  INTENT(in), OPTIONAL        :: restart  !! If "true": restart file will be written
    LOGICAL,  INTENT(in), OPTIONAL        :: finish   !! If "true": routine closes the output file, deallocates memory and exits.

    !! state variables

    logical,save :: initialized = .false.
    integer,save :: fileId   ! netCDF file ID
    integer,save :: box_type_DimIds(1:3)
    integer,save :: tiles_type_DimIds(1:4)
    integer,save :: dim,dimIds(1:4)
    integer,save :: varId_lat,varId_lon,varId_tiles,varId_time
    real(dp),pointer,save ::  array_4D(:,:,:,:)   !! dummy array: lat x lon x tiles x 1 (1 time step)
    real(dp),pointer,save ::  array_3D(:,:,:)     !! dummy array: lat x lon x 1 (1 time step)
    real    ,pointer,save ::  array_sp(:,:,:)     !! dummy array in single precission
    integer,save :: timeStepCounter


    !! local variables

    integer :: iret                                                   ! error status return
    integer :: dimid_lat(1),dimid_lon(1),dimid_tiles(1),dimid_time(1) ! netCDF dimension IDs
    integer :: ivar,itile,dimId
    integer :: tilesArray(1:surface%ntiles)
    real(dp):: timeArray(1:1)
    integer :: start(1:4),count(1:4)
    logical :: close_file
    CHARACTER(32) :: time_string


    !! Check for finishing

    close_file = .false.
    if (present(restart)) then
      close_file = .true.
    else
      if (present(finish)) then
        close_file = finish
      endif
    endif

    !! Initialization

    if(.not. initialized) then

       timeStepCounter = 0

       !! initialize info structure for netcdf output fields

       call initNetcdfInfo(grid, cbalance, cbalance_diag, nbalance, nbalance_offline, dynveg, surface, &
            landcover_change, landuse_transitions, with_nitrogen, with_yasso, run_dynveg, run_disturbance, &
            dynveg_feedback, use_external_landcover_maps, with_landuse_transitions,lcc_scheme)

       ! control output
       if (debug) write(*,*) "writeSingleTimeStep(): open output file ",trim(filename)

       !! allocate memory for future use
       allocate (array_4D(grid%nlon,grid%nlat,surface%ntiles,1))
       allocate (array_3D(grid%nlon,grid%nlat,1))

       ! open the output-file
       if (debug) write (*,*) 'writeSingleTimeStep(): write Cpools: ',filename
       iret = nf_create(filename,nf_clobber,fileId)
       call check_err(iret,"writeSingleTimeStep(): create filename")
       ! define dimensions
       iret = nf_def_dim(fileId,'lon',grid%nlon,dimId)
       call check_err(iret,"writeSingleTimeStep(): def_dim lon")
       dimid_lon(1) = dimId
       tiles_type_DimIds(1) = dimId
       box_type_DimIds(1)   = dimId

       iret = nf_def_dim(fileId,'lat',grid%nlat,dimId)
       call check_err(iret,"writeSingleTimeStep(): def_dim lat")
       dimid_lat(1) = dimId
       tiles_type_DimIds(2) = dimId
       box_type_DimIds(2)   = dimid

       iret = nf_def_dim(fileId,'tiles',surface%ntiles,dimId)
       call check_err(iret,"writeSingleTimeStep(): def_dim tiles")
       dimid_tiles(1) = dimId
       tiles_type_DimIds(3) = dimId

       iret = nf_def_dim(fileId,'time',NF_UNLIMITED,dimId)
       call check_err(iret,"writeSingleTimeStep(): def_dim time")
       dimid_time(1) = dimId
       tiles_type_DimIds(4) = dimId
       box_type_DimIds(3)   = dimId

       !! set dimensions for arrays

       ! define longitudes, latitudes, tiles and time as variables
       iret = nf_def_var(fileId,'lat',NF_DOUBLE,1,dimid_lat(1:1),varId_lat)
       call check_err(iret,"writeSingleTimeStep(): def_var lat")
       iret = nf_put_att_text(fileId, varId_lat, 'units', 9, 'degrees_N')
       call check_err(iret,"writeSingleTimeStep(): def lat unit")

       iret = nf_def_var(fileId,'lon',NF_DOUBLE,1,dimid_lon(1:1),varId_lon)
       call check_err(iret,"writeSingleTimeStep(): def_var lon")
       iret = nf_put_att_text(fileId, varId_lon, 'units', 9, 'degrees_E')
       call check_err(iret,"writeSingleTimeStep(): def lon unit")

       iret = nf_def_var(fileId,'tiles',NF_INT,1,dimid_tiles(1:1),varId_tiles)
       call check_err(iret,"writeSingleTimeStep(): def_var tiles")

       iret = nf_def_var(fileId,'time',NF_DOUBLE,1,dimid_time(1:1),varId_time)
       call check_err(iret,"writeSingleTimeStep(): def_var time")
!       iret = nf_put_att_text(fileId, varId_time, 'units', 16, 'day as %Y%m%d.%f')  ! unit for absolute time
       write (time_string,'(a,I4.4,a)')  'days since ', run_year_first, '-01-01 00:00:00'
       iret = nf_put_att_text(fileId, varId_time, 'units', LEN(TRIM(time_string)), TRIM(time_string))
       call check_err(iret,"writeSingleTimeStep(): def time unit")

       ! define all other output variables

       do ivar =1,numberOfOutputVariables
          select case(netcdfInfo(ivar)%dimType)
          case(BOX_TYPE)
             dim=3
             dimIds(1:3) = box_type_DimIds(1:3)
          case(TILES_TYPE)
             dim=4
             dimIds(1:4) = tiles_type_DimIds(1:4)
          case default
             stop "writeSingleTimeStep(): non-existent dimension type"
          end select
          iret = nf_def_var(fileId,                                &
                            trim(netcdfInfo(ivar)%varName),        &
                            netcdfInfo(ivar)%varType,              &
                            dim,                                   &
                            dimIds(1:dim),                         &
                            netcdfInfo(ivar)%varId                 )
          call check_err(iret,"writeSingleTimeStep(): def_var "//trim(netcdfInfo(ivar)%varName))
          IF (netcdfInfo(ivar)%varType == nf_double) THEN
             iret = nf_put_att_double(fileId, netcdfInfo(ivar)%varId, '_FillValue', &
                 nf_double,1,nf_fill_double)
          END IF
          call check_err(iret,"writeSingleTimeStep(): put_att double for "//trim(netcdfInfo(ivar)%varName))
       end do

       !! end definition mode

       iret = nf_enddef(FILEID)
       call check_err(IRET)

       !! put variables (except time) that derive from dimensions

       iret = nf_put_var_double(fileId,varid_lat,grid%lat(1:grid%nlat))
       call check_err(iret,"writeSingleTimeStep(): put_var_double lat")

       iret = nf_put_var_double(fileId,varid_lon,grid%lon(1:grid%nlon))
       call check_err(iret,"writeSingleTimeStep(): put_var_double lon")

       do itile=1,surface%ntiles
          tilesArray(itile) = itile
       end do
       iret = nf_put_var_int(fileId,varId_tiles,tilesArray(1:surface%ntiles))
       call check_err(iret,"writeSingleTimeStep(): put_var_int tiles")

       !! remember initialization status

       initialized = .true.

    end if

    !! increase counter for time steps

    timeStepCounter = timeStepCounter + 1

    !! put data

    timeArray(1)=timeStep
    start(1)=timeStepCounter
    count(1)=1
    iret = nf_put_vara_double(fileId,varId_time,start(1:1),count(1:1),timeArray(1:1))
    call check_err(iret,"writeSingleTimeStep(): put_vara_double time")

    do ivar =1,numberOfOutputVariables
       select case(netcdfInfo(ivar)%dimType)
       case(BOX_TYPE)
          start(1:3)=(/ 1,1,timeStepCounter /)
          count(1:3)=(/ grid%nlon,grid%nlat,1 /)
          if(debug) write(*,*) "Box-case: Try to write variable "//trim(netcdfInfo(ivar)%varName)
          array_3D(:,:,1) = UNPACK(netcdfInfo(ivar)%array_box(:),mask=grid%mask,field=nf_fill_double)

          select case(netcdfInfo(ivar)%varType)
          case(NF_FLOAT)
             allocate (array_sp(grid%nlon,grid%nlat,1))
             array_sp(:,:,1)=array_3D(1:grid%nlon,1:grid%nlat,1)
             iret = nf_put_vara_real(fileId,netcdfInfo(ivar)%varId,start(1:3),count(1:3),array_sp(:,:,:))
             call check_err(iret,"writeSingleTimeStep(): put_vara_real"//trim(netcdfInfo(ivar)%varName))
             deallocate (array_sp)
          case(NF_DOUBLE)
             iret = nf_put_vara_double(fileId,netcdfInfo(ivar)%varId,start(1:3),count(1:3),array_3D(1:grid%nlon,1:grid%nlat,1:1))
             call check_err(iret,"writeSingleTimeStep(): put_vara_double"//trim(netcdfInfo(ivar)%varName))
          case default
             stop "Non-existent case lhvoe.jv"
          end select

       case(TILES_TYPE)
          start(1:4) = (/ 1,1,1,timeStepCounter /)
          count(1:4) = (/ grid%nlon,grid%nlat,surface%ntiles,1 /)
          if(debug) write(*,*) "Tiles-case: Try to write variable "//trim(netcdfInfo(ivar)%varName)
          do itile = 1,surface%ntiles
             array_4D(:,:,itile,1) = UNPACK(netcdfInfo(ivar)%array_tiles(:,itile),mask=grid%mask,field=nf_fill_double)
          enddo

          select case(netcdfInfo(ivar)%varType)
          case(NF_FLOAT)
             allocate(array_sp(grid%nlon,grid%nlat,surface%ntiles))
             array_sp = array_4D(1:grid%nlon,1:grid%nlat,1:surface%ntiles,1)
             iret = nf_put_vara_real(fileId, netcdfInfo(ivar)%varId, start(1:4), count(1:4), array_sp(:,:,:))
             call check_err(iret,"Tiles-case: Error happened when trying to write "//trim(netcdfInfo(ivar)%varName))
             deallocate(array_sp)
          case(NF_DOUBLE)
             iret = nf_put_vara_double(fileId, netcdfInfo(ivar)%varId, start(1:4), count(1:4), &
                                       array_4D(1:grid%nlon,1:grid%nlat,1:surface%ntiles,1))
             call check_err(iret,"Tiles-case: Error happened when trying to write "//trim(netcdfInfo(ivar)%varName))
          case default
             stop "Non-existent case .kugb-.lqjb"
          end select

       case default
          stop "writeSingleTimeStep(): non-existent dimension type"
       end select
    end do

    if (close_file) then
          
      !! close the output-file
      iret = nf_close(fileId)
      call check_err(iret)

      !! deallocate memory
      deallocate(array_4D,array_3D)
      deallocate(netcdfinfo)

      !! turn into non-initialized status
      initialized = .false.

    end if

  end subroutine writeSingleTimeStep

end module mo_cbalone_io
