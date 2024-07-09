!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_data_scaling
  !!
  !! This module provides routines for scaling of driver data (NPP,LAI,topSoilTemp,alpha). Such a scaling is needed for
  !! transient runs of cbalone-offline when such input data can be provided only for the first and last year of the
  !! transient runs. By the routines in this module intermediate driver data are constructed by using
  !!  (i) the difference between climatologies of the first and last year, and
  !! (ii) a curve for the development of the CO2 concentration between the two years.
  !! (These data have to be provided as external files)
  !! The scaling is performed by adding to the driver data the climatological differences, but scaled according to the
  !! CO2 concentration.
  !!

  USE mo_kind,            ONLY: dp
  USE mo_jsbach_grid,     ONLY: grid_type
  USE mo_cbalone_memory,  ONLY: check_err

  implicit none
  
  private

  public :: initDataScaling  !! Initializes scaling of input data, i.e. allocates memory, reads in CO2 curve and 
                             !! also data for climatological differences
  public :: scaleDailyData   !! Performs the scaling of the data

                                 !! also data for nitrogen deposition differences


  TYPE data_scaling_type     
     !! The following switch indicate whether input data (LAI,NPP,topSoilTemp,alpha) shall be scaled with external CO2-concentration
     !! for performing transient runs. If set to .true. in run.def, two additional files are needed: One containing the yarly
     !! CO2-concentrations from the first year until the last year, and another file containing the difference between the
     !! climatologies of the first and the last year.

     logical      :: do_input_scaling = .false.

     real(dp),pointer :: global_CO2_conc(:)      !! The global CO2-concentration read in from an external data file (for each year),
                                                 !! index is the year (e.g. global_CO2_conc(2000) is the CO2 concentration in 2000)
     real(dp)         :: CO2_conc_past           !! The CO2-concentration used for generating the (climate) driver data
     real(dp)         :: CO2_conc_recent         !! The CO2-concentration used for generating the more recent climate data, from
                                                 !! which, together with the climate driver data the climatological differences are
                                                 !! computed.
     real(dp)         :: CO2_conc_range          !! = CO2_conc_recent - CO2_conc_past

     !! The following fields contain the difference between the climatologies of the first and the last year of the driver data
     !!    (monthly values)
     real(dp),pointer :: diff_LAI_yDayMean(:,:,:)            !! climatology difference of LAI yesterday (nland x ntiles x 12months)
     real(dp),pointer :: diff_NPP_yDayMean(:,:,:)            !! climatology difference of NPP-Rate yesterday 
                                                             !!    [mol(CO2)/(m^2(canopy) s)] (nland x ntiles x 12months)
     real(dp),pointer :: diff_topSoilTemp_yDayMean(:,:,:)    !! climatology difference of upper layer soil temperature yesterday [K]
                                                             !!    (nland x ntiles x 12months)
     real(dp),pointer :: diff_alpha_yDayMean(:,:,:)          !! climatology difference of water stress coefficient alpha yesterday
                                                             !!    (nland x ntiles x 12months)
     real(dp),pointer :: diff_Ndeposition(:)                 !! difference of N Deposition mg per year 
  end TYPE data_scaling_type
  public :: data_scaling_type

  !! parameters

  logical,parameter :: debug=.false.

contains

  !! --- initDataScaling() -----------------------------------------------------------------------
  !!
  !! Initialization for scaling of climate input data
  subroutine initDataScaling(CO2_file,climatology_diff_file,grid,ntiles,ref_year_past,ref_year_recent, &
                                                                        FirstYear,LastYear,data_scaling)
    character(len=*),        intent(in)    :: CO2_file               !! File containing yearly values of CO2-concentration between
                                                                     !! first and last year of run (ASCII)
    character(len=*),        intent(in)    :: climatology_diff_file  !! File containing differences between climatologies of driver
                                                                     !! data between last and first year of run (netcdf)       
    type(grid_type),         intent(in)    :: grid
    integer,                 intent(in)    :: ntiles                 !! The number of tiles used in this run
    integer,                 intent(in)    :: FirstYear,LastYear     !! First and last year of transient run (according to run.def)
    integer,                 intent(in)    :: ref_year_past          !! Reference years for which the more past and more recent
    integer,                 intent(in)    :: ref_year_recent        !! climatologies used for the climatological differences have
                                                                     !! been produecs. The associated CO2-concentrations (from the
                                                                     !! CO2-file) are used to determine the scaling factor (see
                                                                     !! components CO2_conc_past, CO2_conc_recent and CO2_conc_range
                                                                     !! in the structure data_scaling).
    type(data_scaling_type), intent(inout) :: data_scaling           !! Structure that is filled with the input data

    !! locals
    real,pointer :: CO2_hlp(:,:)
    integer :: n, year,count

    !! Read in CO2 data

    write(*,*) "initDataScaling(): Read CO2 data from "//trim(CO2_file)
    call readRealData(CO2_file,1,2,CO2_hlp) !! CO2_hlp is allocated within this routine with appropriate length

    !! Check the time info of CO2 data

    if(size(CO2_hlp,DIM=1) < 1+LastYear-FirstYear) then
       write(*,*) "initDataScaling(): ERROR: CO2-data from "//trim(CO2_file)//&
                            " contains not enough  data points (=",size(CO2_hlp,DIM=1),&
                            ") for the intended run over ",1+LastYear-FirstYear," years."
       stop
    end if

    if(size(CO2_hlp,DIM=1) > 1+LastYear-FirstYear) then
       write(*,*) "initDataScaling(): WARNING: CO2-data from "//trim(CO2_file)//&
                            " contain more data than needed."
    end if

    !! Fill in data-scaling structure

    count=0
    allocate(data_scaling%global_CO2_conc(FirstYear:LastYear))
    data_scaling%global_CO2_conc(FirstYear:LastYear) = -1.0
    data_scaling%CO2_conc_past = -1.
    data_scaling%CO2_conc_recent = -1.
    do n=1,size(CO2_hlp,DIM=1)
       year = nint(CO2_hlp(n,1))

       if(year == ref_year_past)   data_scaling%CO2_conc_past   = CO2_hlp(n,2)
       if(year == ref_year_recent) data_scaling%CO2_conc_recent = CO2_hlp(n,2)
       
       if(FirstYear <= year .and. year <= LastYear) then
          count=count+1
          data_scaling%global_CO2_conc(year) = CO2_hlp(n,2)
       end if

    end do

    !! Check data for more past and more recent CO2-concentrations

    if(data_scaling%CO2_conc_past < 0.0) then
       write(*,*) "initDataScaling(): ERROR: Found no CO2-concentration for more past reference year (", &
                                      ref_year_past," in "//trim(CO2_file)
       stop
    endif

    if(data_scaling%CO2_conc_recent < 0.0) then
       write(*,*) "initDataScaling(): ERROR: Found no CO2-concentration for more recent reference year (", &
                                      ref_year_recent," in "//trim(CO2_file)
       stop
    endif

    !! Set refence CO2-concentration range

    data_scaling%CO2_conc_range = data_scaling%CO2_conc_recent - data_scaling%CO2_conc_past 
    write(*,*) "initDataScaling(): The reference CO2-concentrations for input scaling are:"
    write(*,*) "                   year ",ref_year_past,": ",data_scaling%CO2_conc_past," ppm"
    write(*,*) "                   year ",ref_year_recent,": ",data_scaling%CO2_conc_recent," ppm"

    !! Another check
    if(data_scaling%CO2_conc_range < 1.0_dp) then
       write(*,*) "initDataScaling(): Difference between the two reference CO2 values smaller than 1 ppm. THIS MAKES NO SENSE!"
       stop
    end if

    !! Check completeness of data

    do n=FirstYear,LastYear
       if(data_scaling%global_CO2_conc(n) < 0.0) then
          write(*,*) "initDataScaling(): ERROR: In "//trim(CO2_file)//" CO2-data are missing for years:"
          do year=n,LastYear
             if(data_scaling%global_CO2_conc(n) < 0.0) write(*,*) year
          end do
          write(*,*) "Possible reason: Entries for keywords RUN_YEAR_FIRST and RUN_YEAR_LAST in run.def"
          write(*,*) ".. may not fit to entries in CO2 data file "//trim(CO2_file)
          STOP 
       end if
    end do

    !! Deallocate temporary memory

    deallocate(CO2_hlp)

    !! Allocate memory for climatology differences

    allocate(data_scaling%diff_LAI_yDayMean(grid%nland,ntiles,12))
    allocate(data_scaling%diff_NPP_yDayMean(grid%nland,ntiles,12))
    allocate(data_scaling%diff_topSoilTemp_yDayMean(grid%nland,ntiles,12))
    allocate(data_scaling%diff_alpha_yDayMean(grid%nland,ntiles,12))

    !! Read the climatology difference file

    call getClimatologyDiffs(climatology_diff_file,grid,ntiles,data_scaling)

    !! finish

    write(*,*) "initDataScaling(): Data scaling initialized."

  end subroutine initDataScaling


  subroutine scaleDailyData(nday,data_scaling,year,month,cbalance_diag)
    USE mo_cbalone_memory, only: cbal_offline_type
    integer,                intent(in)    :: nday
    type(data_scaling_type),intent(in)    :: data_scaling
    integer,                intent(in)    :: year
    integer,                intent(in)    :: month
    type(cbal_offline_type),intent(inout) :: cbalance_diag

    !! locals

    integer  :: n
    real(dp) :: scalingFactor

    !! preparations

    scalingFactor = (data_scaling%global_CO2_conc(year) - data_scaling%CO2_conc_past)/data_scaling%CO2_conc_range

    !! scale it

    do n=1,nday
       cbalance_diag%LAI_yDayMean(:,:,n) = &
            cbalance_diag%LAI_yDayMean(:,:,n) + scalingFactor * real(data_scaling%diff_LAI_yDayMean(:,:,month),KIND=dp)
       cbalance_diag%NPP_yDayMean(:,:,n) = &
            cbalance_diag%NPP_yDayMean(:,:,n) + scalingFactor * real(data_scaling%diff_NPP_yDayMean(:,:,month),KIND=dp)
       cbalance_diag%topSoilTemp_yDayMean(:,:,n) = cbalance_diag%topSoilTemp_yDayMean(:,:,n) &
                                              + scalingFactor * real(data_scaling%diff_topSoilTemp_yDayMean(:,:,month),KIND=dp)
       cbalance_diag%alpha_yDayMean(:,:,n) = &
            cbalance_diag%alpha_yDayMean(:,:,n) + scalingFactor * real(data_scaling%diff_alpha_yDayMean(:,:,month),KIND=dp)
    end do

    !! Assure that LAI remains positive

    cbalance_diag%LAI_yDayMean(:,:,1:nday) = max(0.0_dp,cbalance_diag%LAI_yDayMean(:,:,1:nday))

    !! Assure that alpha remains in [0,1]
 
   cbalance_diag%alpha_yDayMean(:,:,1:nday) = max(0.0_dp,min(1._dp,cbalance_diag%alpha_yDayMean(:,:,1:nday)))
    
    
  end subroutine scaleDailyData


  ! --- readRealData() --------------------------------------------------------------------------
  !
  ! Reads columns of real data into an array. Reads all columns from "firstColumn" to "lastColumn".
  ! Reading the data file, two types of lines are distinguished:
  !  (i) comment lines, indicated by "#" in the first column, or
  ! (ii) data lines, where the data are arranged in columns, separated by white spaces (blanks, tabs)
  ! and no other types of lines.
  ! 
  subroutine readRealData(dataFile,firstColumn,lastColumn,data,unitNo)
    implicit none

    character(len=*),intent(in)      :: dataFile    !! Name of the data file
    integer,intent(in)               :: firstColumn !! Number of first column from which the data shall be read
    integer,intent(in)               :: lastColumn  !! Number of last column from which the data shall be read
    real,pointer                     :: data(:,:)   !! Array that on return contains the data. Will
                                                    !! .. be allocated with appropriate length (size(data,1) 
                                                    !! .. will on return be equal to the number of data points and
                                                    !! .. size(data,2) gives the number of columns.)
    integer,intent(in),optional      :: unitNo      !! If present, the unit number to be used for reading the
                                                    !! .. the data. Otherwise a default unit number is used.

    !! --- parameters

    integer,parameter :: maxLineLength = 256      !! maximum length of the lines in the data file
    integer,parameter :: defaultUnitNo = 99       !! unit number to be used for reading data

    !! --- variables

    integer                      :: io_status
    character(len=maxLineLength) :: line
    real                         :: inData(1:lastColumn)
    integer                      :: i,dataCount,inUnit

    !! -- determine unit number

    inUnit = defaultUnitNo
    if(present(unitNo)) inUnit = unitNo

    !! -- open data file

    open(inUnit,FILE=dataFile,ACTION='READ',IOSTAT=io_status)
    if(io_status .ne. 0) then
       write(*,*) "readRealData(): ERROR: Opening of ",trim(dataFile)," failed."
       stop
    endif

    !! --- count number of lines with data

    dataCount =0
    do
       read(inUnit,'(A)',IOSTAT=io_status) line
       if(io_status .ne. 0) exit    !! detect end of file
       if(line(1:1) .eq. '#') cycle !! ignore comment lines
       dataCount = dataCount + 1
    end do

    !! --- allocate memory for data array

    allocate(data(1:dataCount,1:lastColumn-firstColumn+1))

    !! --- read data

    rewind(inUnit) !! go back to begin of file
    i=0
    do
       read(inUnit,'(A)',IOSTAT=io_status) line
       if(io_status .ne. 0) exit       !! detect end of file
       if(line(1:1) .eq. '#') cycle    !! ignore comment lines
       read(line,*) inData(1:lastColumn) !! read all columns until the one we want to know
       i=i+1
       data(i,:) = inData(firstColumn:lastColumn)  !! We want only the data from firstColumn to lastColumn
    end do

    !! close data file

    close(inUnit)

  end subroutine readRealData

  ! --- getClimatologyDiffs()  -----------------------------------------------------------------------------------
  ! This routine reads in the difference between climatologies of the first and last year, for each month one value.
  ! In the file the variables are assumed to have the names LAI_yDayMean, NPP_yDayMean, topSoilTemp_yDayMean and  
  ! alpha_yDayMean. The routine puts the data into the structure data_scaling.
  ! The routine automatiocally detects whether data are packed or in LonLat format.
  !
  subroutine getClimatologyDiffs(filename,grid,ntiles,data_scaling)
    include 'netcdf.inc'
    character(len=*),intent(in)            :: filename
    type(grid_type),intent(in)             :: grid
    integer,intent(in)                     :: ntiles
    type(data_scaling_type),intent(inout)  :: data_scaling 

    ! local variables

    logical :: lai_found,npp_found,temp_found,alpha_found

    real,allocatable :: hlp2D(:,:),hlp3D(:,:,:),hlp4D(:,:,:,:)

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
    integer  ::  i,ii,itile,month,noOfMonths
    ! others
    real(dp) :: min,max

    ! preparations

    lai_found = .false.
    npp_found = .false.
    temp_found = .false.
    alpha_found = .false.
    
    ! open the input-file

    if(debug) write (*,*) 'getClimatologyDiffs(): read from: ',trim(filename)
    iret = nf_open(filename,nf_nowrite,ncid)
    call check_err(iret,"Error in getClimatologyDiffs() when opening file "//trim(filename))

    ! check what is in the input-file

    iret = nf_inq(ncid,ndim,nvar,nattr,nunlimdim)
    call check_err(iret,"Error in getClimatologyDiffs() when inquiring contents of "//trim(filename))
    if (debug) write (*,*) 'getClimatologyDiffs(): ndim,nvar,nattr,nunlimdim',ndim,nvar,nattr,nunlimdim

    ! get the dimension name and length

    allocate (ndimlen(ndim))
    allocate (fdimnam(ndim))
    do i=1,ndim
       iret = nf_inq_dim(ncid,i,fdimnam(i),ndimlen(i))
       call check_err(iret,"Error in getClimatologyDiffs() when inquiring dimension of "//trim(fdimnam(i))//" in "//trim(filename))
       if (debug) write (*,*) 'getClimatologyDiffs(): i,fdimnam,ndimlen',i,trim(fdimnam(i)),ndimlen(i)
    end do

    ! set ndate (number of output intervals of the monthly file) and noutput_per_month

    noOfMonths = 0
    do i=1,ndim
       if (fdimnam(i) == 'time') noOfMonths = ndimlen(i)
    end do
    if (debug) write (*,*) 'getClimatologyDiffs(): noOfMonths: ',noOfMonths
    if (noOfMonths == 0) then
       write(*,*) 'getClimatologyDiffs(): variable time not found in '//trim(filename)
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
       if (debug) write (*,*) 'getClimatologyDiffs(): i,fvarnam,nvartyp,nvardim,nvardimid,nvaratt',&
                                                         i,trim(fvarnam(i)),nvartyp(i),nvardim(i),nvardimid(i,1),nvaratt(i)
    end do

    ! get data
    
    if(debug) write (*,*) "getClimatologyDiffs(): Try to read the data .."
    do i=1,nvar
       select case (trim(fvarnam(i)))
       case('LAI_yDayMean','NPP_yDayMean','alpha_yDayMean')
          if(debug) write (*,*) "getClimatologyDiffs(): Case LAI_yDayMean, NPP_yDayMean or alpha_yDayMean"
          select case (nvardim(i))
          case(3) !! The array has 3 dimensions => We find packed data
             if(debug) write (*,*) "getClimatologyDiffs(): Case nvardim=3 ==> packed data"
             if (ndimlen(nvardimid(i,1)) /= grid%nland .or. ndimlen(nvardimid(i,2)) /= ntiles) then
                write(*,*) "getClimatologyDiffs(): Dimension error for variable "//trim(fvarnam(i))
                write(*,*) "getClimatologyDiffs(): nland=",grid%nland," but ndimlen(nvardimid(i,1)=",ndimlen(nvardimid(i,1))
                write(*,*) "getClimatologyDiffs(): ntiles=",ntiles," but ndimlen(nvardimid(i,2)=",ndimlen(nvardimid(i,2))
                stop
             end if
             if (ndimlen(nvardimid(i,3)) /= noOfMonths) then
                write(*,*) "getClimatologyDiffs(): Number of time steps wrong for input variable "//trim(fvarnam(i))
                stop
             end if
             allocate (hlp3D(grid%nland,ntiles,noOfMonths))
             iret = nf_get_var_real(ncid,i,hlp3D)
             call check_err(iret,"Error(1) in getClimatologyDiffs() when reading variable "//&
                                                                              trim(fvarnam(i))//" from "//trim(filename))
          case(4)!! The array has 4 dimensions => We find latlon data
             if(debug) write (*,*) "getClimatologyDiffs(): Case nvardim=4 ==> latLon-data"
             if(ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat &
                  .or. ndimlen(nvardimid(i,3)) /= ntiles) then
                write(*,*) "getClimatologyDiffs(): Dimension error for variable "//&
                                                                               trim(fvarnam(i))//" in file "//trim(filename)//":"
                write(*,*) "getClimatologyDiffs(): nlon=",grid%nlon," but found nndimlen(nvardimid(i,1))=",ndimlen(nvardimid(i,1))
                write(*,*) "getClimatologyDiffs(): nlat=",grid%nlat," but found nndimlen(nvardimid(i,2))=",ndimlen(nvardimid(i,2))
                write(*,*) "getClimatologyDiffs(): ntiles=",ntiles," but ndimlen(nvardimid(i,3)=",ndimlen(nvardimid(i,3))
                stop
             end if
             if (ndimlen(nvardimid(i,4)) /= noOfMonths) then
                write(*,*) "getClimatologyDiffs(): Number of time steps wrong for input variable "//&
                     trim(fvarnam(i))//" in file "//trim(filename)//":"
                stop
             end if
             allocate (hlp4D(grid%nlon,grid%nlat,ntiles,noOfMonths))
             iret = nf_get_var_real(ncid,i,hlp4D)
             call check_err(iret,"Error(2)in getClimatologyDiffs() when reading variable "//&
                                                                                   trim(fvarnam(i))//" from "//trim(filename))
             ! pack the array
             allocate (hlp3D(grid%nland,ntiles,noOfMonths))
             do itile=1,ntiles
                do month=1,noOfMonths
                   hlp3D(:,itile,month) = pack(hlp4D(:,:,itile,month),mask=grid%mask)
                end do
             end do
             deallocate(hlp4D)
          case default
             write(*,*) "getClimatologyDiffs(): ERROR: Variable "//trim(fvarnam(i))//&
                          " from file "//trim(filename)//" is neither of type 'packed' nor 'LonLat'."
             stop
          end select
       case('topSoilTemp_yDayMean')
          if(debug) write (*,*) "getClimatologyDiffs(): Case topSoilTemp_yDayMean"
          select case (nvardim(i))
          case(2) !! The array has 2 dimensions => We find packed data
             if(debug) write (*,*) "getClimatologyDiffs(): Case nvardim=2 ==> packed data"
             if (ndimlen(nvardimid(i,1)) /= grid%nland) then
                write(*,*) "getClimatologyDiffs(): Dimension error for variable "//&
                                                                                trim(fvarnam(i))//" in file "//trim(filename)//":"
                write(*,*) "getClimatologyDiffs(): nland=",grid%nland," but ndimlen(nvardimid(i,1)=",ndimlen(nvardimid(i,1))
                stop
             end if
             if (ndimlen(nvardimid(i,2)) /= noOfMonths) then
                write(*,*) "getClimatologyDiffs(): Number of time steps wrong for input variable "//trim(fvarnam(i))//&
                     " in file "//trim(filename)//":"
                write(*,*) "getClimatologyDiffs(): noOfMonths=",noOfMonths," but ndimlen(nvardimid(i,2))=",ndimlen(nvardimid(i,2))
                stop
             end if
             allocate (hlp2D(grid%nland,noOfMonths))
             iret = nf_get_var_real(ncid,i,hlp2D)
             call check_err(iret,"Error(3) in getClimatologyDiffs() when reading variable "//&
                                                                                       trim(fvarnam(i))//" from "//trim(filename))
          case(3)!! The array has 3 dimensions => We find latlon data
             if(debug) write (*,*) "getClimatologyDiffs(): Case nvardim=3 ==> latlon data"
             if(ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) then
                write(*,*) "getClimatologyDiffs(): Dimension error for variable "//&
                                                                                 trim(fvarnam(i))//" in file "//trim(filename)//":"
                write(*,*) "getClimatologyDiffs(): nlon=",grid%nlon," but found nndimlen(nvardimid(i,1))=",ndimlen(nvardimid(i,1))
                write(*,*) "getClimatologyDiffs(): nlat=",grid%nlat," but found nndimlen(nvardimid(i,2))=",ndimlen(nvardimid(i,2))
                stop
             end if
             if (ndimlen(nvardimid(i,3)) /= noOfMonths) then
                write(*,*) "getClimatologyDiffs(): Number of time steps wrong for input variable "//&
                     trim(fvarnam(i))//" in file "//trim(filename)//":"
                write(*,*) "getClimatologyDiffs(): noOfMonths=",noOfMonths," but ndimlen(nvardimid(i,3))=",ndimlen(nvardimid(i,3))
                stop
             end if
             allocate (hlp3D(grid%nlon,grid%nlat,noOfMonths))
             iret = nf_get_var_real(ncid,i,hlp3D)
             call check_err(iret,"Error(4) in getClimatologyDiffs() when reading variable "//&
                                                                                       trim(fvarnam(i))//" from "//trim(filename))
             ! pack the array
             allocate (hlp2D(grid%nland,noOfMonths))
             do month=1,noOfMonths
                hlp2D(:,month) = pack(hlp3D(:,:,month),mask=grid%mask)
             end do
             deallocate(hlp3D)
          case(4)!! The array has 4 dimensions => We find latlon data, but tile-dimension is 1
             if(debug) write (*,*) "getClimatologyDiffs(): Case nvardim=4 ==> latlon data"
             if(ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) then
                write(*,*) "getClimatologyDiffs(): Dimension error for variable "//&
                                                                                 trim(fvarnam(i))//" in file "//trim(filename)//":"
                write(*,*) "getClimatologyDiffs(): nlon=",grid%nlon," but found nndimlen(nvardimid(i,1))=",ndimlen(nvardimid(i,1))
                write(*,*) "getClimatologyDiffs(): nlat=",grid%nlat," but found nndimlen(nvardimid(i,2))=",ndimlen(nvardimid(i,2))
                stop
             end if
             if (ndimlen(nvardimid(i,4)) /= noOfMonths) then
                write(*,*) "getClimatologyDiffs(): Number of time steps wrong for input variable "//&
                     trim(fvarnam(i))//" in file "//trim(filename)//":"
                write(*,*) "getClimatologyDiffs(): noOfMonths=",noOfMonths," but ndimlen(nvardimid(i,4))=",ndimlen(nvardimid(i,4))
                stop
             end if
             if (ndimlen(nvardimid(i,3)) /= 1) then
                write(*,*) "getClimatologyDiffs(): Number of tiles wrong for input variable "//&
                     trim(fvarnam(i))//" in file "//trim(filename)//":"
                write(*,*) "getClimatologyDiffs(): should be 1 but is ",ndimlen(nvardimid(i,3))
                stop
             end if
             allocate (hlp4D(grid%nlon,grid%nlat,1,noOfMonths))
             iret = nf_get_var_real(ncid,i,hlp4D)
             call check_err(iret,"Error(4) in getClimatologyDiffs() when reading variable "//&
                                                                                       trim(fvarnam(i))//" from "//trim(filename))
             ! pack the array
             allocate (hlp2D(grid%nland,noOfMonths))
             do month=1,noOfMonths
                hlp2D(:,month) = pack(hlp4D(:,:,1,month),mask=grid%mask)
             end do
             deallocate(hlp4D)
          case default
             write(*,*) "getClimatologyDiffs(): ERROR: Variable "//trim(fvarnam(i))//&
                  " from file "//trim(filename)//" is neither of type 'packed' nor 'LonLat'."
             stop
          end select
       end select

       !! copy data to correct variable

       select case(fvarnam(i))
       case('LAI_yDayMean')
          data_scaling%diff_LAI_yDayMean(:,:,1:noOfMonths) = hlp3D(:,:,1:noOfMonths)
          deallocate(hlp3D)
          lai_found=.true.
       case('NPP_yDayMean')
          data_scaling%diff_NPP_yDayMean(:,:,1:noOfMonths) = hlp3D(:,:,1:noOfMonths)
          deallocate(hlp3D)
          npp_found=.true.
       case('topSoilTemp_yDayMean')
          do itile=1,ntiles
             data_scaling%diff_topSoilTemp_yDayMean(:,itile,1:noOfMonths) = hlp2D(:,1:noOfMonths)
          end do
          deallocate(hlp2D)
          temp_found=.true.
       case('alpha_yDayMean')
          data_scaling%diff_alpha_yDayMean(:,:,1:noOfMonths) = hlp3D(:,:,1:noOfMonths)
          deallocate(hlp3D)
          alpha_found=.true.
       end select

       !! do simple plausibility checks for the data

       select case(fvarnam(i))
       case('LAI_yDayMean')
          min=minval(data_scaling%diff_LAI_yDayMean)
          max=maxval(data_scaling%diff_LAI_yDayMean)
          if(min < -10.) then
             write(*,*) "getClimatologyDiffs(): ERROR: Found nonsense minimum value for diff_LAI_yDayMean ",min,&
                  " at (index,tile,month) ",minloc(data_scaling%diff_LAI_yDayMean)
             write(*,*) "Possible reason: inconsistent land-sea masks between different input files."
             stop
          end if
          if(max > 10.) then
             write(*,*) "getClimatologyDiffs(): ERROR: Found nonsense maximum value for diff_LAI_yDayMean ",max,&
                  " at (index,tile,month) ",maxloc(data_scaling%diff_LAI_yDayMean)
             write(*,*) "Possible reason: inconsistent land-sea masks between different input files."
             stop
          end if
       case('NPP_yDayMean')
          min=minval(data_scaling%diff_NPP_yDayMean)
          max=maxval(data_scaling%diff_NPP_yDayMean)
          if(min < -1E-4) then
             write(*,*) "getClimatologyDiffs(): ERROR: Found nonsense minimum value for diff_NPP_yDayMean ",min,&
                  " at (index,tile,month) ",minloc(data_scaling%diff_NPP_yDayMean)
             write(*,*) "Possible reason: inconsistent land-sea masks between different input files."
             stop
          end if
          if(max > 1E-4) then
             write(*,*) "getClimatologyDiffs(): ERROR: Found nonsense maximum value for diff_NPP_yDayMean ",max,&
                  " at (index,tile,month) ",maxloc(data_scaling%diff_NPP_yDayMean)
             write(*,*) "Possible reason: inconsistent land-sea masks between different input files."
             stop
          end if
       case('topSoilTemp_yDayMean')
          min=minval(data_scaling%diff_topSoilTemp_yDayMean)
          max=maxval(data_scaling%diff_topSoilTemp_yDayMean)
          if(min < -100.) then
             write(*,*) "getClimatologyDiffs(): ERROR: Found nonsense minimum value for diff_topSoilTemp_yDayMean ",min,&
                  " at (index,tile,month) ",minloc(data_scaling%diff_topSoilTemp_yDayMean)
             write(*,*) "Possible reason: inconsistent land-sea masks between different input files."
             stop
          end if
          if(max > 100.) then
             write(*,*) "getClimatologyDiffs(): ERROR: Found nonsense maximum value for diff_topSoilTemp_yDayMean ",max,&
                  " at (index,tile,month) ",maxloc(data_scaling%diff_topSoilTemp_yDayMean)
             write(*,*) "Possible reason: inconsistent land-sea masks between different input files."
             stop
          end if
       case('alpha_yDayMean')
          min=minval(data_scaling%diff_alpha_yDayMean)
          max=maxval(data_scaling%diff_alpha_yDayMean)
          if(min < -1.) then
             write(*,*) "getClimatologyDiffs(): ERROR: Found nonsense minimum value for diff_alpha_yDayMean ",min,&
                  " at (index,tile,month)",minloc(data_scaling%diff_alpha_yDayMean)
             write(*,*) "Possible reason: inconsistent land-sea masks between different input files."
             stop
          end if
          if(max > 1.) then
             write(*,*) "getClimatologyDiffs(): ERROR: Found nonsense maximum value for diff_alpha_yDayMean ",max,&
                  " at (index,tile,month) ",maxloc(data_scaling%diff_alpha_yDayMean)
             write(*,*) "Possible reason: inconsistent land-sea masks between different input files."
             stop
          end if
       end select

    end do

    ! close the file

    iret = nf_close(NCID)
    call check_err(iret,"Error(4) in getClimatologyDiffs() when closing "//trim(filename))

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
       write(*,*) "ERROR: Did not find variable topSoilTemp_yDayMean in "//trim(filename)//"."
       stop
    end if
    if(.not. alpha_found) then
       write(*,*) "ERROR: Did not find variable alpha_yDayMean in "//trim(filename)//"."
       stop
    end if

    ! deallocate local arrays

    deallocate(fdimnam,ndimlen)
    deallocate(fvarnam,nvartyp,nvardim,nvardimid,nvaratt)
  end subroutine getClimatologyDiffs

  ! --- getNDepositionDiffs()  -----------------------------------------------------------------------------------
  ! This routine reads in the difference between ndepositions of the first and last year, for each day one value in molesN/day.
  ! In the file the variables are assumed to have the names ndeposition.  
  ! The routine puts the data into the structure data_scaling.
  ! The routine automatiocally detects whether data are packed or in LonLat format.
  !
  subroutine getNdepositionDiffs(filename, grid, data_scaling_ndep)
    include 'netcdf.inc'
    character(len=*),intent(in)            :: filename
    type(grid_type),        intent(in)     :: grid
    type(data_scaling_type),intent(inout)  :: data_scaling_ndep 

    ! local variables
    real,allocatable :: hlp2D(:,:)

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
    integer  ::  i,ii


    ! open the input-file

    if(debug) write (*,*) 'getNdepositionDiffs(): read from: ',trim(filename)
    iret = nf_open(filename,nf_nowrite,ncid)
    call check_err(iret,"Error in getNdepositionDiffs() when opening file "//trim(filename))

    ! check what is in the input-file

    iret = nf_inq(ncid,ndim,nvar,nattr,nunlimdim)
    call check_err(iret,"Error in getNdepositionDiffs() when inquiring contents of "//trim(filename))
    if (debug) write (*,*) 'getNdepositionDiffs(): ndim,nvar,nattr,nunlimdim',ndim,nvar,nattr,nunlimdim

    ! get the dimension name and length

    allocate (ndimlen(ndim))
    allocate (fdimnam(ndim))
    do i=1,ndim
       iret = nf_inq_dim(ncid,i,fdimnam(i),ndimlen(i))
       call check_err(iret,"Error in getNdepositionDiffs() when inquiring dimension of "//trim(fdimnam(i))//" in "//trim(filename))
       if (debug) write (*,*) 'getNdepositionDiffs(): i,fdimnam,ndimlen',i,trim(fdimnam(i)),ndimlen(i)
    end do

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
       if (debug) write (*,*) 'getNdepositionDiffs(): i,fvarnam,nvartyp,nvardim,nvardimid,nvaratt',&
                                                         i,trim(fvarnam(i)),nvartyp(i),nvardim(i),nvardimid(i,1),nvaratt(i)
    end do

    ! get data
    
    if(debug) write (*,*) "getNdepositionDiffs(): Try to read the data .."
    do i=1,nvar
       select case (trim(fvarnam(i)))
       case('ndeposition')
          if(debug) write (*,*) "getNdepositionDiffs(): Case ndeposition"

             if(ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) then
                write(*,*) "getNdepositionDiffs(): Dimension error for variable "//&
                                                                                 trim(fvarnam(i))//" in file "//trim(filename)//":"
                write(*,*) "getNdepositionDiffs(): nlon=",grid%nlon," but found nndimlen(nvardimid(i,1))=",ndimlen(nvardimid(i,1))
                write(*,*) "getNdepositionDiffs(): nlat=",grid%nlat," but found nndimlen(nvardimid(i,2))=",ndimlen(nvardimid(i,2))
                stop
             end if

             allocate (hlp2D(grid%nlon,grid%nlat))
             iret = nf_get_var_real(ncid,i,hlp2D)
             call check_err(iret,"Error(4) in getNdepositionDiffs() when reading variable "//&
                                                                                       trim(fvarnam(i))//" from "//trim(filename))
             ! pack the array

             data_scaling_ndep%diff_Ndeposition(:) = pack(hlp2D(:,:),mask=grid%mask)
             deallocate(hlp2D)
       end select
    end do

    ! close the file

    iret = nf_close(NCID)
    call check_err(iret,"Error(4) in getNdepositionDiffs() when closing "//trim(filename))

    ! deallocate local arrays

    deallocate(fdimnam,ndimlen)
    deallocate(fvarnam,nvartyp,nvardim,nvardimid,nvaratt)
  end subroutine getNdepositionDiffs

end module mo_data_scaling
