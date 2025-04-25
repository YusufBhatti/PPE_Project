!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!!  mo_moz_jlong
!!
!! \brief
!!  MOZART (WACCM) calculation of photolysis frequencies at wavelength > 200 nm
!!
!! \author Douglas E. Kinnison (NCAR)
!! \author Stacy Walters (NCAR)
!!
!! \responsible_coder
!!  Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!  - DK: original version (2011)
!!
!! \limitations
!!  none
!!
!! \details
!!
!! \bibliographic_references
!!  - D.E. Kinnison et al., JGR 112, D20302, doi:10.1029/2006JD007879, 2007
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!!  Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!!  licencing agreement to be found at:
!!  https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!!  The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      module mo_moz_jlong

      USE mo_kind,    ONLY : dp
!>>SW see feature #415 and feature #469 in HAMMOZ readmine
      USE mo_radiation_parameters, ONLY: isolrad, lyr_perp, yr_perp
!<<SW      

      implicit none

      private
      public :: jlong, jlong_init, numj

!------------------------------------------------------------------------------
!     File names
!------------------------------------------------------------------------------
      character(len=64), parameter  :: tmpprs_file = 'moz_temp_prs_gt200nm.nc'
      character(len=64), parameter  :: rsf_file    = 'moz_rsf_gt200nm.nc'

!------------------------------------------------------------------------------
!     ... Define Wavelength Parameters
!------------------------------------------------------------------------------
      integer           :: nw         ! wavelengths >200nm
!++ost: from moz3 recent version, currently switched off !baustelle!
!     integer           :: nwe_max    ! last wave bin index inside woods grid
!--ost
      integer           :: nt         ! number of temperatures in xsection table
      integer           :: numj       ! number of photorates in xsqy, rsf
!++ost: changed z to p
!     integer           :: numz    		! number of altitudes in rsf
      integer           :: nump       ! number of altitudes in rsf
!--ost
      integer           :: numsza     ! number of zen angles in rsf
      integer           :: numalb     ! number of albedos in rsf
      integer           :: numcolo3   ! number of o3 columns in rsf
!++ost: from moz3: additional pressure dependance
      integer           :: np_xs                ! number of pressure levels in xs    table
      real(dp), allocatable :: prs(:), dprs(:)
!     real(dp), allocatable :: xsqy(:,:,:)
      real(dp), allocatable :: xsqy(:,:,:,:)
!--ost
      real(dp), allocatable :: wc(:)
      real(dp), allocatable :: wlintv(:)
!++ost: from moz3
      real(dp), allocatable :: we(:)
!>>SW see feature #415 in HAMMOZ readmine           
      real(dp), allocatable :: rsf_sclfac(:)
!<<SW      
!--ost
!++ost: changed z to p
!     real(dp), allocatable :: z(:), del_z(:)
      real(dp), allocatable :: p(:), del_p(:)
!--ost
      real(dp), allocatable :: sza(:), del_sza(:)
      real(dp), allocatable :: alb(:), del_alb(:)
      real(dp), allocatable :: o3rat(:), del_o3rat(:)
      real(dp), allocatable :: colo3(:)
      real(dp), allocatable :: rsf_tab(:,:,:,:,:)

      contains

!++ost: use pht_alias
      subroutine jlong_init( lng_indexer )
!--ost

      use mo_exception,   only : finish, message, message_text, em_debug
!++ost: from moz3 recent version, introduce later? !baustelle!
!     use mo_solar_parms, only : rebin
!     use mo_woods,          only : nbins, woods_etf, woods_we => we
!     use mo_neckel,      only : neckel_nw => nw, neckel_we => we, neckel_etf => etf
!--ost

      implicit none

!------------------------------------------------------------------------------
!    ... Dummy arguments
!------------------------------------------------------------------------------
      integer, intent(inout) :: &
                  lng_indexer(:)             ! model to long table xsect mapping

!------------------------------------------------------------------------------
!    ... Local variables
!------------------------------------------------------------------------------
!++ost: from moz3 recent version, currently switched off
!     real(dp), parameter:: woods_limit = 350.          ! limit of woods remapping (nm)
!     logical            :: found
!--ost

!------------------------------------------------------------------------------
!     ... Read Cross Section * QY NetCDF file...
!     ... Find temperature index for given altitude...
!     ... Derive cross*QY results, returns xsqy(nj,nz,nw)
!------------------------------------------------------------------------------
      WRITE(message_text,'(a)') 'before get_xsqy'
      CALL message('moz_jlong', message_text, level=em_debug)
!--ost
!++ost: use pht_alias
!     call get_xsqy( filepath, lpath, mspath )
      call get_xsqy( lng_indexer )
!--ost

!------------------------------------------------------------------------------
!     ... Read Radiative Source Function NetCDF file
!------------------------------------------------------------------------------
!++ost
      WRITE(message_text,'(a)')  'before get_rsf'
      CALL message('moz_jlong', message_text, level=em_debug)
!--ost
      call get_rsf
!>>SW see feature #415 in HAMMOZ readmine      
      ! Set contstant solar minimum/maximum values for rsf_sclfac
      ! if spectrally resolved data are not used
      if ( isolrad /= 1 ) then

        ! solar min etf-Flux (photons/cm2/sec/nm) based on data from J. Lean for Sep 1986.
              rsf_sclfac(:) = (/   &
                7.7132549E+11_dp, 8.7655665E+11_dp, 1.0766029E+12_dp, 1.2111112E+12_dp, 2.2414429E+12_dp, &
                3.4643580E+12_dp, 3.9668566E+12_dp, 4.5926874E+12_dp, 4.9122645E+12_dp, 6.5923402E+12_dp, &
                4.9106314E+12_dp, 5.6781848E+12_dp, 5.0606684E+12_dp, 6.0253012E+12_dp, 4.8162537E+12_dp, &
                7.4213485E+12_dp, 6.4861671E+12_dp, 6.4424955E+12_dp, 5.9519727E+12_dp, 7.8872419E+12_dp, &
                1.5197062E+13_dp, 1.2432981E+13_dp, 3.1368887E+13_dp, 3.3981565E+13_dp, 2.5945438E+13_dp, &
                2.8642267E+13_dp, 2.3444870E+13_dp, 3.7470519E+13_dp, 6.4123065E+13_dp, 7.8811249E+13_dp, &
                7.8341118E+13_dp, 6.7176061E+13_dp, 9.6227249E+13_dp, 9.2894472E+13_dp, 9.9229180E+13_dp, &
                1.1131345E+14_dp, 1.0881218E+14_dp, 1.2176523E+14_dp, 1.3870890E+14_dp, 1.6831589E+14_dp, &
                1.5568178E+14_dp, 1.6814230E+14_dp, 1.6388495E+14_dp, 1.5649723E+14_dp, 1.8620030E+14_dp, &
                1.6546949E+14_dp, 1.8127659E+14_dp, 2.2267347E+14_dp, 1.9913789E+14_dp, 2.2693899E+14_dp, &
                1.8488224E+14_dp, 2.0277141E+14_dp, 2.1040291E+14_dp, 2.5096542E+14_dp, 3.6576163E+14_dp, &
                3.4969785E+14_dp, 3.5263614E+14_dp, 3.7308143E+14_dp, 3.6369968E+14_dp, 3.6416914E+14_dp, &
                4.3881271E+14_dp, 4.8676690E+14_dp, 4.9447140E+14_dp, 5.2070221E+14_dp, 5.2121049E+14_dp, &
                5.0884612E+14_dp, 4.9201706E+14_dp /)

        ! solar max etf-Flux (photons/cm2/sec/nm) based on data from J. Lean for Nov. 1989
        !     rsf_sclfac(:) = (/   &
        !       8.3585427E+11_dp, 9.4828351E+11_dp, 1.1687353E+12_dp, 1.3161200E+12_dp, 2.3580752E+12_dp, &
        !       3.6046048E+12_dp, 4.1235891E+12_dp, 4.7754793E+12_dp, 5.0962812E+12_dp, 6.8141098E+12_dp, &
        !       5.1097224E+12_dp, 5.8712423E+12_dp, 5.2426666E+12_dp, 6.2329822E+12_dp, 5.0043594E+12_dp, &
        !       7.6356298E+12_dp, 6.7051087E+12_dp, 6.6691316E+12_dp, 6.1599719E+12_dp, 8.0930595E+12_dp, &
        !       1.5445745E+13_dp, 1.2727825E+13_dp, 3.1683854E+13_dp, 3.4230254E+13_dp, 2.6223451E+13_dp, &
        !       2.8915368E+13_dp, 2.4068764E+13_dp, 3.7876817E+13_dp, 6.4412539E+13_dp, 7.9006234E+13_dp, &
        !       7.8664641E+13_dp, 6.7578874E+13_dp, 9.6509274E+13_dp, 9.3238346E+13_dp, 9.9627789E+13_dp, &
        !       1.1163716E+14_dp, 1.0914159E+14_dp, 1.2200975E+14_dp, 1.3891006E+14_dp, 1.6844397E+14_dp, &
        !       1.5591028E+14_dp, 1.6835943E+14_dp, 1.6415034E+14_dp, 1.5684640E+14_dp, 1.8649188E+14_dp, &
        !       1.6609434E+14_dp, 1.8165720E+14_dp, 2.2283857E+14_dp, 1.9940575E+14_dp, 2.2735057E+14_dp, &
        !       1.8578948E+14_dp, 2.0404210E+14_dp, 2.1112501E+14_dp, 2.5143359E+14_dp, 3.6614207E+14_dp, &
        !       3.5013114E+14_dp, 3.5310403E+14_dp, 3.7357622E+14_dp, 3.6408079E+14_dp, 3.6466208E+14_dp, &
        !       4.3932747E+14_dp, 4.8722914E+14_dp, 4.9493280E+14_dp, 5.2116901E+14_dp, 5.2167864E+14_dp, &
        !       5.0927910E+14_dp, 4.9235337E+14_dp /)
      endif
!<<SW
      !>>SW see feature #469 in HAMMOZ readmine
      if (lyr_perp .and. yr_perp == 1850) then
          ! for piControl runs use 1850 values based on
          ! reference spectral and total irradiances derived from average over years 1834-1867 (solar cycles 8-10)
          ! values based on spectral_irradiance_Lean_1850_cntl_c100407.nc from CESM 1.0.2 input data repository
          rsf_sclfac(:) = (/   &
               7.75302720E+11_dp, 8.84866286E+11_dp, 1.07605846E+12_dp, 1.21898208E+12_dp, 2.24340845E+12_dp, &
               3.46684869E+12_dp, 3.97125418E+12_dp, 4.59277225E+12_dp, 4.91385402E+12_dp, 6.59374829E+12_dp, &
               4.91158628E+12_dp, 5.67704508E+12_dp, 5.06483344E+12_dp, 6.02533252E+12_dp, 4.81874690E+12_dp, &
               7.41379433E+12_dp, 6.48849492E+12_dp, 6.44285007E+12_dp, 5.95623227E+12_dp, 7.87563915E+12_dp, &
               1.51626300E+13_dp, 1.24202751E+13_dp, 3.12694306E+13_dp, 3.38776959E+13_dp, 2.58999454E+13_dp, &
               2.85668774E+13_dp, 2.33719506E+13_dp, 3.73790387E+13_dp, 6.38498529E+13_dp, 7.85557686E+13_dp, &
               7.75600097E+13_dp, 6.93265942E+13_dp, 9.12200304E+13_dp, 9.51040179E+13_dp, 9.71181964E+13_dp, &
               1.10930608E+14_dp, 1.08861697E+14_dp, 1.20657700E+14_dp, 1.39319244E+14_dp, 1.67877965E+14_dp, &
               1.57017276E+14_dp, 1.65278196E+14_dp, 1.62127557E+14_dp, 1.59017842E+14_dp, 1.85546972E+14_dp, &
               1.64962425E+14_dp, 1.80652325E+14_dp, 2.21887043E+14_dp, 1.98533497E+14_dp, 2.26143670E+14_dp, &
               1.84409320E+14_dp, 2.02135983E+14_dp, 2.09846955E+14_dp, 2.49955281E+14_dp, 3.64490387E+14_dp, &
               3.48502498E+14_dp, 3.51379108E+14_dp, 3.68077407E+14_dp, 3.63122523E+14_dp, 3.62100534E+14_dp, &
               4.36918338E+14_dp, 4.85471342E+14_dp, 4.92892119E+14_dp, 5.19172137E+14_dp, 5.19806977E+14_dp, &
               5.07378572E+14_dp, 4.90769939E+14_dp /)
      endif
!<<SW

!++ost: from moz3, currently switched off, introduce later?
!     found = .false.
!     do nwe_max = nw+1,2,-1
!        if( we(nwe_max) <= woods_limit ) then
!           found = .true.
!           exit 
!        end if
!     end do
!     if( .not. found ) then
!        CALL finish('jlong_init', 'Failed to place long wave spectra in woods')
!     else
!        write(*,*) 'jlong_init: nwe_max = ',nwe_max-1
!     end if
!     nwe_max = nwe_max - 1
!     call rebin( nbins, nwe_max, woods_we, we, woods_etf, etfphot )
!     call rebin( neckel_nw, nw-nwe_max, neckel_we, we(nwe_max+1), neckel_etf, etfphot(nwe_max+1) )
!--ost

      end subroutine jlong_init

      subroutine get_xsqy( lng_indexer )
!=============================================================================!
!   Subroutine GET_XSQY                                                       !
!                                                                             !
!=============================================================================!
!   PURPOSE:                                                                  !
!   Reads a NetCDF file that contains:                                        !
!     cross section * QY temperature dependence, >200nm                       !
!                                                                             !
!=============================================================================!
!   PARAMETERS:                                                               !
!     Input:                                                                  !
!      filepath.... NetCDF filepath that contains the "cross sections"        !
!=============================================================================!
!   EDIT HISTORY:                                                             !
!   Created by Doug Kinnison, 3/14/2002                                       !
!=============================================================================!

!     use NETCDF
      use MO_NETCDF
      use mo_exception,   only : finish, message_text
!++ost: not available anymore !baustelle!
!     use MO_FILE_UTILS, only : open_netcdf_file
!--ost
!++ost: use pht_alias
      use mo_moz_mods, only: phtcnt, pht_alias_lst, rxt_tag_lst
!--ost

      implicit none

!------------------------------------------------------------------------------
!    ... Dummy arguments
!------------------------------------------------------------------------------
      integer, intent(inout) :: &
                  lng_indexer(:)             ! model to long table xsect mapping

!------------------------------------------------------------------------------
!       ... Local variables
!------------------------------------------------------------------------------
      integer :: varid, dimid
      integer :: ncid
      integer :: iret
!++ost: use pht_alias
      integer :: i, m, ndx
      integer :: wrk_ndx(phtcnt)
!--ost
!++ost: from moz3: additional pressure dependance, but use of pht_alias
      integer :: k, n
      real(dp), allocatable :: xsqy_in(:,:,:)
!     real(dp), allocatable :: xsqy_in(:,:,:,:)
!--ost

!------------------------------------------------------------------------------
!       ... Open NetCDF File
!------------------------------------------------------------------------------
!++ost: not available anymore !baustelle!
!     ncid = open_netcdf_file( ncfile, lpath, mspath )
      CALL nf_check(nf__open(TRIM(tmpprs_file), nf_nowrite, chunksize, ncid),   &
                    fname=TRIM(tmpprs_file))
!--ost
!------------------------------------------------------------------------------
!       ... Get dimensions
!------------------------------------------------------------------------------
!+++ost: from moz3
!     iret = nf_inq_dimid( ncid, 'numj', dimid )
!     if( iret /= nf_noerr) then 
!        write(message_text,*) 'get_xsqy : failed to get numj id ; error = ',iret
!        CALL finish('get_xsqy (jlong)', message_text)
!     end if
!     iret = nf_inq_dimlen( ncid, dimid, numj )
!     if( iret /= nf_noerr) then 
!        write(message_text,*) 'get_xsqy : failed to get numj ; error = ',iret
!        CALL finish('get_xsqy (jlong)', message_text)
!     end if
!++ost: from moz3: additional pressure dependance
      iret = nf_inq_dimid( ncid, 'numprs', dimid )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get numprs id ; error = ',iret
         CALL finish('get_xsqy (jlong)', message_text)
      end if
      iret = nf_inq_dimlen( ncid, dimid, np_xs )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get np_xs ; error = ',iret
         CALL finish('get_xsqy (jlong)', message_text)
      end if
!--ost
      iret = nf_inq_dimid( ncid, 'numtemp', dimid )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get numtempj id ; error = ',iret
         CALL finish('get_xsqy (jlong)', message_text)
      end if
      iret = nf_inq_dimlen( ncid, dimid, nt )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get nt ; error = ',iret
         CALL finish('get_xsqy (jlong)', message_text)
      end if
      iret = nf_inq_dimid( ncid, 'numwl', dimid )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get numwl id ; error = ',iret
         CALL finish('get_xsqy (jlong)', message_text)
      end if
      iret = nf_inq_dimlen( ncid, dimid, nw )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get nw ; error = ',iret
         CALL finish('get_xsqy (jlong)', message_text)
      end if
!++ost: use pht_alias
!------------------------------------------------------------------------------
!       ... check for cross section in dataset
!------------------------------------------------------------------------------
      do m = 1,phtcnt
         if( pht_alias_lst(m,2) == ' ' ) then
            iret = nf_inq_varid( ncid, rxt_tag_lst(m), varid )
            if( iret == nf_noerr ) then 
               lng_indexer(m) = varid
            end if
         else if( pht_alias_lst(m,2) == 'userdefined' ) then
            lng_indexer(m) = -1
         else
            iret = nf_inq_varid( ncid, pht_alias_lst(m,2), varid )
            if( iret == nf_noerr ) then 
               lng_indexer(m) = varid
            else        
               write(message_text,'(a,a,a,a)') rxt_tag_lst(m)(:len_trim(rxt_tag_lst(m))),' alias ', &
                          pht_alias_lst(m,2)(:len_trim(pht_alias_lst(m,2))),' not in dataset'       
               CALL finish('get_xsqy (jlong)', message_text)
            end if
         end if
      end do
      numj = 0
      do m = 1,phtcnt
         if( lng_indexer(m) > 0 ) then
            if( any( lng_indexer(:m-1) == lng_indexer(m) ) ) then
               cycle
            end if
            numj = numj + 1
         end if 
      end do


!--ost
!------------------------------------------------------------------------------
!       ... Allocate arrays
!------------------------------------------------------------------------------
!++ost: from moz3: additional pressure dependance and use of pht alias
      allocate( xsqy(numj,nw,nt,np_xs),stat=iret )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to allocate xsqy ; error = ',iret
         CALL finish('get_xsqy (jlong)', message_text)
      end if
!     allocate( xsqy_in(numj,nt,nw),stat=iret )
!     allocate( xsqy_in(numj,nt,np_xs,nw),stat=iret )
      allocate( xsqy_in(nw,nt,np_xs),stat=iret )
!--ost
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to allocate xsqy_in ; error = ',iret
         CALL finish('get_xsqy (jlong)', message_text)
      end if
!++ost: from moz3: additional pressure dependance
      allocate( prs(np_xs),dprs(np_xs-1),stat=iret )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to allocate prs,dprs ; error = ',iret
         CALL finish('get_xsqy (jlong)', message_text)
      end if
!------------------------------------------------------------------------------
!       ... Read variables
!------------------------------------------------------------------------------
      iret = nf_inq_varid( ncid, 'pressure', varid )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get pressure id ; error = ',iret
         CALL finish('get_xsqy (jlong)', message_text)
      end if
      iret = nf_get_var_double( ncid, varid, prs )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get prs ; error = ',iret
         CALL finish('get_xsqy (jlong)', message_text)
      end if
      dprs(:np_xs-1) = 1._dp/(prs(1:np_xs-1) - prs(2:np_xs))
!--ost
!++ost: use pht_alias
!     iret = nf_inq_varid( ncid, 'cross_sections', varid )
!     if( iret /= nf_noerr) then 
!	 write(*,*) 'get_xsqy : failed to get xsqy id ; error = ',iret
!        CALL finish('get_xsqy (jlong)', message_text)
!     end if
!     iret = nf_get_var_double( ncid, varid, xsqy_in )
!     if( iret /= nf_noerr) then 
!	 write(*,*) 'get_xsqy : failed to get xsqy_in ; error = ',iret
!        CALL finish('get_xsqy (jlong)', message_text)
!     end if

      ndx = 0
      do m = 1,phtcnt
         if( lng_indexer(m) > 0 ) then
            if( any( lng_indexer(:m-1) == lng_indexer(m) ) ) then
               cycle
            end if
            iret = nf_get_var_double( ncid, lng_indexer(m), xsqy_in )
            if( iret /= nf_noerr) then 
               write(message_text,'(a,a,a,i0)') 'failed to read cross section for species ',rxt_tag_lst(m),' ; error = ',iret
               CALL finish('get_xsqy (jlong)', message_text)
            end if
            ndx = ndx + 1
            xsqy(ndx,:,:,:) = xsqy_in(:,:,:)
         end if
      end do
      deallocate( xsqy_in )
      if( ndx /= numj ) then
         write(message_text,'(a)') 'ndx count /= cross section count'
         CALL finish('get_xsqy (jlong)', message_text)
      end if
!------------------------------------------------------------------------------
!       ... setup final lng_indexer
!------------------------------------------------------------------------------
      ndx = 0
      wrk_ndx(:) = lng_indexer(:)
      do m = 1,phtcnt
         if( wrk_ndx(m) > 0 ) then
            ndx = ndx + 1
            i = wrk_ndx(m)
            where( wrk_ndx(:) == i )
               lng_indexer(:) = ndx
               wrk_ndx(:)     = -100000
            end where
         end if
      end do

!!------------------------------------------------------------------------------
!!       ... Transfer from input to working xsqy array
!!------------------------------------------------------------------------------
!!++ost: from moz3: additional pressure dependance
!!     allocate( xsqy(numj,nw,nt),stat=iret )
!     allocate( xsqy(numj,nw,nt,np_xs),stat=iret )
!!--ost
!      if( iret /= nf_noerr) then 
!	 write(*,*) 'get_xsqy : failed to allocate xsqy ; error = ',iret
!         CALL finish('get_xsqy (jlong)', message_text)
!     end if
!!++ost: from moz3: additional pressure dependance
!!     do n = 1,nt
!!        xsqy(:,:,n) = xsqy_in(:,n,:)
!!     end do
!     do k = 1,np_xs
!        do n = 1,nt
!           xsqy(:,:,n,k) = xsqy_in(:,n,k,:)
!        end do
!     end do
!!--ost
!
!     deallocate( xsqy_in )

!!      iret = nf_close( ncid )
      CALL nf_check(nf_close(ncid))

      end subroutine get_xsqy

      subroutine get_rsf
!=============================================================================!
!   Subroutine get_rsf                                                        !
!=============================================================================!
!   PURPOSE:                                                                  !
!   Reads a NetCDF file that contains:
!     Radiative Souce function                                                !
!=============================================================================!
!   PARAMETERS:                                                               !
!     Input:                                                                  !
!      filepath.... NetCDF file that contains the RSF                         !
!                                                                             !
!     Output:                                                                 !
!      rsf ........ Radiative Source Function (quanta cm-2 sec-1              !
!                                                                             !
!   EDIT HISTORY:                                                             !
!   Created by Doug Kinnison, 3/14/2002                                       !
!=============================================================================!

!     use netcdf
      use mo_netcdf
      use mo_exception,         only : finish, message_text
!++ost: not available anymore !baustelle!
!     use mo_file_utils, only : open_netcdf_file
!--ost

      implicit none

!------------------------------------------------------------------------------
!       ... Local variables
!------------------------------------------------------------------------------
      integer :: varid, dimid
      integer :: ncid
      integer :: i,j,k,l,wn
      integer :: iret
!++ost: additional variable in moz3 file, interim version
!     real(dp), allocatable :: etfphot(:)
!     real(dp), allocatable :: wrk(:)
!     real(dp), allocatable :: rsf_in(:,:,:,:,:)
      real(dp)  :: wrk
!--ost

!------------------------------------------------------------------------------
!       ... Open NetCDF File
!------------------------------------------------------------------------------
!!      ncid = open_netcdf_file( ncfile, lpath, mspath )
      CALL nf_check(nf__open(TRIM(rsf_file), nf_nowrite, chunksize, ncid),  &
                    fname=TRIM(rsf_file))

!------------------------------------------------------------------------------
!       ... Get dimensions
!------------------------------------------------------------------------------
      iret = nf_inq_dimid( ncid, 'numz', dimid )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get numz id ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
!++ost: changed z to p
!     iret = nf_inq_dimlen( ncid, dimid, numz )
      iret = nf_inq_dimlen( ncid, dimid, nump )
!--ost
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get nump ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
      iret = nf_inq_dimid( ncid, 'numsza', dimid )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get numsza id ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
      iret = nf_inq_dimlen( ncid, dimid, numsza )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get nunsza ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
      iret = nf_inq_dimid( ncid, 'numalb', dimid )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get numalb id ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
      iret = nf_inq_dimlen( ncid, dimid, numalb )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get numalb ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
      iret = nf_inq_dimid( ncid, 'numcolo3fact', dimid )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get numcolo3fact id ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
      iret = nf_inq_dimlen( ncid, dimid, numcolo3 )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get numcolo3 ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if

!------------------------------------------------------------------------------
!       ... Allocate arrays
!------------------------------------------------------------------------------
      allocate( wc(nw),stat=iret )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to allocate wc ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
!>>SW: changed etfphot --> rsf_sclfac
      allocate( wlintv(nw),we(nw+1),rsf_sclfac(nw),stat=iret ) 
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to allocate wlintv, we, rsf_sclfac ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
!<<SW
!++ost: changed z to p
!     allocate( z(numz),del_z(numz-1),stat=iret )
      allocate( p(nump),del_p(nump-1),stat=iret )
!--ost
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to allocate p,del_p ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
      allocate( sza(numsza),del_sza(numsza-1),stat=iret )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to allocate sza,del_sza ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
      allocate( alb(numalb),del_alb(numalb-1),stat=iret )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to allocate alb,del_alb ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
      allocate( o3rat(numcolo3),del_o3rat(numcolo3-1),stat=iret )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to allocate o3rat,del_o3rat ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
!++ost: changed z to p
!     allocate( colo3(numz),stat=iret )
      allocate( colo3(nump),stat=iret )
!--ost
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to allocate colo3; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
!++ost: changed z to p, rsf_in needed only in interim version
!     allocate( rsf_in(numalb,numcolo3,numsza,numz,nw),stat=iret )
!     allocate( rsf_in(numalb,numcolo3,numsza,nump,nw),stat=iret )
!     if( iret /= nf_noerr) then 
!	 write(*,*) 'get_rsf : failed to allocate rsf_in; error = ',iret
!         CALL finish('get_rsf (jlong)', message_text)
!     end if
!--ost
!++ost: changed z to p
!     allocate( rsf_tab(nw,numz,numsza,numcolo3,numalb),stat=iret )
      allocate( rsf_tab(nw,nump,numsza,numcolo3,numalb),stat=iret )
!--ost
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to allocate rsf_tab; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if

!------------------------------------------------------------------------------
!       ... Read variables
!------------------------------------------------------------------------------
      iret = nf_inq_varid( ncid, 'wc', varid )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get wc id ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
      iret = nf_get_var_double( ncid, varid, wc )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get wc ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
      iret = nf_inq_varid( ncid, 'wlintv', varid )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get wc id ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
      iret = nf_get_var_double( ncid, varid, wlintv )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get wlintv ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
!++ost: additional variables from moz3 file z->p, etfphot only in interim version
!     iret = nf_inq_varid( ncid, 'etfphot', varid )
!     if( iret /= nf_noerr) then 
!        write(*,*) 'get_rsf : failed to get etfphot id ; error = ',iret
!        CALL finish('get_rsf (jlong)', message_text)
!     end if
!     iret = nf_get_var_double( ncid, varid, etfphot )
!     if( iret /= nf_noerr) then 
!        write(*,*) 'get_rsf : failed to get etfphot ; error = ',iret
!        CALL finish('get_rsf (jlong)', message_text)
!     end if
      iret = nf_inq_varid( ncid, 'pm', varid )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get pm id ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
      iret = nf_get_var_double( ncid, varid, p )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get pm ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
      del_p(:nump-1) = 1._dp/(p(1:nump-1) - p(2:nump))
!     iret = nf_inq_varid( ncid, 'z', varid )
!     if( iret /= nf_noerr) then 
!	 write(*,*) 'get_rsf : failed to get z id ; error = ',iret
!        CALL finish('get_rsf (jlong)', message_text)
!     end if
!     iret = nf_get_var_double( ncid, varid, z )
!     if( iret /= nf_noerr) then 
!	 write(*,*) 'get_rsf : failed to get z ; error = ',iret
!        CALL finish('get_rsf (jlong)', message_text)
!     end if
!     del_z(:numz-1) = 1._dp/(z(2:numz) - z(1:numz-1))
!--ost
      iret = nf_inq_varid( ncid, 'sza', varid )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get sza id ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
      iret = nf_get_var_double( ncid, varid, sza )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get sza ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
      del_sza(:numsza-1) = 1._dp/(sza(2:numsza) - sza(1:numsza-1))
      iret = nf_inq_varid( ncid, 'alb', varid )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get alb id ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
      iret = nf_get_var_double( ncid, varid, alb )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get alb ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
      del_alb(:numalb-1) = 1._dp/(alb(2:numalb) - alb(1:numalb-1))
      iret = nf_inq_varid( ncid, 'colo3fact', varid )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get colo3fact id ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
      iret = nf_get_var_double( ncid, varid, o3rat )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get colo3fact ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
      del_o3rat(:numcolo3-1) = 1._dp/(o3rat(2:numcolo3) - o3rat(1:numcolo3-1))
      iret = nf_inq_varid( ncid, 'colo3', varid )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get colo3 id ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
      iret = nf_get_var_double( ncid, varid, colo3 )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get colo3 ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
      iret = nf_inq_varid( ncid, 'RSF', varid )
      if( iret /= nf_noerr) then
         write(message_text,'(a,i0)') 'failed to get rsf id ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if
!++ost: from moz3: rsf_in needed only in interim version
      iret = nf_get_var_double( ncid, varid, rsf_tab )
!     iret = nf_get_var_double( ncid, varid, rsf_in )
!--ost
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get rsf ; error = ',iret
         CALL finish('get_rsf (jlong)', message_text)
      end if

!------------------------------------------------------------------------------
!       ... form table rsf field
!------------------------------------------------------------------------------
!++ost: from moz3
      do wn = 1,nw
         wrk = wlintv(wn)
         do l = 1,numalb
            do k = 1,numcolo3
               do j = 1,numsza
                  do i = 1,nump
                     rsf_tab(wn,i,j,k,l) = wrk*rsf_tab(wn,i,j,k,l)
                  end do
               end do
            end do
         end do
      end do
!--ost
!++ost: interim version
!     wrk(:) = wlintv(:) * etfphot(:)
!     do i = 1,nump
!        do j = 1,numsza
!           do k = 1,numcolo3
!              do l = 1,numalb
!                 rsf_tab(:,i,j,k,l) = wrk(:)*rsf_in(l,k,j,i,:)
!              end do
!           end do
!        end do
!     end do
!     deallocate( rsf_in, etfphot, wrk )
!--ost
!++ost: original version
!     do i = 1,numz
!        do j = 1,numsza
!           do k = 1,numcolo3
!              do l = 1,numalb
!	          rsf_tab(:,i,j,k,l) = rsf_in(l,k,j,i,:)
!              end do
!           end do
!        end do
!     end do
!     deallocate( rsf_in )
!--ost

!!      iret = nf_close( ncid )
      CALL nf_check(nf_close(ncid))


!     colo3 for highest level is set to 0. in RSF input file. This
!     can lead to division by zero. Therefor it is set to half the value of the layer below, here.
!     (H. Schmidt, April 2005)

!++ost: changed z to p
!!     if (colo3(numz).eq.0._dp) colo3(numz) = .5_dp * colo3(numz-1)
      if (colo3(nump).eq.0._dp) colo3(nump) = .5_dp * colo3(nump-1)
!--ost

      end subroutine get_rsf

      subroutine jlong( nlev, sza_in, alb_in, p_in, t_in, &
                        colo3_in, j_long, rheat_o2, rheat_o3, & 
                        amu, cp, ip, lat, long )
!--ost

!==============================================================================
!   Subroutine JLONG                                                           
!==============================================================================
!   Purpose:                                                                   
!     To calculate the total J for selective species longward of 200nm.        
!==============================================================================
!                                                                              
!   Approach:
!     1) Reads the Cross Section*QY NetCDF file
!     2) Given a temperature profile, derives the appropriate XS*QY
!
!     3) Reads the Radiative Source function (RSF) NetCDF file
!        Units = quanta cm-2 sec-1
!
!     4) Indices are supplied to select a RSF that is consistent with
!        the reference atmosphere in TUV (for direct comparision of J's).
!        This approach will be replaced in the global model. Here colo3, zenith
!        angle, and altitude will be inputed and the correct entry in the table
!        will be derived.
!==============================================================================
!   EDIT HISTORY:
!   Created by Doug Kinnison, 3/14/2002
!==============================================================================

        use mo_mpi, only : p_pe
!hs heating only used in hammonia: use mo_rad_srbc, only: effhart
        use mo_physical_constants,  only: grav
!>>SW see feature #415 and #469 in HAMMOZ readmine
        USE mo_moz_dailyetf,      ONLY: rsf_sclfac_daily
        USE mo_interpo,       ONLY: ndw1, ndw2, wgtd1, wgtd2
        USE mo_time_control,  ONLY: get_time_step
        USE mo_exception,     ONLY: message, message_text
!<<SW
!hs 27day cycle only used in hammonia: use mo_param_switches,  only: l27day
!hs heating and 27day cycle only used in hammonia: use mo_radiation,       only: sin_27d, ituvheat

      implicit none

!------------------------------------------------------------------------------
!    ... Dummy arguments
!------------------------------------------------------------------------------
      integer, intent (in) :: nlev                             ! number vertical levels
      integer, intent (in) :: ip, lat
      integer, intent(in)  :: long
      real(dp),               intent(in)     :: sza_in         ! solar zenith angle (degrees)
      real(dp), dimension(:), intent(in)     :: alb_in         ! albedo
      real(dp), dimension(:),   intent(in)     :: p_in         ! midpoint pressure (hPa)
      real(dp), dimension(:),   intent(in)     :: t_in         ! Temperature profile (K)
      real(dp), dimension(:),   intent(in)     :: colo3_in     ! o3 column density (molecules/cm^3)
      real(dp), dimension(:,:), intent(out)    :: j_long       ! photo rates (1/s)

      real(dp), dimension(:),   intent (inout) :: rheat_o2     ! radiative O2 heating
      real(dp), dimension(:),   intent (inout) :: rheat_o3     ! radiative O3 heating
      real(dp), dimension(:),   intent (in)    :: amu 
      real(dp), dimension(:),   intent (in)    :: cp 

!----------------------------------------------------------------------
!        ... Local variables
!----------------------------------------------------------------------
      integer, parameter                     :: ituv2vis=117            ! use bands up to 680nm for heating
      real(dp), parameter                    :: ginv=1._dp/grav       ! 1/gravity acceleration (s2/m)
      REAL(dp), PARAMETER                    :: zll = 55._dp    ! lower limit for merging in km
      REAL(dp), PARAMETER                    :: zlu = 65._dp    ! upper limit for merging in km

      integer                                :: wn
      integer                                :: k, m, kk
      integer                                :: mi,mj,mk
!++ost
      integer                                :: iz2, km, pndx
!--ost

      integer                                :: is, iz
      integer                                :: iv, ial
      integer                                :: isp1, ivp1, ialp1
      integer                                :: izl
      integer                                :: t_index             ! Temperature index
      integer                                :: ratindl, ratindu
      integer                                :: pind, albind
      real(dp)                               :: psum_u
      real(dp)                               :: dels2
      real(dp)                               :: ptarget
      real(dp)                               :: wrk0, wrk1, wght1
      real(dp)                               :: v3ratl, v3ratu
      real(dp), dimension(3)                 :: dels
      real(dp), dimension(0:1,0:1,0:1)       :: wghtl, wghtu
      real(dp), dimension(nw)                :: psum_l              ! kk switch to automatic array
      real(dp), dimension(nw,nlev)           :: rsf                 ! Radiative source function
      real(dp), dimension(numj,nw)           :: xswk                ! working xsection array
      real(dp)                               :: delp

      real(dp)                               :: rsf_tmp
!     real(dp), dimension(nr_act_po)         :: jlong_tmp,jlong_tmp2,jlong_tmp3,jlong_tmp4

      real(dp)                               :: hc, be_o3, be_o2, ghm1, dmergei
!     real(dp), dimension(nlev)              :: heatcon


      real(dp) :: fact_27d(122)
!     factors for Jan-Jun 1990
      data fact_27d /0.0225526_dp, 0.0221079_dp, 0.0230692_dp, &
        0.022817_dp, 0.0146334_dp, 0.0113572_dp, 0.0105924_dp, &
        0.0113433_dp, 0.0109363_dp, 0.0106283_dp, 0.00940983_dp, &
        0.0105915_dp, 0.00993771_dp, 0.00966768_dp, 0.00989523_dp, &
        0.00969311_dp, 0.0107902_dp, 0.00807054_dp, 0.00929381_dp, &
        0.00974323_dp, 0.00963514_dp, 0.00688561_dp, 0.00456899_dp, &
        0.00698288_dp, 0.00282427_dp, 0.00203014_dp, 0.00226594_dp, &
        0.00265129_dp, 0.0091942_dp, 0.00298724_dp, 0.00190015_dp, &
        0.000531288_dp, 0.000951969_dp, 0.00156414_dp, 0.000808901_dp, &
        0.00057375_dp, 0.000723464_dp, 0.000941753_dp, 0.000941669_dp, &
        0.000888388_dp, 0.00122112_dp, 0.0011386_dp, 0.000806581_dp, &
        0.000723553_dp, 0.000544715_dp, 0.000835828_dp, 0.000732925_dp, &
        0.000448639_dp, 0.000291086_dp, 0.00010777_dp, 0.00032372_dp, &
        0.000268596_dp, 0.000365961_dp, 0.000393348_dp, 0.000611362_dp, &
        0.000820564_dp, 0.000266503_dp, 8.7437e-05_dp, 0.000598898_dp, &
        0.000556249_dp, 0.00195306_dp, 0.000777345_dp, 0.00116164_dp, &
        0.000180279_dp, 0.00028685_dp, 0.000274347_dp, 0.000263847_dp, &
        0._dp, 0._dp, 0._dp, 0._dp, &
        0._dp, 0._dp, 0._dp, 0._dp, &
        0._dp, 0._dp, 0._dp, 0._dp, &
        0._dp, 0._dp, 0._dp, 0._dp, &
        0._dp, 0._dp, 0._dp, 0._dp, &
        0._dp, 0._dp, 0._dp, 0._dp, &
        0._dp, 0._dp, 0._dp, 0._dp, &
        0._dp, 0._dp, 0._dp, 0._dp, &
        0._dp, 0._dp, 0._dp, 0._dp, &
        0._dp, 0._dp, 0._dp, 0._dp, &
        0._dp, 0._dp, 0._dp, 0._dp, &
        0._dp, 0._dp, 0._dp, 0._dp, &
        0._dp, 0._dp, 0._dp, 0._dp, &
        0._dp, 0._dp, 0._dp/

!kk
      integer                           :: ii

      hc    = 6.626196e-34_dp*2.997925e8_dp*1.e9_dp
      be_o2 = 8.20224e-19_dp                 ! [J], energy needed to break the O2 bond (5.12eV, Strobel et al., 1978)
      be_o3 = 1.6821e-19_dp                  ! [J], energy needed to break the O3 bond (1.05eV, Strobel et al., 1978)

      j_long(:,:) = 0.0_dp

!     heatcon(:) = 6.022e23_dp*1.e3_dp/(amu(:)*cp(:))

!----------------------------------------------------------------------
!        ... Interpolate for rsf
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!        ... Find the zenith angle index ( same for all levels )
!----------------------------------------------------------------------
!     if( sza_in < sza(1) .or. sza_in > sza(numsza) ) then
!       write(*,*) 'jlong: sza = ',sza_in,' is out of bounds at lon,lat,ip = ',long,lat,ip
!     end if
      do is = 1,numsza
         if( sza(is) > sza_in ) then
            exit
         end if
      end do
      is   = max( min( is,numsza ) - 1,1 )
      isp1 = is + 1
      dels(1)  = max( 0._dp,min( 1._dp,(sza_in - sza(is)) * del_sza(is) ) )
      wrk0     = 1._dp - dels(1)

      izl = 2
!>>SW see feature #415 and feature #469 in HAMMOZ readmine
      !----------------------------------------------------------------------
      !        ... Interpolate time varying scale factor
      !----------------------------------------------------------------------
      ! interpolate daily rsf_sclfac onto current time of the day, only use if
      ! spectrally resolved solar radiation is used in radiation as well
      ! if isolard /= 1 then rsf_sclfac(:) from jlong_init will be used
      if ( .not. lyr_perp .and. isolrad  == 1 ) then
         rsf_sclfac(:) = wgtd1*rsf_sclfac_daily(:,ndw1)+wgtd2*rsf_sclfac_daily(:,ndw2)
      endif
      ! Sebastian Wahl DEBUG
      ! write(message_text,'(a,3g15.7)') 'rsf_sclfac(1:3) = ',rsf_sclfac(1:3)
      ! CALL message('jlong', message_text)

      !IF ((lat.eq.1) .AND. ( MOD( get_time_step(),6 ) == 1 )) THEN
      !   WRITE(message_text,'(a20,f9.6)') 'RSF_sclfac @ 201nm: ',rsf_sclfac(1)
      !   CALL message('jlong',message_text)
      !ENDIF
!<<SW
Level_loop : &
      do k = nlev,1,-1
!----------------------------------------------------------------------
!        ... Find albedo indicies
!----------------------------------------------------------------------
!        if( alb_in(k) < alb(1) .or. alb_in(k) > alb(numalb) ) then
!          write(*,*) 'jlong: alb = ',alb_in(k),' is out of bounds at lon,lat,ip,k = ',long,lat,ip,k
!        end if
         do ial = 1,numalb
           if( alb(ial) > alb_in(k) ) then
              exit
           end if
        end do
        albind = max( min( ial,numalb ) - 1,1 )
!----------------------------------------------------------------------
!        ... Find pressure level indicies
!----------------------------------------------------------------------
!        if( p_in(k) > p(1) .or. p_in(k) < p(nump) ) then
!          write(*,*) 'jlong: p = ',p_in(k),' is out of bounds at lon,lat,ip,k = ',long,lat,ip,k
!        end if
         if( p_in(k) > p(1) ) then
            pind  = 2
            wght1 = 1._dp
         else if( p_in(k) <= p(nump) ) then
            pind  = nump
            wght1 = 0._dp
         else
            do iz = izl,nump
               if( p(iz) < p_in(k) ) then
	          izl = iz
	          exit
               end if
            end do
            pind  = max( min( iz,nump ),2 )
            wght1 = max( 0._dp,min( 1._dp,(p_in(k) - p(pind)) * del_p(pind-1) ) )
         end if
!----------------------------------------------------------------------
!        ... Find "o3 ratio" indicies; lower then upper
!----------------------------------------------------------------------
         if( colo3(pind) > 0._dp ) then
            v3ratl = colo3_in(k) / colo3(pind)
         else
            v3ratl = o3rat(numcolo3)
         end if
!        if( v3ratl < o3rat(1) .or. v3ratl > o3rat(numcolo3) ) then
!          write(*,*) 'jlong: v3ratl = ',v3ratl,' is out of bounds at lon,lat,ip,k = ',long,lat,ip,k
!        end if
         do iv = 1,numcolo3
            if( o3rat(iv) > v3ratl ) then
               exit
            end if
         end do
         ratindl = max( min( iv,numcolo3 ) - 1,1 )

         v3ratu  = colo3_in(k) / colo3(pind-1)
#ifdef DEBUG
!        if( v3ratu < o3rat(1) .or. v3ratu > o3rat(numcolo3) ) then
!          write(*,*) 'jlong: v3ratu = ',v3ratu,' is out of bounds at lon,lat,ip,k = ',long,lat,ip,k
!        end if
#endif
         do iv = 1,numcolo3
            if( o3rat(iv) > v3ratu ) then
               exit
            end if
         end do
         ratindu = max( min( iv,numcolo3 ) - 1,1 )

!----------------------------------------------------------------------
!        ... Compute the weigths
!----------------------------------------------------------------------
	 ial   = albind
	 ialp1 = ial + 1
	 iv    = ratindl

         dels(2)  = max( 0._dp,min( 1._dp,(v3ratl - o3rat(iv)) * del_o3rat(iv) ) )
	 dels2    = dels(2)
         dels(3)  = max( 0._dp,min( 1._dp,(alb_in(k) - alb(ial)) * del_alb(ial) ) )

	 wrk1         = (1._dp - dels(2))*(1._dp - dels(3))
	 wghtl(0,0,0) = wrk0*wrk1
	 wghtl(1,0,0) = dels(1)*wrk1
	 wrk1         = (1._dp - dels(2))*dels(3)
	 wghtl(0,0,1) = wrk0*wrk1
	 wghtl(1,0,1) = dels(1)*wrk1
	 wrk1         = dels(2)*(1._dp - dels(3))
	 wghtl(0,1,0) = wrk0*wrk1
	 wghtl(1,1,0) = dels(1)*wrk1
	 wrk1         = dels(2)*dels(3)
	 wghtl(0,1,1) = wrk0*wrk1
	 wghtl(1,1,1) = dels(1)*wrk1

	 iv  = ratindu
         dels(2)  = max( 0._dp,min( 1._dp,(v3ratu - o3rat(iv)) * del_o3rat(iv) ) )

	 wrk1         = (1._dp - dels(2))*(1._dp - dels(3))
	 wghtu(0,0,0) = wrk0*wrk1
	 wghtu(1,0,0) = dels(1)*wrk1
	 wrk1         = (1._dp - dels(2))*dels(3)
	 wghtu(0,0,1) = wrk0*wrk1
	 wghtu(1,0,1) = dels(1)*wrk1
	 wrk1         = dels(2)*(1._dp - dels(3))
	 wghtu(0,1,0) = wrk0*wrk1
	 wghtu(1,1,0) = dels(1)*wrk1
	 wrk1         = dels(2)*dels(3)
	 wghtu(0,1,1) = wrk0*wrk1
	 wghtu(1,1,1) = dels(1)*wrk1

	 iz   = pind
	 iv   = ratindl
	 ivp1 = iv + 1
         do wn = 1,nw
            psum_l(wn) = wghtl(0,0,0) * rsf_tab(wn,iz,is,iv,ial) &
                         + wghtl(0,0,1) * rsf_tab(wn,iz,is,iv,ialp1) &
                         + wghtl(0,1,0) * rsf_tab(wn,iz,is,ivp1,ial) &
                         + wghtl(0,1,1) * rsf_tab(wn,iz,is,ivp1,ialp1) &
                         + wghtl(1,0,0) * rsf_tab(wn,iz,isp1,iv,ial) &
                         + wghtl(1,0,1) * rsf_tab(wn,iz,isp1,iv,ialp1) &
                         + wghtl(1,1,0) * rsf_tab(wn,iz,isp1,ivp1,ial) &
                         + wghtl(1,1,1) * rsf_tab(wn,iz,isp1,ivp1,ialp1)
         end do

	 iz   = iz - 1
	 iv   = ratindu
	 ivp1 = iv + 1
         do wn = 1,nw
            psum_u = wghtu(0,0,0) * rsf_tab(wn,iz,is,iv,ial) &
                     + wghtu(0,0,1) * rsf_tab(wn,iz,is,iv,ialp1) &
                     + wghtu(0,1,0) * rsf_tab(wn,iz,is,ivp1,ial) &
                     + wghtu(0,1,1) * rsf_tab(wn,iz,is,ivp1,ialp1) &
                     + wghtu(1,0,0) * rsf_tab(wn,iz,isp1,iv,ial) &
                     + wghtu(1,0,1) * rsf_tab(wn,iz,isp1,iv,ialp1) &
                     + wghtu(1,1,0) * rsf_tab(wn,iz,isp1,ivp1,ial) &
                     + wghtu(1,1,1) * rsf_tab(wn,iz,isp1,ivp1,ialp1)
            rsf(wn,k) = (psum_l(wn) + wght1*(psum_u - psum_l(wn)))
         end do
!------------------------------------------------------------------------------
!     ... convert to photons/cm^2/s
!------------------------------------------------------------------------------
!>>SW see feature #415 in HAMMOZ readmine
         rsf(:,k) = rsf_sclfac(:) * rsf(:,k)
!<<SW         
      end do Level_loop

!------------------------------------------------------------------------------
!     ... calculate total Jlong for wavelengths >200nm
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!     LLNL LUT approach to finding temperature index...
!     Calculate the temperature index into the cross section
!     data which lists coss sections for temperatures from
!     150 to 350 degrees K.  Make sure the index is a value
!     between 1 and 201.
!------------------------------------------------------------------------------
level_loop_1 : &
      do k = 1,nlev
!----------------------------------------------------------------------
!    	... Get index into xsqy_in
!----------------------------------------------------------------------
        t_index = int( t_in(k) - 148.5_dp )
        t_index = min( 201,max( t_index,1) )
!----------------------------------------------------------------------
! 	... find pressure
!----------------------------------------------------------------------
        ptarget = p_in(k)
        if( ptarget >= prs(1) ) then
           do wn = 1,nw
              xswk(:numj,wn) = xsqy(:numj,wn,t_index,1)
           end do
        else if( ptarget <= prs(np_xs) ) then
           do wn = 1,nw
              xswk(:numj,wn) = xsqy(:numj,wn,t_index,np_xs)
           end do
        else
           do km = 2,np_xs
             if( ptarget >= prs(km) ) then
               pndx = km - 1
               delp = (prs(pndx) - ptarget)*dprs(pndx)
               exit
             end if
           end do
           do wn = 1,nw
              xswk(:numj,wn) = xsqy(:numj,wn,t_index,pndx) &
                               + delp*(xsqy(:numj,wn,t_index,pndx+1) - xsqy(:numj,wn,t_index,pndx))
           end do
        end if
#ifdef USE_ESSL
        call dgemm( 'N', 'N', numj, 1, nw, &
		    1., xswk, numj, rsf(1,k), nw, &
		    0., j_long(1,k), numj )
#else
        j_long(1:numj,k) = matmul( xswk(1:numj,1:nw),rsf(1:nw,k) )
#endif
      end do level_loop_1

  end subroutine jlong

end module mo_moz_jlong

