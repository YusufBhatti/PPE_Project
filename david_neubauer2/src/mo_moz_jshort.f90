!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!!  mo_moz_jshort
!!
!! \brief
!!  Calculation of photolysis frequencies in MOZART for wavelengths < 200 nm
!!
!! \author Douglas E. Kinnison (NCAR)
!! \author Stacy Walters (NCAR)
!!
!! \responsible_coder
!!  Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!  - DK: original version (2011)
!!  - Sebastian Wahl, GEOMAR (06/2015, 03/2016): adjustment to allow daily varying photolysis rates
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

module mo_moz_jshort

      USE mo_kind,   ONLY : dp
!>>SW see feature #415, and feature #469 in HAMMOZ readmine
      USE mo_radiation_parameters, ONLY: isolrad, lyr_perp, yr_perp
!<<SW      

      implicit none

      private
      public :: jshort_init, jshort, nj

!------------------------------------------------------------------------------
!     File names
!------------------------------------------------------------------------------
      character(len=64), parameter  :: effxs_file   = 'moz_effxs.nc'    ! xs_init
      character(len=64), parameter  :: xsshort_file = 'moz_xs_short.nc' ! get_crs

!------------------------------------------------------------------------------
!     ... Define Altitude and Wavelength Parameters
!------------------------------------------------------------------------------
      integer, parameter :: num_ms93tuv = 4       ! wavelength bins for MS, 93
      integer, parameter :: nsrb_tot    = 14      ! total bins <200nm for SRB
      integer, parameter :: nsrbtuv     = 17      ! total SRB bins in TUV
!++ost: from moz3
      integer, parameter :: nsrc_tot    = 19      ! total bins for SRC
!     integer, parameter :: nsrc_tot    = 20      ! total bins for SRC
      integer, parameter :: nw_ms93     = 4       ! wavelength bins for MS, 93
      real(dp), parameter    :: we_ms(nw_ms93+1) = (/ 181.6_dp, 183.1_dp, 184.6_dp, 190.2_dp, 192.5_dp/)
!     integer :: nw_ms93             ! Number of wavelength bins for ?
!--ost
      integer :: nw                  ! Number of wavelength bins <200nm
      integer :: nj                  ! Number of photorates
      real(dp), dimension(6,2) :: &
                 wtno50, csno50, wtno90, csno90, wtno100, csno100
!++ost: not needed
!     real(dp)    :: wave_num(nsrbtuv)
!--ost
      real(dp), allocatable :: wc(:)
!++ost: from moz3
      real(dp), allocatable :: we(:)
!--ost
      real(dp), allocatable :: wlintv(:)
      real(dp), allocatable :: etfphot(:)
      real(dp), allocatable :: etfphot_ms93(:)
!     real(dp), allocatable :: etfphot_sv(:)
!     real(dp), allocatable :: etfphot_ms93_sv(:)
      real(dp), allocatable :: xs_o2src(:)
!++ost: use pht_alias
      real(dp), allocatable :: xs_o3a(:)
      real(dp), allocatable :: xs_o3b(:)
!     real(dp), allocatable :: xs_species(:,:), xs_wl(:,:)
      real(dp), allocatable :: xs_species(:), xs_wl(:,:)
!     real(dp), allocatable :: xs_wl_o3(:)
!     real(dp)    :: ac(20,nsrbtuv)
!     real(dp)    :: bc(20,nsrbtuv)           ! chebyshev polynomial coeffs
      real(dp), allocatable :: ac(:,:)       ! chebyshev polynomial coefficients
      real(dp), allocatable :: bc(:,:)       ! chebyshev polynomial coefficients
!--ost

      contains

      subroutine jshort_init( sht_indexer )
!------------------------------------------------------------------------------
!     initialize the short wavelength photolysis module
!------------------------------------------------------------------------------
!--ost
!++ost: !!baustelle!! from moz3 recent version, currently switched off, introduce later?
!     use mo_solar_parms, only : rebin
!     use mo_woods,          only : nbins, woods_etf, woods_we => we
!--ost


      implicit none

!------------------------------------------------------------------------------
!    ... Dummy arguments
!------------------------------------------------------------------------------            
!++ost: use pht_alias
      integer, intent(inout) :: &
                  sht_indexer(:)             ! model to short table xsect mapping
!--ost

!------------------------------------------------------------------------------
!    ... Local variables
!------------------------------------------------------------------------------            
!++ost: debug
!     integer :: iret
!--ost

!------------------------------------------------------------------------------
!     ... Set the wavelength grid, exoatmospheric flux,
!         and cross sections (for <200nm) - contained in
!         a NetCDF file
!------------------------------------------------------------------------------
!++ost: use pht_alias
!     call get_crs( ncfile, lpath, mspath )
      call get_crs( sht_indexer )
!--ost
!++ost: from moz3
!++ost: !! baustelle !!
!     allocate (we(nw+1),stat=iret)
!--ost
      we(:nw)  = wc(:nw) - .5_dp*wlintv(:nw)
      we(nw+1) = wc(nw) + .5_dp*wlintv(nw)

!>>SW see feature #415 in HAMMOZ readmine
      ! Set contstant solar minimum/maximum values for etfphot, etfphot_ms93 and rsf_sclfac
      ! if spectrally resolved data are not used
      if ( isolrad /= 1) then
        ! solar min etf-Flux (photons/cm2/sec/nm) based on data from J. Lean for Sep 1986.
              etfphot(:) = (/ &
                3.4061225E+11_dp, 2.8171548E+09_dp, 2.6659666E+09_dp, 1.5017036E+09_dp, 2.0961071E+09_dp, &
                2.2475681E+09_dp, 1.3364187E+09_dp, 9.2329380E+09_dp, 4.1266365E+09_dp, 2.7929208E+09_dp, &
                3.7399716E+09_dp, 5.5155220E+09_dp, 9.9781540E+09_dp, 1.3864776E+10_dp, 2.0411505E+10_dp, &
                3.4521383E+10_dp, 4.8794321E+10_dp, 6.0502860E+10_dp, 6.7451658E+10_dp, 9.6101663E+10_dp, &
                1.3012759E+11_dp, 1.4225395E+11_dp, 1.8776085E+11_dp, 2.0622737E+11_dp, 1.9592492E+11_dp, &
                2.3245351E+11_dp, 2.9165466E+11_dp, 3.4027772E+11_dp, 3.9665172E+11_dp, 3.7337534E+11_dp, &
                5.3476999E+11_dp, 5.9912428E+11_dp, 6.3834312E+11_dp /)
              etfphot_ms93(:) = (/ 2.1345835E+11_dp, 2.0718194E+11_dp, 2.7371060E+11_dp, 3.9363851E+11_dp /)
        ! solar max etf-Flux (photons/cm2/sec/nm) based on data from J. Lean for Nov. 1989
        !     etfphot(:) = (/   &
        !       5.2330083E+11_dp, 3.6327593E+09_dp, 3.4755945E+09_dp, 1.9499300E+09_dp, 2.9057126E+09_dp, &
        !       3.3600837E+09_dp, 1.7534728E+09_dp, 1.2137259E+10_dp, 5.7640617E+09_dp, 3.4999340E+09_dp, &
        !       4.4747218E+09_dp, 6.3042381E+09_dp, 1.1855944E+10_dp, 1.5913691E+10_dp, 2.3294069E+10_dp, &
        !       3.9386784E+10_dp, 5.4136459E+10_dp, 6.7957748E+10_dp, 7.5022901E+10_dp, 1.0670839E+11_dp, &
        !       1.4547997E+11_dp, 1.5934888E+11_dp, 2.1573376E+11_dp, 2.3007265E+11_dp, 2.1550454E+11_dp, &
        !       2.5531680E+11_dp, 3.2102783E+11_dp, 3.7440016E+11_dp, 4.3498724E+11_dp, 4.0725417E+11_dp, &
        !       5.8695123E+11_dp, 6.5330971E+11_dp, 6.9177174E+11_dp /)
        !     etfphot_ms93(:) = (/ 2.4030988E+11_dp, 2.2864300E+11_dp, 3.0107890E+11_dp, 4.3160250E+11_dp /)
      endif
!<<SW
!>>SW see feature #469 in HAMMOZ readmine
      if (lyr_perp .and. yr_perp == 1850) then
          ! for piControl runs use 1850 values based on
          ! reference spectral and total irradiances derived from average over years 1834-1867 (solar cycles 8-10)
          ! values based on spectral_irradiance_Lean_1850_cntl_c100407.nc from CESM 1.0.2 input data repository
          etfphot(:) = (/ &
                3.56941437E+11_dp, 3.05022722E+09_dp, 2.15550135E+09_dp, 1.83959780E+09_dp, 2.00579806E+09_dp, &
                2.39609533E+09_dp, 1.52149954E+09_dp, 8.91598660E+09_dp, 4.29327198E+09_dp, 2.84611445E+09_dp, &
                3.79854121E+09_dp, 5.56089641E+09_dp, 1.00862048E+10_dp, 1.39998185E+10_dp, 2.05588278E+10_dp, &
                3.48370579E+10_dp, 4.89969349E+10_dp, 6.08803917E+10_dp, 6.82109038E+10_dp, 9.63390606E+10_dp, &
                1.31838959E+11_dp, 1.44369036E+11_dp, 1.90349554E+11_dp, 2.06143823E+11_dp, 1.99048076E+11_dp, &
                2.31979632E+11_dp, 2.93391263E+11_dp, 3.43669585E+11_dp, 3.91760530E+11_dp, 3.81096044E+11_dp, &
                5.34518732E+11_dp, 5.98253095E+11_dp, 6.39982649E+11_dp /)
          etfphot_ms93(:) = (/ 2.04567998E+11_dp, 2.04126334E+11_dp, 2.76532894E+11_dp, 3.89781838E+11_dp /)
      endif
!<<SW
!++ost: !!baustelle!! from moz3 recent version, currently switched off, introduce later?
!     call rebin( nbins, nw, woods_we, we, woods_etf, etfphot )
!     call rebin( nbins, nw_ms93, woods_we, we_ms, woods_etf, etfphot_ms93 )
!--ost

!------------------------------------------------------------------------------
!     ... Loads Chebyshev polynomial Coeff.
!------------------------------------------------------------------------------
      call xs_init
!------------------------------------------------------------------------------
!     ... Initialize no photolysis
!------------------------------------------------------------------------------
      call jno_init

      end subroutine jshort_init

      subroutine get_crs( sht_indexer )
!=============================================================================!
!   Subroutine get_crs                                                        !
!=============================================================================!
!   PURPOSE:                                                                  !
!   Reads a NetCDF file that contains:                                        !
!     Cross_sections*quantum yield data for species <200nm.                   !
!     Exoatmospheric flux on the wavelength grid of the crs                   !
!                                                                             !
!=============================================================================!
!   PARAMETERS:                                                               !
!     Input:                                                                  !
!      ncfile.... NetCDF file that contains the crs*QY for wavenum species    !
!                                                                             !
!     Output:                                                                 !
!      xs_species.. Cross Sections * QY data for each species                 !
!      etfphot..... Exoatmospheric flux in photons cm-2 s-1 nm-1              !
!      etfphot_ms93.Minshwanner and Siskind JNO SRB etf (on MS93 grid)        !
!      wc.......... wavelength center (nm)                                    !
!      numwl ...... Number of wavelengths < 200nm in NetCDF input file        !
!      wlintv ..... Wavelength inteval for grid, nm                           !
!=============================================================================!
!   EDIT HISTORY:                                                             !
!   Created by Doug Kinnison, 1/14/2002                                       !
!=============================================================================!

!     use NETCDF
      use MO_NETCDF
      use mo_exception,     only : finish, message_text
!!      use MO_FILE_UTILS, only : open_netcdf_file
!++ost: use pht_alias
      use mo_moz_mods, only: phtcnt, pht_alias_lst, rxt_tag_lst
!--ost

      implicit none

!------------------------------------------------------------------------------
!       ... Dummy arguments
!------------------------------------------------------------------------------
!++ost: use pht_alias
      integer, intent(inout) :: &
                  sht_indexer(:)             ! model to short table xsect mapping
!--ost

!------------------------------------------------------------------------------
!       ... Local variables
!------------------------------------------------------------------------------
      integer :: j

      integer :: ncid
      integer :: iret
      integer :: varid, dimid
!++ost: use pht_alias
      integer :: i, m, ndx
      integer :: wrk_ndx(phtcnt)
!--ost

!------------------------------------------------------------------------------
!       ... Open NetCDF File
!------------------------------------------------------------------------------
!!      ncid = open_netcdf_file( ncfile, lpath, mspath )
      CALL nf_check(nf__open(TRIM(xsshort_file), nf_nowrite, chunksize, ncid), &
                    fname=TRIM(xsshort_file))

!------------------------------------------------------------------------------
!       ... Get dimensions
!------------------------------------------------------------------------------
      iret = nf_inq_dimid( ncid, 'numwl', dimid )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get numwl id ; error = ',iret
         CALL finish('get_crs (jshort)', message_text)
      end if
      iret = nf_inq_dimlen( ncid, dimid, nw )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get nw ; error = ',iret
         CALL finish('get_crs (jshort)', message_text)
      end if
!++ost: from moz3
!     iret = nf_inq_dimid( ncid, 'numwl_ms93', dimid )
!     if( iret /= nf_noerr) then 
!	 write(*,*) 'failed to get numwl_ms93 id ; error = ',iret
!        CALL finish('get_crs (jshort)', message_text)
!     end if
!     iret = nf_inq_dimlen( ncid, dimid, nw_ms93 )
!     if( iret /= nf_noerr) then 
!	 write(*,*) 'failed to get nw_ms93 ; error = ',iret
!        CALL finish('get_crs (jshort)', message_text)
!     end if
!     iret = nf_inq_dimid( ncid, 'numj', dimid )
!     if( iret /= nf_noerr) then 
!	 write(*,*) 'failed to get numj id ; error = ',iret
!         CALL finish('get_crs (jshort)', message_text)
!     end if
!     iret = nf_inq_dimlen( ncid, dimid, nj )
!     if( iret /= nf_noerr) then 
!	 write(*,*) 'failed to get nj ; error = ',iret
!        CALL finish('get_crs (jshort)', message_text)
!     end if
!--ost

!++ost: use pht_alias
!------------------------------------------------------------------------------
!       ... check for cross section in dataset
!------------------------------------------------------------------------------
      do m = 1,phtcnt
         if( pht_alias_lst(m,1) == ' ' ) then
            iret = nf_inq_varid( ncid, rxt_tag_lst(m), varid )
            if( iret == nf_noerr ) then 
               sht_indexer(m) = varid
            end if
         else if( pht_alias_lst(m,1) == 'userdefined' ) then
            sht_indexer(m) = -1
         else
            iret = nf_inq_varid( ncid, pht_alias_lst(m,1), varid )
            if( iret == nf_noerr ) then 
               sht_indexer(m) = varid
            else    
               write(message_text,'(a,a,a,a)') rxt_tag_lst(m)(:len_trim(rxt_tag_lst(m))),' alias ', &
                          pht_alias_lst(m,1)(:len_trim(pht_alias_lst(m,1))),' not in dataset' 
               CALL finish('get_crs (jshort)', message_text)
            end if
         end if
      end do
      nj = 0
      do m = 1,phtcnt
         if( sht_indexer(m) > 0 ) then
            if( any( sht_indexer(:m-1) == sht_indexer(m) ) ) then
               cycle
            end if
            nj = nj + 1
         end if 
      end do
!--ost

!------------------------------------------------------------------------------
!       ... Allocate arrays
!------------------------------------------------------------------------------
      allocate( wc(nw),stat=iret )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to allocate wc ; error = ',iret
         CALL finish('get_crs (jshort)', message_text)
      end if
!++ost: from moz3
      allocate( we(nw+1),stat=iret )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to allocate we ; error = ',iret
         CALL finish('get_crs (jshort)', message_text)
      end if
!--ost
      allocate( wlintv(nw),stat=iret )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to allocate wlintv ; error = ',iret
         CALL finish('get_crs (jshort)', message_text)
      end if
      allocate( etfphot(nw),stat=iret )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to allocate etfphot ; error = ',iret
         CALL finish('get_crs (jshort)', message_text)
      end if
      allocate( etfphot_ms93(nw_ms93),stat=iret )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to allocate etfphot_ms93 ; error = ',iret
         CALL finish('get_crs (jshort)', message_text)
      end if
!++ost: this variable is not contained in the most recent file version
!     allocate( etfphot_sv(nw),stat=iret )
!     if( iret /= nf_noerr) then
!        write(*,*) 'failed to allocate etfphot_sv ; error = ',iret
!        CALL finish('get_crs (jshort)', message_text)
!     end if
!     allocate( etfphot_ms93_sv(nw_ms93),stat=iret )
!     if( iret /= nf_noerr) then
!        write(*,*) 'failed to allocate etfphot_ms93_sv ; error = ',iret
!        CALL finish('get_crs (jshort)', message_text)
!     end if
!--ost
      allocate( xs_o2src(nw),stat=iret )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to allocate xs_o2src ; error = ',iret
         CALL finish('get_crs (jshort)', message_text)
      end if
!++ost: use pht_alias
       allocate( xs_o3a(nw),xs_o3b(nw),stat=iret )
       if( iret /= 0 ) then
          write(message_text,'(a,i0)') 'failed to allocate xs_o3a,xs_o3b ; error = ',iret
          CALL finish('get_crs (jshort)', message_text)
       end if
!     allocate( xs_species(nj,nw),xs_wl(nw,nj),stat=iret )
      allocate( xs_species(nw),xs_wl(nw,nj),stat=iret )
!--ost
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to allocate xs_species ; error = ',iret
         CALL finish('get_crs (jshort)', message_text)
      end if
!++ost: from moz3, not contained
!     allocate( xs_wl_o3(nw),stat=iret )
!     if( iret /= nf_noerr) then
!        write(*,*) 'failed to allocate xs_wl_o3 ; error = ',iret
!        CALL finish('get_crs (jshort)', message_text)
!     end if
!--ost

!------------------------------------------------------------------------------
!       ... Read arrays
!------------------------------------------------------------------------------
      iret = nf_inq_varid( ncid, 'wc', varid )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get wc id ; error = ',iret
         CALL finish('get_crs (jshort)', message_text)
      end if
      iret = nf_get_var_double( ncid, varid, wc )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get wc ; error = ',iret
         CALL finish('get_crs (jshort)', message_text)
      end if
      iret = nf_inq_varid( ncid, 'wlintv', varid )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get wlintv id ; error = ',iret
         CALL finish('get_crs (jshort)', message_text)
      end if
      iret = nf_get_var_double( ncid, varid, wlintv )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get wlintv ; error = ',iret
         CALL finish('get_crs (jshort)', message_text)
      end if
!++ost: from moz3, not contained
!     iret = nf_inq_varid( ncid, 'etfphot', varid )
!     if( iret /= nf_noerr) then 
!	 write(*,*) 'failed to get etfphot id ; error = ',iret
!        CALL finish('get_crs (jshort)', message_text)
!     end if
!     iret = nf_get_var_double( ncid, varid, etfphot )
!     if( iret /= nf_noerr) then 
!	 write(*,*) 'failed to get etfphot ; error = ',iret
!        CALL finish('get_crs (jshort)', message_text)
!     end if
!     iret = nf_inq_varid( ncid, 'etfphot_ms93', varid )
!     if( iret /= nf_noerr) then 
!	 write(*,*) 'failed to get etfphot_ms93 id ; error = ',iret
!        CALL finish('get_crs (jshort)', message_text)
!     end if
!     iret = nf_get_var_double( ncid, varid, etfphot_ms93 )
!     if( iret /= nf_noerr) then 
!	 write(*,*) 'failed to get etfphot_ms93 ; error = ',iret
!        CALL finish('get_crs (jshort)', message_text)
!     end if
!--ost
      iret = nf_inq_varid( ncid, 'xs_o2src', varid )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get xs_o2src id ; error = ',iret
         CALL finish('get_crs (jshort)', message_text)
      end if
      iret = nf_get_var_double( ncid, varid, xs_o2src )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get xs_o2src ; error = ',iret
         CALL finish('get_crs (jshort)', message_text)
      end if
!++ost: use pht_alias
!     iret = nf_inq_varid( ncid, 'xs_species', varid )
!     if( iret /= nf_noerr) then 
!	 write(*,*) 'failed to get xs_species id ; error = ',iret
!        CALL finish('get_crs (jshort)', message_text)
!     end if
!     iret = nf_get_var_double( ncid, varid, xs_species )
!     if( iret /= nf_noerr) then 
!	 write(*,*) 'failed to get xs_species ; error = ',iret
!        CALL finish('get_crs (jshort)', message_text)
!     end if
!     do j = 1,nw
!        xs_wl(j,:) = xs_species(:,j)*wlintv(j)
!     end do

      ndx = 0
      do m = 1,phtcnt
         if( sht_indexer(m) > 0 ) then
            if( any( sht_indexer(:m-1) == sht_indexer(m) ) ) then
               cycle
            end if
            iret = nf_get_var_double( ncid, sht_indexer(m), xs_species )
            if( iret /= nf_noerr) then 
               write(message_text,'(a,a,a,i0)') 'failed to read cross section for species ',rxt_tag_lst(m),' ; error = ',iret
               CALL finish('get_crs (jshort)', message_text)
            end if
            ndx = ndx + 1
            xs_wl(:,ndx) = xs_species(:)*wlintv(:)
         end if
      end do
      deallocate( xs_species )
      if( ndx /= nj ) then
         write(message_text,'(a)') 'ndx count /= cross section count'
         CALL finish('get_crs (jshort)', message_text)
      end if
!------------------------------------------------------------------------------
!       ... get jo3 cross sections
!------------------------------------------------------------------------------
      iret = nf_inq_varid( ncid, 'jo3_a', varid )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get jo3_a id ; error = ',iret
         CALL finish('get_crs (jshort)', message_text)
      end if
      iret = nf_get_var_double( ncid, varid, xs_o3a )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get xs_o3a ; error = ',iret
         CALL finish('get_crs (jshort)', message_text)
      end if
      iret = nf_inq_varid( ncid, 'jo3_b', varid )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get jo3_b id ; error = ',iret
         CALL finish('get_crs (jshort)', message_text)
      end if
      iret = nf_get_var_double( ncid, varid, xs_o3b )
      if( iret /= nf_noerr) then 
         write(message_text,'(a,i0)') 'failed to get xs_o3b ; error = ',iret
         CALL finish('get_crs (jshort)', message_text)
      end if
!------------------------------------------------------------------------------
!       ... setup final sht_indexer
!------------------------------------------------------------------------------
      ndx = 0
      wrk_ndx(:) = sht_indexer(:)
      do m = 1,phtcnt
         if( wrk_ndx(m) > 0 ) then
            ndx = ndx + 1
            i = wrk_ndx(m)
            where( wrk_ndx(:) == i )
               sht_indexer(:) = ndx
               wrk_ndx(:)     = -100000
            end where
         end if
      end do

      CALL nf_check(nf_close(ncid))

      end subroutine get_crs

      subroutine xs_init
!-------------------------------------------------------------
!    	... Loads XS_COEFFS containing the Chebyshev
!           polynomial coeffs necessary to calculate O2 effective
!           cross-sections
!-------------------------------------------------------------

      use mo_mpi
      use mo_netcdf
      use mo_exception,        only : finish, message_text

      implicit none

!-------------------------------------------------------------
!     ... Local variables
!-------------------------------------------------------------
!++ost: now read netcdf file
      integer :: ncid
      integer :: iret
      integer :: ncof, nsrb
      integer :: varid, dimid
!--ost
      integer :: unit   ! file unit number
      integer :: istat  ! i/o status
      integer :: i
      logical :: cosb

!++ost: !!baustelle!! how to handle parallel_io?
!     if( masternode ) then
!     if( p_parallel_io ) then
!--ost
!++ost: now read netcdf file
!------------------------------------------------------------------------------
!       ... Open NetCDF File
!------------------------------------------------------------------------------
!!         ncid = open_netcdf_file( filename, lpath, mspath )
         CALL nf_check(nf__open(TRIM(effxs_file), nf_nowrite, chunksize, ncid),   &
                       fname=TRIM(effxs_file))

!------------------------------------------------------------------------------
!       ... get dimensions
!------------------------------------------------------------------------------
         iret = nf_inq_dimid( ncid, 'ncof', dimid )
         if( iret /= nf_noerr) then 
            write(message_text,'(a,i0)') 'xs_init : failed to get ncof id ; error = ',iret
            CALL finish('xs_init (jshort)', message_text)
         end if
         iret = nf_inq_dimlen( ncid, dimid, ncof )
         if( iret /= nf_noerr) then 
            write(message_text,'(a,i0)') 'failed to get ncof ; error = ',iret
            CALL finish('xs_init (jshort)', message_text)
         end if
         iret = nf_inq_dimid( ncid, 'nsrbtuv', dimid )
         if( iret /= nf_noerr) then 
            write(message_text,'(a,i0)') 'xs_init : failed to get nsrbtuv id ; error = ',iret
            CALL finish('xs_init (jshort)', message_text)
         end if
         iret = nf_inq_dimlen( ncid, dimid, nsrb )
         if( iret /= nf_noerr) then 
            write(message_text,'(a,i0)') 'xs_init : failed to get nsrb ; error = ',iret
            CALL finish('xs_init (jshort)', message_text)
         end if
         if( nsrb /= nsrbtuv ) then
            write(message_text,'(a,i0,i0)') 'xs_init : file nsrb != model nsrb ; ',nsrb,nsrbtuv
            CALL finish('xs_init (jshort)', message_text)
         end if
!----------------------------------------------------------------------
!       ... allocate arrays
!----------------------------------------------------------------------
         allocate( ac(ncof,nsrbtuv), bc(ncof,nsrbtuv), stat=iret )
         if( iret /= nf_noerr ) then
            write(message_text,'(a,i0)') 'xs_init : failed to allocate ac,bc ; error = ',iret
            CALL finish('xs_init (jshort)', message_text)
         end if
!----------------------------------------------------------------------
!       ... read coefficient arrays
!----------------------------------------------------------------------
         iret = nf_inq_varid( ncid, 'ac', varid )
         if( iret /= nf_noerr) then 
            write(message_text,'(a,i0)') 'xs_init : failed to get ac id ; error = ',iret
            CALL finish('xs_init (jshort)', message_text)
         end if
         iret = nf_get_var_double( ncid, varid, ac )
         if( iret /= nf_noerr) then 
            write(message_text,'(a,i0)') 'xs_init : failed to get ac ; error = ',iret
            CALL finish('xs_init (jshort)', message_text)
         end if
         iret = nf_inq_varid( ncid, 'bc', varid )
         if( iret /= nf_noerr) then 
            write(message_text,'(a,i0)') 'xs_init : failed to get bc id ; error = ',iret
            CALL finish('xs_init (jshort)', message_text)
         end if
         iret = nf_get_var_double( ncid, varid, bc )
         if( iret /= nf_noerr) then 
            write(message_text,'(a,i0)') 'xs_init : failed to get bc ; error = ',iret
            CALL finish('xs_init (jshort)', message_text)
         end if
         CALL nf_check(nf_close(ncid))

      end subroutine xs_init

      subroutine jno_init
!-------------------------------------------------------------
!    	... Initialization for no photolysis
!-------------------------------------------------------------

      implicit none

!-------------------------------------------------------------
!     ... Local variables
!-------------------------------------------------------------
      real(dp), dimension(24) :: a, b, c

!------------------------------------------------------------------------------
!   	... 6 sub-intervals for O2 5-0 at 265K,
!	    2 sub-sub-intervals for NO 0-0 at 250K
!------------------------------------------------------------------------------
      a(:) = (/ 0._dp, 0._dp, 0._dp, 0.0_dp, &
                   5.12e-02_dp, 5.68e-03_dp, 1.32e-18_dp, 4.41e-17_dp, &
                   1.36e-01_dp, 1.52e-02_dp, 6.35e-19_dp, 4.45e-17_dp, &
                   1.65e-01_dp, 1.83e-02_dp, 7.09e-19_dp, 4.50e-17_dp, &
                   1.41e-01_dp, 1.57e-02_dp, 2.18e-19_dp, 2.94e-17_dp, &
                   4.50e-02_dp, 5.00e-03_dp, 4.67e-19_dp, 4.35e-17_dp /)

!------------------------------------------------------------------------------
!   	... sub-intervals for o2 9-0 band,
!	    2 sub-sub-intervals for no 1-0 at 250 k
!------------------------------------------------------------------------------
      b(:) = (/ 0._dp, 0._dp, 0._dp, 0._dp, &
                 0.00e+00_dp, 0.00e+00_dp, 0.00e+00_dp, 0.00e+00_dp, &
                 1.93e-03_dp, 2.14e-04_dp, 3.05e-21_dp, 3.20e-21_dp, &
                 9.73e-02_dp, 1.08e-02_dp, 5.76e-19_dp, 5.71e-17_dp, &
                 9.75e-02_dp, 1.08e-02_dp, 2.29e-18_dp, 9.09e-17_dp, &
                 3.48e-02_dp, 3.86e-03_dp, 2.21e-18_dp, 6.00e-17_dp /)

!------------------------------------------------------------------------------
! 	... sub-intervals for o2 10-0 band,
!	    2 sub-sub-intervals for no 1-0 at 250 k
!------------------------------------------------------------------------------
      c(:) = (/  4.50e-02_dp, 5.00e-03_dp, 1.80e-18_dp, 1.40e-16_dp, &
                 1.80e-01_dp, 2.00e-02_dp, 1.50e-18_dp, 1.52e-16_dp, &
                 2.25e-01_dp, 2.50e-02_dp, 5.01e-19_dp, 7.00e-17_dp, &
                 2.25e-01_dp, 2.50e-02_dp, 7.20e-20_dp, 2.83e-17_dp, &
                 1.80e-01_dp, 2.00e-02_dp, 6.72e-20_dp, 2.73e-17_dp, &
                 4.50e-02_dp, 5.00e-03_dp, 1.49e-21_dp, 6.57e-18_dp /)

      wtno50 (1:6,1) = a(1:24:4)
      wtno50 (1:6,2) = a(2:24:4)
      csno50 (1:6,1) = a(3:24:4)
      csno50 (1:6,2) = a(4:24:4)
      wtno90 (1:6,1) = b(1:24:4)
      wtno90 (1:6,2) = b(2:24:4)
      csno90 (1:6,1) = b(3:24:4)
      csno90 (1:6,2) = b(4:24:4)
      wtno100(1:6,1) = c(1:24:4)
      wtno100(1:6,2) = c(2:24:4)
      csno100(1:6,1) = c(3:24:4)
      csno100(1:6,2) = c(4:24:4)

      end subroutine jno_init

      subroutine jshort( nlev, zen, n2cc, o2cc, o3cc, &
                         nocc, tlev, zkm, jo2_sht, jno_sht, &
                         jsht, rheat_o2, rheat_o3, amu, cp, & 
                         lat, ip, long, do_diag )
!==============================================================================!
!   Subroutine Jshort                                                          !
!                                                                              !
!==============================================================================!
!   Purpose:                                                                   !
!     To calculate the total J for JO2, JNO, and selective species below 200nm.!
!                                                                              !
!==============================================================================!
!   This routine uses JO2 parameterizations based on:                          !
!        Lyman alpha... Chabrillat and Kockarts, GRL, 25, 2659, 1998           !
!        SRC .......... Brasseur and Solomon, 1986 (from TUV)                  !
!        SRB .......... Koppers and Murtagh, Ann. Geophys., 14, 68-79, 1996    !
!                        (supplied by Dan Marsh, NCAR ACD                      !
!   and JNO:                                                                   !
!        SRB .......... Minschwanner and Siskind, JGR< 98, 20401, 1993.        !
!                                                                              !
!==============================================================================!
!   Input:                                                                     !
!	n2cc....... N2 concentration, molecule cm-3			       !
!     o2cc....... O2 concentration, molecule cm-3			       !
!	o3cc....... O3 concentration, molecule cm-3			       !
!	nocc....... NO concentration, molecule cm-3			       !
!	n2cc....... N2 concentration, molecule cm-3			       !
!     zen........ zenith angle, units = degrees                                !
!     tlev....... Temperature Profile (K)                                      !
!     zkm ....... Altitude, km                                                 !
!                                                                              !
!   Output:                                                                    !
!     jo2_sht ... O2 photolytic rate constant, sec-1, <200nm                   !
!     jno_sht ... NO photolytic rate constant, sec-1, SRB                      !
!     jsht. Photolytic rate constant for other species below 200nm             !
!                                                                              !
!==============================================================================!
!                                                                              !
!   Approach:                                                                  !
!                                                                              !
!    1) Calls subroutine sphers (taken from TUV)                               !
!         -> derives dsdh and nid used in slant column routines                !
!		-> zenith angle dependent                                      !
!                                                                              !
!    2) Calls subroutine's slant_col (taken from TUV)            !
!		-> derives the slant column for each species                   !
!                                                                              !
!    3) Calls get_crs                                                          !
!		-> read a NetCDF file                                          !
!		-> returns cross sections*quantum yields for all species that  !
!		   have absorption below 200nm.                                !
!                                                                              !
!    4) Derives transmission and photolysis rates for selective species        !
!                                                                              !
!==============================================================================!
!   EDIT HISTORY:                                                              !
!   Created by Doug Kinnison, 3/14/2002                                        !
!==============================================================================!

!use mo_mpi, only : thisnode
use mo_mpi,             only: p_pe
use mo_exception,       only: finish, message_text
!>>SW see feature #415 in HAMMOZ readmine     
USE mo_moz_dailyetf,        ONLY: etfphot_daily, etfphot_ms93_daily
USE mo_interpo,         ONLY: ndw1, ndw2, wgtd1, wgtd2
!<<SW
!hs only used in hammonia: use mo_rad_srbc,        only: effsrc, effhart
!hs only used in hammonia: use mo_param_switches,  only: l27day
!hs only used in hammonia: use mo_radiation,       only: sin_27d

  implicit none

  integer, parameter :: branch    = 2       ! two photolysis pathways for JO2
!++ost: from moz3
        integer, parameter :: lya_ndx     = 1       ! lyman-alpha wavelength index
!--ost
	real(dp), parameter    :: km2cm    = 1.e5_dp

!------------------------------------------------------------------------------
!     ... Dummy arguments
!------------------------------------------------------------------------------
  integer, intent(in) :: nlev                     ! model vertical levels
  integer, intent(in) :: lat, ip
  integer, intent(in) :: long

  real(dp), intent (in)    :: zen           ! Zenith angle (degrees)
  real(dp), dimension(:),   intent (in)    :: n2cc      ! Molecular Nitrogen conc.
  real(dp), dimension(:),   intent (in)    :: o2cc      ! Molecular Oxygen conc.
  real(dp), dimension(:),   intent (in)    :: o3cc      ! Ozone concentration
  real(dp), dimension(:),   intent (in)    :: nocc      ! Nitric Oxide conc.
  real(dp), dimension(:),   intent (in)    :: tlev      ! Temperature profile
  real(dp), dimension(:),   intent (in)    :: zkm       ! Altitude, km
  real(dp), dimension(:),   intent (in)    :: amu
  real(dp), dimension(:),   intent (in)    :: cp

  real(dp), dimension(:,:), intent (out)   :: jo2_sht       ! JO2, sec-1, <200nm
  real(dp), dimension(:,:), intent (out)   :: jsht          ! Additional J's
  real(dp), dimension(:),   intent (out)   :: jno_sht       ! JNO, sec-1, SRB
 
  real(dp), dimension(:),   intent (inout) :: rheat_o2      ! radiative O2 heating
  real(dp), dimension(:),   intent (inout) :: rheat_o3      ! radiative O3 heating

  logical, intent(in)                    :: do_diag

!------------------------------------------------------------------------------
!     ... Local variables
!------------------------------------------------------------------------------
        integer             :: mi, mj,mk
        integer             :: ii, k, n

  integer :: iz, izz              ! Altitude index
  integer :: iw                   ! Wavelength index
  integer :: astat
  integer, dimension(0:nlev)   :: nid   ! Number of layers crossed by the direct
                                        ! beam when travelling from the top of the
                ! atmosphere to layer i; NID(i), i = 0..NZ-1
  real(dp), dimension(0:nlev,nlev) :: dsdh  ! Slant path of direct beam through each
                                            ! layer crossed  when travelling from the top of
                ! the atmosphere to layer i; DSDH(i,j), i = 0.
                ! NZ-1, j = 1..NZ-1   (see sphers.f)
        real(dp), allocatable :: fnorm(:,:)      ! Normalized ETF
	real(dp) :: jo2_lya(nlev)      ! Total photolytic rate constant for Ly alpha
	real(dp) :: jo2_srb(nlev)      ! Total JO2  for SRB
	real(dp) :: jo2_src(nlev)      ! Total JO2 for SRC
	real(dp) :: delz(nlev)         ! layer thickness (cm)
	real(dp) :: o2scol(nlev)       ! O2 Slant Column
        real(dp) :: o3scol(nlev)       ! O3 Slant Column
        real(dp) :: noscol(nlev)       ! NO Slant Column
        real(dp) :: rmla(nlev)         ! Transmission, Lyman Alpha (other species)
	real(dp) :: ro2la(nlev)        ! Transmission, Lyman Alpha (for JO2)
	real(dp) :: tlevmin(nlev)
	real(dp) :: abs_col(nlev)
        real(dp) :: tsrb(nlev,nsrbtuv) ! Transmission in the SRB
	real(dp) :: xs_o2srb(nlev,nsrbtuv)       ! Cross section * QY for O2 in SRB

        real(dp) :: hc, be_o2, be_o3

        real(dp), dimension(nlev)       :: heatcon
        real(dp), dimension(nlev)       :: xh
!       real(dp), dimension(nw)         :: wrk,wrk2               ! wrk arrays
        real(dp), dimension(nw)         :: wrk                    ! wrk array
        real(dp), dimension(nsrb_tot)   :: wrk3,wrk4
        real(dp), dimension(nlev,nw)    :: trans_o3        ! Transmission, ozone
	real(dp), dimension(nlev,nw)    :: trans_o2        ! Transmission o2 (total)

!>>SW see feature #415 and feature #469 in HAMMOZ readmine
      ! interpolate daily etf onto current time of the day, only use if
      ! spectrally resolved solar radiation is used in radiation as well
      ! if isolard /= 1 then etfphot*(:) from jshort_init will be used
      if (.not. lyr_perp .and. isolrad  == 1) then
         etfphot(:) = wgtd1*etfphot_daily(:,ndw1)+wgtd2*etfphot_daily(:,ndw2)
         etfphot_ms93(:) =  wgtd1*etfphot_ms93_daily(:,ndw1)+wgtd2*etfphot_ms93_daily(:,ndw2)
      endif
!<<SW

!hs only used in hammonia
!       real(dp) :: fact_27d(34), fact_27d_ms93(4)

!     factors for Jan-Jun 1990
!       data fact_27d /0.111117_dp, 0.10801_dp, 0.0922031_dp, &
!         0.067548_dp, 0.0715905_dp, 0.0709129_dp, 0.0697145_dp, &
!         0.0811673_dp, 0.103972_dp, 0.0687287_dp, 0.0722705_dp, &
!         0.0753679_dp, 0.0686696_dp, 0.0912396_dp, 0.046475_dp, &
!         0.0556049_dp, 0.0438403_dp, 0.0371127_dp, 0.0326213_dp, &
!         0.0308256_dp, 0.0290922_dp, 0.0309952_dp, 0.0315118_dp, &
!         0.0382804_dp, 0.0304522_dp, 0.0266075_dp, 0.0262314_dp, &
!         0.0268217_dp, 0.0266829_dp, 0.0258015_dp, 0.0243137_dp, &
!         0.0260429_dp, 0.0242706_dp, 0.0225591_dp/

!       data fact_27d_ms93 /0.0328889_dp, 0.0275014_dp, 0.0266294_dp, 0.0257492_dp/

!       if (l27day) then
!         do iw = 1,nw
!           etfphot_sv(iw) = etfphot(iw) * (1._dp + fact_27d(iw)*sin_27d)
!         enddo
!         do iw = 1,nw_ms93
!           etfphot_ms93_sv(iw) = etfphot_ms93(iw) * (1._dp + fact_27d_ms93(iw)*sin_27d)
!         enddo
!       else
!         etfphot_sv = etfphot
!         etfphot_ms93_sv = etfphot_ms93
!       endif

!       heatcon=6.022e23_dp*1.e3_dp/(amu*cp)
!       hc=6.626196e-34_dp*2.997925e8_dp*1.e9_dp
!       be_o2=8.20224e-19_dp
!       be_o3=1.6821e-19_dp

      if( .not. allocated (fnorm) ) then
        allocate( fnorm(nlev,nw),stat=astat )
        if( astat /= 0 ) then
          write(message_text,'(a,i0)') 'failed to allocate fnorm; error = ',astat
          CALL finish('jshort', message_text)
        end if
      end if

      call sphers( nlev, zkm, zen, dsdh, nid )

!------------------------------------------------------------------------------
!     ... Derive O2, O3, and NO Slant Column
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!     ... Derive Slant Path for Spherical Atmosphere
!------------------------------------------------------------------------------
      delz(1:nlev-1) = km2cm*(zkm(1:nlev-1) - zkm(2:nlev))
      call slant_col( nlev, delz, dsdh, nid, o2cc, o2scol, lat, ip )
      call slant_col( nlev, delz, dsdh, nid, o3cc, o3scol, lat, ip )
      call slant_col( nlev, delz, dsdh, nid, nocc, noscol, lat, ip )

!------------------------------------------------------------------------------
!     ... Transmission due to ozone
!------------------------------------------------------------------------------

        do iw = 1,nw
           abs_col(:)     = min( (xs_o3a(iw) + xs_o3b(iw))*o3scol(:),100._dp )
           trans_o3(:,iw) = exp( -abs_col(:) )
        end do

!------------------------------------------------------------------------------
!    ... Derive the cross section and transmission for
!        molecular oxygen Lya, SRC, and SRB's
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!     ... Transmission due to molecular oxygen in the SRC
!------------------------------------------------------------------------------
        do iw = 1,nsrc_tot
           abs_col(:) = min( xs_o2src(iw)*o2scol(:),100._dp )
           trans_o2(:,iw) = exp( -abs_col(:) )
        end do

!------------------------------------------------------------------------------
!     ... Transmission and cross sections due to O2 at lyman alpha
!------------------------------------------------------------------------------
        call lymana( nlev, o2scol, rmla, ro2la )

!------------------------------------------------------------------------------
!    ... Place lya reduction faction in transmission array
!    ... This must follow the SRC placement (above)
!------------------------------------------------------------------------------

!++ost: from moz3
!	trans_o2(:,:,2) = rmla(:,:)  ! old setting
        trans_o2(:,lya_ndx) = rmla(:)
!--ost

!------------------------------------------------------------------------------
!     ... Molecular Oxygen, SRB
!------------------------------------------------------------------------------
!     ... Koppers Grid (see Koppers and Murtagh, Ann. Geophys., 14, 68-79, 1996)
!        #    wl(i)       wl(i+1)
!        1   174.4        177.0
!        2   177.0        178.6
!        3   178.6        180.2
!        4   180.2        181.8
!        5   181.8        183.5
!        6   183.5        185.2
!        7   185.2        186.9
!        8   186.9        188.7
!        9   188.7        190.5
!        10  190.5        192.3
!        11  192.3        194.2
!        12  194.2        196.1
!        13  196.1        198.0
!        14  198.0        200.0  <<last wl bin for <200nm
!        ----------------------
!        15  200.0        202.0
!        16  202.0        204.1
!        17  204.1        205.8
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!     ... Kopper SRB parameterization is used here
!------------------------------------------------------------------------------
      tlevmin(nlev:1:-1) = min( max( tlev(:),150._dp ),350._dp )
      call calc_o2srb( nlev, nid, o2scol, tlevmin, tsrb, xs_o2srb, lat, ip, long )

!------------------------------------------------------------------------------
!     ... Place Koppers SRB transmission (1-14) on user
!         grid #21 (wc=175.7nm)
!------------------------------------------------------------------------------
      do iw = 1,nsrb_tot
         trans_o2(:,iw+nsrc_tot) = tsrb(:,iw)
      end do

!------------------------------------------------------------------------------
!     ... Derive the normalize flux at each altitude,
!         corrected for O2 and O3 absorption
!------------------------------------------------------------------------------
      do iw = 1,nw     ! nw = 34 (nsrb_tot+nsrc_tot)
!       fnorm(:,:,iw) = etfphot_sv(iw)*trans_o2(:,:,iw)*trans_o3(:,:,iw)
        fnorm(:,iw) = etfphot(iw)*trans_o2(:,iw)*trans_o3(:,iw)
      end do

!------------------------------------------------------------------------------
!     ... Derive the O2 rate constant and apply branching ratio (QY)
!------------------------------------------------------------------------------
!     ... SRC and SRB QY
!         Longward  of 174.65 the product is O2 + hv => O(3P) + O(3P)
!         Shortward of 174.65 the product is O2 + hv => O(3P) + O(1D)
!         The QY is assumed to be unity in both wavelength ranges.
!
!     ... Lyman Alpha QY
!         O2 + hv -> O(3P) + O(3P) at Lyman Alpha has a QY = 0.47
!         O2 + hv -> O(3P) + O(1D) at Lyman Alpha has a QY = 0.53
!         Lacoursiere et al., J. Chem. Phys. 110., 1949-1958, 1999.
!------------------------------------------------------------------------------
!     ... Lyman Alpha
!----------------------------
!++ost: from moz3
!     jo2_lya(:,:) = etfphot_sv(2)*ro2la(:,:)*wlintv(2)
      jo2_lya(:) = etfphot(lya_ndx)*ro2la(:)*wlintv(lya_ndx)
!--ost

!hs heating computation only used in hammonia
!     rheat_o2(:,nlev:1:-1) = rheat_o2(:,nlev:1:-1) + &
!        jo2_lya(:,:) * ( hc/wc(2) - be_o2 ) * heatcon(:, nlev:1:-1) * 0.95_dp


      wrk(1:nsrc_tot) = xs_o2src(1:nsrc_tot)*wlintv(1:nsrc_tot)
!++ost: from moz3
!     wrk(2)          = 0._dp
      wrk(lya_ndx)    = 0._dp
!--ost

!     wrk2(1:nsrc_tot) = wrk(1:nsrc_tot)/wc(1:nsrc_tot)

!----------------------------
!     ... SRC
!----------------------------
      
      jo2_src = 0.0_dp
      jo2_src(:) = matmul( fnorm(:,1:nsrc_tot),wrk(1:nsrc_tot) )
!--ost

!hs heating computation only used in hammonia
!     do k=1,nlev
!       do ii=1,nr_act_po
!         rheat_o2(ii,nlev-k+1) = rheat_o2(ii,nlev-k+1) + ( xh(ii,k) * hc - &
!            jo2_src(ii,k) * be_o2 ) * heatcon(ii, nlev-k+1) * effsrc(k)
!       end do
!     end do
         
!----------------------------
!     ... SRB
!----------------------------
!++ost: from moz3
      jo2_srb = 0.0_dp

      do iz = 1,nlev
        wrk3(1:nsrb_tot) = xs_o2srb(iz,1:nsrb_tot)*wlintv(nsrc_tot+1:nsrc_tot+nsrb_tot)
        jo2_srb(iz) = dot_product( fnorm(iz,nsrc_tot+1:nsrc_tot+nsrb_tot),wrk3(1:nsrb_tot) )
      end do

!      do iz = 1,nlev
!        izz=nlev-iz+1
!        do n=1,nsrb_tot
!         do ii=1,nr_act_po
!            wrk3(ii,n) = xs_o2srb(ii,iz,n)*wlintv(nsrc_tot+n)
!!            wrk4(ii,n) = wrk3(ii,n)/wc(nsrc_tot+n)
!          end do
!        end do
!        jo2_srb(:,iz)   = 0._dp
!        dot2            = 0._dp
!        do n=1,nsrb_tot
!          do ii=1,nr_act_po
!            jo2_srb(ii,iz) = jo2_srb(ii,iz)+fnorm(ii,iz,nsrc_tot+n)*wrk3(ii,n)
!!            dot2(ii)       = dot2(ii)      +fnorm(ii,iz,nsrc_tot+n)*wrk4(ii,n)
!          end do
!        end do
!hs heating computation only used in hammonia
!!        rheat_o2(:,izz) = rheat_o2(:,izz) + ( dot2 * hc - jo2_srb(:,iz) * be_o2 ) * heatcon(:, izz)
!     end do
!--ost

!------------------------------------------------------------------------------
!     ... Total JO2
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!     ... Branch 1, O2 + hv => O(3P) + O(3P); wavelengths >175nm
!------------------------------------------------------------------------------
      jo2_sht(:,1) = jo2_lya(:)*.47_dp + jo2_srb(:)
!------------------------------------------------------------------------------
!     ... Branch 2, O2 + hv => O(3P) + O(1D);  wavelengths <175nm
!------------------------------------------------------------------------------
      jo2_sht(:,2) = jo2_lya(:)*.53_dp + jo2_src(:)

!------------------------------------------------------------------------------
!     ... Derive the NO rate constant Minsch. and Siskind, JGR, 98, 20401, 1993
!------------------------------------------------------------------------------
!     call calc_jno( nr_act_po, nlev, etfphot_ms93_sv, n2cc, o2scol, o3scol, &
!                    noscol, jno_sht )
      call calc_jno( nlev, etfphot_ms93, n2cc, o2scol, o3scol, &
                     noscol, jno_sht )

!------------------------------------------------------------------------------
!    ... Derive addtional rate constants for species with wl < 200 nm.!
!        Temperature dependence of the cross sections are not included in this
!        version.
!------------------------------------------------------------------------------
      jsht(:,:) = matmul( fnorm,xs_wl )

!     do mj=1,nj
!        do mk=1,nlev
!           do mi = 1,nw
!              do ii=1,nr_act_po
!                jsht(ii,mk,mj) = jsht(ii,mk,mj)+fnorm(ii,mk,mi)*xs_wl(mi,mj)
!              end do
!           end do
!        end do
!     end do


!  Computation of heating through O3 absorption at wavelengths below 200 nm
!  Note: The contribution is very small (less than 0.005 K/day at 1hP and less than 0.015 K/day at 0.001 hPa).
!       Therefore, the computation is currently ignored.
!     jsht_o3(:,:) = jsht(:,:,1) + jsht(:,:,2)
!     xs_wl_o3(:)=(xs_wl(:,1)+xs_wl(:,2))/wc(:)

!     fnorm_xs_o3(:,:)=0.
!     do mk=1,nlev
!       do mi = 1,nw
!          do ii=1,nr_act_po
!             fnorm_xs_o3(ii,mk) = fnorm_xs_o3(ii,mk) + fnorm(ii,mk,mi) * xs_wl_o3(mi)
!          end do
!       end do
!     end do

!     do mk=1,nlev
!       mkk=nlev-mk+1
!       do ii=1,nr_act_po
!          rheat_o3(ii,mkk) = rheat_o3(ii,mkk) + &
!                   ( fnorm_xs_o3(ii,mk) * hc - &
!                   jsht_o3(ii,mk) * be_o3 ) * heatcon(ii,mkk) * effhart(mk)
!       end do
!     end do

    deallocate( fnorm )

    end subroutine jshort

    subroutine sphers( nlev, z, zenith_angle, dsdh, nid )
!=============================================================================!
!   Subroutine sphers                                                         !
!=============================================================================!
!   PURPOSE:                                                                  !
!   Calculate slant path over vertical depth ds/dh in spherical geometry.     !
!   Calculation is based on:  A.Dahlback, and K.Stamnes, A new spheric model  !
!   for computing the radiation field available for photolysis and heating    !
!   at twilight, Planet.Space Sci., v39, n5, pp. 671-683, 1991 (Appendix B)   !
!=============================================================================!
!   PARAMETERS:                                                               !
!   NZ      - INTEGER, number of specified altitude levels in the working (I) !
!             grid                                                            !
!   Z       - REAL, specified altitude working grid (km)                  (I) !
!   ZEN     - REAL, solar zenith angle (degrees)                          (I) !
!   DSDH    - REAL, slant path of direct beam through each layer crossed  (O) !
!             when travelling from the top of the atmosphere to layer i;      !
!             DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                             !
!   NID     - INTEGER, number of layers crossed by the direct beam when   (O) !
!             travelling from the top of the atmosphere to layer i;           !
!             NID(i), i = 0..NZ-1                                             !
!=============================================================================!
!   EDIT HISTORY:                                                             !
!   Original: Taken By Doug Kinnison from Sasha Madronich, TUV Code, V4.1a,   !
!             on 1/1/02                                                       !
!=============================================================================!

      use mo_moz,   only : pi, d2r, rearth

!------------------------------------------------------------------------------
!       ... Dummy arguments
!------------------------------------------------------------------------------
      integer, intent(in)    :: nlev                ! number model vertical levels
      integer, intent(out)   :: nid(0:nlev)         ! see above
      real(dp), intent (in)  :: zenith_angle        ! zenith_angle
      real(dp), intent (in)  :: z(nlev)             ! geometric altitude (km)
      real(dp), intent (out) :: dsdh(0:nlev,nlev)   ! see above


!------------------------------------------------------------------------------
!       ... Local variables
!------------------------------------------------------------------------------
      real(dp), parameter ::  radius = rearth*1.e-3_dp   ! radius earth (km)

      real(dp) :: re
      real(dp) :: zenrad
      real(dp) :: rpsinz
      real(dp) :: const0
      real(dp) :: rj
      real(dp) :: rjp1
      real(dp) :: dsj
      real(dp) :: dhj
      real(dp) :: ga
      real(dp) :: gb
      real(dp) :: sm
      real(dp) :: zd(0:nlev-1)

      integer :: i
      integer :: j
      integer :: k
      integer :: id
      integer :: nlayer

!------------------------------------------------------------------------------
!       ... set zenith angle in radians
!------------------------------------------------------------------------------
      zenrad = zenith_angle*d2r
      const0 = SIN( zenrad )

!------------------------------------------------------------------------------
!       ... set number of layers:
!------------------------------------------------------------------------------
      nlayer = nlev - 1

!------------------------------------------------------------------------------
!       ... include the elevation above sea level to the radius of the earth:
!------------------------------------------------------------------------------
      re = radius + z(nlev)

!------------------------------------------------------------------------------
!       ... inverse coordinate of z
!------------------------------------------------------------------------------
      do k = 0,nlayer
        zd(k) = z(k+1) - z(nlev)
      end do

!------------------------------------------------------------------------------
!       ... initialize dsdh(i,j), nid(i)
!------------------------------------------------------------------------------
      nid(:) = 0
      do j = 1,nlev
        dsdh(:,j) = 0._dp
      end do

!------------------------------------------------------------------------------
!       ... calculate ds/dh of every layer
!------------------------------------------------------------------------------
      do i = 0,nlayer
        rpsinz = (re + zd(i)) * const0
        if( zenith_angle <= 90._dp .or. rpsinz >= re ) then
!------------------------------------------------------------------------------
! Find index of layer in which the screening height lies
!------------------------------------------------------------------------------
           id = i 
           if( zenith_angle > 90._dp ) then
              do j = 1,nlayer
                 if( rpsinz < (zd(j-1) + re) .and.  rpsinz >= (zd(j) + re) ) then
		    id = j
		    exit
		 end if
              end do
           end if
 
           do j = 1,id
             sm = 1._dp
             if( j == id .and. id == i .and. zenith_angle > 90._dp ) then
                sm = -1._dp
             end if
             rj   = re + zd(j-1)
             rjp1 = re + zd(j)
             dhj  = zd(j-1) - zd(j)
             ga   = max( rj*rj - rpsinz*rpsinz,0._dp )
             gb   = max( rjp1*rjp1 - rpsinz*rpsinz,0._dp )
             if( id > i .and. j == id ) then
                dsj = SQRT( ga )
             else
                dsj = SQRT( ga ) - sm*SQRT( gb )
             end if
             dsdh(i,j) = dsj / dhj
           end do
           nid(i) = id
        else
           nid(i) = -1
        end if
      end do

      end subroutine sphers

      subroutine slant_col( nlev, delz, dsdh, nid, absden, scol, lat, ip )
!=============================================================================!
!   PURPOSE:                                                                  !
!   Derive Column
!=============================================================================!
!   PARAMETERS:                                                               !
!   NLEV   - INTEGER, number of specified altitude levels in the working  (I) !
!            grid                                                             !
!   DELZ   - REAL, specified altitude working grid (km)                   (I) !
!   DSDH   - REAL, slant path of direct beam through each layer crossed  (O)  !
!             when travelling from the top of the atmosphere to layer i;      !
!             DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                             !
!   NID    - INTEGER, number of layers crossed by the direct beam when   (O)  !
!             travelling from the top of the atmosphere to layer i;           !
!             NID(i), i = 0..NZ-1                                             !
!            specified altitude at each specified wavelength                  !
!   absden - REAL, absorber concentration, molecules cm-3                       !
!   SCOL   - REAL, absorber Slant Column, molecules cm-2                        !
!=============================================================================!
!   EDIT HISTORY:                                                             !
!   09/01  Read in profile from an input file, DEK                            !
!   01/02  Taken from Sasha Madronich's TUV code                              !
!=============================================================================!

      use mo_submodel,     ONLY: lhammonia

!------------------------------------------------------------------------------
!       ... Dummy arguments
!------------------------------------------------------------------------------
      integer, intent(in) :: nlev
      integer, intent(in) :: lat, ip
      integer, intent(in) :: nid(0:nlev)           ! see above
      real(dp), intent(in)    :: delz(nlev)        ! layer thickness (cm)
      real(dp), intent(in)    :: dsdh(0:nlev,nlev) ! see above
      real(dp), intent(in)    :: absden(nlev)      ! absorber concentration (molec. cm-3)
      real(dp), intent(out)   :: scol(nlev)        ! absorber Slant Column (molec. cm-2)

!------------------------------------------------------------------------------
!       ... Local variables
!------------------------------------------------------------------------------
      real, parameter :: largest = 1.e+36_dp

      real(dp) :: sum
      real(dp) :: hscale
      real(dp) :: numer, denom
      real(dp) :: cz(nlev)

      integer :: id
      integer :: j
      integer :: k

!------------------------------------------------------------------------------
!     ... compute column increments (logarithmic integrals)
!------------------------------------------------------------------------------
      do k = 1,nlev-1
        if( absden(k) > 0._dp .and. absden(k+1) > 0._dp ) then
          cz(nlev-k) = (absden(k) - absden(k+1))/LOG( absden(k)/absden(k+1) ) * delz(k)
        else
          cz(nlev-k) = .5_dp*(absden(k) + absden(k+1)) * delz(k)
        end if
      end do

!------------------------------------------------------------------------------
!     ... Include exponential tail integral from infinity to model top
!         specify scale height near top of data.
!------------------------------------------------------------------------------
      if( .not. lhammonia ) then
        hscale  = 10.0e5_dp
      else
        hscale  = 50.0e5_dp
      endif
      cz(nlev-1) = cz(nlev-1) + hscale * absden(1)

!------------------------------------------------------------------------------
!       ...  Calculate vertical and slant column from each level:
!            work downward
!------------------------------------------------------------------------------
      do id = 0,nlev-1
         sum = 0.
         if( nid(id) >= 0 ) then
!------------------------------------------------------------------------------
!       ...  Single pass layers:
!------------------------------------------------------------------------------
            do j = 1, MIN(nid(id), id)
               sum = sum + cz(nlev-j)*dsdh(id,j)
            end do
!------------------------------------------------------------------------------
!       ...  Double pass layers:
!------------------------------------------------------------------------------
            do j = MIN(nid(id),id)+1, nid(id)
               sum = sum + 2._dp*cz(nlev-j)*dsdh(id,j)
            end do
         else
            sum = largest
         end if
         scol(nlev-id) = sum
      end do
      scol(nlev) = .95_dp*scol(nlev-1)

      end subroutine slant_col

      subroutine lymana( nlev, o2scol, rm, ro2 )
!-----------------------------------------------------------------------------!
!   PURPOSE:                                                                  !
!   Calculate the effective absorption cross section of O2 in the Lyman-Alpha !
!   bands and an effective O2 optical depth at all altitudes.  Parameterized  !
!   after:  Chabrillat, S., and G. Kockarts, Simple parameterization of the   !
!   absorption of the solar Lyman-Alpha line, Geophysical Research Letters,   !
!   Vol.24, No.21, pp 2659-2662, 1997.                                        !
!-----------------------------------------------------------------------------!
!   PARAMETERS:                                                               !
!   nz      - INTEGER, number of specified altitude levels in the working (I) !
!             grid                                                            !
!   o2scol  - REAL, slant overhead O2 column (molec/cc) at each specified (I) !
!             altitude                                                        !
!   dto2la  - REAL, optical depth due to O2 absorption at each specified  (O) !
!             vertical layer                                                  !
!   xso2la  - REAL, molecular absorption cross section in LA bands        (O) !
!-----------------------------------------------------------------------------!
!   EDIT HISTORY:                                                             !
!   01/15/2002 Taken from Sasha Madronich's TUV Version 4.1a, Doug Kinnison   !                  !
!   01/15/2002 Upgraded to F90, DK                                            !
!-----------------------------------------------------------------------------!

!------------------------------------------------------------------------------
!       ... Dummy arguments
!------------------------------------------------------------------------------
      integer, intent(in)     :: nlev
      real(dp), intent(in)    :: o2scol(nlev)
      real(dp), intent(out)   :: ro2(nlev)
      real(dp), intent(out)   :: rm(nlev)

!------------------------------------------------------------------------------
!     ... Local variables
!------------------------------------------------------------------------------
      real(dp),save :: b(3)
      real(dp),save :: c(3)
      real(dp),save :: d(3)
      real(dp),save :: e(3)

      data b / 6.8431e-01_dp, 2.29841e-01_dp,  8.65412e-02_dp /, &
           c / 8.22114e-21_dp, 1.77556e-20_dp,  8.22112e-21_dp /, &
           d / 6.0073e-21_dp, 4.28569e-21_dp,  1.28059e-20_dp /, &
           e / 8.21666e-21_dp, 1.63296e-20_dp,  4.85121e-17_dp /

      integer  :: i, k
      real(dp) :: wrk, term

!------------------------------------------------------------------------------
!     ... Calculate reduction factors at every altitude
!------------------------------------------------------------------------------
      do k = 1,nlev
        wrk = 0._dp
        do i = 1,2 ! pc Dan Marsh
          term = e(i)*o2scol(k)
          if( term < 100._dp ) then
            wrk = wrk + d(i) * exp( -term )
          end if
        end do
        ro2(k) = wrk
        wrk = 0._dp
        do i = 1,3
          term = c(i)*o2scol(k)
          if( term < 100._dp ) then
            wrk = wrk + b(i) * exp( -term )
          end if
        end do
        rm(k) = wrk
      end do

      end subroutine lymana

      subroutine calc_o2srb( nlev, nid, o2col, tlev, tsrb, xscho2, lat, ip, long )

      USE mo_exception,       ONLY: message, message_text, em_debug

!-----------------------------------------------------------------------------!
!   PURPOSE:                                                                  !
!   Calculate the equivalent absorption cross section of O2 in the SR bands.  !
!   The algorithm is based on parameterization of G.A. Koppers, and           !
!   D.P. Murtagh [ref. Ann.Geophys., 14 68-79, 1996]                          !
!   Final values do include effects from the Herzberg continuum.              !
!-----------------------------------------------------------------------------!
!   PARAMETERS:                                                               !
!   NZ      - INTEGER, number of specified altitude levels in the working (I) !
!             grid                                                            !
!   O2COL   - REAL, slant overhead O2 column (molec/cc) at each specified (I) !
!             altitude                                                        !
!   TLEV    - tmeperature at each level                                   (I) !
!   TSRB    - REAL, transmission for the SRB                                  !
!   XSCHO2  - REAL, molecular absorption cross section in SR bands at     (O) !
!             each specified wavelength.  Includes Herzberg continuum         !
!-----------------------------------------------------------------------------!
!   EDIT HISTORY: Taken from TUV, 1/17/2002                                   !
!   This code was supplied to TUV by Dan Marsh.                               !
!-----------------------------------------------------------------------------!

!------------------------------------------------------------------------------
!       ... Dummy arguments
!------------------------------------------------------------------------------
      integer, intent(in) :: nlev
      integer, intent(in) :: lat, ip, long
      integer, intent(in) :: nid(0:nlev)
      real(dp), intent (in)   :: o2col(nlev)
      real(dp), intent (in)   :: tlev(nlev)
      real(dp), intent (out)  :: tsrb(nlev,nsrbtuv)
      real(dp), intent (out)  :: xscho2(nlev,nsrbtuv)

!------------------------------------------------------------------------------
!     ... Local variables
!------------------------------------------------------------------------------
      real(dp), parameter :: min_o2col = 38._dp
      real(dp), parameter :: max_o2col = 56._dp

      integer  :: i, k
      real(dp) :: dto2
      real(dp) :: den, num
      real(dp) :: term1, term2
      real(dp) :: dtsrb(nlev)
      real(dp) :: xs(nsrbtuv)

!------------------------------------------------------------------------------
!     ... Calculate cross sections
!------------------------------------------------------------------------------
      do k = 1,nlev
        call effxs( min( max( log( o2col(k) ),min_o2col ),max_o2col ), tlev(k), xs )
        xscho2(k,:) = xs(:)
      end do

!-------------------------------------------------------
!     ... Calculate incremental optical depths
!-------------------------------------------------------
wave_loop : &
      do i = 1,nsrbtuv 
        do k = 1,nlev-1
          if( nid(nlev-k) /= -1 ) then
!-------------------------------------------------------
!     ... Calculate an optical depth weighted by density
!-------------------------------------------------------
            num   = xscho2(k+1,i)*o2col(k+1) - xscho2(k,i)*o2col(k)
            if( num == 0. ) then
              WRITE(message_text,'(a,4e25.15,a,5i0)') 'calc_o2srb : o2col(k:k+1),xscho2(k:k+1,i) = ',       &
                                                      o2col(k:k+1),xscho2(k:k+1,i),' @ i,k,lat,ip,long = ', &
                                                      i,k,lat,ip,long
              CALL message('moz_jshort', message_text, level=em_debug)
            end if
            term1 = log( xscho2(k+1,i)/xscho2(k,i) )
            term2 = log( o2col(k+1)/o2col(k) )
            if( term2 == 0. ) then
              WRITE(message_text,'(a,4e25.15,a,5i0)') 'calc_o2srb : o2col(k:k+1),xscho2(k:k+1,i) = ',       &
                                                    o2col(k:k+1),xscho2(k:k+1,i),' @ i,k,lat,ip,long = ', & 
                                                    i,k,lat,ip,long
             CALL message('moz_jshort', message_text, level=em_debug)
            end if
            den = 1._dp + log( xscho2(k+1,i)/xscho2(k,i) )/log( o2col(k+1)/o2col(k) )
            dto2 = -num/den
            if( abs( dto2 ) < 100._dp ) then
              dtsrb(k) = exp( -dto2 )
            else
              dtsrb(k) = 0._dp
            end if
          else
            dtsrb(k) = 0._dp
          end if
        end do
!-----------------------------------------------
!     ... Calculate Transmission for SRB
!-----------------------------------------------
        tsrb(nlev,i) = 1._dp
        do k = nlev-1,1,-1
          tsrb(k,i) = tsrb(k+1,i)*dtsrb(k)
        end do
      end do wave_loop

      end subroutine calc_o2srb

      subroutine effxs( x, t, xs )
!-------------------------------------------------------------
!     Subroutine for evaluating the effective cross section
!     of O2 in the Schumann-Runge bands using parameterization
!     of G.A. Koppers, and D.P. Murtagh [ref. Ann.Geophys., 14
!     68-79, 1996]
!
!     method:
!     ln(xs) = A(X)[T-220]+B(X)
!     X = log of slant column of O2
!     A,B calculated from chebyshev polynomial coeffs
!     AC and BC using NR routine chebev.  Assume interval
!     is 38<ln(NO2)<56.
!
!     Revision History:
!
!     drm 2/97  initial coding
!-------------------------------------------------------------

!-------------------------------------------------------------
!	... Dummy arguments
!-------------------------------------------------------------
      real(dp), intent(in)  :: x
      real(dp), intent(in)  :: t
      real(dp), intent(out) :: xs(nsrbtuv)

!-------------------------------------------------------------
!	... Local variables
!-------------------------------------------------------------
      real(dp) :: a(nsrbtuv), b(nsrbtuv)

      call calc_params( x, a, b )

      xs(:) = exp( a(:)*(t - 220._dp) + b(:) )

      end subroutine effxs

      subroutine calc_params( x, a, b )
!-------------------------------------------------------------
!       calculates coefficients (A,B), used in calculating the
!       effective cross section, for 17 wavelength intervals
!       as a function of log O2 column density (X)
!       Wavelength intervals are defined in WMO1985
!-------------------------------------------------------------

!-------------------------------------------------------------
!	... Dummy arguments
!-------------------------------------------------------------
      real(dp), intent(in)  :: x
      real(dp), intent(out) :: a(nsrbtuv), b(nsrbtuv)

!-------------------------------------------------------------
!	... Local variables
!-------------------------------------------------------------
      real(dp), parameter :: min_o2col = 38._dp
      real(dp), parameter :: max_o2col = 56._dp
      integer :: i

!-------------------------------------------------------------
!     ... call chebyshev evaluation routine to calc a and b from
!	    set of 20 coeficients for each wavelength
!-------------------------------------------------------------
      do i = 1,nsrbtuv
        a(i) = jchebev( min_o2col, max_o2col, ac(1,i), 20, x )
        b(i) = jchebev( min_o2col, max_o2col, bc(1,i), 20, x )
      end do

      contains

      function jchebev( a, b, c, m, x )

      USE mo_exception,        ONLY : message, message_text, em_debug

!-------------------------------------------------------------
!     Chebyshev evaluation algorithm
!     See Numerical recipes p193
!-------------------------------------------------------------

!-------------------------------------------------------------
!	... Dummy arguments
!-------------------------------------------------------------
      integer, intent(in) :: m
      real(dp), intent(in)    :: a, b, x
      real(dp), intent(in)    :: c(m)

      real :: jchebev
!-------------------------------------------------------------
!	... Local variables
!-------------------------------------------------------------
      integer  :: j
      real(dp) :: d, dd, sv, y, y2

      if( (x-a)*(x-b) > 0._dp ) then
        WRITE(message_text,'(a,e25.15)') 'x not in range in chebev', x
        CALL message('moz_jchebev', message_text, level=em_debug)
        jchebev = 0._dp
        return
      end if

      d  = 0._dp
      dd = 0._dp
      y  = (2._dp*x - (a + b))/(b-a)
      y2 = 2._dp*y
      do j = m,2,-1
        sv = d
        d  = y2*d - dd + c(j)
        dd = sv
      end do

      jchebev = y*d - dd + .5_dp*c(1)

      end function jchebev

      end subroutine calc_params

      subroutine calc_jno( nlev, etfphot_ms93, n2cc, o2scol, o3scol, &
                           noscol, jno )
!-----------------------------------------------------------------------------!
!   PURPOSE:                                                                  !
!   Compute the total photolytic rate constant for NO in the SR bands         !
!     - following the approach of Minshwanner and Siskind, JGR,               !
!       98, D11, 20401-20412, 1993.                                           !
!                                                                             !
!-----------------------------------------------------------------------------!
!   PARAMETERS:                                                               !
!   NZ           - INTEGER, number of specified altitude levels               !
!                                                                             !
!   etfphot_ms93 - Extraterrestrial Flux, within the MS 1993 Grid             !
!                  units of photons cm-2 sec-1 nm-1                           !
!   n2cc         - N2 conc (molecules cm-3)                                   !
!   o3scol       - Ozone Slant Column (molecules cm-2)                        !
!   o2scol       - Oxygen Slant Column (molecules cm-2)                       !
!   noscol       - Nitric Oxide Slant Column(molecules cm-2)                  !
!                                                                             !
!   LOCAL VARIABLES:                                                          !
!   tauo3        - Transmission factor in the Hartley Band of O3              !
!   etfphot_ms93 - Solar Irr. on Minschwaner and Siskind 1993 (MS93) Grid     !
!   xs_o3ms93    - O3 cross section on the MS93 Grid                          !
!                                                                             !
!   OUTPUT VARIABLES:                                                         !
!   jno          - photolytic rate constant                                   !
!                  each specified altitude                                    !
!                                                                             !
!-----------------------------------------------------------------------------!
!   EDIT HISTORY:                                                             !
!   08/01  Created, Doug Kinnison, NCAR, ACD                                  !
!-----------------------------------------------------------------------------!

!------------------------------------------------------------------------------
!       ... Dummy arguments
!------------------------------------------------------------------------------
      integer, intent(in)     :: nlev
      real(dp), intent(in)    :: etfphot_ms93(num_ms93tuv)
      real(dp), intent(in)    :: n2cc(nlev)
      real(dp), intent(in)    :: o3scol(nlev)
      real(dp), intent(in)    :: o2scol(nlev)
      real(dp), intent(in)    :: noscol(nlev)
      real(dp), intent(out)   :: jno(nlev)

!------------------------------------------------------------------------------
!	... Local variables
!------------------------------------------------------------------------------
      integer     :: i, wn, lev
      real(dp)    :: jno50
      real(dp)    :: jno90
      real(dp)    :: jno100
      real(dp)    :: tauo3(nlev,num_ms93tuv)

!------------------------------------------------------------------------------
!   	... O3 SRB Cross Sections from WMO 1985, interpolated onto MS, 1993 grid
!------------------------------------------------------------------------------
      real(dp), save :: xso3_ms93(num_ms93tuv) = (/ 7.3307600e-19_dp, &
                                                    6.9660105E-19_dp, 5.9257699E-19_dp, 4.8372219E-19_dp /)

!------------------------------------------------------------------------------
!   	... delta wavelength of the MS, 1993 grid
!------------------------------------------------------------------------------
      real(dp), save :: delta_wc_ms93(num_ms93tuv) = (/ 1.50_dp, 1.50_dp, 6.0_dp, 2.3_dp /)

!------------------------------------------------------------------------------
!   	... O2 SRB Cross Sections for the six ODF regions, MS, 1993
!------------------------------------------------------------------------------
      real(dp), save :: cs250(6)  = (/ 1.117e-23_dp, 2.447e-23_dp, 7.188e-23_dp, &
                                       3.042e-22_dp, 1.748e-21_dp, 1.112e-20_dp /)
      real(dp), save :: cs290(6)  = (/ 1.350e-22_dp, 2.991e-22_dp, 7.334e-22_dp, &
                                       3.074e-21_dp, 1.689e-20_dp, 1.658e-19_dp /)
      real(dp), save :: cs2100(6) = (/ 2.968e-22_dp, 5.831e-22_dp, 2.053e-21_dp, &
                                       8.192e-21_dp, 4.802e-20_dp, 2.655e-19_dp /)

!------------------------------------------------------------------------------
!     ... derive tauo3 for the three o2 srb
!     ... wn = 1,2, and 4 are used below for jno
!------------------------------------------------------------------------------
      do wn = 1,num_ms93tuv
         tauo3(:,wn) = exp( -xso3_ms93(wn)*o3scol(:) )
      end do

!------------------------------------------------------------------------------
!   	... Call PJNO Function to derive SR Band JNO contributions
!         Called in order of wavelength interval (shortest firs)
!------------------------------------------------------------------------------
      do lev = 1,nlev
         jno100   = pjno( 1, cs2100, wtno100, csno100 )
         jno90    = pjno( 2, cs290,  wtno90,  csno90 )
         jno50    = pjno( 4, cs250,  wtno50,  csno50 )
         jno(lev) = jno50 + jno90 + jno100
      end do

      contains

      real function pjno( w, cso2, wtno, csno )
!------------------------------------------------------------------------------
!   	... uses xsec at center of g subinterval for o2
!           uses mean values for no
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!	... parameters
!------------------------------------------------------------------------------
      integer, parameter :: ngint = 6
      integer, parameter :: nno = 2

!----------------------------------------------------------------
!	... Dummy arguments
!----------------------------------------------------------------
      integer, intent(in)     :: w
      real(dp),    intent(in) :: cso2(ngint)
      real(dp),    intent(in) :: csno(ngint,nno)
      real(dp),    intent(in) :: wtno(ngint,nno)


!----------------------------------------------------------------
!	... Local variables
!----------------------------------------------------------------
      integer ::  jj, i, k
      real(dp) :: tauno
      real(dp) :: transno
      real(dp) :: transo2
      real(dp) :: tauo2
      real(dp) :: jno
      real(dp) :: jno1

!----------------------------------------------------------------
!	... derive the photolysis frequency for no within a given
!         srb (i.e., 5-0, 9-0, 10-0)
!----------------------------------------------------------------
      jno = 0._dp

      do k = 1,ngint
        tauo2 = o2scol(lev) * cso2(k)
        if( tauo2 < 50._dp ) then
          transo2 = exp( -tauo2 )
        else
          transo2 = 0._dp
        end if
        jno1 = 0._dp
        do jj = 1,nno
          tauno = noscol(lev)*csno(k,jj)
          if( tauno < 50._dp ) then
            transno = exp( -tauno )
          else
            transno = 0._dp
          end if
          jno1 = jno1 + csno(k,jj) * wtno(k,jj) * transno
        end do
        jno = jno + jno1*transo2
      end do

      pjno = delta_wc_ms93(w)*etfphot_ms93(w)*tauo3(lev,w)*jno

!----------------------------------------------------------------
!	... correct for the predissociation of the deltq 1-0
!         transition in the srb (5-0)
!----------------------------------------------------------------
      if( w == 4 ) then
        pjno = 1.65e9_dp/(5.1e7_dp + 1.65e9_dp + (1.5e-9_dp*n2cc(nlev-lev+1)))*pjno
      end if

      end function pjno

      end subroutine calc_jno

end module mo_moz_jshort

