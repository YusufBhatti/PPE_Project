!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!!  mo_moz_waccm_photo
!!
!! \brief
!!  This module calculates photolysis frequencies based on a tabulated approach
!!
!! \author Douglas E. Kinnison (NCAR)
!!
!! \responsible_coder
!!  Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!  - DK: original version (1999-11-08)
!!  - Stacy Walters and Martin Schultz (2011-04-08): revision for HAMMOZ
!!  - MGS (2013-07-25): clean-up, decoupling of ndx_ and rid_
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

module mo_moz_waccm_photo
!----------------------------------------------------------------------
!	... Photolysis interp table and related arrays for stratosphere
!----------------------------------------------------------------------

      USE mo_kind,    ONLY : dp
!++mgs 20130228: implemented fastj as "replacement" for jlong only
      USE mo_moz,     ONLY : lfastj
!--mgs

      implicit none

      private
      public :: waccm_prate_inti, waccm_photo, set_ub_col, setcol, diurnal_geom, sundis

      save

!+++sschr20140311: changed
      real(dp), parameter :: max_zen_angle = 97.01_dp
! comment H.Schmidt: A test with max_zen_angle = 105 was crashing, with
!     max_zen_angle = 100 running. However, the table in jlong is only valid
!     until 94, so I stay with this for the moment.

      integer, allocatable :: sht_indexer(:), lng_indexer(:)
      real(dp) :: o2_exo_col, o3_exo_col

      ! species indices
      integer :: ndx_n2, ndx_no, ndx_o2, ndx_o3
      ! reaction indices
      integer :: rid_jo2_a, rid_jo2_b, rid_jno
      logical :: lo_inv_o2, lo_inv_n2

      contains

      subroutine waccm_prate_inti
!----------------------------------------------------------------------
!     ... Read in the stratospheric photorate tables and arrays
!----------------------------------------------------------------------

      use mo_exception,  only : finish, message, message_text, em_debug
      use mo_moz_util,   only : get_spc_ndx, get_inv_ndx, get_rxt_ndx
      use mo_moz_mods,   only : phtcnt, pht_alias_lst, pht_alias_mult, rxt_tag_lst
      use mo_moz_jshort, only : jshort_init, nsht => nj
      use mo_moz_jlong,  only : jlong_init, nlng => numj ! added nlng for fastj
      use mo_moz_photo,  only : jdefined   
!++mgs
      use mo_moz_fastj,  only : setfjx, jvn_
!--mgs
      use mo_submodel,   only : lhammonia

      use mo_mpi,        only : p_parallel_io

      implicit none

!----------------------------------------------------------------------
!        ... Local variables
!----------------------------------------------------------------------
      integer :: i, astat, ndx

!----------------------------------------------------------------------
!        ... check for short wavelength photorates
!----------------------------------------------------------------------
      rid_jo2_a    = get_rxt_ndx('jo2_a')
      rid_jo2_b    = get_rxt_ndx('jo2_b')
      rid_jno      = get_rxt_ndx('jno')  
      if( rid_jo2_a <= 0 .or. rid_jo2_b <= 0 .or. rid_jno <= 0 ) then
        write(message_text,'(a)') 'Chemistry scheme must have O2 and NO photolysis'
        CALL finish('waccm_prate_inti', message_text)
      endif
!----------------------------------------------------------------------
!        ... check for absorber species
!----------------------------------------------------------------------
      ! check presence of species as solution or invariant species
      lo_inv_n2 = .false.
      ndx_n2 = get_spc_ndx( 'N2' )
      if( ndx_n2 <= 0 ) then
        ndx_n2 = get_inv_ndx( 'N2' )
        lo_inv_n2 = .true.
      end if
      lo_inv_o2 = .false.
      ndx_o2 = get_spc_ndx( 'O2' )
      if( ndx_o2 <= 0 ) then
        ndx_o2 = get_inv_ndx( 'O2' )
        lo_inv_o2 = .true.
      end if
      ndx_no = get_spc_ndx( 'NO' )
      ndx_o3 = get_spc_ndx( 'O3' )

      if( ndx_o2 <= 0 .or. ndx_o3 <= 0  .or. &
          ndx_no <= 0 .or. ndx_n2 <= 0 ) then
        write(message_text,'(a)') 'Chemistry must have O2, O3, NO, and N2 columns'
        CALL finish('waccm_prate_inti', message_text)
      endif

!----------------------------------------------------------------------
!        ... Initialize photolysis modules
!----------------------------------------------------------------------
      !-- oxygen and ozone column above model top
      ! values for 0.1 hPa model top
!+++sschr20140311: !!!!!!!! commented out because hammoz model top is 0.01 hPa
!     if( lhammonia ) then
!        o2_exo_col = 0._dp
!     else
!        o2_exo_col = 4.203e20_dp
!     end if

!     if( lhammonia ) then
!        o3_exo_col = 0._dp
!     else
!        o3_exo_col = 3.588e14_dp
!     end if
      o2_exo_col = 0._dp
      o3_exo_col = 0._dp
!---sschr
      allocate( sht_indexer(phtcnt),stat=astat )
      if( astat /= 0 ) then
        write(message_text,'(a,i0)') 'Failed to allocate sht_indexer; error = ',astat
        CALL finish('waccm_prate_inti', message_text)
      end if
      sht_indexer(:) = 0
      call jshort_init( sht_indexer )

!     WRITE(message_text,'(a)') 'done with jshort_init'
!     CALL message('waccm_prate_inti', message_text, level=em_debug)

      allocate( lng_indexer(phtcnt),stat=astat )
      if( astat /= 0 ) then
        write(message_text,'(a,i0)') 'Failed to allocate lng_indexer; error = ',astat
        CALL finish('waccm_prate_inti', message_text)
      end if
      lng_indexer(:) = 0
!++mgs
      IF (lfastj) THEN
        call setfjx( lng_indexer )
        nlng = jvn_
      ELSE
        call jlong_init( lng_indexer )
      END IF
!--mgs

!     WRITE(message_text,'(a)') 'done with jlong_init'
!     CALL message('waccm_prate_inti', message_text, level=em_debug)

!----------------------------------------------------------------------
!        ... check that each photorate is in short or long datasets
!----------------------------------------------------------------------

      IF( any( abs(sht_indexer(:)) + abs(lng_indexer(:)) == 0 ) ) THEN
        WRITE(message_text,'(2a)') 'waccm_prate_inti: the following photorate(s) are not in',&
                                  ' either the short or long datasets'
        CALL message('waccm_prate_inti',message_text, level=em_debug)
        DO ndx = 1,phtcnt
           IF( abs(sht_indexer(ndx)) + abs(lng_indexer(ndx)) == 0 ) THEN
              WRITE(message_text,'(a,a)') '           ',TRIM( rxt_tag_lst(ndx) )
              CALL message('waccm_prate_inti',message_text, level=em_debug)
           END IF
        END DO
        CALL finish('waccm_prate_inti', 'Missing photolysis rate parameters')
      END IF

!----------------------------------------------------------------------
!        ... output any aliased photorates
!----------------------------------------------------------------------
!     IF( any( pht_alias_lst(:,1) /= ' ' ) ) THEN

!        WRITE(message_text,'(a)') 'waccm_prate_inti: the following short photorate(s) are aliased'
!        CALL message('waccm_prate_inti', message_text, level=em_debug)

!          DO ndx = 1,phtcnt
!             IF( pht_alias_lst(ndx,1) /= ' ' ) THEN
!                IF( pht_alias_mult(ndx,1) == 1._dp ) THEN
!                   WRITE(message_text,'(4a)') ' ',TRIM(rxt_tag_lst(ndx)),' -> ',TRIM(pht_alias_lst(ndx,1))
!                   CALL message('waccm_prate_inti', message_text, level=em_debug)
!                ELSE
!                   WRITE(message_text,'(3a,e25.15,2a)') ' ',TRIM(rxt_tag_lst(ndx)),' -> ', &
!                                                        pht_alias_mult(ndx,1),'*',TRIM(pht_alias_lst(ndx,1))
!                   CALL message('waccm_prate_inti', message_text, level=em_debug)
!                END IF
!             END IF
!          END DO
!     END IF

!     IF( any( pht_alias_lst(:,2) /= ' ' ) ) THEN

!        WRITE(message_text,'(a)') 'waccm_prate_inti: the following long photorate(s) are aliased'
!        CALL message('waccm_prate_inti', message_text, level=em_debug)

!          DO ndx = 1,phtcnt
!             IF( pht_alias_lst(ndx,2) /= ' ' ) THEN
!                IF( pht_alias_mult(ndx,2) == 1. ) THEN
!                   WRITE(message_text,'(4a)') ' ',TRIM(rxt_tag_lst(ndx)),' -> ',TRIM(pht_alias_lst(ndx,2))
!                   CALL message('waccm_prate_inti', message_text, level=em_debug)
!                ELSE
!                   WRITE(message_text,'(3a,e25.15,2a)') ' ',TRIM(rxt_tag_lst(ndx)),' -> ', pht_alias_mult(ndx,2), &
!                                                        '*',TRIM(pht_alias_lst(ndx,2))
!                   CALL message('waccm_prate_inti', message_text, level=em_debug)
!                END IF
!             END IF
!          END DO
!     END IF

! tell ECHAM which photorates are defined
      jdefined(:) = .true.

!     WRITE(message_text,'(a)') 'end of waccm_prate_inti'
!     CALL message('waccm_prate_inti', message_text, level=em_debug)

      end subroutine waccm_prate_inti


!!    subroutine waccm_photo( lat, ip, photos, pmid, pdel, &
!!                            temper, zmid, col_dens, zen_angle, srf_alb, &
!!                            lwc, clouds, sunon, sunoff, esfact, &
!!                            vmr, invariants, amu, cp, rheat_o2, rheat_o3)
!++mgs 20130228: added arguments for fastj
      subroutine waccm_photo( lat, ip, caldayn, photos, pmid, pfull, pdel, &
                              temper, tsurf, zmid, col_dens, zen_angle, srf_alb, &
                              lwc, iwc, clouds, cdnc, sunon, sunoff, esfact, &
                              vmr, invariants, amu, cp, rheat_o2, rheat_o3)
!--mgs
!-----------------------------------------------------------------
!	Calculate the photolysis rate constants
!-----------------------------------------------------------------

      use mo_moz_mods, only   : plonl, plev, pcnstm1, ncol_abs, phtcnt, indexm
      use mo_moz_mods, only   : pht_alias_mult
      use mo_moz_jshort, only : nsht => nj, jshort
      use mo_moz_jlong,  only : nlng => numj, jlong
      use mo_moz,        only : r2d
      use mo_moz_fastj,  only : fjx_interface
!     use mo_timer,      only : time_diff, elapsed, cdate, ctime
      use mo_exception,  only : finish, message_text

      use mo_mpi,        only : p_parallel_io

      implicit none

!-----------------------------------------------------------------
!   	... Dummy arguments
!-----------------------------------------------------------------
      integer, intent(in)   ::  lat, ip
      real(dp), intent(out) ::  photos(plonl,plev,phtcnt)                ! photodissociation rates (1/s)
      real(dp), intent(in)  ::  caldayn, &                               ! calendar day
                                pmid(plonl,plev), &                      ! midpoint pressure (Pa)
                                pfull(plonl,plev+1), &                   ! halflevel pressure (Pa)
                                pdel(plonl,plev), &                      ! del pressure about midpoint (Pa)
                                temper(plonl,plev), &                    ! midpoint temperature (K)
                                tsurf(plonl), &                          ! surface temperature (K)
                                zmid(plonl,plev), &                      ! midpoint height (km)
                                col_dens(plonl,plev,max(1,ncol_abs)), &  ! column densities (molecules/cm^2)
                                zen_angle(plonl), &                      ! solar zenith angle (radians)
                                srf_alb(plonl), &                        ! surface albedo
                                lwc(plonl,plev), &                       ! liquid water content (kg/kg)     
                                iwc(plonl,plev), &                       ! ice water content (kg/kg)
                                clouds(plonl,plev), &                    ! cloud fraction                  
                                cdnc(plonl,plev), &                      ! condensation nuclei
                                sunon(plonl), &                          ! angle for sunrise (radians)
                                sunoff(plonl), &                         ! angle for sunset (radians)
                                esfact, &                                ! earth-sun distance factor
                                vmr(plonl,plev,max(1,pcnstm1)), &        ! tracer vmr
!++ost: allow for use of more invariants
                                invariants(plonl,plev,10), &     ! invariant densities (molecules/cm^3)
!--ost
                                amu(plonl,plev), &
                                cp(plonl,plev)
      real(dp), intent(inout) ::rheat_o2(plonl,plev), &                  ! temperature tendency due to O2 heating
                                rheat_o3(plonl,plev)                     ! temperature tendency due to O3 heating

!-----------------------------------------------------------------
!    	... Local variables
!-----------------------------------------------------------------
      integer  ::  i,  &                  ! indicies
                   ii, &
                   k,  &
                   m,  &
                   n,  &
                   astat

      real(dp)                  :: sza    ! solar zenith angle (degrees)
      real(dp), dimension(plev) :: &
                   eff_alb, &             ! effective albedo from cloud modifications
                   cld_mult, &            ! clould multiplier
                   tline, &               ! vertical temperature array
                   zarg, &                ! vertical height array
!++ost: additional pressure dependance for jlong
                   parg, &                ! vertical pressure array (hPa)
!--ost
                   colo3, &               ! vertical o3 column density
                   fac1, &                ! work array
                   cld_line, &            ! vertical cloud array
                   lwc_line, &            ! vertical lwc array
                   o2cc, &                ! o2 density (molecules/cm^3)
                   o3cc, &                ! o3 density (molecules/cm^3)
                   nocc, &                ! no density (molecules/cm^3)
                   n2cc, &                ! n2 density (molecules/cm^3)
                   heat_o2, &             ! work version of rheat_o2
                   heat_o3, &             ! work version of rheat_o2
                   amu_c, &               ! work version of amu
                   cp_c                   ! work version of cp

      real(dp), allocatable :: &
                   sht_prates(:,:), &     ! short photorates
                   lng_prates(:,:), &     ! long  photorates
                   fjx_prates(:,:,:)      ! long  photorates for fastj (extra dimension)

      logical :: do_diag

      real(dp), dimension(plev,2) :: jo2_sht            ! o2 short photorate
      real(dp), dimension(plev)   :: jno_sht            ! no short photorate
!++ost: use pht_alias
      real(dp) ::  alias_factor           ! aliasing multiplier
!--ost


      do_diag = .false.
!-----------------------------------------------------------------
!	... Zero all photorates
!-----------------------------------------------------------------
      photos(:,:,:) = 0.0_dp

!>>>##DEBUG##
!sschr: this produces a lot of output!!!
!sschr: this will slow down the run
!sschr: therefore commented out
!     if( lat == 6 ) then
!       write(0,*) '--------------------------------'
!       write(0,*) 'waccm_photo: surface albedo'
!       write(0,'(1p5g15.7)') srf_alb(:)
!       write(0,*) '--------------------------------'
!     endif
!<<<
!-----------------------------------------------------------------
!	... Allocate short, long work arrays
!-----------------------------------------------------------------
      allocate( sht_prates(plev,nsht),stat=astat )
      if( astat /= 0 ) then
        write(message_text,'(a,i0)') 'Failed to allocate sht_prates; error = ',astat
        CALL finish('waccm_photo', message_text)
      end if
      if ( lfastj ) then
        allocate( fjx_prates(plonl,plev,nlng),stat=astat )
        if( astat /= 0 ) then
          write(message_text,'(a,i0)') 'Failed to allocate fjx_prates; error = ',astat
          CALL finish('waccm_photo', message_text)
        end if
      else
        allocate( lng_prates(nlng,plev),stat=astat )
        if( astat /= 0 ) then
          write(message_text,'(a,i0)') 'Failed to allocate lng_prates; error = ',astat
          CALL finish('waccm_photo', message_text)
        end if
      end if

column_loop : &
      do i = 1,plonl
        sza = zen_angle(i)*r2d
has_light : &
        if( sza < max_zen_angle ) then
           zarg(:)  = zmid(i,:)
           parg(:)  = 1.e-2_dp*pmid(i,:)
           colo3(:) = col_dens(i,:,1)
           fac1(:)  = pdel(i,:)
           tline(:) = temper(i,:)
           amu_c(:) = amu(i,:)
           cp_c(:)  = cp(i,:)
           heat_o2(:)  = 0._dp
           heat_o3(:)  = 0._dp
           IF (lo_inv_o2) THEN
             o2cc(:)  = invariants(i,:,ndx_o2)
           ELSE
             o2cc(:)  = vmr(i,:,ndx_o2) * invariants(i,:,indexm)
           ENDIF
           o3cc(:)  = vmr(i,:,ndx_o3) * invariants(i,:,indexm)
           nocc(:)  = vmr(i,:,ndx_no) * invariants(i,:,indexm)
           IF (lo_inv_n2) THEN
             n2cc(:)  = invariants(i,:,ndx_n2)
           ELSE
             n2cc(:)  = vmr(i,:,ndx_n2) * invariants(i,:,indexm)
           ENDIF
           lwc_line(:) = lwc(i,:)
           cld_line(:) = clouds(i,:)
           call cloud_mod( zen_angle(i), cld_line, lwc_line, fac1, srf_alb(i), &
                            eff_alb, cld_mult )
           cld_mult(:) = esfact * cld_mult(:)
            call jshort( plev, sza, n2cc, o2cc, o3cc, &
                         nocc, tline, zarg, jo2_sht, jno_sht, &
                         sht_prates, heat_o2, heat_o3, amu_c, cp_c, &
                         lat, ip, i, do_diag )
!++mgs 20130228
           if (.not. lfastj) then
              call jlong( plev, sza, eff_alb, parg, tline, &
                          colo3, lng_prates, heat_o2, heat_o3, amu_c, &
                          cp_c, ip, lat, i )
           end if
!--mgs
!++lp just use longwave j values to compare with fastj
            do m = 1,phtcnt
               if( sht_indexer(m) > 0 ) then
                  alias_factor = pht_alias_mult(m,1)
                  if( alias_factor == 1._dp ) then
                     photos(i,plev:1:-1,m) = sht_prates(:,sht_indexer(m))
                  else
                     photos(i,plev:1:-1,m) = alias_factor * sht_prates(:,sht_indexer(m))
                  end if
               end if
            end do
            if( rid_jno > 0 ) then
               photos(i,plev:1:-1,rid_jno) = jno_sht(:)
            end if
            if( rid_jo2_a > 0 ) then
               photos(i,plev:1:-1,rid_jo2_a) = jo2_sht(:,2)
            end if
            if( rid_jo2_b > 0 ) then
               photos(i,plev:1:-1,rid_jo2_b) = jo2_sht(:,1)
            end if
!--lp finish test
           ! sort lng photorates into photos array (only for MOZ routine - fastj see below)
           if (.not. lfastj) then
              do m = 1,phtcnt
                 if( lng_indexer(m) > 0 ) then
                    alias_factor = pht_alias_mult(m,2)
                    if( alias_factor == 1._dp ) then
                       photos(i,:,m) = photos(i,:,m) + lng_prates(lng_indexer(m),:)
                    else
                       photos(i,:,m) = photos(i,:,m) + alias_factor * lng_prates(lng_indexer(m),:)
                    end if
                 end if
              end do
           end if
!++mgs 20130228
           ! cloud correction not needed for fastj, because clouds are treated explicitly there
           if (.not. lfastj) then
              do m = 1,phtcnt
                 photos(i,:,m) = photos(i,:,m)*cld_mult(:)
              end do
           end if
!--mgs
            rheat_o2(i,:) = heat_o2(:)
            rheat_o3(i,:) = heat_o3(:)
        end if has_light
      end do column_loop

      ! need to call fjx_interface outside column loop, because FJX has its own loop included.
      if ( lfastj ) then
         call fjx_interface(plonl,plev,lat,    & !MOZART indices mapped onto ECHAM indices kproma,klev,krow
                            caldayn,           & !calendar day
                            srf_alb,           & !surface albedo
                            pmid,              & !half level pressure
                            pfull,             & !full level pressure
                            pdel,              & !layer thickness
                            temper,            & !full level temperature
                            tsurf,             & !surface temperature
                            lwc,               & !totel water content (water+ice)
                            iwc,               & !ice content
                            clouds,            & !cloud fraction
                            cdnc,              & !cdnc field
                            vmr(:,:,ndx_o3),   & !ozone vmr 
                            fjx_prates         ) !photolysis values

           do m = 1,phtcnt
              if( lng_indexer(m) > 0 ) then
                 alias_factor = pht_alias_mult(m,2)
                 if( alias_factor == 1._dp ) then
                    photos(:,:,m) = photos(:,:,m) + fjx_prates(:,:,lng_indexer(m))
                 else
                    photos(:,:,m) = photos(:,:,m) + alias_factor * fjx_prates(:,:,lng_indexer(m))
                 end if
              end if
           end do
        end if

      deallocate( sht_prates )
      if ( lfastj ) then
         deallocate( fjx_prates )
      else
         deallocate( lng_prates )
      endif

!++mgs: moz2_photo has a lot of analogy tests here (jmpan, jmacr, jglyald, etc. ...)

      end subroutine waccm_photo

      subroutine cloud_mod( zen_angle, clouds, lwc, delp, srf_alb, &
                            eff_alb, cld_mult )
!-----------------------------------------------------------------------
! 	... Cloud alteration factors for photorates and albedo
!-----------------------------------------------------------------------

!++mgs: replaced mo_grid with mo_moz_mods
      use mo_moz_mods, only : plev, plevm

      implicit none

      real(dp), parameter :: gi = 1._dp/9.80616_dp

!-----------------------------------------------------------------------
! 	... Dummy arguments
!-----------------------------------------------------------------------
      real(dp), intent(in)    ::  zen_angle         ! zenith angle
      real(dp), intent(in)    ::  srf_alb           ! surface albedo
      real(dp), intent(in)    ::  clouds(plev)       ! cloud fraction
      real(dp), intent(in)    ::  lwc(plev)          ! liquid water content (mass mr)
      real(dp), intent(in)    ::  delp(plev)         ! del press about midpoint in pascals
      real(dp), intent(out)   ::  eff_alb(plev)      ! effective albedo
      real(dp), intent(out)   ::  cld_mult(plev)     ! photolysis mult factor

!-----------------------------------------------------------------------
! 	... Local variables
!-----------------------------------------------------------------------
      integer :: k
      real(dp)    :: coschi
      real(dp)    :: del_lwp(plev)
      real(dp)    :: del_tau(plev)
      real(dp)    :: above_tau(plev)
      real(dp)    :: below_tau(plev)
      real(dp)    :: above_cld(plev)
      real(dp)    :: below_cld(plev)
      real(dp)    :: above_tra(plev)
      real(dp)    :: below_tra(plev)
      real(dp)    :: fac1(plev)
      real(dp)    :: fac2(plev)

!---------------------------------------------------------
!	... Modify lwc for cloud fraction and form
!	    liquid water path for each layer
!---------------------------------------------------------
      where( clouds(:) /= 0._dp )
         del_lwp(:) = gi * lwc(:) * delp(:) * 1.e3_dp / clouds(:)
      elsewhere
         del_lwp(:) = 0._dp
      endwhere
!---------------------------------------------------------
!    	... Form tau for each model layer
!---------------------------------------------------------
      where( clouds(:) /= 0._dp )
         del_tau(:) = del_lwp(:) *.155_dp * clouds(:)**1.5_dp
      elsewhere
         del_tau(:) = 0._dp
      end where
!---------------------------------------------------------
!    	... Form integrated tau from top down
!---------------------------------------------------------
      above_tau(1) = 0._dp
      do k = 1,plevm
         above_tau(k+1) = del_tau(k) + above_tau(k)
      end do
!---------------------------------------------------------
!    	... Form integrated tau from bottom up
!---------------------------------------------------------
      below_tau(plev) = 0._dp
      do k = plevm,1,-1
         below_tau(k) = del_tau(k+1) + below_tau(k+1)
      end do
!---------------------------------------------------------
!	... Form vertically averaged cloud cover above and below
!---------------------------------------------------------
      above_cld(1) = 0._dp
      do k = 1,plevm
         above_cld(k+1) = clouds(k) * del_tau(k) + above_cld(k)
      end do
      do k = 2,plev
         if( above_tau(k) /= 0._dp ) then
            above_cld(k) = above_cld(k) / above_tau(k)
         else
            above_cld(k) = above_cld(k-1)
         end if
      end do
      below_cld(plev) = 0._dp
      do k = plevm,1,-1
         below_cld(k) = clouds(k+1) * del_tau(k+1) + below_cld(k+1)
      end do
      do k = plevm,1,-1
         if( below_tau(k) /= 0._dp ) then
            below_cld(k) = below_cld(k) / below_tau(k)
         else
            below_cld(k) = below_cld(k+1)
         end if
      end do
!---------------------------------------------------------
!	... Modify above_tau and below_tau via jfm
!---------------------------------------------------------
      where( above_cld(2:plev) /= 0._dp )
         above_tau(2:plev) = above_tau(2:plev) / above_cld(2:plev)
      end where
      where( below_cld(:plevm) /= 0._dp )
         below_tau(:plevm) = below_tau(:plevm) / below_cld(:plevm)
      end where
      where( above_tau(2:plev) < 5._dp )
            above_cld(2:plev) = 0._dp
      end where
      where( below_tau(:plevm) < 5._dp )
         below_cld(:plevm) = 0._dp
      end where
!---------------------------------------------------------
!	... Form transmission factors
!---------------------------------------------------------
      above_tra(:) = 11.905_dp / (9.524_dp + above_tau(:))
      below_tra(:) = 11.905_dp / (9.524_dp + below_tau(:))
!---------------------------------------------------------
!	... Form effective albedo
!---------------------------------------------------------
      where( below_cld(:) /= 0._dp )
         eff_alb(:) = srf_alb + below_cld(:) * (1._dp - below_tra(:)) &
                                             * (1._dp - srf_alb)
      elsewhere
         eff_alb(:) = srf_alb
      end where
      coschi = MAX( COS( zen_angle ),.5_dp )
      where( del_lwp(:)*.155_dp < 5._dp )
         fac1(:) = 0._dp
      elsewhere
         fac1(:) = 1.4_dp * coschi - 1._dp
      end where
      fac2(:)     = MIN( 0._dp,1.6_dp*coschi*above_tra(:) - 1._dp )
      cld_mult(:) = 1._dp + fac1(:) * clouds(:) + fac2(:) * above_cld(:)
      cld_mult(:) = MAX( .05_dp,cld_mult(:) )

      end subroutine cloud_mod


      subroutine set_ub_col( col_delta, vmr, pdel, pmid )
!---------------------------------------------------------------
!        ... Set the column densities at the upper boundary
!---------------------------------------------------------------

      use mo_moz_mods,     only : plonl, plev, pcnstm1, ncol_abs

      implicit none

!---------------------------------------------------------------
!        ... Dummy args
!---------------------------------------------------------------
      real(dp), intent(out) ::    col_delta(plonl,0:plev,MAX(1,ncol_abs))  ! /cm**2
      real(dp), intent(in)  ::    vmr(plonl,plev,pcnstm1), &               ! xported species vmr
                                  pdel(plonl,plev),        &
                                  pmid(plonl,plev)         ! ++mgs

!---------------------------------------------------------------
!        NOTE: xfactor = 10.*R/(K*g) in cgs units.
!              The factor 10. is to convert pdel
!              from pascals to dyne/cm**2.
!---------------------------------------------------------------
      real(dp), parameter :: xfactor = 2.8704e21_dp/(9.80616_dp*1.38044_dp)
      integer :: k

!---------------------------------------------------------------
!        ... Assign ozone column density at the upper boundary (index = 1)
!---------------------------------------------------------------
      col_delta(:,0,1) = o3_exo_col
      if (ndx_o3 > 0) then
         do k = 1,plev
            col_delta(:,k,1) = xfactor * pdel(:,k) * vmr(:,k,ndx_o3)
         enddo
      endif
!---------------------------------------------------------------
!        ... Assign oxygen column density at the upper boundary (index = 2)
!---------------------------------------------------------------
         col_delta(:,0,2) = o2_exo_col
      if (ndx_o2 > 0) then
         do k = 1,plev
            col_delta(:,k,2) = xfactor * pdel(:,k) * vmr(:,k,ndx_o2)
         end do
      else
         do k = 1,plev
            col_delta(:,k,2) = xfactor * pdel(:,k) * .21_dp
         end do
      endif
         

      end subroutine set_ub_col

!++ost: vmr and pdel not needed here. Perhaps for HAMMONIA??? Leave in for now.
      subroutine setcol( col_delta, col_dens, vmr, pdel )
!---------------------------------------------------------------
!     	... Set the column densities
!---------------------------------------------------------------

      use mo_moz_mods, only : plonl, plev, ncol_abs, pcnstm1

      implicit none

!---------------------------------------------------------------
!     	... Dummy arguments
!---------------------------------------------------------------
      real(dp), intent(in)  ::   col_delta(plonl,0:plev,MAX(1,ncol_abs))      ! layer column densities (molecules/cm^2)
      real(dp), intent(out) ::   col_dens(plonl,plev,MAX(1,ncol_abs))         ! column densities ( /cm**2 )
      real(dp), intent(in)  ::   vmr(plonl,plev,pcnstm1)                      ! xported species vmr
      real(dp), intent(in)  ::   pdel(plonl,plev)                             ! delta about midpoints
!---------------------------------------------------------------
!        The local variables
!---------------------------------------------------------------
!! !---------------------------------------------------------------
!! !        NOTE: xfactor = 10.*R/(K*g) in cgs units.
!! !              The factor 10. is to convert pdel
!! !              from pascals to dyne/cm**2.
!! !---------------------------------------------------------------
!!       real(dp), parameter :: xfactor = 2.8704e21_dp/(9.80616_dp*1.38044_dp)

      integer  ::   k, km1      ! alt indices

!---------------------------------------------------------------
!   	... Compute column densities down to the
!           current eta index in the calling routine.
!           The first column is O3 and the second is O2.
!---------------------------------------------------------------
      col_dens(:,1,1) = col_delta(:,0,1) + .5_dp * col_delta(:,1,1)
      col_dens(:,1,2) = col_delta(:,0,2) + .5_dp * col_delta(:,1,2)
      do k = 2,plev
         km1 = k - 1
         col_dens(:,k,1) = col_dens(:,km1,1) + .5_dp * (col_delta(:,km1,1) + col_delta(:,k,1))
         col_dens(:,k,2) = col_dens(:,km1,2) + .5_dp * (col_delta(:,km1,2) + col_delta(:,k,2))
      end do

!++mgs: note - HAMMOZ (moz2) code has query for ox/o3 species index. Not sure why??
!--mgs
      end subroutine setcol

!++mgs: only parameter needed lateron in waccm_photo is zen_angle.
!       therefore calculation of other quantities disabled
      subroutine diurnal_geom( ip, lat, time_of_year, polar_night, polar_day, &
                               sunon, sunoff, loc_angle, zen_angle )
!------------------------------------------------------------------
!    	... Diurnal geometry factors
!------------------------------------------------------------------

      use mo_moz,       only : pi, twopi, pid2, dayspy, d2r
!++mgs
!     use mo_grid,      only : plong => plon, plonl
      use mo_moz_mods,  only : plonl
!--mgs
!     use mo_mpi,       only : base_lat
      use mo_geoloc,    only : philat_2d, philon_2d

!sschr: commented out (just to be sure that nothing from mo_debugs is used)
!use mo_debugs

      implicit none

!------------------------------------------------------------------
!    	... Dummy arguments
!------------------------------------------------------------------
      integer, intent(in)  ::     ip                 ! longitude index
      integer, intent(in)  ::     lat                ! latitude index
      real(dp), intent(in)     ::     time_of_year       ! time of year
      real(dp), intent(out)    ::     sunon(plonl)       ! sunrise angle in radians
      real(dp), intent(out)    ::     sunoff(plonl)      ! sunset angle in radians
      real(dp), intent(out)    ::     zen_angle(plonl)   ! solar zenith angle
      real(dp), intent(out)    ::     loc_angle(plonl)   ! "local" time angle
      logical, intent(out) ::     polar_day(plonl)   ! continuous daylight flag
      logical, intent(out) ::     polar_night(plonl) ! continuous night flag

!------------------------------------------------------------------
!        ... Local variables
!------------------------------------------------------------------
      integer ::  i
      real(dp)::  dec_max
      real(dp)::  declination
      real(dp)::  latitude
      real(dp)::  doy_loc            ! day of year
      real(dp)::  tod                ! time of day
      real(dp)::  sin_dec, cos_dec   ! sin, cos declination
      real(dp)::  cosphi             ! cos latitude
      real(dp)::  sinphi             ! sin latitude

      dec_max     = 23.45_dp * d2r
      sunon(:) = 0._dp
      sunoff(:) = 0._dp

!------------------------------------------------------------------
!        Note: this formula assumes a 365 day year !
!------------------------------------------------------------------
      doy_loc     = AINT( time_of_year )
      declination = dec_max * COS((doy_loc - 172._dp)*twopi/dayspy)
      sin_dec     = SIN( declination )
      cos_dec     = COS( declination )

!------------------------------------------------------------------
!	... Compute base for zenith angle
!------------------------------------------------------------------
      tod = (time_of_year - doy_loc) + .5_dp

      DO i=1,plonl

        latitude = philat_2d(i,lat)*d2r
        sinphi = SIN( latitude )
        cosphi = COS( latitude )
        polar_day(i) = .false.
        polar_night(i) = .false.

!-------------------------------------------------------------------
!        Note: Longitude 0 (Greenwich) at 0:00 hrs
!              maps to local angle = pi
!-------------------------------------------------------------------
!      loc_angle(:) = (/ ((tod + REAL(i+(ip-1)*plonl-1)/REAL(plong))*twopi,i = 1,plonl) /)
!      loc_angle(:) = MOD( loc_angle(:),twopi )
!      zen_angle(:) = ACOS( sinphi*sin_dec + cosphi*cos_dec*COS(loc_angle(:)) )

        loc_angle(i) = tod*twopi + philon_2d(i,lat)*d2r
        loc_angle(i) = MOD( loc_angle(i),twopi )
        zen_angle(i) = ACOS( sinphi*sin_dec + cosphi*cos_dec*COS(loc_angle(i)) )

!------------------------------------------------------------------
!        Determine if in polar day or night
!        If NOT in polar day or night then
!        calculate terminator longitudes
!------------------------------------------------------------------
!       if( ABS(latitude) >= (pid2 - ABS(declination)) ) then
!         if( SIGN(1._dp,declination) == SIGN(1._dp,latitude) ) then
!           polar_day(i) = .true.
!           sunoff(i) = 2._dp*twopi
!           sunon(i)  = -twopi
!         else
!           polar_night(i) = .true.
!         end if
!       else
!         sunoff(i) = ACOS( -TAN(declination)*TAN(latitude) )
!         sunon(i)  = twopi - sunoff(i)
!       end if

      enddo
!sschr: this will slow down the run
!sschr: therefore commented out
!     zdf3(1:plonl,lat) = philon_2d(1:plonl,lat)
!     zdf4(1:plonl,lat) = philat_2d(1:plonl,lat)
!     zdf5(1:plonl,lat) = loc_angle(1:plonl)

      end subroutine diurnal_geom

      real(dp) function sundis( idate )
!!++mgs : this routine should be obsolete in ECHAM (see mo_radiation)
!-----------------------------------------------------------------------------
!=  PURPOSE:                                                                 =*
!=  Calculate Earth-Sun distance variation for a given date.  Based on       =*
!=  Fourier coefficients originally from:  Spencer, J.W., 1971, Fourier      =*
!=  series representation of the position of the sun, Search, 2:172          =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  IDATE  - INTEGER, specification of the date, from YYMMDD              (I)=*
!=  ESRM2  - REAL, variation of the Earth-sun distance                    (O)=*
!=           ESRM2 = (average e/s dist)^2 / (e/s dist on day IDATE)^2        =*
!-----------------------------------------------------------------------------*
!=  EDIT HISTORY:                                                            =*
!=  01/95  Changed computation of trig function values                       =*
!-----------------------------------------------------------------------------*
!= This program is free software;  you can redistribute it and/or modify     =*
!= it under the terms of the GNU General Public License as published by the  =*
!= Free Software Foundation;  either version 2 of the license, or (at your   =*
!= option) any later version.                                                =*
!= The TUV package is distributed in the hope that it will be useful, but    =*
!= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
!= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
!= License for more details.                                                 =*
!= To obtain a copy of the GNU General Public License, write to:             =*
!= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
!-----------------------------------------------------------------------------*
!= To contact the authors, please mail to:                                   =*
!= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
!= send email to:  sasha@ucar.edu                                            =*
!-----------------------------------------------------------------------------*
!= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
!-----------------------------------------------------------------------------

      use mo_math_constants, only : pi
      use mo_exception,   only : finish, message_text

      implicit none

!-----------------------------------------------------------------------------
!	... Dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in) :: idate

!-----------------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------------
      integer :: iyear, imonth, iday, mday, month, jday
      integer, save :: imn(12) = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
      real(dp)    :: dayn, thet0
      real(dp)    :: sinth, costh, sin2th, cos2th

!-----------------------------------------------------------------------------
! 	... Parse date to find day number (Julian day)
!-----------------------------------------------------------------------------
      iyear  = INT( idate/10000 )
      imonth = INT( (idate - 10000*iyear)/100 )
      iday   = idate - (10000*iyear + 100*imonth)

      if( imonth > 12 ) then
         write(message_text,'(a,i0,a,i0)') 'Month in date exceeds 12: date = ', idate,  &
                               ', month = ',imonth
         CALL finish('sundis (waccm_photo)', message_text)
      end if

      if( MOD(iyear,4) == 0 ) then
         imn(2) = 29
      else
         imn(2) = 28
      end if

      if( iday > imn(imonth) ) then
         write(message_text,'(a,i0,a,i0)') 'Day in date exceeds days in month: date = ', &
                               idate, ', day = ', iday
         CALL finish('sundis (waccm_photo)', message_text)
      end if

      mday = 0
      do month = 1,imonth-1
         mday = mday + imn(month)
      end do
      jday = mday + iday
      dayn = REAL(jday - 1) + .5_dp

!-----------------------------------------------------------------------------
! 	... Define angular day number and compute esrm2:
!-----------------------------------------------------------------------------
      thet0 = 2._dp*pi*dayn/365._dp

!-----------------------------------------------------------------------------
! 	... Calculate SIN(2*thet0), COS(2*thet0) 
!-----------------------------------------------------------------------------
      sinth   = SIN( thet0 )
      costh   = COS( thet0 )
      sin2th  = 2._dp*sinth*costh
      cos2th  = costh*costh - sinth*sinth
      sundis  = 1.000110_dp + .034221_dp*costh  +  .001280_dp*sinth &
                + .000719_dp*cos2th +  .000077_dp*sin2th

      end function sundis

end module mo_moz_waccm_photo
