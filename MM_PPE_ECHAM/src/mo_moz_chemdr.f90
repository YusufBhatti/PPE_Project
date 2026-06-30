!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!!  mo_moz_chemdr
!!
!! \brief
!!  This module is the actual interface to the MOZART chemistry solver
!!
!! \author S. Walters (NCAR)
!!
!! \responsible_coder
!!  Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!  - SW: original version (pre 2004)
!!  - SR, HS: adaptation to ECHAM5 (2004)
!!  - MGS, OS: harmonized merge between ECHAM5-MOZ and HAMMONIA (2009-01-27)
!!  - MGS: cleanups (2011)
!!  - MGS: cleaned and added diagnostics (2014-03-19)
!!
!! \limitations
!!  This version only supports the implicit solver. The explicit and Rodas solvers
!!  are not implemented.
!!
!! \details
!!
!! \bibliographic_references
!!  none
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

#if defined (NAG)
#define ARGCHECK 1
#endif

module mo_moz_chemdr


      use mo_kind,   only: dp
      
      implicit none

      private
      public :: chemdr_inti
      public :: chemdr

      save
      !-- private module variables
      real(dp) :: esfact                   ! relative earth sun distance factor (1/R^2)
      !-- MOZART species indices
      integer  :: ndx_o2, ndx_o3, ndx_h2o, ndx_oh, ndx_hno3

      contains

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

      subroutine chemdr_inti

      use mo_exception,                 ONLY: finish, message, em_debug, em_error
      use mo_moz,                       ONLY: lstrathet, ltrophet
      use mo_moz_util,                  ONLY: get_spc_ndx, get_rxt_ndx
      use mo_moz_mods,                  ONLY: grpcnt, clscnt4, clscnt1, clscnt5
      use mo_moz_imp_sol,               ONLY: imp_slv_inti
      use mo_moz_strato_sad,            ONLY: sad_inti
      use mo_moz_het_strato_rates,      ONLY: ratecon_sfstrat_inti
      use mo_hammoz_het_tropo_rates,    ONLY: ratecon_sftropo_inti


  !-----------------------------------------------------------------------
  !       ... Obtain species and reaction indices
  !-----------------------------------------------------------------------
      ndx_o2   = get_spc_ndx('O2')
      ndx_o3   = get_spc_ndx('O3')
      ndx_h2o  = get_spc_ndx('H2O')
      IF (ndx_h2o <= 0)                     &
        CALL message('chemdr_inti',         &
                     'No species H2O in MOZ solution matrix. This code requires H2O!', level=em_error)
      ndx_oh   = get_spc_ndx('OH')
      ndx_hno3 = get_spc_ndx('HNO3')

  !-----------------------------------------------------------------------
  !       ... Initialize the Euler implicit solver
  !-----------------------------------------------------------------------
  ! note: the explicit and Rodas solvers are not supported in this version
      if ( clscnt4 > 0 ) then
        call imp_slv_inti
        CALL message('moz_init','finished imp_slv_inti', level=em_debug)
      end if
      if (clscnt1 > 0 .OR. clscnt5 > 0) CALL finish('moz_chemdr', 'Solver classes 1 and 5 are not implemented (explicit and Rodas)')

  !-----------------------------------------------------------------------
  !       ... Initialize stratospheric aerosols
  !       ratecon_sfstrat_inti must always be called to avoid undefined rates
  !-----------------------------------------------------------------------
      CALL ratecon_sfstrat_inti 
      IF (lstrathet) THEN
      ! initialize
        CALL sad_inti
        CALL message('moz_init','finished sad_inti', level=em_debug)
      ENDIF

  !-----------------------------------------------------------------------
  !       ... Initialize tropospheric heterogeneous reactions
  !-----------------------------------------------------------------------

      IF (ltrophet) CALL ratecon_sftropo_inti()

  !-----------------------------------------------------------------------
  !       ... More consistency checking
  !-----------------------------------------------------------------------
      IF (grpcnt > 0) CALL finish('moz_chemdr', 'Grpcnt > 0. Wrong mozpp file.')

      end subroutine chemdr_inti




!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                        MOZ names:                          ECHAM names:
      subroutine chemdr( ncdate, ncsec, kproma, lat,   &   ! yymmdd, sec, kproma, krow
                         mmr,                          &   ! tracer mass mixing ratios (inout)
                         sh, cwat, cice,               &   ! water content             (inout)
                         cheat, rheat,                 &   ! zhch, zhrad               (inout)
                         delt,                         &   ! time_step_len
                         pfull, pmid, pdel, pgeom1,    &   ! paphp1, papp1, zpdel, pgeom1
                         tsurf, zma, zi,               &   ! tsurfm, zma, zmah
                         cldfr,                        &   ! aclc, xlm1+xim1, xim1
                         tfld,                         &   ! ptm1
                         prho,                         &   ! rhoam1
                         sad_sage,                     &   ! stratospheric aerosol density from SAGE
                         albs,                         &   ! zalbedo
                         amu, cp,                      &   ! amu, cp
                         cdnc )                            ! CDNC field (for photolysis) 

!-----------------------------------------------------------------------
!     ... Chemdr advances the volumetric mixing ratio
!         forward one time step via a combination of explicit,
!         ebi, hov, fully implicit, and/or rodas algorithms.
!-----------------------------------------------------------------------

!-- ECHAM use statements
      use mo_exception,            only : finish, message, message_text,        &
                                          em_warn, em_info, em_debug, em_error
      use mo_submodel,             only : lhammonia
      USE mo_control,              only : ltimer
      USE mo_timer,                only : timer_start,                          &
                                          timer_stop
      use mo_decomposition,        only : ldc => local_decomposition
      use mo_radiation_parameters, only : flx_ratio_cur

!-- MOZART use statements
      use mo_moz,                  only : lchemsolv, lphotolysis, lstrathet
      use mo_moz_mods,             only : plonl, plev, plevp, pcnstm1, plnplv,  &
                                          indexm, phtcnt, gascnt, rxntot,       &
                                          rxt_tag_lst, clscnt4,                 &
                                          ncol_abs, nfs, extcnt, hetcnt
      use mo_moz_util,             only : inti_mr_xform, caldayr, negtrc,       &
                                          mmr2vmr, vmr2mmr,                     &
                                          h2o_to_vmr
      use mo_moz_photo,            only : jdefined, photo_mem
      use mo_moz_waccm_photo,      only : set_ub_col, setcol, waccm_photo,      &
                                          diurnal_geom
      use mo_moz_subs,             only : setrxt, adjrxt, phtadj
      use mo_moz_tropo_rates,      only : tropo_rates
      use mo_moz_setinv,           only : setinv
      use mo_moz_imp_sol,          only : imp_sol
      use mo_moz_strato_sad,       only : sad_strat_calc, sad_top
      use mo_moz_het_strato_rates, only : ratecon_sfstrat, zero_sfstrat
      use mo_moz_aero_settling,    only : strat_aer_settling

!-- HAMMOZ use statements
      use mo_hammoz_timer,         only : timer_moz_chem_prep,                  &
                                          timer_moz_chem_photo,                 &
                                          timer_moz_chem_impsol
!
!-- HAMMONIA use statements
!     use mo_setext_moz3,          only : setext_moz3
!     use mo_noyion,               only : noyion
!     use mo_noycosmic,            only : noycosmic
!-- Diagnostics
      use mo_moz_diag,             only : dpohconc,                             &
                                          dpsad1, dpsad2, dpsad3,               &
                                          dpratech4, dpratemcl,                 &
                                          dphno3_cond, moz_family_diag,         &
                                          moz_rates_diag, mozh_diag
!-- Debugging 
      use mo_debugs
      use mo_mpi,                  only: p_pe, p_parallel_io
      use mo_time_control,         only: lstart, lresume

      implicit none

!-----------------------------------------------------------------------
!        ... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) ::  ncdate,        &        ! date at end present time step
                              ncsec,         &        ! seconds relative to ncdate
                              kproma,        &        ! ECHAM vector length (for diagnostics)
                              lat                     ! latitude index
      real(dp), intent(inout) ::  mmr(plonl,plev,pcnstm1) ! xported species ( mmr )
      real(dp), intent(inout) ::  sh(plonl,plev),     &   ! specific humidity ( kg/kg )
                                  cwat(plonl,plev),   &   ! total cloud water (kg/kg)
                                  cice(plonl,plev),   &   ! ice cloud water (kg/kg )
                                  cheat(plonl,plev),  &   ! chemical heating rate
                                  rheat(plonl,plev)       ! radiative heating rate
      real(dp), intent(in)    ::  delt,               &   ! timestep in seconds
                                  pfull(plonl,plevp), &   ! interface (full level) press ( pascals )
                                  pmid(plonl,plev),   &   ! midpoint press ( pascals )
                                  pdel(plonl,plev),   &   ! delta press across midpoints
                                  pgeom1(plonl,plev), &   ! geopotential height
                                  tsurf(plonl),       &   ! surface temperature                    
                                  zma(plonl,plev),    &   ! abs geopot height at midpoints ( m )
                                  zi(plonl,plevp),    &   ! abs geopot height at interfaces ( m )
                                  cldfr(plonl,plev),  &   ! cloud fraction                   
                                  tfld(plonl,plev),   &   ! midpoint temperature
                                  prho(plonl,plev),   &   ! air density                
                                  sad_sage(plonl,plev), &   ! sageii surface area density (cm^2/cm^3)
                                  albs(plonl),        &   ! surface albedo
                                  amu(plonl,plev),    &   ! molecular mass of air
                                  cp(plonl,plev)          ! specific heat
      real(dp), intent(in)    ::  cdnc(plonl,plev)

!-----------------------------------------------------------------------
!     	... Local variables
!-----------------------------------------------------------------------
      integer, parameter  ::  ip = 1                   ! longitude tile index
      integer, parameter  :: nmaxinv = 10              ! maximum number of "invariants"
      real(dp), parameter ::  m2km = 1.e-3_dp
      real(dp), parameter ::  Pa2hPa = 1.e-2_dp
      real(dp), parameter :: secperyear = 365._dp*24._dp*3600._dp  ! seconds per year

      integer             ::  i, k, m, n, jt
      integer             ::  o3_ndx, o_ndx, o1d_ndx
      integer             ::  n_ndx, no_ndx, no2_ndx, no3_ndx, n2o5_ndx,  &
                              hno3_ndx, ho2no2_ndx, clono2_ndx, brono2_ndx
      real(dp)            ::  caldayn                  ! day of year at end of time step


      real(dp)     ::  ps(plonl)               ! surface press ( pascals )
      real(dp)     ::  invariants(plonl,plev,nmaxinv)
      real(dp)     ::  col_dens(plonl,plev,max(1,ncol_abs))     ! column densities (molecules/cm^2)
      real(dp)     ::  col_delta(plonl,0:plev,max(1,ncol_abs))  ! layer column densities (molecules/cm^2)
      real(dp)     ::  het_rates(plonl,plev,max(1,hetcnt))
      real(dp)     ::  extfrc(plonl,plev,max(1,extcnt))
      real(dp)     ::  vmr(plonl,plev,pcnstm1), &               ! xported species ( vmr )
                       reaction_rates(plonl,plev,rxntot)
! Variables for stratospheric heterogeneous chemistry
      real(dp), dimension(plonl,plev) :: &
                   hno3_gas, &           ! gas phase hno3 (vmr)
                   h2o_gas, &            ! gas phase h2o (vmr)
                   h2o_cond, &           ! condensed phase h2o (vmr)
                   pmb                   ! midpoint pressure (hpa)
      real(dp), dimension(plonl,plev,3) :: &
                   hno3_cond             ! condensed phase hno3 (vmr; sulfate,sm NAT, lg NAT)
      real(dp), dimension(plonl,plev,4) :: &
                   radius_strat, &       ! radius of aerosols (cm; sulfate, sm NAT, lg NAT, water-ice)
                   sad_strat             ! surface area density (cm-1; sulfate, sm NAT, lg NAT, water-ice)
! Other stuff...
      real(dp), dimension(plonl,plev) :: &
                   h2ovmr, &             ! water vapor volume mixing ratio
                   mbar, &               ! mean dry atmospheric mass ( amu ) - was wet mass in ori MOZART
                   zmid                  ! midpoint geopotential in km
      real(dp), dimension(plonl,plevp) :: &
                   zint                  ! interface geopotential in km
! Variables for diurnal_geometry. Not used except for zen_angle in waccm_photo.
      real(dp)           :: zen_angle(plonl), loc_angle(plonl)
      real(dp)           :: sunon(plonl), sunoff(plonl)
      logical            :: polar_night(plonl), polar_day(plonl)
      REAL(dp), POINTER  :: photo_mem_p(:,:,:)
      real(dp), dimension(plonl,plev) :: rheat_o2, rheat_o3

      real(dp), dimension(plonl,plevp) :: &
                   pio_hv, pio_e, pio_c  ! external no production
!++mgs20140319: added these (moved here from moz_chemistry)
! Diagnostics
      real(dp)::   zkch4(plonl,plev),             &  ! reaction rate constant OH+CH4
                   zkmcl(plonl,plev),             &  ! reaction rate constant OH+CH3CCl3
                   zrate(plonl,plev)                 ! buffer for reaction turnover
!--mgs




! #debug# >>>>
      INTEGER :: itest
      LOGICAL :: ldebug
      LOGICAL, SAVE :: lfirst = .true.
! #debug# <<<<
!
!-----------------------------------------------------------------------      
!        ... Initialize
!-----------------------------------------------------------------------      
      IF (ltimer) CALL timer_start( timer_moz_chem_prep )
      ldebug = ( (lstart .OR. lresume) .AND. ( lat == 1 .OR. lat == DNINT(0.5_dp * ldc% ngpblks) ) )
      IF ( ldebug ) THEN
        write(message_text,'(a,i0,a,i0)') 'Starting on CPU ', p_pe, ' for lat = ',lat
        CALL message('moz_chemdr', message_text, level=em_debug)
      END IF

      reaction_rates(:,:,:) = -999._dp

!--- extract surface pressure for convenience and compatibility with MOZART
      ps(:) = pfull(:,plevp)

!--- Date and time
      caldayn = caldayr( ncdate, ncsec )

!--- Calculate parameters for diurnal geometry
!        The only parameter which is really used lateron is zen_angle
!        in waccm_photo.
!        We keep the other parameters from diurnal_geometry to preserve
!        the calling sequence to waccm_photo. Note that they are no longer
!        calculated in the diurnal_geometry subroutine.

      if( phtcnt /= 0 ) then
         call diurnal_geom( ip, lat, caldayn, polar_night, polar_day, &
                            sunon, sunoff, loc_angle, zen_angle )
      end if

!--- Set water vapor: always use sh (q from ECHAM) to initialize h2ommr variable
      mmr(:,:,ndx_h2o) = sh(:,:)
!--- Xform between mass and volume mixing ratios
      call inti_mr_xform( amu, mbar )
      call mmr2vmr( vmr, mmr, mbar )
      h2ovmr(:,:) = vmr(:,:,ndx_h2o)    ! convenience variable

!--- Xform geopotential height from m to km 
      do k = 1,plev
         zmid(:,k) = m2km * zma(:,k)
      end do
      do k = 1,plevp
         zint(:,k) = m2km * zi(:,k)
      end do

!--- Set the "invariants"
      if( nfs > 0 ) then
         call setinv( invariants, tfld, h2ovmr, pmid, vmr )
      end if

     IF (ldebug) THEN
        WRITE(message_text,'(a)') 'after setinv'
        CALL message('moz_chemdr', message_text, level=em_debug)
     ENDIF

!-----------------------------------------------------------------------
!        ... Stratospheric heterogeneous chemistry
!-----------------------------------------------------------------------
!--- Initialize hno3 condensed, gas phase: all to gas
!    zero_sfstrat must always be called to avoid undefined rates
      call zero_sfstrat( reaction_rates )
inti_strathet : &
      if (lstrathet) then
        IF (ldebug) THEN
           WRITE(message_text,'(a)') 'initializing hno3 for lstrathet'
           CALL message('moz_chemdr', message_text, level=em_debug)
        ENDIF

        do k = 1,plev
           hno3_gas(:,k)  = vmr(:,k,ndx_hno3)
           hno3_cond(:,k,:) = 0._dp
           h2o_gas(:,k)=h2ovmr(:,k)
        end do

!--- Xform condensed water from mmr to mole fraction (vmr equivalent)
!    cice is mmr (xi in ECHAM) -> h2o_cond
        call h2o_to_vmr( h2o_cond, cice, mbar, plonl )
!--- Convert mid level pressure from Pa to hPa
        pmb(:,:)=pmid(:,:)*Pa2hPa
        IF (ldebug) THEN
           WRITE(message_text,'(a)') 'before call to sad_strat_calc'
           CALL message('moz_chemdr', message_text, level=em_debug)
        ENDIF

!--- calculation is only done on levels for which 2 hPa < P < 300 hPa is valid
!--- initialize sad_strat

        sad_strat = 0._dp

!--- SAD calculation
#ifdef ARGCHECK
        call sad_strat_calc( lat, invariants(:,:,indexm), pmb, tfld, &
                             hno3_gas, hno3_cond, h2o_gas, h2o_cond, sad_sage, &
                             radius_strat, sad_strat, plonl )
#else
        call sad_strat_calc( lat, invariants(1,1,indexm), pmb, tfld, &
                             hno3_gas, hno3_cond, h2o_gas, h2o_cond, sad_sage, &
                             radius_strat, sad_strat, plonl )
#endif

        do k = 1,plev
           vmr(:,k,ndx_hno3) = hno3_gas(:,k)
        end do

!--- Set aerosol reaction rates
        do k = 1,plev
!--- ratecon_sfstrat needs sum of small and large NAT particles as input
           sad_strat(:,k,2) = sad_strat(:,k,2) + sad_strat(:,k,3)
        end do
#ifdef ARGCHECK
        call ratecon_sfstrat( invariants(:,:,indexm), pmid, tfld, radius_strat, sad_strat, &
                            sad_strat(:,:,2), sad_strat(:,:,4), h2ovmr, vmr, reaction_rates )
#else
        call ratecon_sfstrat( invariants(1,1,indexm), pmid, tfld, radius_strat, sad_strat, &
                            sad_strat(1,1,2), sad_strat(1,1,4), h2ovmr, vmr, reaction_rates )
#endif
      endif inti_strathet

!-----------------------------------------------------------------------      
!       ...  Set rates for "tabular" and user specified reactions
!-----------------------------------------------------------------------  
      if( gascnt > 0 ) then

#ifdef ARGCHECK
         call setrxt( reaction_rates, tfld, invariants(:,:,indexm), plonl )
         call tropo_rates( reaction_rates, tfld, prho, h2ovmr, invariants, &
                           pmid, invariants(:,:,indexm), lat, kproma )
         call adjrxt( reaction_rates, invariants, invariants(:,:,indexm), plnplv )
#else
! ### remove this branch? Probably obsolete (Stacy)
         call setrxt( reaction_rates, tfld, invariants(1,1,indexm), plonl )
         call tropo_rates( reaction_rates, tfld, prho, h2ovmr, invariants, &
                           pmid, invariants(1,1,indexm), lat, kproma )
         call adjrxt( reaction_rates, invariants, invariants(1,1,indexm), plnplv )
#endif

         IF (ldebug) THEN
            WRITE(message_text,'(a)') 'after setrxt, tropo_rates, adjrxt'
            CALL message('moz_chemdr', message_text, level=em_debug)
         ENDIF
      end if

!++mgs: HAMMOZ hetwet coupling: used to be call to het_n2o5_no3, now included in tropo_rates

!--- Set MOZART field for het_rates to zero (field kept for compatibility)
      het_rates(:,:,:) = 0.0_dp
      IF (ltimer) CALL timer_stop( timer_moz_chem_prep )


!-----------------------------------------------------------------------
!        ... Compute the photolysis rates at time = t(n+1)
!        Note that some accounting and diagnostics for photolysis frequencies is handled 
!        in the general mo_moz_photo module.
!-----------------------------------------------------------------------      

      IF (lfirst) THEN
        IF (phtcnt <= 0) THEN
          lphotolysis = .FALSE.
          CALL message('moz_chemdr', 'lphotolysis set to .FALSE. because phtcnt == 0!', &
                       level=em_warn)
        END IF
      END IF

      IF (ltimer) CALL timer_start( timer_moz_chem_photo )
has_photolysis : &
      if( lphotolysis ) then

!--- Set the column densities at the upper boundary
        if( ncol_abs > 0 .and. phtcnt > 0 ) then
          IF (ldebug) THEN
             WRITE(message_text,'(a)') 'before call to set_ub_col'
             CALL message('moz_chemdr', message_text, level=em_debug)
          ENDIF
          call set_ub_col( col_delta, vmr, pdel, pmid )
        end if

!--- Calculate chemical heating  (HAMMONIA)
!hs   if( lchemheat ) then
!hs      call CHEMHEAT( vmr, reaction_rates, amu, cp, cheat)
!hs   endif
! --- end of HAMMONIA code ---
      
!!!!     reaction_rates(:,:,1:phtcnt) = 0._dp      !### for debugging...

        esfact = flx_ratio_cur        ! correction factor for solar constant (mo_radiation_parameters)

!--- Compute the photolysis rates with WACMM photo module
!--- Set the column densities
        if( ncol_abs > 0 ) then
          call setcol( col_delta, col_dens, vmr, pdel )
        end if

!--- Calculate the photodissociation and heating rates
        rheat_o2(:,:)=0._dp
        rheat_o3(:,:)=0._dp

!--- Call WACMM photo routine
!++mgs 20130228: added arguments for fastj
          call waccm_photo( lat,ip,caldayn,reaction_rates, pmid, pfull, pdel, tfld, &
                            tsurf, zmid, col_dens, zen_angle, albs, &
                            cwat, cice, cldfr, cdnc, sunon, sunoff, esfact, &
                            vmr, invariants, amu, cp, rheat_o2, rheat_o3)
!--mgs

        IF (ndx_o2 > 0) THEN
          rheat_o2(:,:)=rheat_o2(:,:)*vmr(:,:,ndx_o2)
        ENDIF
        IF (ndx_o3 > 0) THEN
          rheat_o3(:,:)=rheat_o3(:,:)*vmr(:,:,ndx_o3)
        ENDIF
        rheat(:,:)=rheat_o2(:,:)+rheat_o3(:,:)

!--- Adjust the photodissociation rates if O2 is an invariant species (see mozpp)
        call phtadj( reaction_rates, invariants, invariants(:,:,indexm), plnplv )

! Photolysis rate diagnostics
        DO jt=1,phtcnt
           photo_mem_p => photo_mem(jt)% ptr
           photo_mem_p(:,:,lat) = reaction_rates(:,:,jt)
        END DO

      ELSE has_photolysis

        DO itest=1, phtcnt
          reaction_rates(:,:,itest) = 1.e-9_dp        ! set to dummy value
        END DO

      END IF has_photolysis
      IF (ltimer) CALL timer_stop( timer_moz_chem_photo )

!--- Make sure all reaction rates are properly defined
      IF (ldebug .AND. p_parallel_io) THEN
         IF ( any(reaction_rates < 0._dp) ) THEN
            CALL message('moz_chemdr', 'Undefined reaction rates!', &
                         level=em_error)
            CALL message('', 'The following reaction rates were undefined:', &
                         level=em_error)
            DO itest=1,rxntot
              IF ( minval(reaction_rates(:,:,itest)) < 0._dp ) THEN
                WRITE(message_text,'(a,i4,2a)') '    reaction ', itest, ': ', &
                      TRIM(rxt_tag_lst(itest))
                CALL message('', message_text, level=em_error)
              END IF
            END DO
            CALL finish('moz_chemdr', 'Undefined reaction rates detected!' )    
         END IF
      END IF

!-----------------------------------------------------------------------
!        ... Compute the extraneous forcing at time = t(n+1)
!-----------------------------------------------------------------------      
! ECHAM-HAMMOZ doesn't use the extfrc fields and setext routine. Keep field for now because of HAMMONIA
      extfrc(:,:,:) = 0.0_dp
      if( extcnt > 0 ) then
         IF (ldebug) THEN
            WRITE(message_text,'(a)') 'Detected extcnt > 0. Note that setext routine has been removed.'
            CALL finish('moz_chemdr', message_text)
         ENDIF
!        call setext( pgeom1, extfrc, lat )
      end if
!++hs: upper atmosphere boundary conditions for HAMMONIA
!     if (lhammonia) THEN
!IF (ldebug .AND. p_parallel_io) write(0,*) '#debug# mo_moz_chemdr: before call to noyion, noycosmis and setext_moz3'
!        call noyion( zmid, lat, zen_angle, pio_hv, pio_e )
!        call noycosmic( lat, invariants(:,:,indexm), pio_c )
!        call setext_moz3( pgeom1, extfrc, lat, ip, zint, cldtop, pio_hv, pio_e, pio_c )
!     endif
!     do m = 1,extcnt
!        do k = 1,plev
!           extfrc(:,k,m) = extfrc(:,k,m) / invariants(:,k,indexm)
!        end do
!     end do
!--hs


!=======================================================================
!        ... Call the class solution algorithms
! This version only supports the implicit solver
!=======================================================================
has_chemsolv:  &
  if ( lchemsolv ) then
      IF (ldebug) THEN
         WRITE(message_text,'(a)') 'before call to imp_sol'
         CALL message('moz_chemdr', message_text, level=em_debug)
      ENDIF

!!### DEBUGGING
!!    IF (ldebug .AND. lat == 3) THEN
!!       DO i = 1, pcnstm1
!!          WRITE(message_text,'(a,i0,a,1p5g15.7)')         &
!!                             'i = ', i, ' vmr', vmr(5,34,i)
!!          CALL message('moz_chemdr', message_text, level=em_debug)
!!       ENDDO
!!       DO i = 1,rxntot
!!          WRITE(message_text,'(a,i0,a,1p5g15.7)')         &
!!                             'i = ', i, 'rxt rates', reaction_rates(5,34,i)
!!          CALL message('moz_chemdr', message_text, level=em_debug)
!!       ENDDO
!!    ENDIF
!!### END DEBUGGING

      if( clscnt4 > 0 ) then
!-----------------------------------------------------------------------
!	... Solve for "Implicit" species
!-----------------------------------------------------------------------
            IF (ltimer) CALL timer_start( timer_moz_chem_impsol )
            call imp_sol( vmr, reaction_rates, het_rates, extfrc, &
                          delt, lat, plonl, plnplv )
            IF (ltimer) CALL timer_stop( timer_moz_chem_impsol )
      end if
      IF (ldebug) THEN
         WRITE(message_text,'(a)') 'after call to imp_sol'
         CALL message('moz_chemdr', message_text, level=em_debug)
      ENDIF

  end if has_chemsolv


!-----------------------------------------------------------------------
!         ... aerosol settling
!-----------------------------------------------------------------------
has_strathet:  &
      if (lstrathet) then
         IF (ldebug) THEN
            WRITE(message_text,'(a)') 'before call to strat_aer_settling'
            CALL message('moz_chemdr', message_text, level=em_debug)
         ENDIF

#ifdef ARGCHECK
        call strat_aer_settling( invariants(:,:,indexm), pmid, delt, zmid, tfld, &
                                hno3_cond(:,:,3), radius_strat(:,:,3), &
                                plonl, lat )
#else
        call strat_aer_settling( invariants(1,1,indexm), pmid, delt, zmid, tfld, &
                                hno3_cond(1,1,3), radius_strat(1,1,3), &
                                plonl, lat )
#endif
!       add condensed and gas phase hno3 to transport them together
        if( ndx_hno3 > 0 ) then
           do k = sad_top+1,plev
              vmr(:,k,ndx_hno3) = vmr(:,k,ndx_hno3) + hno3_cond(:,k,1) + &
                                 hno3_cond(:,k,2) + hno3_cond(:,k,3)
           end do
        end if
      endif has_strathet

!-----------------------------------------------------------------------      
!         ... Diagnostics and Cleanup
!-----------------------------------------------------------------------      

!--- Check for negative values and reset to zero
! ### TODO: test if necessary (should actually be obsolete now)
      call negtrc( lat, 'After chemistry ', vmr )

!--- Xform from vmr to mmr
      call vmr2mmr( vmr, mmr, mbar )
!--- Store updated H2O mixing ratio in sh field in order to return it to ECHAM 
      sh(:,:) = mmr(:,:,ndx_h2o)

!--- Diagnose OH concentration and lifetimes (for global loss rates and output)
!    Conversion constant: NA / R * 1.e-6 from ideal gas law
!++mgs20140319: bug fix - conversion constant was 1.233468e17_dp which is only correct for MMR
      IF (ndx_oh > 0 .AND. associated(dpohconc)) THEN
         dpohconc(1:kproma,1:plev,lat) = 7.242963582053075e+16                  &
                                         * pmid(1:kproma,1:plev)                &
                                         * vmr(1:kproma,1:plev,ndx_oh)          &
                                         / tfld(1:kproma,1:plev)
!--mgs
!++mgs20140319: moved these diagnostics here from moz_chemistry
         IF (associated(dpratech4)) THEN
            ! Compute rate coefficient manually here in case we would run without CH4 (although this doesn't make much sense)
            zkch4(1:kproma, 1:plev) = 2.45e-12_dp * exp(-1775._dp / tfld(1:kproma, 1:plev))  ! JPL17_10-6 (2011)
            dpratech4(1:kproma, 1:plev, lat) = zkch4(1:kproma, 1:plev)*dpohconc(1:kproma,1:plev,lat)*secperyear
         END IF
         IF (associated(dpratemcl)) THEN
            ! Compute rate coefficient manually here in case we would run without CH3CCl3
            zkmcl(1:kproma, 1:plev) = 1.64e-12_dp * exp(-1520._dp / tfld(1:kproma, 1:plev))  ! JPL17_10-6 (2011)
            dpratemcl(1:kproma, 1:plev, lat) = zkmcl(1:kproma, 1:plev)*dpohconc(1:kproma,1:plev,lat)*secperyear
         END IF
!--mgs
      END IF

!++mgs20140319: added stratospheric aerosol diagnostics
      ! Stratospheric aerosol
      ! (calculation of sad_strat and hno3_cond is only done if lstrathet = .true.)
      IF (lstrathet) THEN
         IF (associated(dpsad1)) THEN
            dpsad1(1:kproma, 1:plev, lat) = sad_strat(1:kproma, 1:plev, 1)    ! sulfate
         END IF
         IF (associated(dpsad2)) THEN
            dpsad2(1:kproma, 1:plev, lat) = sad_strat(1:kproma, 1:plev, 2) + sad_strat(1:kproma, 1:plev, 3)    ! small and large NAT
         END IF
         IF (associated(dpsad3)) THEN
            dpsad3(1:kproma, 1:plev, lat) = sad_strat(1:kproma, 1:plev, 4)    ! ice
         END IF
         IF (associated(dphno3_cond)) THEN
            dphno3_cond(1:kproma, 1:plev, lat) = hno3_cond(1:kproma, 1:plev, 1)    &
                                               + hno3_cond(1:kproma, 1:plev, 2)    &
                                               + hno3_cond(1:kproma, 1:plev, 3)
         END IF
      ELSE
         IF (associated(dpsad1)) THEN
            dpsad1(1:kproma, 1:plev, lat) = 0._dp
         END IF
         IF (associated(dpsad2)) THEN
            dpsad2(1:kproma, 1:plev, lat) = 0._dp
         END IF
         IF (associated(dpsad3)) THEN
            dpsad3(1:kproma, 1:plev, lat) = 0._dp
         END IF
         IF (associated(dphno3_cond)) THEN
            dphno3_cond(1:kproma, 1:plev, lat) = 0._dp
         END IF
      END IF
!--mgs

      CALL moz_family_diag(kproma, plonl, plev, lat, pcnstm1, vmr)
      CALL mozh_diag(kproma, plonl, plev, lat, pcnstm1, vmr) ! s.stadtler
      CALL moz_rates_diag(kproma, plonl, plev, lat, pcnstm1, rxntot, vmr, reaction_rates)

      IF (ldebug) THEN
        write(message_text,*) 'Done with subroutine on CPU ', p_pe
        CALL message('moz_chemdr', message_text, level=em_debug)
      END IF

      lfirst = .FALSE.

      end subroutine chemdr

end module mo_moz_chemdr
