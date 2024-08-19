!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!!  mo_moz_het_strato_rates
!!
!! \brief
!!  Derivation of the rate constant for reactions on sulfate, NAD, and ICE aerosols.
!!
!! \author Douglas E. Kinnison (NCAR)
!!
!! \responsible_coder
!!  Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!  - DK: original version (2002-08-15)
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


      module mo_moz_het_strato_rates
!=======================================================================
! ROUTINE
!   ratecon_sfstrat.f
!
! DESCRIPTION
!
! Sulfate Aerosol Reactions               Rxn#   Gamma
!   N2O5   + H2O(l)     =>  2HNO3         (1)    f(wt%)
!   ClONO2 + H2O(l)     =>  HOCl + HNO3   (2)    f(T,P,HCl,H2O,r)
!   BrONO2 + H2O(l)     =>  HOBr + HNO3   (3)    f(T,P,H2O,r)
!   ClONO2 + HCl(l)     =>  Cl2  + HNO3   (4)    f(T,P,HCl,H2O,r)
!   HOCl   + HCl(l)     =>  Cl2  + H2O    (5)    f(T,P,HCl,HOCl,H2O,r)
!   HOBr   + HCl(l)     =>  BrCl + H2O    (6)    f(T,P,HCl,HOBr,H2O,r)
!
! Nitric Acid Di-hydrate Reactions        Rxn#    Gamma   Reference
!   N2O5   + H2O(s)     =>  2HNO3         (7)     4e-4   JPL00-03
!   ClONO2 + H2O(s)     =>  HOCl + HNO3   (8)     4e-3   JPL00-03
!   ClONO2 + HCl(s)     =>  Cl2  + HNO3   (9)     0.2    JPL00-03
!   HOCl   + HCl(s)     =>  Cl2  + H2O    (10)    0.1    JPL00-03
!   BrONO2 + H2O(s)     =>  HOBr + HNO3   (11)    0.3    David Hanson PC
!
! ICE Aersol Reactions                    Rxn#    Gamma
!   N2O5   + H2O(s)     =>  2HNO3         (12)     0.02   JPL00-03
!   ClONO2 + H2O(s)     =>  HOCl + HNO3   (13)     0.3    JPL00-03
!   BrONO2 + H2O(s)     =>  HOBr + HNO3   (14)     0.3    JPL00-03
!   ClONO2 + HCl(s)     =>  Cl2  + HNO3   (15)     0.3    JPL00-03
!   HOCl   + HCl(s)     =>  Cl2  + H2O    (16)     0.2    JPL00-03
!   HOBr   + HCl(s)     =>  BrCl + H2O    (17)     0.3    JPL00-03
!
! NOTE: The rate constants derived from species reacting with H2O are
!       first order (i.e., sec-1 units) - an example is N2O5 + H2O = 2HNO3.
!       Other reactions, e.g., ClONO2 + HCl have rate constants that
!       are second order (i.e., cm+3 molecules-1 sec-1 units). In all
!       of these types of reactions the derived first order rate constant
!       {0.25*(mean Velocity)*SAD*gamma} is divided by the HCl abundance
!       to derive the correct second order units.
!
! NOTE: Liquid Sulfate Aerosols...
!       See coding for references on how the Sulfate Aerosols were handled.
!       Data was used that was more recent than JPL00.
!
!
! INPUT:
!  ad      .    .... air density, molec. cm-3
!  pmid        ..... pressures, hPa
!  temp        ..... temperatures, K
!  rad_sulfate ..... Surface area density, cm2 cm-3
!  sad_sulfate ..... Surface area density, cm2 cm-3
!  sad_nad     ..... Surface area density, cm2 cm-3
!  sad_ice     ..... Surface area density, cm2 cm-3
!  brono2mv    ..... BrONO2 Volume Mixing Ratio
!  clono2mvr   ..... ClONO2 Volume Mixing Ratio
!  h2omvr      ..... H2O Volume Mixing Ratio
!  hclmvr      ..... HCl Volume Mixing Ratio
!  hobrmvr     ..... HOBr Volume Mixing Ratio
!  hoclmvr     ..... HOCl Volume Mixing Ratio
!  n2o5mvr     ..... N2O5 Volume Mixing Ratio
!
! OUTPUT:
!
!  rxt         ..... Rate constant (s-1 and cm3 sec-1 molec-1)
!=======================================================================

      USE mo_kind,     ONLY : dp

      implicit none

      private
      public :: ratecon_sfstrat_inti, ratecon_sfstrat, zero_sfstrat

      integer, parameter :: nrid_het = 17
      integer            :: ndx_brono2, ndx_clono2, ndx_hcl, ndx_hocl, &
                            ndx_hobr, ndx_n2o5
      integer            :: rid_het(17)


      contains


      subroutine ratecon_sfstrat_inti ()

      USE mo_exception, ONLY : message, em_error
      USE mo_moz_util,  ONLY : get_spc_ndx, get_rxt_ndx

         ! obtain necessary species ids
         ndx_brono2   = get_spc_ndx('BRONO2')
         ndx_clono2   = get_spc_ndx('CLONO2')
         ndx_hcl      = get_spc_ndx('HCL')
         ndx_hocl     = get_spc_ndx('HOCL')
         ndx_hobr     = get_spc_ndx('HOBR')
         ndx_n2o5     = get_spc_ndx('N2O5')

         IF (ndx_brono2 == -1) CALL message('moz_het_strato_rates', 'No BRONO2. Need this in ratecon_sfstrat!', level=em_error)
         IF (ndx_clono2 == -1) CALL message('moz_het_strato_rates', 'No CLONO2. Need this in ratecon_sfstrat!', level=em_error)
         IF (ndx_hcl == -1) CALL message('moz_het_strato_rates', 'No HCL. Need this in ratecon_sfstrat!', level=em_error)
         IF (ndx_hocl == -1) CALL message('moz_het_strato_rates', 'No HOCL. Need this in ratecon_sfstrat!', level=em_error)
         IF (ndx_hobr == -1) CALL message('moz_het_strato_rates', 'No HOBR. Need this in ratecon_sfstrat!', level=em_error)
         IF (ndx_n2o5 == -1) CALL message('moz_het_strato_rates', 'No N2O5. Need this in ratecon_sfstrat!', level=em_error)

         ! obtain necessary reaction indices
         rid_het( 1)   = get_rxt_ndx('het_N2O5_strat_sulfate')
         rid_het( 2)   = get_rxt_ndx('het_CLONO2_strat_sulfate')
         rid_het( 3)   = get_rxt_ndx('het_BRONO2_strat_sulfate')
         rid_het( 4)   = get_rxt_ndx('het_CLONO2_HCL_strat_sulfate')
         rid_het( 5)   = get_rxt_ndx('het_HOCL_HCL_strat_sulfate')
         rid_het( 6)   = get_rxt_ndx('het_HOBR_HCL_strat_sulfate')
         rid_het( 7)   = get_rxt_ndx('het_N2O5_strat_NAD')
         rid_het( 8)   = get_rxt_ndx('het_CLONO2_strat_NAD')
         rid_het( 9)   = get_rxt_ndx('het_CLONO2_HCL_strat_NAD')
         rid_het(10)   = get_rxt_ndx('het_HOCL_HCL_strat_NAD')
         rid_het(11)   = get_rxt_ndx('het_BRONO2_strat_NAD')
         rid_het(12)   = get_rxt_ndx('het_N2O5_strat_ice')
         rid_het(13)   = get_rxt_ndx('het_CLONO2_strat_ice')
         rid_het(14)   = get_rxt_ndx('het_BRONO2_strat_ice')
         rid_het(15)   = get_rxt_ndx('het_CLONO2_HCL_strat_ice')
         rid_het(16)   = get_rxt_ndx('het_HOCL_HCL_strat_ice')
         rid_het(17)   = get_rxt_ndx('het_HOBR_HCL_strat_ice')

         IF (ANY(rid_het(:) == -1)) CALL message('moz_het_strato_rates',    &
                                                'One of rid_het1..rid_het17 is 0. Need these for ratecon_sfstrat!',  &
                                                level=em_error)

      end subroutine ratecon_sfstrat_inti


      subroutine ratecon_sfstrat( ad, pmid, temp, rad_sulfate, sad_sulfate, &
                                  sad_nad, sad_ice, h2ovmr, vmr, rxt )

!++mgs: replaced mo_grid and chem_mods with mo_moz_mods
      use mo_moz_mods, only : ncol => plonl, pcols => plonl, pver => plev
      use mo_moz_mods, only : adv_mass, rxntot, gas_pcnst => pcnstm1
!--mgs
      use mo_moz_strato_sad, only : sad_top


!-----------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------
!!    integer, intent(in) :: ncol                                 ! columns in chunk
      real(dp), dimension(ncol,pver,gas_pcnst), intent(in) :: &   ! species concentrations (mol/mol)
        vmr
      real(dp), dimension(ncol,pver), intent(in) :: &
        ad, &                                                     ! Air Density (molec. cm-3)
        rad_sulfate, &                                            ! Radius of Sulfate Aerosol (cm)
        sad_ice, &                                                ! ICE Surface Area Density (cm-1)
        sad_nad, &                                                ! NAD Surface Area Density (cm-1)
        sad_sulfate, &                                            ! Sulfate Surface Area Density (cm-1)
        h2ovmr                                                    ! water vapor volume mixing ratio( gas phase )
      real(dp), dimension(pcols,pver), intent(in) :: &
        pmid, &                                                   ! pressure (Pa)
        temp                                                      ! temperature (K)

      real(dp), intent(inout) :: &
        rxt(ncol,pver,rxntot)                                     ! rate constants

!-----------------------------------------------------------------------
!  	... local variables
!-----------------------------------------------------------------------
	real(dp), parameter :: small_conc = 1.e-30_dp
	real(dp), parameter :: av_const   = 2.117265e4_dp  ! (8*8.31448*1000 / PI)
	real(dp), parameter :: pa2mb      = 1.e-2_dp       ! Pa to mb
	real(dp), parameter :: m2cm       = 100._dp        ! meters to cms

      integer :: &
        i, &                      ! altitude loop index
        k, &                      ! level loop index
        m                         ! species index

!-----------------------------------------------------------------------
!   	... variables for gamma calculations
!-----------------------------------------------------------------------
      real(dp) :: &
        brono2vmr, &                            ! BrONO2 Volume Mixing Ratio
        clono2vmr, &                            ! ClONO2 Volume Mixing Ratio
        hclvmr, &                               ! HCl Volume Mixing Ratio
        hcldeni, &                              ! inverse of HCl density
        hoclvmr, &                              ! HOCl Volume Mixing Ratio
        hobrvmr, &                              ! HOBr Volume Mixing Ratio
        n2o5vmr                                 ! N2O5 Volume Mixing Ratio

        real(dp) :: &
        av_n2o5, &                              ! N2O5 Mean Velocity (cm s-1)
        av_clono2, &                            ! ClONO2 Mean Velocity (cm s-1)
        av_brono2, &                            ! BrONO2Mean Velocity (cm s-1)
        av_hocl, &                              ! HOCl Mean Velocity (cm s-1)
        av_hobr                                 ! HOBr Mean Velocity (cm s-1)

      real(dp) :: &
        pzero_h2o, &                            ! H2O sat vapor press (mbar)
        e0, e1, e2, e3, &                       ! coefficients for H2O sat vapor press.
        aw, &                                   ! Water activity
        m_h2so4, &                              ! H2SO4 molality (mol/kg)
        wt, &                                   ! wt % H2SO4
        y1, y2, &                               ! used in H2SO4 molality
        &  a1, b1, c1, d1, &                    ! used in H2SO4 molality
        a2, b2, c2, d2                          ! used in H2SO4 molality

        real(dp) :: &
        z1, z2, z3, &                           ! used in H2SO4 soln density
        den_h2so4, &                            ! H2SO4 soln density, g/cm3
        mol_h2so4, &                            ! Molality of H2SO4, mol / kg
        molar_h2so4, &                          ! Molarity of H2SO4, mol / l
        x_h2so4, &                              ! H2SO4 mole fraction
        aconst, tzero, &                        ! used in viscosity of H2SO4
        vis_h2so4, &                            ! H2SO4 viscosity
        ah, &                                   ! Acid activity, molarity units
        term1,term2,term3,term4, &              ! used in ah
        term5,term6,term7,term0, &
        T_limit, &                              ! temporary variable for temp (185-260K range)
        T_limiti, &                             ! 1./T_limit
        T_limitsq, &                            ! sqrt( T_limit )
        rad_sulf, &                             ! temporary variable for sulfate radius (cm)
        rad_sad_sulf, &                         ! temporary variable for sulfate radius (cm)
        sadice, &                               ! temporary variable for sulfate radius (cm)
        rad_nad_sulf                            ! temporary variable for sulfate radius (cm)

      real(dp) :: &
        C_cnt, S_cnt, &                         ! used in H_cnt
        H_cnt, &                                ! Henry's law coeff. for ClONO2
        H_hcl, &                                ! Henry's law coeff. for HCl
        D_cnt, &
        k_hydr, &
        k_h2o, &
        k_h, &
        k_hcl, &
        rdl_cnt, &
        f_cnt, &
        M_hcl, &
        atmos

      real(dp) :: &
        Gamma_b_h2o, &
        Gamma_cnt_rxn, &
        Gamma_b_hcl, &
        Gamma_s, &
        Fhcl, &
        Gamma_s_prime, &
        Gamma_b_hcl_prime, &
        Gamma_b, &
        gprob_n2o5, &
        gprob_rxn, &
        gprob_tot, &
        gprob_cnt, &
        gprob_cnt_hcl, &
        gprob_cnt_h2o

        real(dp) :: &
        D_hocl, &
        k_hocl_hcl, &
        C_hocl, &
        S_hocl, &
        H_hocl, &
        Gamma_hocl_rxn, &
        rdl_hocl, &
        f_hocl, &
        gprob_hocl_hcl

        real(dp) :: &
        h1, h2, h3, &
        alpha, &
        gprob_bnt_h2o

      real(dp) :: &
        C_hobr, &
        D_hobr, &
        aa, bb, cc, dd, &
        k_hobr_hcl, &
        k_dl, &
        k_wasch, &
        H_hobr, &
        rdl_hobr, &
        Gamma_hobr_rxn, &
        f_hobr, &
        gprob_hobr_hcl

      real(dp) :: &
        pmb,&					! Pressure, mbar (hPa)
        pH2O_atm,&				! Partial press. H2O (atm)
        pH2O_hPa,&				! Partial press. H2O (hPa)
        pHCl_atm,&				! Partial press. HCl (atm)
        pCNT_atm                                ! Partial press. ClONO2 (atm)

!-----------------------------------------------------------------------
!     	... Used in pzero h2o calculation
!-----------------------------------------------------------------------
      real(dp), parameter :: wt_e0 = 18.452406985_dp
      real(dp), parameter :: wt_e1 = 3505.1578807_dp
      real(dp), parameter :: wt_e2 = 330918.55082_dp
      real(dp), parameter :: wt_e3 = 12725068.262_dp

      real(dp) :: &
        wrk, tmp
      REAL(dp):: zeps

      zeps=EPSILON(1._dp)
!-----------------------------------------------------------------------
!     	... set rate constants
!-----------------------------------------------------------------------
Level_loop : &
      do k = sad_top+1,pver
column_loop : &
         do i = 1,ncol
!-----------------------------------------------------------------------
!	... set species, pmb, and atmos
!-----------------------------------------------------------------------
	    brono2vmr = vmr(i,k,ndx_brono2)
	    clono2vmr = vmr(i,k,ndx_clono2)
	    hclvmr    = vmr(i,k,ndx_hcl)
!           if( hclvmr > 0._dp ) then
!              hcldeni  = 1._dp/(hclvmr*ad(i,k))
!           end if
            hcldeni      = 1._dp/(MAX(1.e-30_dp,hclvmr)*ad(i,k)) ! changed by H. Schmidt, MPI, to enable vectorization
	    hoclvmr      = vmr(i,k,ndx_hocl)
	    hobrvmr      = vmr(i,k,ndx_hobr)
	    n2o5vmr      = vmr(i,k,ndx_n2o5)
            rad_sad_sulf = sad_sulfate(i,k)
            rad_nad_sulf = sad_nad(i,k)
	    sadice       = sad_ice(i,k)
            pmb          = pa2mb*pmid(i,k)
            atmos        = pmb/1013.25_dp

!-----------------------------------------------------------------------
!  	... setup for stratospheric aerosols
!           data range set: 185K - 240K;    GRL, 24, 1931, 1997
!-----------------------------------------------------------------------
            T_limit   = max( temp(i,k),185._dp )
            T_limit   = min( T_limit,240._dp )
            T_limiti  = 1._dp/T_limit
            T_limitsq = sqrt( T_limit )

!-----------------------------------------------------------------------
!     .... Average velocity (8RT*1000/(PI*MW))**1/2 * 100.(units cm s-1)
!     .... or (av_const*T/M2)**1/2
!-----------------------------------------------------------------------
	    wrk       = av_const*T_limit
            av_n2o5   = sqrt( wrk/adv_mass(ndx_n2o5) )*m2cm
            av_clono2 = sqrt( wrk/adv_mass(ndx_clono2) )*m2cm
            av_brono2 = sqrt( wrk/adv_mass(ndx_brono2) )*m2cm
            av_hocl   = sqrt( wrk/adv_mass(ndx_hocl) )*m2cm
            av_hobr   = sqrt( wrk/adv_mass(ndx_hobr) )*m2cm
has_rad_sad_sulf : &
            if( rad_sad_sulf > 0._dp ) then
!-----------------------------------------------------------------------
!     .... Partial Pressure of H2O, ClONO2, and HCl in atmospheres
!-----------------------------------------------------------------------

!+++sschr: just to test
!write(*,*) "moz_het_strato_rates: reaction_rates set for ", i, k 
!---sschr 

               if( hclvmr > 0._dp ) then
                  pHCl_atm  = hclvmr*atmos
               else
                  pHCl_atm  = 0._dp
               end if

               if( clono2vmr > 0._dp ) then
                  pCNT_atm  = clono2vmr*atmos
               else
                  pCNT_atm  = 0._dp
               end if

               if( h2ovmr(i,k) > 0._dp ) then
                  pH2O_atm  = h2ovmr(i,k)*atmos
               else
                  pH2O_atm  = 0._dp
               end if

!-----------------------------------------------------------------------
!     .... Partial Pressure of H2O in hPa
!-----------------------------------------------------------------------
               pH2O_hpa = h2ovmr(i,k)*pmb
!-----------------------------------------------------------------------
!     .... Calculate the h2so4 Wt% and Activity of H2O - JPL-00
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     ... Saturation Water Vapor Pressure (mbar)
!-----------------------------------------------------------------------
               pzero_h2o = exp( wt_e0 - T_limiti*(wt_e1 + T_limiti*(wt_e2 - T_limiti*wt_e3)) )
!-----------------------------------------------------------------------
!     ... H2O activity
!     ... if the activity of H2O goes above 1.0, wt% can go negative
!-----------------------------------------------------------------------
               aw = ph2o_hpa / pzero_h2o
               aw = min( aw,1._dp )
               aw = max( aw,.001_dp )
!-----------------------------------------------------------------------
!     ... h2so4 Molality (mol/kg)
!-----------------------------------------------------------------------
               if( aw <= .05_dp ) then
                  a1 = 12.37208932_dp
                  b1 = -0.16125516114_dp
                  c1 = -30.490657554_dp
                  d1 = -2.1133114241_dp
                  a2 = 13.455394705_dp
                  b2 = -0.1921312255_dp
                  c2 = -34.285174607_dp
                  d2 = -1.7620073078_dp
               else if( aw > .05_dp .and. aw < .85_dp ) then
                  a1 = 11.820654354_dp
                  b1 = -0.20786404244_dp
                  c1 = -4.807306373_dp
                  d1 = -5.1727540348_dp
                  a2 = 12.891938068_dp
                  b2 = -0.23233847708_dp
                  c2 = -6.4261237757_dp
                  d2 = -4.9005471319_dp
               else
                  a1 = -180.06541028_dp
                  b1 = -0.38601102592_dp
                  c1 = -93.317846778_dp
                  d1 = 273.88132245_dp
                  a2 = -176.95814097_dp
                  b2 = -0.36257048154_dp
                  c2 = -90.469744201_dp
                  d2 = 267.45509988_dp
               end if
!-----------------------------------------------------------------------
!     ... h2so4 mole fraction
!-----------------------------------------------------------------------
               y1       = a1*(aw**b1) + c1*aw + d1
               y2       = a2*(aw**b2) + c2*aw + d2
               m_h2so4  = y1 + ((T_limit - 190._dp)*(y2 - y1)) / 70._dp
!-----------------------------------------------------------------------
!     ... h2so4 Weight Percent
!-----------------------------------------------------------------------
               wt = 9800._dp*m_h2so4 / (98._dp*m_h2so4  + 1000._dp)
!-----------------------------------------------------------------------
!     .... Parameters for h2so4 Solution, JPL-00
!-----------------------------------------------------------------------
!     ... h2so4 Solution Density (g/cm3)
!-----------------------------------------------------------------------
	       wrk = T_limit*T_limit
               z1 =  .12364_dp  - 5.6e-7_dp*wrk
               z2 = -.02954_dp  + 1.814e-7_dp*wrk
               z3 =  2.343e-3_dp - T_limit*1.487e-6_dp - 1.324e-8_dp*wrk
!-----------------------------------------------------------------------
!     ... where mol_h2so4 is molality in mol/kg
!-----------------------------------------------------------------------
               den_h2so4 = 1._dp + m_h2so4*(z1 + z2*sqrt(m_h2so4) + z3*m_h2so4)
!-----------------------------------------------------------------------
!     ... h2so4 Molarity, mol / l
!-----------------------------------------------------------------------
               molar_h2so4 = den_h2so4*wt/9.8_dp
!-----------------------------------------------------------------------
!     ... h2so4 Mole fraction
!-----------------------------------------------------------------------
               x_h2so4   = wt / (wt + ((100._dp - wt)*98._dp/18._dp))
               term1     = .094_dp - x_h2so4*(.61_dp - 1.2_dp*x_h2so4)
               term2     = (8515._dp - 10718._dp*(x_h2so4**.7_dp))*T_limiti
               H_hcl     = term1 * exp( -8.68_dp + term2 )
               M_hcl     = H_hcl*pHCl_atm

!-----------------------------------------------------------------------
!     ... h2so4 solution viscosity
!-----------------------------------------------------------------------
               aconst    = 169.5_dp + wt*(5.18_dp - wt*(.0825_dp - 3.27e-3_dp*wt))
               tzero     = 144.11_dp + wt*(.166_dp - wt*(.015_dp - 2.18e-4_dp*wt))
               vis_h2so4 = aconst/(T_limit**1.43_dp) * exp( 448._dp/(T_limit - tzero) )

!-----------------------------------------------------------------------
!     ... Acid activity in molarity
!-----------------------------------------------------------------------
               term1 = 60.51_dp
               term2 = .095_dp*wt
	       wrk   = wt*wt
               term3 = .0077_dp*wrk
               term4 = 1.61e-5_dp*wt*wrk
               term5 = (1.76_dp + 2.52e-4_dp*wrk) * T_limitsq
               term6 = -805.89_dp + (253.05_dp*(wt**.076_dp))
               term7 = T_limitsq
               ah    = exp( term1 - term2 + term3 - term4 - term5 + term6/term7 )
!              removed by H.Schmidt, MPI, to enable vectorization
!	       if( ah <= 0._dp ) then
!	          write(*,*) 'ratecon: ah <= 0 at i,k, = ',i,k
!	          write(*,*) 'ratecon: term1,term2,term3,term4,term5,term6,term7,wt,T_limit,ah = ', &
!	                               term1,term2,term3,term4,term5,term6,term7,wt,T_limit,ah 
!	       end if

	       wrk      = .25_dp*rad_sad_sulf
               rad_sulf = max( rad_sulfate(i,k),1.e-7_dp )
!-----------------------------------------------------------------------
!     N2O5 + H2O(liq) =>  2.00*HNO3  Sulfate Aerosol Reaction
!-----------------------------------------------------------------------
               if( n2o5vmr > 0._dp ) then
                  term0 = -25.5265_dp - wt*(.133188_dp - wt*(.00930846_dp - 9.0194e-5_dp*wt))
                  term1 = 9283.76_dp + wt*(115.345_dp - wt*(5.19258_dp - .0483464_dp*wt))
                  term2 = -851801._dp - wt*(22191.2_dp - wt*(766.916_dp - 6.85427_dp*wt))
                  gprob_n2o5 = exp( term0 + T_limiti*(term1 + term2*T_limiti) )
                  rxt(i,k,rid_het(1)) = max( 0._dp,wrk*av_n2o5*gprob_n2o5 )
               end if

!-----------------------------------------------------------------------
!     ClONO2 + H2O(liq) =  HOCl + HNO3   Sulfate Aerosol Reaction
!-----------------------------------------------------------------------
! 	... NOTE: Aerosol radius in units of cm.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     	... Radius sulfate set (from sad module)
!           Set min radius to 0.001 microns (1e-7 cm)
!           Typical radius is 0.1 microns (1e-5 cm)
!           f_cnt may go negative under if not set.
!-----------------------------------------------------------------------
                  C_cnt         = 1474._dp*T_limitsq
                  S_cnt         = .306_dp + 24._dp*T_limiti
                  term1         = exp( -S_cnt*molar_h2so4 )
                  H_cnt         = 1.6e-6_dp * exp( 4710._dp*T_limiti )*term1
                  D_cnt         = 5.e-8_dp*T_limit / vis_h2so4
                  k_h           = 1.22e12_dp*exp( -6200._dp*T_limiti )
                  k_h2o         = 1.95e10_dp*exp( -2800._dp*T_limiti )
                  k_hydr        = (k_h2o + k_h*ah)*aw
                  k_hcl         = 7.9e11_dp*ah*D_cnt*M_hcl
                  rdl_cnt       = sqrt( D_cnt/(k_hydr + k_hcl) )
                  term1         = 1._dp/tanh( rad_sulf/rdl_cnt )
                  term2         = rdl_cnt/rad_sulf
                  f_cnt         = term1 - term2
                  if( f_cnt > 0._dp ) then
                     term1         = 4._dp*H_cnt*.082_dp*T_limit
                     term2         = sqrt( D_cnt*k_hydr )
                     Gamma_b_h2o   = term1*term2/C_cnt
                     term1         = sqrt( 1._dp + k_hcl/k_hydr )
                     Gamma_cnt_rxn = f_cnt*Gamma_b_h2o*term1
                     Gamma_b_hcl   = Gamma_cnt_rxn*k_hcl/(k_hcl + k_hydr)
                     term1         = exp( -1374._dp*T_limiti )
                     Gamma_s       = 66.12_dp*H_cnt*M_hcl*term1
		     if( pHCl_atm > 0._dp ) then
                        term1      = .612_dp*(Gamma_s*Gamma_b_hcl)* pCNT_atm/pHCl_atm
                        Fhcl       = 1._dp/(1._dp + term1)
		     else
                        Fhcl       = 1._dp
		     end if
                     Gamma_s_prime     = Fhcl*Gamma_s
                     Gamma_b_hcl_prime = Fhcl*Gamma_b_hcl
                     term1         = Gamma_cnt_rxn*k_hydr
                     term2         = k_hcl + k_hydr
                     Gamma_b       = Gamma_b_hcl_prime + (term1/term2)
                     term1         = 1._dp / (Gamma_s_prime + Gamma_b)
                     gprob_cnt     = 1._dp / (1._dp + term1)
                     term1         = Gamma_s_prime + Gamma_b_hcl_prime
                     term2         = Gamma_s_prime + Gamma_b
                     gprob_cnt_hcl = gprob_cnt * term1/term2
                     gprob_cnt_h2o = gprob_cnt - gprob_cnt_hcl
                  else
                     gprob_cnt_h2o = 0._dp
                     gprob_cnt_hcl = 0._dp
                     Fhcl          = 1._dp
                  end if
                  if( clono2vmr > 0._dp ) then
                     rxt(i,k,rid_het(2)) = max( 0._dp,wrk*av_clono2*gprob_cnt_h2o )
                  end if

!-----------------------------------------------------------------------
!  	... BrONO2 + H2O(liq) =  HOBr + HNO3   Sulfate Aerosol Reaction
!-----------------------------------------------------------------------
               if( brono2vmr > 0._dp ) then
                  h1    = 29.24_dp
                  h2    = -.396_dp
                  h3    = .114_dp
                  alpha = .805_dp
                  gprob_rxn = exp( h1 + h2*wt ) + h3
                  term1     = 1._dp/alpha
                  term2     = 1._dp/gprob_rxn
                  gprob_bnt_h2o = 1._dp / (term1 + term2)
                  rxt(i,k,rid_het(3)) = max( 0._dp,wrk*av_brono2*gprob_bnt_h2o )
               end if

!-----------------------------------------------------------------------
!     	... ClONO2 + HCl(liq) =  Cl2  + HNO3  Sulfate Aerosol Reaction
!-----------------------------------------------------------------------
               if( hclvmr > 0._dp .and. clono2vmr > 0._dp ) then
                  rxt(i,k,rid_het(4)) = max( 0._dp,wrk*av_clono2*gprob_cnt_hcl )*hcldeni
               end if

!-----------------------------------------------------------------------
!     	... HOCl + HCl(liq) =  Cl2 + H2O   Sulfate Aerosol Reaction
!-----------------------------------------------------------------------
               if( hclvmr > 0._dp .and. hoclvmr > 0._dp ) then
!-----------------------------------------------------------------------
!     	... Radius sulfate set (from sad module)
!           Set min radius to 0.001 microns (1e-7 cm)
!           Typical radius is 0.1 microns (1e-5 cm)
!           f_hocl may go negative under if not set.
!-----------------------------------------------------------------------
	          if( pCNT_atm > 0._dp ) then
                     D_hocl          = 6.4e-8_dp*T_limit/vis_h2so4
                     k_hocl_hcl      = 1.25e9_dp*ah*D_hocl*M_hcl
                     C_hocl          = 2009._dp*T_limitsq
                     S_hocl          = .0776_dp + 59.18_dp*T_limiti
                     term1           = exp( -S_hocl*molar_h2so4 )
                     H_hocl          = 1.91e-6_dp * exp( 5862.4_dp*T_limiti )*term1
                     term1           = 4._dp*H_hocl*.082_dp*T_limit
                     term2           = sqrt( D_hocl*k_hocl_hcl )
                     Gamma_hocl_rxn  = term1*term2/C_hocl
                     rdl_hocl        = sqrt( D_hocl/k_hocl_hcl )
                     term1           = 1._dp/tanh( rad_sulf/rdl_hocl )
                     term2           = rdl_hocl/rad_sulf
                     f_hocl          = term1 - term2
                     if( f_hocl > 0._dp ) then
                        term1           = 1._dp / (f_hocl*Gamma_hocl_rxn*Fhcl)
                        gprob_hocl_hcl  = 1._dp / (1._dp + term1)
                     else
                        gprob_hocl_hcl  = 0._dp
                     end if
                     rxt(i,k,rid_het(5)) = max( 0._dp,wrk*av_hocl*gprob_hocl_hcl )*hcldeni
	          end if
               end if

!-----------------------------------------------------------------------
!     	... HOBr + HCl(liq) =  BrCl + H2O  Sulfate Aerosol Reaction
!-----------------------------------------------------------------------
               if( hclvmr > 0._dp .and. hobrvmr > 0._dp ) then
!-----------------------------------------------------------------------
!   	... Radius sulfate set (from sad module)
!           Set min radius to 0.001 microns (1e-7 cm)
!           Typical radius is 0.1 microns (1e-5 cm)
!           f_hobr may go negative under if not set.
!-----------------------------------------------------------------------
                  C_hobr          = 1477._dp*T_limitsq
                  D_hobr          = 9.e-9_dp
!-----------------------------------------------------------------------
!     	...  Taken from Waschewsky and Abbat
!            Dave Hanson (PC) suggested we divide this rc by eight to agee
!            with his data (Hanson, in press, 2002).
!            k1=k2*Mhcl for gamma(HOBr)
!-----------------------------------------------------------------------
                  k_wasch         = .125_dp * exp( .542_dp*wt - 6440._dp*T_limiti + 10.3_dp)
!-----------------------------------------------------------------------
!     	... Taken from Hanson 2002.
!-----------------------------------------------------------------------
                  H_hobr          = exp( -9.86_dp + 5427._dp*T_limiti )
                  k_dl            = 7.5e14_dp*D_hobr*2._dp                        ! or  7.5e14*D *(2nm)
!-----------------------------------------------------------------------
!  	... If k_wasch is GE than the diffusion limit...
!-----------------------------------------------------------------------
		  if( M_hcl > 0._dp ) then
                     if( k_wasch >= k_dl ) then
                        k_hobr_hcl   = k_dl * M_hcl
                     else
                        k_hobr_hcl   = k_wasch * M_hcl
                     end if
                     term1           = 4._dp*H_hobr*.082_dp*T_limit
                     term2           = sqrt( D_hobr*k_hobr_hcl )
                     tmp = 2.e2_dp
                     if (term2 .gt. zeps) tmp = rad_sulf/term2
                     Gamma_hobr_rxn  = term1*term2/C_hobr
                     rdl_hobr        = sqrt( D_hobr/k_hobr_hcl )
		     if( tmp < 1.e2_dp ) then
                        term1           = 1._dp/tanh( rad_sulf/rdl_hobr )
		     else
                        term1           = 1._dp
		     end if
                     term2           = rdl_hobr/rad_sulf
                     f_hobr          = term1 - term2
                     if( f_hobr > 0._dp .and. Gamma_hobr_rxn > 0._dp ) then
                        term1            = 1._dp / (f_hobr*Gamma_hobr_rxn)
                        gprob_hobr_hcl   = 1._dp / (1._dp + term1)
                     else
                         gprob_hobr_hcl  = 0._dp
                     end if
                     rxt(i,k,rid_het(6)) = max( 0._dp,wrk*av_hobr*gprob_hobr_hcl )*hcldeni
		  end if
               end if
            end if has_rad_sad_sulf

has_rad_nad_sulf : &
	    if( rad_nad_sulf > 0._dp ) then
	       wrk = .25_dp*rad_nad_sulf
!-----------------------------------------------------------------------
!     	... N2O5 + H2O(s) => 2HNO3  NAD Aerosol Reaction
!-----------------------------------------------------------------------
               if( n2o5vmr > 0._dp ) then
!-----------------------------------------------------------------------
!     ... gprob based on JPL00-03 for NAT.
!         also see Hanson and Ravi, JPC, 97, 2802-2803, 1993.
!                 gprob_tot     = 4.e-4
!-----------------------------------------------------------------------
                  rxt(i,k,rid_het(7))  = wrk*av_n2o5*4.e-4_dp
	       end if

!-----------------------------------------------------------------------
!     ClONO2 + H2O(s) => HNO3 + HOCl  NAD Aerosol Reaction
!-----------------------------------------------------------------------
               if( clono2vmr > 0._dp ) then
!-----------------------------------------------------------------------
!     ... gprob based on JPL00-03 for NAT.
!         also see Hanson and Ravi, JPC, 97, 2802-2803, 1993.
!                 gprob_tot    = 0.004
!-----------------------------------------------------------------------
                  rxt(i,k,rid_het(8)) = wrk*av_clono2*4.e-3_dp
   	       end if

!-----------------------------------------------------------------------
!     	... ClONO2 + HCl(s) => HNO3 + Cl2, NAD Aerosol Reaction
!-----------------------------------------------------------------------
               if( hclvmr > 0._dp ) then
                  if( clono2vmr > 0._dp ) then
!-----------------------------------------------------------------------
!     ... gprob based on JPL00-03 for NAT.
!         also see Hanson and Ravi, JPC, 96, 2682-2691, 1992.
!                 gprob_tot   = 0.2
!-----------------------------------------------------------------------
                     rxt(i,k,rid_het(9)) = wrk*av_clono2*.2_dp*hcldeni
                  end if

!-----------------------------------------------------------------------
!     	... HOCl + HCl(s) => H2O + Cl2  NAD Aerosol Reaction
!-----------------------------------------------------------------------
                  if( hoclvmr > 0._dp ) then
!-----------------------------------------------------------------------
!     ... gprob based on JPL00-03 for NAT.
!         see Hanson and Ravi, JPC, 96, 2682-2691, 1992.
!         and      Abbatt and Molina, GRL, 19, 461-464, 1992.
!                 gprob_tot   = 0.1
!-----------------------------------------------------------------------
                     rxt(i,k,rid_het(10)) = wrk*av_hocl*.1_dp*hcldeni
                  end if
               end if

!-----------------------------------------------------------------------
!     	... BrONO2 + H2O(s) => HOBr + HNO3  NAD Aerosol Reaction
!-----------------------------------------------------------------------
               if( brono2vmr > 0._dp ) then
!-----------------------------------------------------------------------
!       ... Personel Communication, 11/4/99, David Hanson
!                 gprob_tot   = 0.3
!-----------------------------------------------------------------------
                  rxt(i,k,rid_het(11)) = wrk*av_brono2*.3_dp
               end if
            end if has_rad_nad_sulf

has_sadice : &
	    if( sadice > 0._dp ) then
	       wrk = .25_dp*sadice
!-----------------------------------------------------------------------
!     N2O5 + H2O(s) => 2HNO3  ICE Aerosol Reaction
!-----------------------------------------------------------------------
               if( n2o5vmr > 0._dp ) then
!-----------------------------------------------------------------------
!       ... gprob based on JPL00-03 for ICE.
!           also see Hanson and Ravi, JPC, 97, 2802-2803, 1993.
!                 gprob_tot    = .02
!-----------------------------------------------------------------------
                  rxt(i,k,rid_het(12)) = wrk*av_n2o5*.02_dp
 	       end if
!-----------------------------------------------------------------------
!     	... ClONO2 + H2O(s) => HNO3 + HOCl  ICE Aerosol Reaction
!-----------------------------------------------------------------------
               if( clono2vmr > 0._dp ) then
!-----------------------------------------------------------------------
!     	... gprob based on JPL00-03 for ICE.
!     	    also see Hanson and Ravi, JGR, 96, 17307-17314, 1991.
!                 gprob_tot    = .3
!-----------------------------------------------------------------------
                  rxt(i,k,rid_het(13)) = wrk*av_clono2*.3_dp
	       end if

!-----------------------------------------------------------------------
!     	... BrONO2 + H2O(s) => HNO3 + HOBr  ICE Aerosol Reaction
!-----------------------------------------------------------------------
               if( brono2vmr > 0._dp ) then
!-----------------------------------------------------------------------
!     	... gprob based on JPL00-03 for ICE.
!           also see Hanson and Ravi, JPC, 97, 2802-2803, 1993.
!           could be as high as 1.0
!                 gprob_tot    = .3
!-----------------------------------------------------------------------
                  rxt(i,k,rid_het(14)) = wrk*av_brono2*.3_dp
      	       end if

!-----------------------------------------------------------------------
!     ClONO2 + HCl(s) => HNO3 + Cl2, ICE Aerosol Reaction
!-----------------------------------------------------------------------
               if( hclvmr > 0._dp ) then
                  if( clono2vmr > 0._dp ) then
!-----------------------------------------------------------------------
!       ... gprob based on JPL00-03 for ICE.
!           also see Hanson and Ravi, GRL, 15, 17-20, 1988.
!           also see Lue et al.,
!                 gprob_tot    = .3
!-----------------------------------------------------------------------
                     rxt(i,k,rid_het(15)) = wrk*av_clono2*.3_dp*hcldeni
                  end if
!
!-----------------------------------------------------------------------
!     	... HOCl + HCl(s) => H2O + Cl2, ICE Aerosol Reaction
!-----------------------------------------------------------------------
                  if( hoclvmr > 0._dp .and. hclvmr > 0._dp ) then
!-----------------------------------------------------------------------
!       ... gprob based on JPL00-003 for ICE.
!           also see Hanson and Ravi, JPC, 96, 2682-2691, 1992.
!           also see Abbatt and Molina, GRL, 19, 461-464, 1992.
!                 gprob_tot   = .2
!-----------------------------------------------------------------------
                     rxt(i,k,rid_het(16)) = wrk*av_hocl*.2_dp*hcldeni
                  end if

!-----------------------------------------------------------------------
!     HOBr + HCl(s) => H2O + BrCl, ICE Aerosol Reaction
!-----------------------------------------------------------------------
                  if( hobrvmr > 0._dp .and. hclvmr > 0._dp ) then
!-----------------------------------------------------------------------
!       ... gprob based on JPL00-03 for ICE.
!           Abbatt GRL, 21, 665-668, 1994.
!                    gprob_tot   = .3
!-----------------------------------------------------------------------
                    rxt(i,k,rid_het(17)) = wrk*av_hobr*.3_dp*hcldeni
                  end if
	       end if
            end if has_sadice
         end do column_loop
      end do Level_loop

      end subroutine ratecon_sfstrat


      subroutine zero_sfstrat ( rxt )
     
      use mo_moz_mods, only : ncol => plonl, pver => plev, rxntot

      real(dp), intent(inout) :: rxt(ncol,pver,rxntot)  ! rate constants

      integer :: i

!-----------------------------------------------------------------------
!      .. set het. reaction rates to zero (lstrathet = false)
!-----------------------------------------------------------------------

      do i = 1, nrid_het
         if ( rid_het(i) > 0 )  rxt(:,:,rid_het(i))  = 0._dp
      end do

      end subroutine zero_sfstrat

      end module mo_moz_het_strato_rates
