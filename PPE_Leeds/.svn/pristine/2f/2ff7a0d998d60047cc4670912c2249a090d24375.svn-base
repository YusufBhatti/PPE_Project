!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!!  mo_moz_tropo_rates
!!
!! \brief
!!  Define "user-defined" reaction rates for MOZART chemistry scheme
!!
!! \author Stacy Walters (NCAR)
!!
!! \responsible_coder
!!  Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!  - SW: original version (pre 2004)
!!  - MGS (2013-07-25): revised rate expressions from JPL17(2011) and removed IBM vexp commands
!!
!! \limitations
!!  none
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


MODULE mo_moz_tropo_rates

  USE mo_kind,     ONLY : dp
  USE mo_moz_util, ONLY : get_rxt_ndx

  IMPLICIT NONE
      
  PRIVATE      

  PUBLIC :: tropo_rates_inti, tropo_rates

  ! reaction indices
  INTEGER ::  rid_O_O2, rid_O_O, rid_HO2_HO2, rid_NO2_NO3, rid_N2O5, &
              rid_HNO3_OH, rid_NO2_HO2_a, rid_NO2_HO2_b, rid_NO2_HO2, rid_HO2NO2, &
              rid_CLO_CLO_d, rid_CL2O2, &
              rid_CO_OH_a, rid_CH3CO3_NO2, rid_PAN, rid_CH3COCH3_OH, rid_MCO3_NO2, &
              rid_MPAN, rid_C59OOH_OH, rid_LC5PAN1719, rid_PBZNIT, rid_DMS_OH_b, &
              rid_CH3O2NO2, rid_CH3O2_NO2, rid_C2H5O2NO2, rid_C2H5O2_NO2, &
              rid_CH3CO_O2_a, rid_CH3CO_O2_b, rid_HOCH2CO_O2_a, rid_HOCH2CO_O2_b
  INTEGER ::  rid_het_NO3_tropo,  rid_het_N2O5_tropo, rid_het_HO2_tropo, &
              rid_het_HNO3_tropo, rid_het_O3_tropo, rid_het_NO2_tropo
  !Scarlet and Domenico
  INTEGER ::  rid_HO2_O3


     CONTAINS


!++mgs: the following cde was in mo_moz_chemdr before (didn't work properly though). 
!++mgs: Must re-think the use of sulfate_fld!!baustelle!!
!       Reason for failure might have been that mo_moz_sulf was nowhere used (read in etc.)
!! !-----------------------------------------------------------------------
!! ! get pointer to sulfate (climatology or from other model)  (MOZART2)
!! !-----------------------------------------------------------------------
!!       IF ( ichem == 2 ) THEN
!!          CALL get_tracer('sulfate', pxt=sulfate_fld, ierr=ierr)
!!          IF (ierr == 1) THEN
!!             CALL message ('moz_chemistry:', 'no sulfate tracer found, set sulfate vmr to 0')
!!             IF (ASSOCIATED(sulfate_fld)) DEALLOCATE(sulfate_fld)
!!             ALLOCATE(sulfate_fld(1:plonl,1:plev,lat:lat))
!!             sulfate_fld(:,:,:)=0._dp
!!          END IF
!!       ELSE      ! make sure sulfate_fld is defined (!!baustelle!!)
!! !!!         IF (ASSOCIATED(sulfate_fld)) DEALLOCATE(sulfate_fld)
!!          IF (.NOT. lallocsulfate) THEN
!!             ALLOCATE(sulfate_fld(1:plonl,1:plev,lat:lat))
!!             lallocsulfate = .true.
!!          END IF
!!          sulfate_fld(:,:,:)=0._dp
!!       ENDIF
!--mgs


  SUBROUTINE tropo_rates_inti ()

    rid_O_O2           = get_rxt_ndx('O_O2')
    rid_O_O            = get_rxt_ndx('O_O')
    rid_HO2_HO2        = get_rxt_ndx('HO2_HO2')
    rid_NO2_NO3        = get_rxt_ndx('NO2_NO3')
    rid_N2O5           = get_rxt_ndx('N2O5')
    rid_HNO3_OH        = get_rxt_ndx('HNO3_OH')
    rid_NO2_HO2_a      = get_rxt_ndx('NO2_HO2_a')
    rid_NO2_HO2_b      = get_rxt_ndx('NO2_HO2_b')
    rid_NO2_HO2        = get_rxt_ndx('NO2_HO2')
    rid_HO2NO2         = get_rxt_ndx('HO2NO2')
    rid_CH3O2NO2       = get_rxt_ndx('CH3O2NO2')
    rid_CH3O2_NO2      = get_rxt_ndx('CH3O2_NO2')
    rid_C2H5O2NO2      = get_rxt_ndx('C2H5O2NO2')
    rid_C2H5O2_NO2     = get_rxt_ndx('C2H5O2_NO2')
    rid_CLO_CLO_d      = get_rxt_ndx('CLO_CLO_d')
    rid_CL2O2          = get_rxt_ndx('CL2O2')
    rid_CO_OH_a        = get_rxt_ndx('CO_OH_a')
    rid_CH3CO3_NO2     = get_rxt_ndx('CH3CO3_NO2')
    rid_PAN            = get_rxt_ndx('PAN')
    rid_CH3COCH3_OH    = get_rxt_ndx('CH3COCH3_OH')
    rid_MCO3_NO2       = get_rxt_ndx('MCO3_NO2')
    rid_MPAN           = get_rxt_ndx('MPAN')
    rid_C59OOH_OH      = get_rxt_ndx('C59OOH_OH')
    rid_DMS_OH_b       = get_rxt_ndx('DMS_OH_b')
    rid_LC5PAN1719     = get_rxt_ndx('LC5PAN1719')
    rid_PBZNIT         = get_rxt_ndx('PBZNIT')
    rid_CH3CO_O2_a     = get_rxt_ndx('CH3CO_O2_a')
    rid_CH3CO_O2_b     = get_rxt_ndx('CH3CO_O2_b')
    rid_HOCH2CO_O2_a   = get_rxt_ndx('HOCH2CO_O2_a')
    rid_HOCH2CO_O2_b   = get_rxt_ndx('HOCH2CO_O2_b')

    rid_het_NO3_tropo  = get_rxt_ndx('het_NO3_tropo')
    rid_het_NO2_tropo  = get_rxt_ndx('het_NO2_tropo')
    rid_het_N2O5_tropo = get_rxt_ndx('het_N2O5_tropo')
    rid_het_HO2_tropo  = get_rxt_ndx('het_HO2_tropo')
    rid_het_HNO3_tropo = get_rxt_ndx('het_HNO3_tropo')
    rid_het_O3_tropo   = get_rxt_ndx('het_O3_tropo')

    rid_HO2_O3         = get_rxt_ndx('HO2_O3') !Scarlet and Domenico
  END SUBROUTINE tropo_rates_inti


  SUBROUTINE tropo_rates( rxt, temp, prho, h2ovmr, invariants, pmid, m, lat, kproma )

!-----------------------------------------------------------------
!        ... set the user specified reaction rates
!-----------------------------------------------------------------

      USE mo_moz_mods,       ONLY : rxntot, plonl, plev
      USE mo_moz,            ONLY : ltrophet
      USE mo_math_constants, ONLY : pi
      USE mo_vphysc,         ONLY : vphysc
      USE mo_moz_eslookup,   ONLY : aqsat
!++mgs: added reference to sulfate_fld
!++ost: currently not available, introduce later?
!     USE mo_moz,           ONLY : sulfate => sulfate_fld
      USE mo_hammoz_het_tropo_rates, ONLY : ratecon_sftropo

      IMPLICIT NONE

!-----------------------------------------------------------------
!        ... dummy arguments
!-----------------------------------------------------------------
      REAL(dp), INTENT(in) :: temp(plonl,plev), &           ! temperature (K)
                              prho(plonl, plev), &          ! air density         
                              m(plonl,plev), &              ! total atm density (1/cm^3)
                              h2ovmr(plonl,plev), &         ! water vapor (mol/mol)
                              pmid(plonl,plev), &           ! midpoint pressure (Pa)
                              invariants(plonl,plev,10)     ! invariants density (1/cm^3)

!!baustelle!! sulfate set to zero in mo_moz_chemdr -- clean up and eliminate??
      REAL(dp), INTENT(inout) :: rxt(plonl,plev,rxntot)     ! gas phase rates
      INTEGER, INTENT(in)     :: lat, kproma                ! latitude index, kproma
      
!-----------------------------------------------------------------
!        ... local variables
!-----------------------------------------------------------------
      REAL(dp), PARAMETER :: boltz = 1.38044e-16_dp      ! erg / K
      REAL(dp), PARAMETER :: avo   = 6.023e23_dp         ! molecules/mole
!-----------------------------------------------------------------
!	... Density of sulfate aerosol
!-----------------------------------------------------------------
      REAL(dp), PARAMETER :: gam1 = 0.04_dp              ! N2O5+SUL ->2HNO3
      REAL(dp), PARAMETER :: wso4 = 98._dp
      REAL(dp), PARAMETER :: den  = 1.15_dp              ! each molecule of SO4(aer) density g/cm3
!-------------------------------------------------
! 	... Volume of sulfate particles
!           assuming mean rm 
!           continent  0.05_dpum  0.07_dpum  0.09_dpum
!           ocean      0.09_dpum  0.25_dpum  0.37_dpum
!                      0.16_dpum                  Blake JGR,7195, 1995
!-------------------------------------------------
      REAL(dp), PARAMETER :: rm1  = 0.16_dp*1.e-4_dp          ! mean radii in cm

      INTEGER  ::  i, k, ktop
      INTEGER  ::  nltrop(plonl)
      REAL(dp)     ::  amas, fare
      REAL(dp), DIMENSION(plonl) :: &
                   tp, &
                   tinv, &
                   ko, &
                   kinf, &
                   fc, &
                   xr, &                    ! factor to increase particle radii depending on rel hum
                   sur, &                   ! sulfate particle surface area (cm^2/cm^3)
                   term1, term2, term3, term4, &
                   exp_fac, &               ! for vector exponential (ibm only)
                   kwho2                    ! Equilibrium constant for HO2-water complex
      REAL(dp), DIMENSION(plonl) :: &
                   satq, &                  ! saturation specific humidity
                   satv                     ! saturation vapor pressure
      REAL(dp), DIMENSION(plonl,plev) :: &
                   wrk,          &          ! holding array
                   zk_het_no3,   &          ! HAMMOZ reaction rates
                   zk_het_no2,   &          ! HAMMOZ reaction rates
                   zk_het_n2o5,  &
                   zk_het_ho2,   &
                   zk_het_hno3,  &
                   zk_het_o3

      fare = 4._dp*pi*rm1*rm1              ! particle mean area(r=0.1_dpu) (cm2/cm3)
level_loop : &
      do k = 1,plev
         tinv(:) = 1._dp / temp(:,k)
         tp(:)   = 300._dp * tinv(:)

         kwho2(:) = 2.4e-25_dp * exp(4350._dp * tinv(:))   ! JPL 2011 Equilibrium constant for HO2-water complex

!-----------------------------------------------------------------
!	... o + o2 + m --> o3 + m
!-----------------------------------------------------------------
         if( rid_O_O2 > 0 ) then
            rxt(:,k,rid_O_O2) = 6.e-34_dp * tp(:)**2.4_dp
         end if

!-----------------------------------------------------------------
!	... o + o + m -> o2 + m
!-----------------------------------------------------------------
         if( rid_O_O > 0 ) then
            rxt(:,k,rid_O_O) = 4.23e-28_dp * temp(:,k)**(-2)
         end if

!-----------------------------------------------------------------
!	... ho2 + ho2 --> h2o2
!	Note: this rate involves the water vapor number density
!-----------------------------------------------------------------
         if( rid_HO2_HO2 > 0 ) then
            ko(:)   = 3.0e-13_dp * exp( 460._dp*tinv(:) )
            kinf(:) = 2.1e-33_dp * m(:,k) * exp( 920._dp*tinv(:) )
            fc(:)   = 1._dp + 1.4e-21_dp * m(:,k) * h2ovmr(:,k) * exp( 2200._dp*tinv(:) )
            rxt(:,k,rid_HO2_HO2) = (ko(:) + kinf(:)) * fc(:)
         end if

!-----------------------------------------------------------------
!       ... ho2 + no2 + m --> ho2no2 + m  a
!           ho2-h2o + no2 --> ho2no2 + m  b1  
!                             hono + o2   b2 (only to parameterize the results of Li et al. 2014)
!-----------------------------------------------------------------
         if( rid_NO2_HO2_a > 0 .AND. rid_NO2_HO2_b > 0) then
!           NO2 + HO2 -> HO2NO2          
            kinf(:)  = 2.9e-12_dp * (temp(:,k)/ 300._dp)**(-1.1_dp)
            ko  (:)  = 2.0e-31_dp * (temp(:,k)/ 300._dp)**(-3.4_dp)

            term1(:) = ko(:) / ( (kinf(:) / m(:,k)) )
            term2(:) = ko(:) / (1._dp + term1(:)) ! * m(:,k) is done by the solver because M appears explicitely in the reaction

            term1(:) = log10( term1(:) )
            term1(:) = 1._dp / (1._dp + term1(:)*term1(:))


            term3(:) = 1._dp /(kwho2(:) * m(:,k) * h2ovmr(:,k) + 1._dp) ! fraction of bare HO2 (1./(K[H2O]+1.))
            
            term4(:) = 1._dp  ! scaling factor with respect to JPL 2011 for newer/more accurate k determinations like the one by M. Rolletter (FZJ-IEK-8) being 1.6_dp/1.06_dp 

            rxt(:,k,rid_NO2_HO2_a) = term2(:) * (0.6_dp)**term1(:) * term3(:) * term4(:) 

!           NO2 + HO2-H2O -> HO2NO2 + H2O or HONO + O2 + H2O          

            term4(:) = 2.5_dp !! enhancement factor of k with respect to k_dry. Used value intermediate from data by Sander and Peterson 1984 being in the 2-3.77 range. Another estimate is ~ 2 (3.3_dp/1.6_dp) by M. Rolletter (FZJ-IEK-8)
            rxt(:,k,rid_NO2_HO2_b) = term2(:) * (0.6_dp)**term1(:) * (1_dp - term3(:)) * term4(:)
         end if

!-----------------------------------------------------------------
!	... n2o5 + m --> no2 + no3 + m
!-----------------------------------------------------------------
         if( rid_NO2_NO3 > 0 .AND. rid_N2O5 > 0 ) then
            rxt(:,k,rid_N2O5) = rxt(:,k,rid_NO2_NO3) * 3.333e26_dp * exp( -10991._dp*tinv(:) )
         end if

!-----------------------------------------------------------------
! 	... hno3 + oh --> no3 + h2o
!-----------------------------------------------------------------
         if( rid_HNO3_OH > 0 ) then
            ko(:) = m(:,k) * 6.5e-34_dp * exp( 1335._dp*tinv(:) )
            ko(:) = ko(:) / (1._dp + ko(:)/(2.7e-17_dp*exp( 2199._dp*tinv(:) )))
            rxt(:,k,rid_HNO3_OH)  = ko(:) + 2.4e-14_dp*exp( 460._dp*tinv(:) )
         end if

!-----------------------------------------------------------------
!       ... ho2no2 + m --> ho2 + no2 + m
!-----------------------------------------------------------------
         if( rid_NO2_HO2 > 0 .AND. rid_HO2NO2 > 0 ) then
            rxt(:,k,rid_HO2NO2)  = rxt(:,k,rid_NO2_HO2) * exp( -10900._dp*tinv(:) ) / 2.1e-27_dp
         end if
         if( rid_NO2_HO2_a > 0 .AND. rid_NO2_HO2_b > 0 .AND. rid_HO2NO2 > 0 ) then
            rxt(:,k,rid_HO2NO2)  = (rxt(:,k,rid_NO2_HO2_a) + rxt(:,k,rid_NO2_HO2_b)) * exp( -10900._dp*tinv(:) ) / 2.1e-27_dp
         end if


!-----------------------------------------------------------------
!       ... ch3o2no2 + m --> ch3o2 + no2 + m
!-----------------------------------------------------------------
         if( rid_CH3O2_NO2 > 0 .AND. rid_CH3O2NO2 > 0 ) then
            rxt(:,k,rid_CH3O2NO2)  = rxt(:,k,rid_CH3O2_NO2) * exp( -11234._dp*tinv(:) ) / 9.5e-29_dp
         end if

!-----------------------------------------------------------------
!       ... c2h5o2no2 + m --> c2h5o2 + no2 + m
!-----------------------------------------------------------------
         if( rid_C2H5O2_NO2 > 0 .AND. rid_C2H5O2NO2 > 0 ) then
            rxt(:,k,rid_C2H5O2NO2)  = rxt(:,k,rid_C2H5O2_NO2) * exp( -11234._dp*tinv(:) ) / 9.5e-29_dp *5. !5 times the one for CH3O2NO2 as estimated by Zabel et al(1989) at the tropopause 
         end if

!-----------------------------------------------------------------
! 	... cl2o2 + m -> 2*clo + m
!-----------------------------------------------------------------
         if( rid_CLO_CLO_d > 0 .AND. rid_CL2O2 > 0 ) then
            ko(:)            = 1.3e-27_dp * exp( 8744._dp* tinv(:) )
            rxt(:,k,rid_CL2O2) = rxt(:,k,rid_CLO_CLO_d)/ko(:)
         end if

!-----------------------------------------------------------------
!       ... co + oh + m --> co2 + h + m (second branch JPL06; pg2.2_dp; 2.10_dp)
!-----------------------------------------------------------------
         if( rid_CO_OH_a > 0 ) then
            kinf(:)  = 2.1e+09_dp * (temp(:,k)/ 300._dp)**(6.1_dp)
            ko  (:)  = 1.5e-13_dp * (temp(:,k)/ 300._dp)**(0.6_dp)

            term1(:) = ko(:) / ( (kinf(:) / m(:,k)) )
            term2(:) = ko(:) / (1._dp + term1(:))

            term1(:) = log10( term1(:) )
            term1(:) = 1._dp / (1._dp + term1(:)*term1(:))

            rxt(:,k,rid_CO_OH_a) = term2(:) * (0.6_dp)**term1(:)
         end if
!! old JPL97 expression:     rxt(:,k,rid_CO_OH_a) = 1.5e-13_dp * (1._dp + 6.e-7_dp*boltz*m(:,k)*temp(:,k))

!-----------------------------------------------------------------
!       ... mco3 + no2 -> mpan
!-----------------------------------------------------------------
! relict from previous MOZART version
!        if( rid_MCO3_NO2 > 0 ) then
!           rxt(:,k,rid_MCO3_NO2) = 1.1e-11_dp * tp(:) / m(:,k)
!        end if

!-----------------------------------------------------------------
!       ... pan + m --> ch3co3 + no2 + m
!       ... mpan + m --> mco3 + no2 + m
!-----------------------------------------------------------------
         if( rid_PAN > 0 .or. rid_MPAN > 0 ) then
            wrk(:,1) = 1.111e28_dp * exp( -14000._dp*tinv(:) )
            if( rid_PAN > 0 .AND. rid_CH3CO3_NO2 > 0 ) then
               rxt(:,k,rid_PAN) = rxt(:,k,rid_CH3CO3_NO2) * wrk(:,1)
            end if
            if( rid_MPAN > 0 .AND. rid_MCO3_NO2 > 0 ) then
               rxt(:,k,rid_MPAN) = rxt(:,k,rid_MCO3_NO2) * wrk(:,1)
            end if
         end if

!-----------------------------------------------------------------
!       ... ch3co + o2 -> ch3co3 
!                         oh + ch2o + co 
!-----------------------------------------------------------------
         if( rid_CH3CO_O2_a > 0 ) then
            rxt(:,k,rid_CH3CO_O2_a) = 5.1e-12_dp * (1._dp - 1._dp/(1._dp + 9.4e-18_dp * m(:,k))) ! IUPAC and  Groß et al. 2014. Rate constant is the high-pressure limit as recommended by IUPAC. 
         end if

         if( rid_CH3CO_O2_b > 0 ) then
            rxt(:,k,rid_CH3CO_O2_b) = 5.1e-12_dp * 1._dp/(1._dp + 9.4e-18_dp * m(:,k)) ! IUPAC and  Groß et al. 2014. Rate constant is the high-pressure limit as recommended by IUPAC. 
         end if

!-----------------------------------------------------------------
!       ... hoch2co + o2 -> hoch2co3 
!                           oh + ch2o + co2 
!-----------------------------------------------------------------
         if( rid_HOCH2CO_O2_a > 0 ) then
            rxt(:,k,rid_HOCH2CO_O2_a) = 5.1e-12_dp * (1._dp - 1._dp/(1._dp + 1.85e-18_dp * m(:,k))) ! IUPAC and  Groß et al. 2014. Rate constant is the high-pressure limit as recommended by IUPAC. 
         end if

         if( rid_HOCH2CO_O2_b > 0 ) then
            rxt(:,k,rid_HOCH2CO_O2_b) = 5.1e-12_dp * 1._dp/(1._dp + 1.85e-18_dp * m(:,k)) ! IUPAC and  Groß et al. 2014. Rate constant is the high-pressure limit as recommended by IUPAC. 
         end if

!-----------------------------------------------------------------
!       ... ch3coch3 + oh -> ro2 + h2o
!-----------------------------------------------------------------
         if( rid_CH3COCH3_OH > 0 ) then
            rxt(:,k,rid_CH3COCH3_OH) = 1.33e-13_dp + 3.82e-11_dp * exp( -2000._dp*tinv(:) )
         end if

!-----------------------------------------------------------------
!       ... xooh + oh -> h2o + oh
!-----------------------------------------------------------------
! relict from previous MOZART version
!        if( rid_C59OOH_OH > 0 ) then
!           rxt(:,k,rid_C59OOH_OH) = temp(:,k)**2 * 7.69e-17_dp * exp( 253._dp*tinv(:) )
!        end if

!-----------------------------------------------------------------
!       ... old reaction: DMS + OH  --> .5_dp * SO2
!       ... now:          DMS + OH  --> .5_dp * SO2 + .5_dp * HO2
!-----------------------------------------------------------------
! old code:
!        if( rid_DMS_OH_b > 0 ) then
!           ko(:) = 1._dp + 5.5e-31_dp * exp( 7460._dp*tinv(:) ) * m(:,k) * 0.21_dp
!           rxt(:,k,rid_DMS_OH_b) = 1.7e-42_dp * exp( 7810._dp*tinv(:) ) * m(:,k) * 0.21_dp / ko(:)
!        end if
! now:
         if( rid_DMS_OH_b > 0 ) then
             rxt(:,k,rid_DMS_OH_b) =  8.2e-39_dp * exp(5376._dp*tinv(:)) * 0.21_dp * m(:,k) &
                                      / (1._dp + 1.05e-5_dp * exp(3644._dp*tinv(:)) * 0.21_dp)
         end if
  
!-----------------------------------------------------------------
!       ... NH3 --> NH4
!-----------------------------------------------------------------
!         if( rid_usr25 > 0 )  then
!             rxt(:,k,rid_usr25) = 0._dp      !!baustelle!!  ???
!         end if

!-----------------------------------------------------------------
!       ... LC5PAN1719 -> LHC4ACCO3 + NO2
!-----------------------------------------------------------------
         if( rid_LC5PAN1719 > 0 ) then
            rxt(:,k,rid_LC5PAN1719) =  rxt(:,k,rid_PAN)
         end if
 
!-----------------------------------------------------------------
!       ... PBZNIT + M -> ACBZO2 + NO2 + M
!-----------------------------------------------------------------
         if( rid_PBZNIT > 0 ) then
            rxt(:,k,rid_PBZNIT) =  rxt(:,k,rid_PAN)
         end if

 !-----------------------------------------------------------------
 !      ...  HO2 + O3 -> OH + 2*O2 !Scarlet and Domenico
 !-----------------------------------------------------------------
         if( rid_HO2_O3 > 0 ) then
             rxt(:,k,rid_HO2_O3) =  2.03e-16_dp*exp(693.0_dp*tinv(:))*(temp(:,k)/300.0_dp)**4.57_dp
         end if

      end do level_loop


 !-----------------------------------------------------------------
 !  .. heterogeneous reactions ..
 !-----------------------------------------------------------------

      ! Call heterogeneous tropospheric chemistry (if ltrophet == .true)
      IF ( ltrophet ) THEN
        CALL ratecon_sftropo(temp, prho, h2ovmr, pmid, m, lat, kproma, &
                             zk_het_no3,zk_het_no2, zk_het_n2o5, &
                             zk_het_ho2, zk_het_hno3, zk_het_o3)

      ELSE
        zk_het_n2o5(:,:) = 0._dp
        zk_het_no3(:,:)  = 0._dp
        zk_het_no2(:,:)  = 0._dp
        zk_het_ho2(:,:)  = 0._dp
        zk_het_hno3(:,:)  = 0._dp
        zk_het_o3(:,:)  = 0._dp
      END IF

      IF (rid_het_N2O5_tropo > 0) rxt(:,:,rid_het_N2O5_tropo) = zk_het_n2o5(:,:)
      IF (rid_het_NO3_tropo > 0)  rxt(:,:,rid_het_NO3_tropo)  = zk_het_no3(:,:)
      IF (rid_het_NO2_tropo > 0)  rxt(:,:,rid_het_NO2_tropo)  = zk_het_no2(:,:)
      IF (rid_het_HO2_tropo > 0)  rxt(:,:,rid_het_HO2_tropo)  = zk_het_ho2(:,:)
      IF (rid_het_HNO3_tropo > 0) rxt(:,:,rid_het_HNO3_tropo) = zk_het_hno3(:,:)
      IF (rid_het_O3_tropo > 0)   rxt(:,:,rid_het_O3_tropo)   = zk_het_o3(:,:)

  END SUBROUTINE tropo_rates

END MODULE mo_moz_tropo_rates
