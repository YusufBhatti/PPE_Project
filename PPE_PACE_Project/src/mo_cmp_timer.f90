!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_cmp_timer.f90
!!
!! \brief
!! implementation of timers specific to the cloud microphysics
!!
!! \author Ulrike Proske (ETHZ)
!!
!! \responsible_coder
!!
!! \revision_history
!! This code is based on a cp of mo_hammoz_timer.f90
!!
!! \limitations
!! None
!!
!! \details
!! None
!!
!! \bibliographic_references
!! None
!!
!! \belongs_to
!!  
!!
!! \copyright
!! Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!! licencing agreement to be found at:
!! https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!! The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_cmp_timer

  !---inherited functions, types and data
  USE mo_real_timer, ONLY: new_timer, timer_report, &
                           timer_start, timer_stop, &
                           timer_reset_all, timer_reset_top

  IMPLICIT NONE

  !---public member functions
  PUBLIC :: init_cmp_timers 

  !---allow users of this module to USE these functions from mo_real_timer:
  PUBLIC :: timer_start, timer_stop

  !---high-level (overall process) timers 
  INTEGER, PUBLIC :: timer_cmp_bulk
  INTEGER, PUBLIC :: timer_cmp_0inits
  INTEGER, PUBLIC :: timer_cmp_1airdens
  INTEGER, PUBLIC :: timer_cmp_1cbasetop
  INTEGER, PUBLIC :: timer_cmp_2zerolocaltend
  INTEGER, PUBLIC :: timer_cmp_3melting
  INTEGER, PUBLIC :: timer_cmp_3subandevap
  INTEGER, PUBLIC :: timer_cmp_4cicesed
  INTEGER, PUBLIC :: timer_cmp_4incloudwi
  INTEGER, PUBLIC :: timer_cmp_5conddepevapsubl
  INTEGER, PUBLIC :: timer_cmp_5clearairevap
  INTEGER, PUBLIC :: timer_cmp_5incloudconddepevapsubl
  INTEGER, PUBLIC :: timer_cmp_6freezingTs238
  INTEGER, PUBLIC :: timer_cmp_6freezingTl238
  INTEGER, PUBLIC :: timer_cmp_7cloudphyscandsurfprecip
  INTEGER, PUBLIC :: timer_cmp_7warmclouds
  INTEGER, PUBLIC :: timer_cmp_7coldclouds
  !>>UP #777
  ! cold clouds sub-timers (within precipitation_formation_cold)
  ! aggregation
  INTEGER, PUBLIC :: timer_cmp_7ccpfcaggr
  ! riming
  INTEGER, PUBLIC :: timer_cmp_7ccpfcrime
  ! accretion between snow and ice
  INTEGER, PUBLIC :: timer_cmp_7ccpfcaccr
  ! secondary ice production
  INTEGER, PUBLIC :: timer_cmp_7ccpfcsepr
  ! self-collection of ice
  INTEGER, PUBLIC :: timer_cmp_7ccpfcsci
  ! finish up, belongs to multiple proceses
  INTEGER, PUBLIC :: timer_cmp_7ccpfcfinish
  ! whole subroutine inside
  INTEGER, PUBLIC :: timer_cmp_7ccpfcinside
  ! whole subroutine outside
  INTEGER, PUBLIC :: timer_cmp_7ccpfcoutside
  ! whole CMP subroutine inside
  INTEGER, PUBLIC :: timer_cmp_inside
  ! CMP subroutine inside, but only counting what we consider "real" CMPs
  INTEGER, PUBLIC :: timer_cmp_cmp
  !<<UP #777
  INTEGER, PUBLIC :: timer_cmp_7precipflux
  INTEGER, PUBLIC :: timer_cmp_8updating
  INTEGER, PUBLIC :: timer_cmp_9wetchem
  INTEGER, PUBLIC :: timer_cmp_diagnostics
  INTEGER, PUBLIC :: timer_cmp_icenucl
  INTEGER, PUBLIC :: timer_ccn_activation
  INTEGER, PUBLIC :: timer_wetdep_strat
  INTEGER, PUBLIC :: timer_wetdep_conv

CONTAINS
  
  SUBROUTINE init_cmp_timers

    ! submodel timers (process level)
    timer_cmp_bulk                     = new_timer('cmp_bulk')
    timer_cmp_0inits                   = new_timer('cmp_0initialisations')
    timer_cmp_1airdens                 = new_timer('cmp_1airdensity')
    timer_cmp_1cbasetop                = new_timer('cmp_1calc-cloud-base-and-top')
    timer_cmp_2zerolocaltend           = new_timer('cmp_2set-to-zero-local-tendencies')
    timer_cmp_3melting                 = new_timer('cmp_3melting')
    timer_cmp_3subandevap              = new_timer('cmp_3snow-ice-sub_rain-evap')
    timer_cmp_4cicesed                 = new_timer('cmp_4cloudice-sed')
    timer_cmp_4incloudwi               = new_timer('cmp_4incloud-water-ice')
    timer_cmp_5conddepevapsubl         = new_timer('cmp_5cond-dep-evap-subl')
    timer_cmp_5clearairevap            = new_timer('cmp_5clearair-evap')
    timer_cmp_5incloudconddepevapsubl  = new_timer('cmp_5incloud-cond-dep-evap-subl')
    timer_cmp_6freezingTs238           = new_timer('cmp_freezingTs238')
    timer_cmp_6freezingTl238           = new_timer('cmp_freezingTl238')
    timer_cmp_7cloudphyscandsurfprecip = new_timer('cmp_7tot')
    timer_cmp_7warmclouds              = new_timer('cmp_7warmclouds')
    timer_cmp_7coldclouds              = new_timer('cmp_7coldclouds')
    !>>UP #777
    timer_cmp_7ccpfcaggr               = new_timer('cmp_7ccpfcaggr')
    timer_cmp_7ccpfcrime               = new_timer('cmp_7ccpfcrime')
    timer_cmp_7ccpfcaccr               = new_timer('cmp_7ccpfcaccr')
    timer_cmp_7ccpfcsci                = new_timer('cmp_7ccpfcsci')
    timer_cmp_7ccpfcsepr               = new_timer('cmp_7ccpfcsepr')
    timer_cmp_7ccpfcfinish             = new_timer('cmp_7ccpfcfinish')
    timer_cmp_7ccpfcinside             = new_timer('cmp_7ccpfcinside')
    timer_cmp_7ccpfcoutside            = new_timer('cmp_7ccpfcoutside')
    timer_cmp_inside                   = new_timer('cmp_inside')
    timer_cmp_cmp                      = new_timer('cmp_cmp')
    !<<UP #777
    timer_cmp_7precipflux              = new_timer('cmp_7precipfluxes')
    timer_cmp_8updating                = new_timer('cmp_8updatetendencies')
    timer_cmp_9wetchem                 = new_timer('cmp_9wetchem')
    timer_cmp_diagnostics              = new_timer('cmp_diagnostics')
    timer_cmp_icenucl                  = new_timer('cmp_icenucleation')
    timer_ccn_activation               = new_timer('cmp_ccnactiv')
    timer_wetdep_strat                 = new_timer('cmp_wetdepstrat')
    timer_wetdep_conv                  = new_timer('cmp_wetdepconv')

  END SUBROUTINE init_cmp_timers
END MODULE mo_cmp_timer
