!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_timer.f90
!!
!! \brief
!! implementation of timers specific to the HAMMOZ submodel
!!
!! \author Declan O'Donnell (FMI)
!!
!! \responsible_coder
!! Declan O'Donnell, Declan.Odonnell@fmi.fi
!!
!! \revision_history
!!   -# Declan O'Donnell (FMI) - original code (2011-07)
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
!!  HAMMOZ
!!
!! \copyright
!! Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!! licencing agreement to be found at:
!! https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!! The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_hammoz_timer

  !---inherited functions, types and data
  USE mo_real_timer, ONLY: new_timer, timer_report, &
                           timer_start, timer_stop, &
                           timer_reset_all, timer_reset_top

  IMPLICIT NONE

  !---public member functions
  PUBLIC :: init_hammoz_timers 

  !---allow users of this module to USE these functions from mo_real_timer:
  PUBLIC :: timer_start, timer_stop

  !---high-level (overall process) timers 
  INTEGER, PUBLIC :: timer_ham_bulk !sf
  INTEGER, PUBLIC :: timer_ham_m7
  INTEGER, PUBLIC :: timer_ham_salsa !sf
  INTEGER, PUBLIC :: timer_ham_wetchem
  INTEGER, PUBLIC :: timer_ham_gaschem
  INTEGER, PUBLIC :: timer_ham_activation !dod 
  INTEGER, PUBLIC :: timer_moz_chem
  INTEGER, PUBLIC :: timer_moz_chem_prep    !new
  INTEGER, PUBLIC :: timer_moz_chem_diag    !new
  INTEGER, PUBLIC :: timer_moz_chem_photo   !new
  INTEGER, PUBLIC :: timer_moz_chem_impsol  !new
  INTEGER, PUBLIC :: timer_hammoz_emissions
  INTEGER, PUBLIC :: timer_hammoz_wetdep
  INTEGER, PUBLIC :: timer_hammoz_drydep
  INTEGER, PUBLIC :: timer_hammoz_sedimentation
  INTEGER, PUBLIC :: timer_hammoz_burden

  ! aerosol radiation effects timer
  INTEGER, PUBLIC :: timer_ham_rad
  INTEGER, PUBLIC :: timer_ham_rad_diag
  INTEGER, PUBLIC :: timer_ham_rad_fitplus
  INTEGER, PUBLIC :: timer_ham_rad_refrac

  !---detailed M7 timers
  INTEGER, PUBLIC :: timer_ham_m7_main
  INTEGER, PUBLIC :: timer_ham_m7_avg
  INTEGER, PUBLIC :: timer_ham_m7_hygro
  INTEGER, PUBLIC :: timer_ham_m7_cs
  INTEGER, PUBLIC :: timer_ham_m7_cond
  INTEGER, PUBLIC :: timer_ham_m7_nucl
  INTEGER, PUBLIC :: timer_ham_m7_dnum
  INTEGER, PUBLIC :: timer_ham_m7_coaset
  INTEGER, PUBLIC :: timer_ham_m7_delcoa
  INTEGER, PUBLIC :: timer_ham_m7_concoag
  INTEGER, PUBLIC :: timer_ham_m7_redistr

CONTAINS
  
  SUBROUTINE init_hammoz_timers

    ! submodel timers (process level)
    timer_ham_bulk          = new_timer('ham_bulk')
    timer_ham_m7            = new_timer('ham_m7')
    timer_ham_salsa         = new_timer('ham_salsa')
    timer_ham_wetchem       = new_timer('ham_wetchem')
    timer_ham_gaschem       = new_timer('ham_gaschem')
    timer_ham_activation    = new_timer('ham_activation') !dod
    timer_moz_chem          = new_timer('moz_chem')
    timer_moz_chem_prep     = new_timer('moz_chem_prep')
    timer_moz_chem_diag     = new_timer('moz_chem_diag')
    timer_moz_chem_photo    = new_timer('moz_chem_photo')
    timer_moz_chem_impsol   = new_timer('moz_chem_impsol')
    timer_hammoz_emissions  = new_timer('hammoz_emissions')
    timer_hammoz_wetdep     = new_timer('hammoz_wetdep')
    timer_hammoz_drydep     = new_timer('hammoz_drydep')
    timer_hammoz_sedimentation  = new_timer('hammoz_sedimentation')
    timer_hammoz_burden     = new_timer('hammoz_burden')
    timer_ham_rad           = new_timer('ham_rad')

    ! Aerosol direct radiative effect timers
    timer_ham_rad_diag      = new_timer('ham_rad_diag')
    timer_ham_rad_fitplus   = new_timer(' ham_rad_fitplus')
    timer_ham_rad_refrac    = new_timer(' ham_rad_refrac')

    !---M7 timers
    timer_ham_m7_main       = new_timer('ham_m7_main')
    timer_ham_m7_avg        = new_timer('ham_m7_averages')
    timer_ham_m7_hygro      = new_timer('ham_m7_hygroscopicity')
    timer_ham_m7_cs         = new_timer('ham_m7_cond_sink')
    timer_ham_m7_cond       = new_timer('ham_m7_H2SO4_cond')
    timer_ham_m7_nucl       = new_timer('ham_m7_nucleation')
    timer_ham_m7_dnum       = new_timer('ham_m7_dnum')
    timer_ham_m7_coaset     = new_timer('ham_m7_coaset')
    timer_ham_m7_delcoa     = new_timer('ham_m7_delcoa')
    timer_ham_m7_concoag    = new_timer('ham_m7_concoag')
    timer_ham_m7_redistr    = new_timer('ham_m7_redistr')

  END SUBROUTINE init_hammoz_timers
END MODULE mo_hammoz_timer
