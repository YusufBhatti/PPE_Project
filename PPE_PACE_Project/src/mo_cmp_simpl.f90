MODULE mo_cmp_simpl
! Based on mo_ccnclim
! Author: UP
! Date: 07 2020

  USE mo_kind,                     ONLY: dp
  USE mo_external_field_processor, ONLY: EF_FILE, EF_LATLEV, EF_IGNOREYEAR, &
                                         EF_NOINTER, EF_CONSTANT, EF_VALUE
  USE mo_boundary_condition,       ONLY: bc_nml, bc_define

  IMPLICIT NONE

  PUBLIC :: init_cmp_simpl

  INTEGER, PUBLIC :: ibc_cmpsimpl_rime
  INTEGER, PUBLIC :: ibc_cmpsimpl_icnucl
  INTEGER, PUBLIC :: ibc_cmpsimpl_icaccr
  INTEGER, PUBLIC :: ibc_cmpsimpl_sci
  INTEGER, PUBLIC :: ibc_cmpsimpl_subfi
  INTEGER, PUBLIC :: ibc_cmpsimpl_subfs
  INTEGER, PUBLIC :: ibc_cmpsimpl_mlt

  TYPE(bc_nml) :: bc_cmpsimpl_rime
  TYPE(bc_nml) :: bc_cmpsimpl_icnucl
  TYPE(bc_nml) :: bc_cmpsimpl_icaccr
  TYPE(bc_nml) :: bc_cmpsimpl_sci
  TYPE(bc_nml) :: bc_cmpsimpl_subfi
  TYPE(bc_nml) :: bc_cmpsimpl_subfs
  TYPE(bc_nml) :: bc_cmpsimpl_mlt

  REAL(dp), PARAMETER :: cmpsimpl_rime_const = 1.1e-5_dp ! [kg kg-1] Globally constant rime value
                                                     ! when relevant (ncmpsimpl_prescr_rime = 1)
  REAL(dp), PARAMETER :: cmpsimpl_icnucl_const = 3.7e4_dp ! [m-3] Globally constant ICnucl value
                                                     ! when relevant (ncmpsimpl_prescr_icnucl = 1)
  REAL(dp), PARAMETER :: cmpsimpl_icaccr_const = 1.7e-6_dp ! [kg kg-1] Globally constant ICaccr value
                                                     ! when relevant (ncmpsimpl_prescr_icaccr = 1)
  REAL(dp), PARAMETER :: cmpsimpl_sci_const = 71._dp ! [m-3] Globally constant sci value
                                                     ! when relevant (ncmpsimpl_prescr_sci = 1)
  REAL(dp), PARAMETER :: cmpsimpl_subfi_const = 1.1e-6_dp ! [kg kg-1] Globally constant subfi value
                                                     ! when relevant (ncmpsimpl_prescr_subfi = 1)
  REAL(dp), PARAMETER :: cmpsimpl_subfs_const = 2.6e-6_dp ! [kg kg-1] Globally constant subfs value
                                                     ! when relevant (ncmpsimpl_prescr_subfs = 1)
  REAL(dp), PARAMETER :: cmpsimpl_mlt_const = 8e-3_dp ! [kg m-2 s-1] Globally constant mlt value
                                                     ! when relevant (ncmpsimpl_prescr_mlt = 1)

CONTAINS

  SUBROUTINE init_cmp_simpl
      ! #821
      ! This subroutine sets the boundary conditions when using prescribed
      ! quantities for riming, ice crystal nucleation and ice crystal accretion

      USE mo_param_switches,  ONLY: ncmpsimpl_prescr_rime,   &
                                    ncmpsimpl_prescr_icnucl, &
                                    ncmpsimpl_prescr_icaccr, &
                                    ncmpsimpl_prescr_sci,    &
                                    ncmpsimpl_prescr_subfis, &
                                    ncmpsimpl_prescr_mlt

      !-- Riming
      SELECT CASE(ncmpsimpl_prescr_rime)
          CASE(1)
              bc_cmpsimpl_rime%ef_type        = EF_VALUE
              bc_cmpsimpl_rime%ef_value       = cmpsimpl_rime_const
              ibc_cmpsimpl_rime = bc_define('Riming prescribed input', &
                                            bc_cmpsimpl_rime, &
                                            ! ndims
                                            3)!, &
                                            !.TRUE. ) !lverbose
          CASE(2)
              bc_cmpsimpl_rime%ef_type        = EF_FILE
              bc_cmpsimpl_rime%ef_template    = 'rime.nc'
              bc_cmpsimpl_rime%ef_timedef     = EF_CONSTANT
              bc_cmpsimpl_rime%ef_timeindex   = 1
              bc_cmpsimpl_rime%ef_interpolate = EF_NOINTER
              bc_cmpsimpl_rime%ef_geometry    = EF_LATLEV
        !EF_LATLEV
              bc_cmpsimpl_rime%ef_actual_unit = 'kg kg-1'
              bc_cmpsimpl_rime%ef_varname        = 'rime'

              ibc_cmpsimpl_rime = bc_define('Riming prescribed input', &
                                            bc_cmpsimpl_rime, &
                                            ! ndims
                                            3)!, &
                                            !.TRUE. ) !lverbose
      END SELECT

      !-- IC nucleation
      SELECT CASE(ncmpsimpl_prescr_icnucl)
          CASE(1)
              bc_cmpsimpl_icnucl%ef_type        = EF_VALUE
              bc_cmpsimpl_icnucl%ef_value       = cmpsimpl_icnucl_const
              ibc_cmpsimpl_icnucl = bc_define('ICnucl prescribed input', &
                                            bc_cmpsimpl_icnucl, &
                                            ! ndims
                                            3)!, &
                                            !.TRUE. ) !lverbose
          CASE(2)
              bc_cmpsimpl_icnucl%ef_type        = EF_FILE
              bc_cmpsimpl_icnucl%ef_template    = 'icnucl.nc'
              bc_cmpsimpl_icnucl%ef_timedef     = EF_CONSTANT
              bc_cmpsimpl_icnucl%ef_timeindex   = 1
              bc_cmpsimpl_icnucl%ef_interpolate = EF_NOINTER
              bc_cmpsimpl_icnucl%ef_geometry    = EF_LATLEV
        !EF_LATLEV
              bc_cmpsimpl_icnucl%ef_actual_unit = 'm-3'
              bc_cmpsimpl_icnucl%ef_varname        = 'icnucl'

              ibc_cmpsimpl_icnucl = bc_define('ICnucl prescribed input', &
                                            bc_cmpsimpl_icnucl, &
                                            ! ndims
                                            3)!, &
                                            !.TRUE. ) !lverbose
      END SELECT

      !-- IC accretion
      SELECT CASE(ncmpsimpl_prescr_icaccr)
          CASE(1)
              bc_cmpsimpl_icaccr%ef_type        = EF_VALUE
              bc_cmpsimpl_icaccr%ef_value       = cmpsimpl_icaccr_const
              ibc_cmpsimpl_icaccr = bc_define('ICaccr prescribed input', &
                                            bc_cmpsimpl_icaccr, &
                                            ! ndims
                                            3)!, &
                                            !.TRUE. ) !lverbose
          CASE(2)
              bc_cmpsimpl_icaccr%ef_type        = EF_FILE
              bc_cmpsimpl_icaccr%ef_template    = 'icaccr.nc'
              bc_cmpsimpl_icaccr%ef_timedef     = EF_CONSTANT
              bc_cmpsimpl_icaccr%ef_timeindex   = 1
              bc_cmpsimpl_icaccr%ef_interpolate = EF_NOINTER
              bc_cmpsimpl_icaccr%ef_geometry    = EF_LATLEV
        !EF_LATLEV
              bc_cmpsimpl_icaccr%ef_actual_unit = 'kg kg-1'
              bc_cmpsimpl_icaccr%ef_varname        = 'icaccr'

              ibc_cmpsimpl_icaccr = bc_define('ICaccr prescribed input', &
                                            bc_cmpsimpl_icaccr, &
                                            ! ndims
                                            3)!, &
                                            !.TRUE. ) !lverbose
      END SELECT

      !-- Self-collection of ice
      SELECT CASE(ncmpsimpl_prescr_sci)
          CASE(1)
              bc_cmpsimpl_sci%ef_type        = EF_VALUE
              bc_cmpsimpl_sci%ef_value       = cmpsimpl_sci_const
              ibc_cmpsimpl_sci = bc_define('SCI prescribed input', &
                                            bc_cmpsimpl_sci, &
                                            ! ndims
                                            3)!, &
                                            !.TRUE. ) !lverbose
          CASE(2)
              bc_cmpsimpl_sci%ef_type        = EF_FILE
              bc_cmpsimpl_sci%ef_template    = 'sci.nc'
              bc_cmpsimpl_sci%ef_timedef     = EF_CONSTANT
              bc_cmpsimpl_sci%ef_timeindex   = 1
              bc_cmpsimpl_sci%ef_interpolate = EF_NOINTER
              bc_cmpsimpl_sci%ef_geometry    = EF_LATLEV
        !EF_LATLEV
              bc_cmpsimpl_sci%ef_actual_unit = 'kg kg-1'
              bc_cmpsimpl_sci%ef_varname        = 'sci'

              ibc_cmpsimpl_sci = bc_define('SCI prescribed input', &
                                            bc_cmpsimpl_sci, &
                                            ! ndims
                                            3)!, &
                                            !.TRUE. ) !lverbose
      END SELECT

      !-- Sublimation of falling/sedimenting ice
      SELECT CASE(ncmpsimpl_prescr_subfis)
          CASE(1)
              bc_cmpsimpl_subfi%ef_type        = EF_VALUE
              bc_cmpsimpl_subfi%ef_value       = cmpsimpl_subfi_const
              ibc_cmpsimpl_subfi = bc_define('Subfi prescribed input', &
                                            bc_cmpsimpl_subfi, &
                                            ! ndims
                                            3)!, &
                                            !.TRUE. ) !lverbose
          CASE(2)
              bc_cmpsimpl_subfi%ef_type        = EF_FILE
              bc_cmpsimpl_subfi%ef_template    = 'subfi.nc'
              bc_cmpsimpl_subfi%ef_timedef     = EF_CONSTANT
              bc_cmpsimpl_subfi%ef_timeindex   = 1
              bc_cmpsimpl_subfi%ef_interpolate = EF_NOINTER
              bc_cmpsimpl_subfi%ef_geometry    = EF_LATLEV
        !EF_LATLEV
              bc_cmpsimpl_subfi%ef_actual_unit = 'kg kg-1'
              bc_cmpsimpl_subfi%ef_varname        = 'subfi'

              ibc_cmpsimpl_subfi = bc_define('Subfi prescribed input', &
                                            bc_cmpsimpl_subfi, &
                                            ! ndims
                                            3)!, &
                                            !.TRUE. ) !lverbose
      END SELECT

      !-- Sublimation of falling/sedimenting snow
      SELECT CASE(ncmpsimpl_prescr_subfis)
          CASE(1)
              bc_cmpsimpl_subfs%ef_type        = EF_VALUE
              bc_cmpsimpl_subfs%ef_value       = cmpsimpl_subfs_const
              ibc_cmpsimpl_subfs = bc_define('Subfs prescribed input', &
                                            bc_cmpsimpl_subfs, &
                                            ! ndims
                                            3)!, &
                                            !.TRUE. ) !lverbose
          CASE(2)
              bc_cmpsimpl_subfs%ef_type        = EF_FILE
              bc_cmpsimpl_subfs%ef_template    = 'subfs.nc'
              bc_cmpsimpl_subfs%ef_timedef     = EF_CONSTANT
              bc_cmpsimpl_subfs%ef_timeindex   = 1
              bc_cmpsimpl_subfs%ef_interpolate = EF_NOINTER
              bc_cmpsimpl_subfs%ef_geometry    = EF_LATLEV
        !EF_LATLEV
              bc_cmpsimpl_subfs%ef_actual_unit = 'kg kg-1'
              bc_cmpsimpl_subfs%ef_varname        = 'subfs'

              ibc_cmpsimpl_subfs = bc_define('Subfs prescribed input', &
                                            bc_cmpsimpl_subfs, &
                                            ! ndims
                                            3)!, &
                                            !.TRUE. ) !lverbose
      END SELECT

      ! Melting of snow
      SELECT CASE(ncmpsimpl_prescr_mlt)
          CASE(1)
              bc_cmpsimpl_mlt%ef_type        = EF_VALUE
              bc_cmpsimpl_mlt%ef_value       = cmpsimpl_mlt_const
              ibc_cmpsimpl_mlt = bc_define('Mlt prescribed input', &
                                            bc_cmpsimpl_mlt, &
                                            ! ndims
                                            3)!, &
                                            !.TRUE. ) !lverbose
          CASE(2)
              bc_cmpsimpl_mlt%ef_type        = EF_FILE
              bc_cmpsimpl_mlt%ef_template    = 'mlt.nc'
              bc_cmpsimpl_mlt%ef_timedef     = EF_CONSTANT
              bc_cmpsimpl_mlt%ef_timeindex   = 1
              bc_cmpsimpl_mlt%ef_interpolate = EF_NOINTER
              bc_cmpsimpl_mlt%ef_geometry    = EF_LATLEV
        !EF_LATLEV
              bc_cmpsimpl_mlt%ef_actual_unit = 'kg m-2 s-1'
              bc_cmpsimpl_mlt%ef_varname        = 'mlt'

              ibc_cmpsimpl_mlt = bc_define('Mlt prescribed input', &
                                            bc_cmpsimpl_mlt, &
                                            ! ndims
                                            3)!, &
                                            !.TRUE. ) !lverbose
      END SELECT
  END SUBROUTINE init_cmp_simpl

END MODULE mo_cmp_simpl
