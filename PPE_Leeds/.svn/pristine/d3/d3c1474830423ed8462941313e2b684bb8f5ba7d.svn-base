#ifdef __xlC__
@PROCESS STRICT
#endif
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE iniphy

  ! Description:
  !
  ! Initialises physical constants of uncertain value.
  !
  ! Method:
  !
  ! This routine sets the values for the physical constants used
  ! in the parameterization routines (except for the radiation
  ! black-box) whenever these values are not well enough known to
  ! forbid any tuning or whenever they are subject to an arbitrary
  ! choice of the modeller. These constants will be in *comph2*.
  !
  ! *iniphy* is called from *physc*.
  !
  ! Authors:
  !
  ! J. F. Geleyn, ECMWF, December 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! M. Esch, MPI, June 1999, ECHAM5-modifications
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind,               ONLY: dp
  USE mo_control,            ONLY: nlev, nn, vct
  USE mo_physical_constants, ONLY: grav
  USE mo_physc2,             ONLY: clam, ckap, cb, cc, cd, cchar, cfreec     &
                                 , cgam, cvdifts, ctfreez, cz0ice, cqsncr    &
                                 , csncri, cevapcu, cwlmax
  USE mo_hyb,                ONLY: ceta
  USE mo_echam_conv_constants,       ONLY: cuparam
  USE mo_echam_cloud_params, ONLY: sucloud
  USE mo_surface_ice,        ONLY: init_surface_ice
  
  IMPLICIT NONE

  !  Local scalars: 
  INTEGER :: jk

  !  Executable statements 

!-- 1. Setting of constants

!-- 1.1 Constants for *vdiff*

  clam = 150._dp
  ckap = 0.4_dp
  cb = 5._dp
  cc = 5._dp
  cd = 5._dp
  cchar = 0.018_dp
  cfreec = 0.001_dp
  cgam = 1.25_dp
  cvdifts = 1.5_dp

!-- 1.2 Constants for *vdiff*, *clsst* and *atmice*

  ctfreez = 271.38_dp
  cz0ice = 0.001_dp

!-- 1.3 Constants for *vdiff* and *surf*

  cqsncr = 0.95_dp
  csncri = 5.85036E-3_dp

!-- 1.4 Constant for massflux convection scheme

  CALL cuparam

!-- 1.5 Highest level ncctop where condensation is allowed

  CALL sucloud ( nlev, vct )
!
  DO jk = 1, nlev
   cevapcu(jk) = 1.93E-6_dp*261._dp*SQRT(1.E3_dp/(38.3_dp*0.293_dp)*   &
                                              SQRT(ceta(jk)))*0.5_dp/grav
  END DO

!-- 1.6 Constants for *physc*

  cwlmax = 2.E-4_dp

!-- 1.7 Resolution dependend constants for surface_ice  

  CALL init_surface_ice(nn)

END SUBROUTINE iniphy
